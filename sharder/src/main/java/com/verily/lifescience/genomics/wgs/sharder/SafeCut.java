/*
 * Copyright 2018 Verily Life Sciences LLC
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package com.verily.lifescience.genomics.wgs.sharder;

import static com.google.common.base.Preconditions.checkArgument;
import static com.google.common.collect.ImmutableList.toImmutableList;

import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.verily.lifescience.genomics.wgs.mindex.Position;
import java.io.IOException;
import java.nio.file.Path;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import org.checkerframework.checker.nullness.qual.EnsuresNonNull;
import org.checkerframework.checker.nullness.qual.MonotonicNonNull;
import org.checkerframework.checker.nullness.qual.Nullable;

/**
 * Finds a safe cut.
 *
 * <p>A "safe cut" is a genomic position such that cutting just before that position would not split
 * any of the deletions listed in any of the input VCFs. So for example if there's a deletion from
 * 10-20 then 10 is a safe cut, but 11 is not (and neither is 20).
 *
 * <p>It's safe to cut insertions because they have a single reference base, so all the inserted
 * bases will still always end up in matching shards. Consider for example two insertions that start
 * at chr1:1000. Those inserted bases end before chr:1001, so it's safe to include them in the "up
 * to chr:1001, excluded" shard. None of the information from the "chr:1001 and up" shard is
 * relevant to those inserts, it's OK for those to be separate.
 *
 * <p>See go/joint-genotyping-splitting-plan for details.
 */
public class SafeCut implements AutoCloseable {

  private final ImmutableList<Path> variantsPath;
  /** A reader for each input VCF. Set by init. */
  private VcfReader[] vcfs;
  /** Contents of the shards file. */
  private final ImmutableList<Position> positions;
  /** The contigs listed in the shards file. */
  private final ImmutableMap<String, Integer> contigs;
  /** Each input VCF's corresponding mindex. */
  private final ImmutableList<Mindex> mindexes;
  /** The reference for the vcfs. */
  private final Reference reference;
  /** True iff init() was called. */
  private boolean initialized;
  /** Turn on for extra info only relevant for debugging. */
  private boolean verbose;
  /** Where we're currently thinking of cutting. only null until init is called. */
  @MonotonicNonNull private Position tentativePos;

  private final int threads;

  /** Loads the shards file, stores other paths. */
  public SafeCut(Path shards, List<Path> mindexPaths, ImmutableList<Path> variantsPaths, int threads,
      Reference reference) throws IOException, ParseException {
    checkArgument(
        mindexPaths.size() == variantsPaths.size(),
        "Arrays must have the same length.");
    this.threads = threads;
    this.variantsPath = variantsPaths;
    this.reference = reference;
    contigs = Position.contigsFromFile(shards);
    positions = Position.fromFile(shards, contigs);
    ImmutableList.Builder<Mindex> mindexBuilder = ImmutableList.builder();
    for (Path p : mindexPaths) {
      mindexBuilder.add(new Mindex(p));
    }
    mindexes = mindexBuilder.build();
    initialized = false;
  }

  /** set verbose flag **/
  public SafeCut setVerbose(boolean verbose) {
    this.verbose = verbose;
    return this;
  }

  public int numShards() {
    return positions.size();
  }

  /** Returns the contigs in the shards file. */
  public ImmutableMap<String, Integer> contigs() {
    return this.contigs;
  }

  /**
   * Opens all the files, seek to the tentative cut position. Call before {@link #findSafeCut}.
   *
   * <p>You can call it multiple times, to find the start point of different shards (this way you
   * don't have to reload the shards file and mindex).
   *
   * <p>Can run multithreaded.
   *
   * @throws IllegalArgumentException if we detect that the mindex may be broken because it sends us
   *     beyond the target position instead of at or before.
   */
  @EnsuresNonNull("this.tentativePos")
  public void init(int shardNo) throws IOException, InterruptedException {
    if (null == this.vcfs) {
      vcfs = new VcfReader[variantsPath.size()];
    } else {
      if (verbose) {
        System.out.println("Reusing already-open VCFs");
      }
    }
    tentativePos = positions.get(shardNo);
    ExecutorService exec;
    // only do a limited number of tasks on the executor before we make a new one:
    // this works around a bug where we'd run out of thread-local space otherwise
    // (at around task #1000).
    int step = 250;
    for (int start=0; start < vcfs.length; start += step) {
      if (threads > 1) {
        // These threads mostly just block, waiting on IO - so allocate more
        // than we have processors. The workStealingPool works too, but not as well.
        // We're using daemon threads because we want the program to quit if the main thread exits,
        // as may happen in case of error.
        exec = Executors.newFixedThreadPool(32, new ThreadFactoryBuilder().setDaemon(true).build());
      } else {
        exec = Executors.newSingleThreadExecutor(new ThreadFactoryBuilder().setDaemon(true).build());
      }
      int end = Math.min(start + step, vcfs.length);
      List<Future<Boolean>> done = new ArrayList<>();
      for (int i = start; i < end; i++) {
        final int j = i;
        done.add(exec.submit(() -> initializeVcf(j, shardNo)));
      }
      // wait for all tasks to be done, check for exceptions.
      for (Future<Boolean> d : done) {
        try {
          d.get();
        } catch (ExecutionException e) {
          Throwable inner = e.getCause();
          if (inner instanceof IOException) {
            throw (IOException) inner;
          }
          throw new RuntimeException("Unexpected exception while initializing", e);
        }
      }
      exec.shutdown();
      exec.awaitTermination(10, TimeUnit.SECONDS);
    }
    initialized = true;
  }

  private boolean initializeVcf(int index, int shardNo) throws IOException {
    // TentativePos is always set in init, so non-null by the time we get here.
    Preconditions.<@Nullable Position>checkNotNull(tentativePos);
    checkArgument(
      index<variantsPath.size(), 
      "Cannot initialize file " + index + ", we only have " + variantsPath + " vcfs.");
    checkArgument(
      index<variantsPath.size(), 
      "Cannot initialize file " + index + ", we only have " + mindexes + " mindexes.");
    VcfReader reader = vcfs[index];
    if (null==reader) {
      reader = new VcfReader(variantsPath.get(index), contigs, reference);
      vcfs[index] = reader;
    }
    long offset = mindexes.get(index).get(shardNo);
    reader.offset(offset);
    reader.advanceTo(tentativePos);
    return true;
  }

  /**
   * Parses the files, advance until a safe cut position is found (will stay in place if that's
   * already a safe cut). Uses the number of threads specified in the ctor.
   */
  public Position findSafeCut() throws IOException, InterruptedException {
    ensureInitialized();
    int parallelism = Math.min(1, threads);
    int perWorker = (int) Math.ceil(vcfs.length / (double) parallelism);
    ImmutableList.Builder<SubsetCutter> builder = ImmutableList.builder();
    for (int start = 0; start < vcfs.length; start += perWorker) {
      builder.add(new SubsetCutter(start, start + perWorker));
    }
    ImmutableList<SubsetCutter> cutters = builder.build();
    Position considering = tentativePos;
    boolean change;
    ExecutorService exec =
        Executors.newFixedThreadPool(
            parallelism, new ThreadFactoryBuilder().setDaemon(true).build());
    do {
      change = false;
      Position initialPos = considering;

      ImmutableList.Builder<Future<Position>> cutBuilder = ImmutableList.builder();
      for (SubsetCutter cutter : cutters) {
        cutBuilder.add(exec.submit(() -> cutter.findSafeCut(initialPos)));
      }
      for (Future<Position> cut : cutBuilder.build()) {
        try {
          Position safe = cut.get();
          if (!safe.equals(initialPos)) {
            change = true;
          }
          if (considering.before(safe)) {
            considering = safe;
          }
        } catch (ExecutionException e) {
          Throwable inner = e.getCause();
          if (inner instanceof IOException) {
            throw (IOException) inner;
          }
          throw new RuntimeException("Unexpected exception while computing cut", e);
        }
      }
    } while (change);
    exec.shutdown();
    exec.awaitTermination(10, TimeUnit.SECONDS);
    return considering;
  }

  /** Returns the offset to the line before the safe cut. */
  public ImmutableList<Long> getPreviousOffsets() {
    ensureInitialized();
    return Arrays.asList(vcfs).stream().map(VcfReader::getPreviousOffset).collect(toImmutableList());
  }

  @Override
  public void close() throws IOException {
    for (VcfReader reader : vcfs) {
      if (null != reader) {
        reader.close();
      }
    }
  }

  /**
   * Returns true iff the record is a deletion. We identify them by the length of the "REF" column.
   */
  static boolean isDeletion(String vcfRow) {
    String[] elements = vcfRow.split("\t", 5);
    String ref = elements[3];
    // Only deletion records list more than one reference base.
    return ref.length() > 1;
  }

  /**
   * Throws if we're not initialized, and helps the nullness checked understand
   * that tentativePos is not null (it needs a bit of help).
   */
  @EnsuresNonNull("this.tentativePos")
  private void ensureInitialized() {
    checkArgument(initialized, "Call init first");
    // this should never fail, but it's here for the nullness checker.
    Preconditions.<@Nullable Position>checkNotNull(tentativePos);
  }

  /** Finds a safe cut for a subset of the inputs. */
  private class SubsetCutter {
    private final int start;
    private final int endExcl;

    /** Specifies the range of inputs to process. */
    SubsetCutter(int start, int endExcl) {
      this.start = start;
      this.endExcl = Math.min(endExcl, vcfs.length);
    }

    /** Finds a safe cut for a subset of the inputs. The safe cut is at or after tentativePos. */
    Position findSafeCut(Position tentativePos) throws IOException {
      Position initialPos;
      do {
        initialPos = tentativePos;
        for (int i=start; i<endExcl; i++) {
          VcfReader reader = vcfs[i];
          reader.advanceToAtLeast(tentativePos);
          if (!reader.isEof()) {
            Position actual = reader.getPosition();
            if (actual.pos() > tentativePos.pos() && isDeletion(reader.getPreviousRecord().get())) {
              // We had to move forward a bit, and the record before us isn't splittable.
              // tentativePos must advance. Our actual position is fine
              tentativePos = actual;
            }
          } else {
            Optional<String> previousRecord = reader.getPreviousRecord();
            // it's possible for no previous record to be present, for example if the start
            // position is already past EOF.
            if (previousRecord.isPresent() && isDeletion(previousRecord.get())) {
              // TODO: implement this case. Need to figure out the position
              // by parsing the deleted record, and update tentativePos.
              throw new RuntimeException("Unimplemented edge case: last record is a deletion");
            }
          }
        }
        // Repeat until fixed point.
      } while (!tentativePos.equals(initialPos));
      return tentativePos;
    }
  }
}
