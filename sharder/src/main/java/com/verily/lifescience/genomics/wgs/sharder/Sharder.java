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

import com.google.common.base.Stopwatch;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.verily.lifescience.genomics.wgs.mindex.Position;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.channels.SeekableByteChannel;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.text.ParseException;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import javax.annotation.Nullable;

/** Writes out one shard of a collection of GVCF files */
class Sharder {

  private final Path shards;
  private final ImmutableList<Path> mindexes;
  private final ImmutableList<Path> vcfPaths;
  private final int totalShards;
  private final ImmutableList<Path> outPaths;
  private final Path metricsPath;
  private final Stopwatch totalTime = Stopwatch.createUnstarted();
  private final Stopwatch initTime = Stopwatch.createUnstarted();
  private final Stopwatch writeTime = Stopwatch.createUnstarted();
  private final Reference reference;
  private ImmutableList<Long> beginOffsets;
  private ImmutableList<Long> endOffsets;
  private Position beginCut;
  @Nullable private Position endCut;
  private ImmutableMap<String, Integer> contigs;
  private boolean verbose;
  private boolean skipWriting;

  /** Sets up the sharder. */
  Sharder(
      Path shards,
      ImmutableList<Path> mindexes,

      ImmutableList<Path> vcfPaths,
      int totalShards,
      ImmutableList<Path> outPaths,
      Path metricsPath,
      Reference reference) {
    checkArgument(mindexes.size() == vcfPaths.size(), "Must have as many mindexes as vcfs");
    checkArgument(mindexes.size() == outPaths.size(), "Must have as many mindexes as output files");
    this.shards = shards;
    this.mindexes = mindexes;
    this.vcfPaths = vcfPaths;
    this.totalShards = totalShards;
    this.outPaths = outPaths;
    this.metricsPath = metricsPath;
    this.reference = reference;
  }

  public Sharder setVerbose(boolean verbose) {
    this.verbose = verbose;
    return this;
  }

  public Sharder skipWriting(boolean skipWriting) {
    this.skipWriting = skipWriting;
    return this;
  }

  /** Copies the specified shard, optionally multithreaded. */
  public void shard(int shard, int threads)
      throws IOException, ParseException, InterruptedException {
    checkArgument(
        shard >= 0 && shard < totalShards, "Shard number must be between 0 and shards_total-1");
    checkArgument(threads > 0, "Need at least 1 thread");
    checkOutputIsWriteable();
    initTime.reset();
    totalTime.reset();
    totalTime.start();
    System.out.printf("Starting on shard %d.\n", shard);

    try (SafeCut safeCut = new SafeCut(shards, mindexes, vcfPaths, threads, reference).setVerbose(verbose)) {

      int numShardsInFile = safeCut.numShards();
      checkArgument(
          totalShards <= numShardsInFile,
          "'shards_total' must be at most the number of rows in the shards file");
      checkArgument(
          numShardsInFile % totalShards == 0,
          "'shards_total' must be a divisor of the number of rows in the shards file. Got " + totalShards + ", there are " + numShardsInFile + " shards in the file.");
      int shardsAtATime = numShardsInFile / totalShards;

      // 1. Find safe begin cut point.
      initTime.start();
      if (verbose) {
        System.out.println("Computing first cut.");
      }
      safeCut.init(shard * shardsAtATime);
      initTime.stop();
      beginCut = safeCut.findSafeCut();
      if (verbose) {
        System.out.println("  First cut: " + beginCut);
      }
      beginOffsets = safeCut.getPreviousOffsets();

      // 2. Find safe end cut point.
      int endShard = (shard + 1) * shardsAtATime;
      if (endShard < totalShards * shardsAtATime) {
        // Normal case: find the safe cut.
        initTime.start();
        if (verbose) {
          System.out.println("computing second cut.");
        }
        safeCut.init(endShard);
        initTime.stop();
        endCut = safeCut.findSafeCut();
        endOffsets = safeCut.getPreviousOffsets();
      } else {
        if (verbose) {
          System.out.println("No second cut, we go until end of file.");
        }
        // Final shard: we'll copy until the end.
        ImmutableList.Builder<Long> builder = ImmutableList.builder();
        for (Path element : vcfPaths) {
          builder.add(Files.size(element));
        }
        endOffsets = builder.build();
        endCut = null;
      }

      System.out.printf(
          "Safe cut points found (%s - %s), cutting.\n", beginCut, (endCut == null ? "EOF" : endCut));

      if (!skipWriting) {
        // 3. Copy shard, between the cut points.
        writeTime.reset();
        writeTime.start();
        contigs = safeCut.contigs();
        safeCut.close();
        int perWorker = (int) Math.ceil(vcfPaths.size() / (double) threads);
        ExecutorService exec =
            Executors
                .newFixedThreadPool(threads, new ThreadFactoryBuilder().setDaemon(true).build());
        ImmutableList.Builder<Future<Boolean>> doneBuilder = ImmutableList.builder();
        for (int start = 0; start < vcfPaths.size(); start += perWorker) {
          int sharderStart = start;
          SubsetSharder worker = new SubsetSharder(sharderStart, sharderStart + perWorker);
          doneBuilder.add(exec.submit(worker::copy));
        }
        ImmutableList<Future<Boolean>> done = doneBuilder.build();
        for (Future<Boolean> doneYet : done) {
          try {
            doneYet.get();
          } catch (ExecutionException x) {
            if (x.getCause() instanceof IOException) {
              throw (IOException) x.getCause();
            }
            if (x.getCause() instanceof ParseException) {
              throw (ParseException) x.getCause();
            }
            throw new RuntimeException("Unexpected error", x);
          }
        }
        exec.shutdown();
        exec.awaitTermination(10, TimeUnit.SECONDS);
        writeTime.stop();
        System.out.printf("Done writing %d shards.\n", vcfPaths.size());
      } else {
        System.out.printf("Done. Not writing the %d shards to disk, as requested.\n", vcfPaths.size());
      }
      totalTime.stop();

      // empty metrics file name -> don't save metrics
      if (!metricsPath.toString().isEmpty()) {
        saveMetrics(shard, threads);
      } else {
        System.out.println("Metrics path not specified, skipping metrics.");
      }
    }
  }

  /**
   * Checks the outputs are writeable. The test files are deleted afterwards, so if there's an error
   * then no output files will be written.
   *
   * @throws IOException if the destination isn't writeable.
   */
  private void checkOutputIsWriteable() throws IOException {
    tryWrite(outPaths.get(0));
    if (!metricsPath.toString().isEmpty()) {
      tryWrite(metricsPath);
    }
  }

  private void tryWrite(Path dest) throws IOException {
    try (BufferedWriter writer =
        Files.newBufferedWriter(
            dest,
            StandardOpenOption.WRITE,
            StandardOpenOption.CREATE,
            StandardOpenOption.TRUNCATE_EXISTING)) {
      writer.write("Testing\n");
    }
    Files.delete(dest);
  }

  private void saveMetrics(int shard, int threads) throws IOException {
    try (BufferedWriter writer =
        Files.newBufferedWriter(
            metricsPath,
            StandardOpenOption.WRITE,
            StandardOpenOption.CREATE,
            StandardOpenOption.TRUNCATE_EXISTING)) {
      writer.write("{\n");
      writer.write(String.format("  \"shard_number\": \"%s\"\n", shard));
      writer.write(String.format("  \"shards_total\": \"%s\"\n", totalShards));
      writer.write(String.format("  \"vcf_count\": \"%s\"\n", vcfPaths.size()));
      writer.write(String.format("  \"threads\": \"%s\"\n", threads));
      writer.write(
          String.format("  \"begin_cut\": \"%s:%s\"\n", beginCut.contig(), beginCut.pos()));
      saveOffsetMetric(writer, "begin_offset", beginOffsets);
      if (endCut == null) {
        writer.write("  \"end_cut\": \"null\"\n");
      } else {
        writer.write(String.format("  \"end_cut\": \"%s:%s\"\n", endCut.contig(), endCut.pos()));
      }
      saveOffsetMetric(writer, "end_offset", endOffsets);
      // "init" includes the time to open the file and seek to the index-specified offset.
      writer.write(String.format("  \"init_s\": \"%s\"\n", initTime.elapsed(TimeUnit.SECONDS)));
      // "write" is the time to copy the shard to disk, including possibly splitting the
      // first or last line.
      writer.write(String.format("  \"write_s\": \"%s\"\n", writeTime.elapsed(TimeUnit.SECONDS)));
      // "total" measures the time spend in the "shard" method.
      writer.write(String.format("  \"total_s\": \"%s\"\n", totalTime.elapsed(TimeUnit.SECONDS)));
      // seek-safe-cut time not shown, but can be computed from the rest.

      if (!skipWriting) {
        ImmutableList.Builder<Long> sizesBuilder = ImmutableList.builder();
        for (Path p : outPaths) {
          sizesBuilder.add(Files.size(p));
        }
        saveOffsetMetric(writer, "shard_size", sizesBuilder.build());
        writer.write(String.format("  \"ref_queried\": \"%d\"\n", reference.getQueryCount()));
      } else {
        writer.write("  \"write_skipped\": \"true\"\n");
      }
      writer.write("}\n");
    }
  }

  private void saveOffsetMetric(BufferedWriter writer, String name, ImmutableList<Long> offsets)
      throws IOException {
    double total = 0;
    long min = offsets.get(0);
    long max = offsets.get(0);
    for (long o : offsets) {
      total += o;
      if (o < min) {
        min = o;
      } else if (o > max) {
        max = o;
      }
    }
    writer.write(String.format("  \"%s_min\": \"%s\"\n", name, min));
    writer.write(
        String.format(
            "  \"%s_avg\": \"%s\"\n", name, (long) Math.round(total / (double) offsets.size())));
    writer.write(String.format("  \"%s_max\": \"%s\"\n", name, max));
    // one QC check we want to do is make sure that for each input file, one shard ends exactly
    // where the next one begins. That's why we print one of the values, so we can check (we can't
    // print the value for every input, there are too many of them).
    writer.write(String.format("  \"%s_first\": \"%s\"\n", name, offsets.get(0)));
  }

  /** Copies a subset of the inputs. */
  private class SubsetSharder {
    private final int begin;
    private final int endExcl;

    /** Sets which inputs to shard: from "begin" (included) to "endExcl" (excluded). */
    SubsetSharder(int begin, int endExcl) {
      checkArgument(begin >= 0 && begin <= endExcl, "begin must be 0<=begin<=endExcl");
      this.begin = begin;
      this.endExcl = Math.min(endExcl, vcfPaths.size());
    }

    /** Copies the shard from the specified inputs to disk. */
    boolean copy() throws IOException, ParseException {
      for (int i = begin; i < endExcl; i++) {
        try (SeekableByteChannel outChan =
            Files.newByteChannel(
                outPaths.get(i),
                StandardOpenOption.WRITE,
                StandardOpenOption.CREATE,
                StandardOpenOption.TRUNCATE_EXISTING)) {
          VcfReader vcfReader = new VcfReader(vcfPaths.get(i), contigs, reference);
          vcfReader.copy(beginOffsets.get(i), beginCut, endOffsets.get(i), endCut, outChan);
        }
      }
      return true;
    }
  }

}
