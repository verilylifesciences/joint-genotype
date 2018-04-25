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

import static com.google.common.truth.Truth.assertThat;
import static com.google.common.truth.Truth.assertWithMessage;
import static java.nio.charset.StandardCharsets.UTF_8;

import com.google.common.collect.ImmutableList;
import com.google.common.jimfs.Configuration;
import com.google.common.jimfs.Jimfs;
import com.verily.lifescience.genomics.wgs.mindex.Position;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.LongBuffer;
import java.nio.file.FileSystem;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.text.ParseException;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

@RunWith(JUnit4.class)
public class SafeCutTest {

  private static final String TEST_DATA =
      "src/test/resources/";

  private static final Path SHARDS = Paths.get(TEST_DATA, "shards.tsv");

  private static final ImmutableList<Path> VCFS =
      ImmutableList.of(Paths.get(TEST_DATA, "test1.vcf"), Paths.get(TEST_DATA, "test2.vcf"));

  private static final int THREADS = 4;

  @Test
  public void testIsDeletion() {
    String input =
        "chr1\t465782\t.\tATAT\tA,<NON_REF>\t0\t.\tDP=33;ExcessHet=3.0103;"
            + "MLEAC=0,0;MLEAF=0.00,0.00;RAW_MQ=121400.00\tGT:AD:DP:GQ:PL:SB"
            + "\t0/0:27,0,0:27:84:0,0,0,0,0,0:11,16,0,0";
    boolean got = SafeCut.isDeletion(input);
    assertThat(got).isTrue();
  }

  @Test
  public void testFirstShardBeginsAtOne() throws IOException, ParseException, InterruptedException {
    try (FileSystem fs = Jimfs.newFileSystem(Configuration.unix())) {
      // Set up.
      SafeCut cutter = safeCutFromFiles(fs);
      // Compute very first safe cut.
      cutter.init(0);
      Position safeCut1 = cutter.findSafeCut();
      // There's a deletion from 0 onwards, but that doesn't mean we can't cut BEFORE 0.
      // The first safe cut is at chr1:0 (doing anything else would mean this record
      // would never be included in the output).
      assertThat(safeCut1.contig()).isEqualTo("chr1");
      assertThat(safeCut1.pos()).isEqualTo(1);
    }
  }

  @Test
  public void testAdvancePastDeletion() throws IOException, ParseException, InterruptedException {
    try (FileSystem fs = Jimfs.newFileSystem(Configuration.unix())) {
      // Set up.
      SafeCut cutter = safeCutFromFiles(fs);
      cutter.init(1);
      Position safeCut1 = cutter.findSafeCut();
      // The shard position is inside a deletion in test1.vcf, so we move forward until it's over.
      assertThat(safeCut1.contig()).isEqualTo("chr1");
      assertThat(safeCut1.pos()).isEqualTo(379);
    }
  }

  @Test
  public void testAdvancePastDeletionInSecondFile()
      throws IOException, ParseException, InterruptedException {
    try (FileSystem fs = Jimfs.newFileSystem(Configuration.unix())) {
      // Set up.
      SafeCut cutter = safeCutFromReverseFiles(fs);
      cutter.init(1);
      Position safeCut1 = cutter.findSafeCut();
      // The shard position is inside a deletion in test1.vcf, so we move forward until it's over.
      assertThat(safeCut1.contig()).isEqualTo("chr1");
      assertThat(safeCut1.pos()).isEqualTo(379);
    }
  }

  @Test
  public void testAdvancePastMultipleDeletions()
      throws IOException, ParseException, InterruptedException {
    try (FileSystem fs = Jimfs.newFileSystem(Configuration.unix())) {
      // Set up.
      SafeCut cutter = safeCutFromFiles(fs);
      cutter.init(2);
      Position safeCut1 = cutter.findSafeCut();
      // The shard position is inside a deletion in test1.vcf, so we move forward until it's over.
      // Then that puts us in a deletion in test2.vcf. And that puts us in a deletion in test1.
      assertThat(safeCut1.contig()).isEqualTo("chr2");
      assertThat(safeCut1.pos()).isEqualTo(190);
    }
  }

  @Test
  public void testCutAtContigStart() throws IOException, ParseException, InterruptedException {
    try (FileSystem fs = Jimfs.newFileSystem(Configuration.unix())) {
      // Set up.
      SafeCut cutter = safeCutFromFiles(fs);
      cutter.init(3);
      Position safeCut1 = cutter.findSafeCut();
      // The shard position is inside a deletion in test1.vcf, so we move forward until it's over.
      // Then that puts us in a deletion in test2.vcf. And that puts us in a deletion in test1.
      assertThat(safeCut1.contig()).isEqualTo("chr3");
      assertThat(safeCut1.pos()).isEqualTo(1);
    }
  }

  @Test
  public void testManyInputs() throws IOException, ParseException, InterruptedException {
    int inputCount = 63;
    try (FileSystem fs = Jimfs.newFileSystem(Configuration.unix())) {
      // Set up.
      ImmutableList<Path> mindexPaths = createEmptyMindexes(fs, inputCount);
      ImmutableList.Builder<Path> vcfBuilder = ImmutableList.builder();
      // The VCFs have staggered deletions.
      // So if the splitter doesn't look at all the VCFs, it'll end up at the wrong spot.
      for (int i = 0; i < inputCount; i++) {
        Path vcf = fs.getPath("vcf-" + i);
        BufferedWriter writer = Files.newBufferedWriter(vcf, UTF_8);
        int pos = i * 4 + 349;
        // deletion at the very beginning
        writer.write("chr1\t1\t.\tAAAAA\tA,<NON_REF>\t0\t.\t\t\t\n");
        // splittable zone before
        if (pos > 5) {
          writer.write("chr1\t5\t.\tC\t<NON_REF>\t.\t.\tEND=" + (pos - 1) + "\t\t\n");
        }
        // unsplittable deletion
        writer.write("chr1\t" + pos + "\t.\tAAAAA\tA,<NON_REF>\t0\t.\t\t\t\n");
        // splittable zone after
        writer.write("chr1\t" + (pos + 5) + "\t.\tA\t<NON_REF>\t.\t.\tEND=999000\t\t\n");
        writer.close();
        vcfBuilder.add(vcf);
      }
      ImmutableList<Path> vcfPaths = vcfBuilder.build();
      SafeCut cutter1 = new SafeCut(SHARDS, mindexPaths, vcfPaths, THREADS, new TestReference());
      // shard 1: at chr1:350
      cutter1.init(1);
      Position safeCut1 = cutter1.findSafeCut();
      assertWithMessage("Forward direction").that(safeCut1.contig()).isEqualTo("chr1");
      assertWithMessage("Forward direction")
          .that(safeCut1.pos())
          .isEqualTo(4 * (inputCount - 1) + 349 + 5);

      ImmutableList<Path> vcfPathsReverse = vcfPaths.reverse();
      SafeCut cutter2 = new SafeCut(SHARDS, mindexPaths, vcfPathsReverse, THREADS,
          new TestReference());
      // shard 1: at chr1:350
      cutter2.init(1);
      Position safeCut2 = cutter2.findSafeCut();
      assertWithMessage("Reverse direction").that(safeCut2.contig()).isEqualTo("chr1");
      assertWithMessage("Reverse direction")
          .that(safeCut2.pos())
          .isEqualTo(4 * (inputCount - 1) + 349 + 5);
    }
  }

  private static SafeCut safeCutFromFiles(FileSystem fs) throws IOException, ParseException {
    ImmutableList<Path> mindexPaths = createEmptyMindexes(fs, 2);
    return new SafeCut(SHARDS, mindexPaths, VCFS, THREADS, new TestReference());
  }

  private static SafeCut safeCutFromReverseFiles(FileSystem fs) throws IOException, ParseException {
    ImmutableList<Path> mindexPaths = createEmptyMindexes(fs, 2);
    return new SafeCut(SHARDS, mindexPaths, ImmutableList.of(VCFS.get(1), VCFS.get(0)), THREADS,
        new TestReference());
  }

  private static ImmutableList<Path> createEmptyMindexes(FileSystem fs, int howMany)
      throws IOException {
    Path mindex = fs.getPath("mindex");
    ByteBuffer fourzeroes = ByteBuffer.allocate(4 * 8);
    LongBuffer longBuffer = fourzeroes.asLongBuffer();
    for (int i = 0; i < 4; i++) {
      longBuffer.put(0);
    }
    Files.write(mindex, fourzeroes.array(), StandardOpenOption.WRITE, StandardOpenOption.CREATE);
    ImmutableList.Builder<Path> builder = ImmutableList.builder();
    for (int i = 0; i < howMany; i++) {
      builder.add(mindex);
    }
    return builder.build();
  }


}
