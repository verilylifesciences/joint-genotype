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

import static com.google.common.truth.Truth.assertWithMessage;
import static java.nio.charset.StandardCharsets.UTF_8;

import com.google.common.collect.ImmutableList;
import com.google.common.jimfs.Configuration;
import com.google.common.jimfs.Jimfs;
import java.io.BufferedReader;
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
public class SharderTest {

  private static final String TEST_DATA =
      "src/test/resources/";

  private static final Path SHARDS = Paths.get(TEST_DATA, "shards.tsv");
  private static final Path SHARDS_10 = Paths.get(TEST_DATA, "shards-10.tsv");

  private static final ImmutableList<Path> vcfs =
      ImmutableList.of(Paths.get(TEST_DATA, "test1.vcf"), Paths.get(TEST_DATA, "test2.vcf"));
  private static final ImmutableList<Path> short_vcfs =
      ImmutableList.of(Paths.get(TEST_DATA, "test1.vcf"), Paths.get(TEST_DATA, "test2b.vcf"));

  @Test
  public void testFirstShardBeginsAtOne() throws IOException, ParseException, InterruptedException {
    try (FileSystem fs = Jimfs.newFileSystem(Configuration.unix())) {
      // Set up.
      ImmutableList<Path> mindexPaths = makeEmptyMindexes(fs, 5);

      for (int i = 0; i < 5; i++) {
        Path metricsPath = fs.getPath("metrics.json");
        ImmutableList<Path> outPaths =
            ImmutableList.of(
                fs.getPath("gold1_" + i + "-of-5.vcf"), fs.getPath("gold2_" + i + "-of-5.vcf"));
        Path[] goldPaths =
            new Path[] {
              Paths.get(TEST_DATA, "gold1_" + i + "-of-5.vcf"),
              Paths.get(TEST_DATA, "gold2_" + i + "-of-5.vcf"),
            };
        // Shard.
        new Sharder(SHARDS, mindexPaths, vcfs, 5, outPaths, metricsPath, new TestReference())
            .shard(i, 4);

        // Compare vs known-good output.
        for (int j = 0; j < 2; j++) {
          String got = new String(Files.readAllBytes(outPaths.get(j)), UTF_8);
          String expected = new String(Files.readAllBytes(goldPaths[j]), UTF_8);
          assertWithMessage("mismatch in gold" + (j + 1) + "_" + i + "-of-5")
              .that(got).isEqualTo(expected);
        }
        // Compare metrics vs known-good
        try (BufferedReader metrics = Files.newBufferedReader(metricsPath);
            BufferedReader gold = Files.newBufferedReader(
                Paths.get(TEST_DATA, "gold_metrics-" + i + ".json"))) {
          while (true) {
            String expected = gold.readLine();
            if (expected == null) {
              break;
            }
            String got = metrics.readLine();
            // do not compare the timings
            if (expected.contains("_s\"")) {
              continue;
            }
            assertWithMessage("mismatch in gold_metrics-" + i + ".json")
                .that(got).isEqualTo(expected);
          }
        }
      }
    }
  }

  /**
   * One of the inputs doesn't include any data for the shard, so we hit end of file.
   * Check we don't crash (and return the correct output).
   */
  @Test
  public void testLastShardWhenPastEOF() throws IOException, ParseException, InterruptedException {
    try (FileSystem fs = Jimfs.newFileSystem(Configuration.unix())) {
      // Set up.
      ImmutableList<Path> mindexPaths = makeShortMindexes(fs, 5);
      Path metricsPath = Paths.get("");
      ImmutableList<Path> outPaths =
          ImmutableList.of(fs.getPath("test1_4-of-5.vcf"), fs.getPath("test2b_4-of-5.vcf"));
      ImmutableList<Path> goldPaths =
          ImmutableList.of(
              Paths.get(TEST_DATA, "gold1_4-of-5.vcf"), Paths.get(TEST_DATA, "gold2b_4-of-5.vcf"));
      // Shard.
      new Sharder(SHARDS, mindexPaths, short_vcfs, 5, outPaths, metricsPath, new TestReference())
          .shard(4, 4);
      // Compare vs known-good output.
      for (int j = 0; j < 2; j++) {
        String got = new String(Files.readAllBytes(outPaths.get(j)), UTF_8);
        String expected = new String(Files.readAllBytes(goldPaths.get(j)), UTF_8);
        assertWithMessage("mismatch in file " + j).that(got).isEqualTo(expected);
      }
    }
  }

  /**
   * In this test we have twice as many shards (half the size), but we ask the program to cut shards
   * two at a time. The resulting cut points should be the same as the initial 5-shard case.
   */
  @Test
  public void testMultipleShardsAtATime() throws IOException, ParseException, InterruptedException {
    try (FileSystem fs = Jimfs.newFileSystem(Configuration.unix())) {
      // Set up.
      ImmutableList<Path> mindexPaths = makeEmptyMindexes(fs, 10);

      for (int i = 0; i < 5; i++) {
        Path metricsPath = fs.getPath("metrics.json");
        ImmutableList<Path> outPaths =
            ImmutableList.of(
                fs.getPath("gold1_" + i + "-of-5.vcf"), fs.getPath("gold2_" + i + "-of-5.vcf"));
        // Shard.
        new Sharder(SHARDS_10, mindexPaths, vcfs, 5, outPaths, metricsPath, new TestReference())
            .shard(i, 4);

        // Compare metrics vs known-good (are the cut points correct?)
        try (BufferedReader metrics = Files.newBufferedReader(metricsPath);
            BufferedReader gold =
                Files.newBufferedReader(Paths.get(TEST_DATA, "gold_metrics-" + i + ".json"))) {
          while (true) {
            String expected = gold.readLine();
            String actual = metrics.readLine();
            // do not compare the timings
            if (expected != null && expected.contains("_s\"")) {
              continue;
            }
            assertWithMessage("mismatch in gold_metrics-" + i + ".json")
                .that(actual)
                .isEqualTo(expected);
            if (expected == null) {
              // we check the last line, too.
              break;
            }
          }
        }
      }
    }
  }

  private static ImmutableList<Path> makeEmptyMindexes(FileSystem fs, int nRecords)
      throws IOException {
    Path mindex = makeEmptyMindex(fs, nRecords);
    return ImmutableList.of(mindex, mindex);
  }

  private static Path makeEmptyMindex(FileSystem fs, int nRecords) throws IOException {
    Path mindex = fs.getPath("mindex");
    ByteBuffer fourzeroes = ByteBuffer.allocate(nRecords * 8);
    LongBuffer longBuffer = fourzeroes.asLongBuffer();
    for (int i = 0; i < nRecords; i++) {
      longBuffer.put(0);
    }
    Files.write(mindex, fourzeroes.array(), StandardOpenOption.WRITE, StandardOpenOption.CREATE);
    return mindex;
  }

  private static ImmutableList<Path> makeShortMindexes(FileSystem fs, int nRecords)
      throws IOException {
    Path mindex = fs.getPath("mindex");
    ByteBuffer mindexContents = ByteBuffer.allocate(nRecords * 8);
    LongBuffer mindexContentsAsLongs = mindexContents.asLongBuffer();
    for (int i = 0; i < nRecords - 1; i++) {
      mindexContentsAsLongs.put(0);
    }
    mindexContentsAsLongs.put(99999999);
    Files.write(
        mindex, mindexContents.array(), StandardOpenOption.WRITE, StandardOpenOption.CREATE);
    return ImmutableList.of(makeEmptyMindex(fs, nRecords), mindex);
  }
}
