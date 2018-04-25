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
package com.verily.lifescience.genomics.wgs.mindex;

import static com.google.common.truth.Truth.assertThat;
import static org.junit.jupiter.api.Assertions.assertThrows;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import java.io.IOException;
import java.nio.file.Paths;
import java.text.ParseException;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

@RunWith(JUnit4.class)
public class PositionTest {

  static final String TEST_DATA = "src/test/resources/";

  @Test
  public void loadContigs() throws Exception {
    String fname = TEST_DATA + "positions.tsv";
    ImmutableMap<String, Integer> contigs = Position.contigsFromFile(Paths.get(fname));
    ImmutableMap<String, Integer> expected = ImmutableMap.of("chr1", 0, "chr2", 1, "chrM", 2);
    assertThat(contigs).containsExactlyEntriesIn(expected);
  }

  @Test
  public void loadFromFile() throws IOException, ParseException {
    String fname = TEST_DATA + "positions.tsv";
    ImmutableMap<String, Integer> contigs = Position.contigsFromFile(Paths.get(fname));
    Position[] expected =
        new Position[] {
          Position.of("chr1", 10, contigs),
          Position.of("chr2", 500, contigs),
          Position.of("chrM", 1, contigs)
        };
    ImmutableList<Position> pl = Position.fromFile(Paths.get(fname), contigs);
    assertThat(pl.size()).isEqualTo(expected.length);
    for (int i = 0; i < pl.size(); i++) {
      assertThat(pl.get(i)).isEqualTo(expected[i]);
    }
  }

  @Test
  public void noNulls() throws IOException, ParseException {
    ImmutableMap<String, Integer> contigs = ImmutableMap.of("chr1", 0, "chr10", 1, "chr11", 2);
    assertThrows(NullPointerException.class, () -> Position.of(null, 0, contigs));
    assertThrows(NullPointerException.class, () -> Position.of("chr1", 0, null));
  }

  @Test
  public void loadFromBadFile() throws IOException, ParseException {
    String fname = TEST_DATA + "bad_positions.tsv";
    assertThrows(
        ParseException.class,
        () -> Position.fromFile(Paths.get(fname), Position.contigsFromFile(Paths.get(fname))));
  }

  @Test
  public void ordersCorrectly() throws IOException, ParseException {
    ImmutableMap<String, Integer> contigs = ImmutableMap.of("chr1", 0, "chr10", 1, "chr11", 2);
    assertThat(Position.of("chr10", 10, contigs)).isLessThan(Position.of("chr10", 20, contigs));
    assertThat(Position.of("chr1", 20, contigs)).isLessThan(Position.of("chr10", 10, contigs));
    assertThat(Position.of("chr10", 10, contigs)).isLessThan(Position.of("chr11", 5, contigs));
  }

  @Test
  public void checkEqual() {
    ImmutableMap<String, Integer> contigs = ImmutableMap.of("chr1", 0, "chr10", 1, "chr11", 2);
    ImmutableMap<String, Integer> contigs2 = ImmutableMap.of("chr0", 0, "chr1", 1, "chr4", 2);
    Position apos = Position.of("chr1", 10, contigs);
    Position[][] equalPairs =
        new Position[][] {
          // pointer equal
          new Position[] {apos, apos},
          // all values equal
          new Position[] {Position.of("chr1", 10, contigs), Position.of("chr1", 10, contigs)},
          // contigs different, but rest equals
          new Position[] {Position.of("chr1", 10, contigs), Position.of("chr1", 10, contigs2)},
        };
    for (Position[] pair : equalPairs) {
      assertThat(pair[0].equals(pair[1])).isTrue();
      assertThat(pair[0].compareTo(pair[1])).isEqualTo(0);
      assertThat(pair[0].hashCode()).isEqualTo(pair[1].hashCode());
    }
  }

  @Test
  public void notComparableIfContigOrdersNotPointerEqual() {
    ImmutableMap<String, Integer> contigs1 = ImmutableMap.of("chr1", 0, "chr10", 1, "chr11", 2);
    ImmutableMap<String, Integer> contigs2 = ImmutableMap.of("chr1", 0, "chr10", 1, "chr11", 2);
    Position pos1 = Position.of("chrX", 10, contigs1);
    Position pos2 = Position.of("chrM", 10, contigs2);
    assertThrows(
        IllegalArgumentException.class,
        () -> pos1.compareTo(pos2));
  }
}
