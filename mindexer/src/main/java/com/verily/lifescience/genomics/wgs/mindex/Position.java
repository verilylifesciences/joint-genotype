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

import com.google.auto.value.AutoValue;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.ParseException;
import java.util.Objects;
import java.util.stream.Stream;
import javax.annotation.Nullable;

/**
 * Genomic position, 1-based.
 *
 * <p>Each position is described by a contig name, and position inside that contig.
 */
@AutoValue
public abstract class Position implements Comparable<Position> {

  /**
   * Create a new Position with the given values.
   *
   * <p>The contigs map is used for ordering. Pass what you got from {@link #contigsFromFile}.
   */
  public static Position of(String contig, int pos, ImmutableMap<String, Integer> contigs) {
    return new AutoValue_Position(contig, pos, contigs);
  }

  /**
   * Returns a list of start positions from a file.
   *
   * <p>The file must have a tab-separated list of genomic intervals. Each line in the file is
   * either a comment (starting with "#"), or a tab-separated sequence of interval triplets. Each
   * interval triplet is CONTIG\tSTART_POSITION\tEND_POSITION.
   *
   * <p>The returned list has one entry per non-comment line in the input file: the start position
   * of the first interval in that line. All the entries share the same contigs map, naturally.
   */
  public static ImmutableList<Position> fromFile(Path file) throws IOException, ParseException {
    return fromFile(file, contigsFromFile(file));
  }

  /**
   * Returns a list of start positions from a file.
   *
   * <p>The file must have a tab-separated list of genomic intervals. Each line in the file is
   * either a comment (starting with "#"), or a tab-separated sequence of interval triplets. Each
   * interval triplet is CONTIG\tSTART_POSITION\tEND_POSITION.
   *
   * <p>"contigs" is the value returned by {@link #contigsFromFile}.
   *
   * <p>The returned list has one entry per non-comment line in the input file: the start position
   * of the first interval in that line. All the entries share the same contigs map, naturally.
   */
  public static ImmutableList<Position> fromFile(Path file, ImmutableMap<String, Integer> contigs)
      throws IOException, ParseException {
    try (final Stream<String> lines = Files.lines(file, java.nio.charset.StandardCharsets.UTF_8)) {
      ImmutableList.Builder<Position> positions = ImmutableList.builder();
      for (String s : (Iterable<String>) lines::iterator) {
        Position p = parseLine(s, contigs);
        if (null != p) {
          positions.add(p);
        }
      }
      return positions.build();
    }
  }

  /**
   * Returns a list of contigs mentioned in a file.
   *
   * <p>The file must have a tab-separated list of genomic intervals. Each line in the file is
   * either a comment (starting with "#"), or a tab-separated sequence of interval triplets. Each
   * interval triplet is CONTIG\tSTART_POSITION\tEND_POSITION.
   *
   * <p>The returned map assigns an index number to each contig mentioned in the file. They are
   * assigned sequentially, starting at zero.
   */
  public static ImmutableMap<String, Integer> contigsFromFile(Path file)
      throws IOException, ParseException {
    try (final Stream<String> lines = Files.lines(file, java.nio.charset.StandardCharsets.UTF_8)) {
      ImmutableMap.Builder<String, Integer> contigs = ImmutableMap.builder();
      String lastContig = "";
      int count = 0;
      for (String s : (Iterable<String>) lines::iterator) {
        for (String contig : parseContigs(s)) {
          if (!contig.equals(lastContig)) {
            contigs.put(contig, count);
            count++;
            lastContig = contig;
          }
        }
      }
      return contigs.build();
    }
  }

  private static @Nullable Position parseLine(String line, ImmutableMap<String, Integer> contigs)
      throws ParseException {
    // expected: CONTIG\tSTART_POSITION\tEND_POSITION\tCONTIG2 ...
    // or # comment
    // (or empty line at the end)
    if (line.length()==0 || line.startsWith("#")) {
      return null;
    }
    String[] elems = checkAndSplit(line);
    return Position.of(elems[0], Integer.valueOf(elems[1]), contigs);
  }

  private static ImmutableList<String> parseContigs(String line) throws ParseException {
    // expected: CONTIG\tSTART_POSITION\tEND_POSITION\tCONTIG2 ...
    // or # comment
    if (line.length()==0 || line.startsWith("#")) {
      return ImmutableList.of();
    }
    String[] elems = checkAndSplit(line);
    ImmutableList.Builder<String> builder = ImmutableList.builder();
    for (int i = 0; i < elems.length; i += 3) {
      builder.add(elems[i]);
    }
    return builder.build();
  }

  private static String[] checkAndSplit(String line) throws ParseException {
    String[] elems = line.split("\t");
    if (elems.length % 3 != 0) {
      throw new ParseException("Expected records in 3s, got " + elems.length + ":\n" + line, 0);
    }
    if (elems.length < 3) {
      throw new ParseException(
          "Expected at least 3 records, got " + elems.length + ":\n" + line, line.length() - 1);
    }
    return elems;
  }

  /** Name of the contig. */
  public abstract String contig();

  /** Genomic position within the contig, 1-based. */
  public abstract int pos();

  /**
   * The list of contigs, in order. Encoded as a mapping from contig name to its position in the
   * list (so 0 for the earliest contig).
   */
  public abstract ImmutableMap<String, Integer> contigs();

  /**
   * Two positions are equal if they're the same (contig name, pos in contig), even if not for the
   * same file. That is, we don't compare the contigs pointers.
   */
  @Override
  public boolean equals(Object obj) {
    if (!(obj instanceof Position)) {
      return false;
    }
    Position rhs = (Position) obj;
    return (pos() == rhs.pos() && contig().equals(rhs.contig()));
  }

  /**
   * Two positions that are equal always have the same hash code.
   */
  @Override
  public int hashCode() {
    return Objects.hashCode(contig()) + pos();
  }

  /**
   * Equal positions always return 0.
   *
   * Non-equal positions must have the same contigs map (pointer-equal) in order to be comparable
   * (this is the case if they were loaded from the same file).
   */
  @SuppressWarnings("ReferenceEquality")  // explained below
  @Override
  public int compareTo(Position rhs) {
    if (null == rhs) {
      return 1;
    }
    if (contig().equals(rhs.contig())) {
      return Integer.compare(pos(), rhs.pos());
    }
    // We want pointer equality here.
    // Value equality would be way too slow for something we do this frequently.
    // The code explicitly assigns the same (pointer-wise) contigs value to all positions within
    // the same file. Note that the test code does the same. This is the "no performance surprises"
    // principle: if the client code is correct then it's fast. If the client code is wrong then
    // there's no slow fallback: it'll fail with an exception and you know you did wrong.
    if (contigs() != rhs.contigs()) {
      throw new IllegalArgumentException(
          "Can only compare positions if they use the same contig ordering.");
    }
    return Integer.compare(contigs().get(contig()), contigs().get(rhs.contig()));
  }

  /** Returns true iff this is before the argument (earlier in the genome). */
  public boolean before(Position rhs) {
    return this.compareTo(rhs) < 0;
  }

  @Override
  public String toString() {
    return contig() + ":" + pos();
  }
}
