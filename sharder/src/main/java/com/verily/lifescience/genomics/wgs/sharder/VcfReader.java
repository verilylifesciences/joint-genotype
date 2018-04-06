package com.verily.lifescience.genomics.wgs.sharder;

import static java.nio.charset.StandardCharsets.UTF_8;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ImmutableMap;
import com.verily.lifescience.genomics.wgs.mindex.Position;
import java.io.BufferedReader;
import java.io.EOFException;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.Channels;
import java.nio.channels.SeekableByteChannel;
import java.nio.channels.WritableByteChannel;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.ParseException;
import java.util.Optional;
import javax.annotation.Nullable;

/**
 * Lets you seek through a VCF file and read records.
 *
 * <p>Conceptually, it keeps a cursor on one of the lines, the "current record." You can ask for
 * that record, or the genomic position it represents. You can also get a file offset either just at
 * the beginning of this record, or just after. The former is called "currentOffset", the latter
 * "nextOffset". We also keep the previous record and its offset.
 *
 * <p>Note that this class is not so much about parsing VCF, it's really about having the ability to
 * both read lines and seek to arbitrary byte offsets, as well as containing the logic of seeking to
 * a genomic position and copying shards of VCF.
 */
public class VcfReader implements AutoCloseable {

  /** Returned when the offset is not known. */
  public static final long UNKNOWN_OFFSET = -1;

  private final SeekableByteChannel channel;
  private final ImmutableMap<String, Integer> contigs;
  @Nullable private BufferedReader input = null;
  private final Reference reference;
  /** Offset where we'll next read from. */
  private long offset = 0;
  /**
   * Offset where we read the previous record (or UNKNOWN_OFFSET if we've only read one record since
   * the last seek).
   */
  private long prevOffset = UNKNOWN_OFFSET;
  /** Latest line we read (if any). (null at EOF) */
  @Nullable private String latestLine = null;
  @Nullable private String prevLine = null;
  private boolean primed;

  public VcfReader(Path path, ImmutableMap<String, Integer> contigs, Reference reference)
      throws IOException {
    this.channel = Files.newByteChannel(path);
    this.contigs = contigs;
    this.reference = reference;
  }

  public void semiClose() {
    input = null;
  }


  /**
   * Seek into the file.
   *
   * @param offset where to go (in bytes).
   * @throws IOException
   */
  public void offset(long offset) throws IOException {
    if (null != input) {
      // Closing would close the channel too.
      input = null;
    }
    channel.position(offset);
    this.prevOffset = UNKNOWN_OFFSET;
    this.offset = offset;
    this.prevLine = null;
    this.latestLine = null;
    this.primed = false;
    prime();
  }

  /**
   * Returns the previous record (if any).
   *
   * <p>Does not advance the cursor.
   */
  public Optional<String> getPreviousRecord() throws IOException {
    prime();
    if (prevLine == null) {
      return Optional.empty();
    }
    return Optional.of(prevLine);
  }

  /**
   * Returns the current record.
   *
   * <p>Does not advance the cursor.
   *
   * @throws EOFException at the end of the file.
   */
  public String getRecord() throws IOException {
    prime();
    if (null == latestLine) {
      throw new EOFException();
    }
    return latestLine;
  }

  /** Returns true if we are at the end of the file. */
  public boolean isEof() throws IOException {
    prime();
    return (latestLine == null);
  }

  /**
   * Returns the current genomic position.
   *
   * <p>Does not advance the cursor.
   *
   * @throws EOFException at the end of the file.
   */
  public Position getPosition() throws IOException {
    return parsePosition(getRecord());
  }

  /**
   * Returns the genomic position of the previous record (if any).
   *
   * <p>Does not advance the cursor.
   */
  public Optional<Position> getPreviousPosition() throws IOException {
    Optional<String> previousRecord = getPreviousRecord();
    if (!previousRecord.isPresent()) {
      return Optional.empty();
    }
    return Optional.of(parsePosition(previousRecord.get()));
  }

  /** Parses the record, returns its position. */
  public Position parsePosition(String line) {
    String[] elements = line.split("\t", 3);
    if (elements.length<2) {
      // You may get this error if the index or mindex was wrong, so you jump partway into a record.
      // This would happen if e.g. you change the header after the fact.
      throw new IllegalArgumentException("Expected at least 2 tab-separated fields, got " + elements.length + " in '" + line + "'.");
    }
    return Position.of(elements[0], Integer.valueOf(elements[1]), contigs);
  }

  /**
   * Returns the end genomic position for the given current record, if any.
   *
   * @return empty if record doesn't have an "END" tag.
   */
  public Optional<Position> parseEndPosition(String line) {
    String[] elements = line.split("\t", 9);
    if (elements.length < 8) {
      // no "END" tag because that column isn't there at all.
      return Optional.empty();
    }
    String end = elements[7];
    if (!end.startsWith("END=")) {
      // no "END" tag (probably just a normal record)
      return Optional.empty();
    }
    return Optional.of(
        Position.of(elements[0], Integer.valueOf(elements[7].substring("END=".length())), contigs));
  }

  /**
   * Advances to the next record (skipping comment lines).
   *
   * <p>Note that the reader starts at the first record, so calling next() on a fresh reader will
   * move you to the second record.
   */
  public void next() throws IOException {
    prime();
    inner_next();
  }

  private void inner_next() throws IOException {
    readline();
    while (this.latestLine != null && this.latestLine.startsWith("#")) {
      readline();
    }
  }

  /** Returns the offset for the record before the current record. */
  public long getPreviousOffset() {
    return this.prevOffset;
  }

  /** Returns the offset for the current record. */
  public long getCurrentOffset() {
    if (null == latestLine) {
      return this.offset;
    }
    return this.offset - latestLine.length() - 1;
  }

  /** Returns the offset for the record after the current record. */
  public long getNextOffset() {
    return this.offset;
  }

  /**
   * Moves forward in the file until this position or later.
   *
   * @throws IllegalArgumentException if throwIfPast and the current position is already beyond the
   *     target.
   */
  public void advanceTo(Position target, boolean throwIfPast) throws IOException {
    if (isEof()) {
      return;
    }
    Position current = getPosition();
    if (throwIfPast && target.before(current)) {
      throw new IllegalArgumentException("Current position already beyond the target.");
    }
    while (current.before(target)) {
      next();
      if (isEof()) {
        break;
      }
      current = getPosition();
    }
  }

  /**
   * Moves forward in the file until this position or later.
   *
   * @throws IllegalArgumentException if the current position is already beyond the target.
   */
  public void advanceTo(Position pos) throws IOException {
    advanceTo(pos, true);
  }

  /**
   * Moves forward in the file until this position or later. Does not report an error if we start
   * already after the position.
   */
  public void advanceToAtLeast(Position pos) throws IOException {
    if (isEof()) {
      return;
    }
    advanceTo(pos, false);
  }

  /** "Primes the pump": reads the current record if we haven't yet. */
  private void prime() throws IOException {
    if (!primed) {
      // fetch the first record, skipping over comments.
      inner_next();
      primed = true;
    }
  }

  /** Reads one line from the file, updates our state. If at end of file, then return null. */
  @Nullable
  private String readline() throws IOException {
    BufferedReader myInput;
    synchronized (this) {
      if (null == input) {
        this.input = new BufferedReader(Channels.newReader(channel, "UTF-8"));
      }
      myInput = this.input;
    }
    prevOffset = getCurrentOffset();
    @Nullable String line = myInput.readLine();
    prevLine = latestLine;
    latestLine = line;
    if (null == line) {
      // EOF
      return null;
    }
    // We assume UNIX-style lines, split by a single newline character (LF).
    // We also assume that we only have single-byte characters. Otherwise, it all goes south.
    offset += line.length() + 1;
    return latestLine;
  }

  @Override
  public void close() throws IOException {
    if (null != input) {
      input.close();
    }
    channel.close();
  }

  /**
   * Writes the first record, splitting an existing one if necessary. Returns the offset after the
   * record.
   */
  @VisibleForTesting
  long saveFirstRecord(Position start, WritableByteChannel dest)
      throws IOException, ParseException {
    ByteBuffer newline = ByteBuffer.wrap("\n".getBytes(UTF_8));
    advanceTo(start);
    if (isEof()) {
      return channel.size();
    }
    Position afterFirstCut = getPosition();
    if (!afterFirstCut.equals(start)) {
      // the previous record starts before the cut. We may have to split it.
      Optional<String> lineBeforeMaybe = getPreviousRecord();
      if (!lineBeforeMaybe.isPresent()) {
        throw new IllegalArgumentException("copy() given offset too close to safe cut.");
      }
      String lineBefore = lineBeforeMaybe.get();

      // We wanted "start". The previous record must begin before "start" (can't start there
      // exactly,
      // or we should have been there).
      Position beforePos = parsePosition(lineBefore);
      if (!beforePos.before(start)) {
        throw new RuntimeException(
            "Unexpected (bug?): line before the cut should have been included:\n"
                + lineBefore
                + "\nIts position is "
                + beforePos
                + " and the cut starts at "
                + start);
      }

      // afterFirstCut is where we ended up. There are two possible scenarios where
      // afterFirstCut is after start:
      // the previous record goes beyond start, or there is a gap.
      // The ASCII art below illustrates these scenarios.
      //
      // -[-----------]  +                   +------------
      //  ^beforePos     ^start (in the gap) ^afterFirstCut
      //
      // -[----------+------]+-------
      //  ^beforePos ^start  ^afterFirstCut
      //
      // In this second scenario, we have to split the previous record in two.
      // We include only the second half in our cut.
      //
      // result: (double line means included in the cut).
      // -[---------][======]+=====
      //  ^beforePos ^start  ^afterFirstCut

      Optional<Position> endOfBefore = parseEndPosition(lineBefore);
      if (endOfBefore.isPresent() && endOfBefore.get().compareTo(start) >= 0) {
        String[] parts = lineBefore.split("\t");
        parts[1] = "" + start.pos();
        String refBase = reference.getBase(start.contig(), start.pos());
        if (refBase.length() != 1) {
          throw new RuntimeException(
              String.format(
                  "Bug: reference base at %s should be 1 letter, but is: '%s'", start.toString(),
                  refBase));
        }
        parts[3] = refBase;
        String secondHalf = String.join("\t", parts);
        dest.write(ByteBuffer.wrap(secondHalf.getBytes(UTF_8)));
        newline.position(0);
        newline.mark();
        newline.limit(1);
        dest.write(newline);
      }
    }
    if (!isEof()) {
      String firstIncludedLine = getRecord();
      dest.write(ByteBuffer.wrap(firstIncludedLine.getBytes(UTF_8)));
      newline.position(0);
      newline.mark();
      newline.limit(1);
      dest.write(newline);
    }
    return getNextOffset();
  }

  /**
   * Using 1-based, include-include position because that's what everyone else does. File offsets
   * are 0-based as normal.
   *
   * @param startOffset - byte offset into the file at the beginning of a line, and before start.
   * @param startPosition - the first genomic position we want included.
   * @param nextShardPosition - the first genomic position we do not want included (null if want
   *     until the end).
   * @param dest - file to write the output to.
   * @return total bytes read
   * @throws IOException
   */
  public long copy(
      long startOffset,
      Position startPosition,
      long endOffset,
      @Nullable Position nextShardPosition,
      WritableByteChannel dest)
      throws IOException, ParseException {
    if (startOffset > channel.size()) {
      // No data to copy from this file.
      return 0;
    }
    offset(startOffset);
    // Find and write the first record (with possible split)
    long offset = saveFirstRecord(startPosition, dest);

    // Bulk copy in the middle
    long total = offset - startOffset;
    final int bufSize = 1024 * 1024;
    ByteBuffer buf = ByteBuffer.allocate(bufSize);
    channel.position(offset);
    while (offset < endOffset) {
      if (offset + bufSize > endOffset) {
        // Read only as much as we need to end up at the end.
        int newLimit = (int) (endOffset - offset);
        buf.limit(newLimit);
      }
      int count = channel.read(buf);
      if (count < 0) {
        // EOF
        break;
      }
      offset += count;
      total += count;
      buf.flip();
      dest.write(buf);
      buf.clear();
    }

    if (null != nextShardPosition) {
      // At this point we have read until just before endOffset.
      // One of the upcoming lines will be the boundary.
      offset(offset);
      total += VcfReaderHelper.saveLastRecord(this, nextShardPosition, dest);
    }
    return total;
  }
}
