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

import static java.nio.charset.StandardCharsets.UTF_8;

import com.verily.lifescience.genomics.wgs.mindex.Position;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.WritableByteChannel;
import java.util.Optional;

/** Contains static helper methods for {@link VcfReader}. */
class VcfReaderHelper {

  // Static methods only, don't instantiate.
  private VcfReaderHelper() {}

  /**
   * Writes every line (starting at the previous record), up to and excluding "excluded". If the
   * last record extends past "excluded", it'll be trimmed.
   */
  static long saveLastRecord(VcfReader reader, Position excluded, WritableByteChannel dest)
      throws IOException {
    long total = 0;
    ByteBuffer newline = ByteBuffer.wrap("\n".getBytes(UTF_8));
    String oldLine = null;
    Position pos;
    do {
      if (null != oldLine) {
        // Write the previous line, it's included.
        total += (oldLine.length() + 1);
        dest.write(ByteBuffer.wrap(oldLine.getBytes(UTF_8)));
        newline.position(0);
        newline.mark();
        newline.limit(1);
        dest.write(newline);
      }
      oldLine = reader.getPreviousRecord().orElse(null);
      pos = reader.getPosition();
      reader.next();
      // Exit when we reach or pass the final position.
    } while (pos.before(excluded));
    // The last line we read was past the end, so we don't write it.
    // The one before that, we include. We may have to trim the end though.
    if (null != oldLine) {
      // Write the previous line, it's included.
      String[] parts = oldLine.split("\t");
      // Split the previous line if necessary.
      Optional<Position> endPos = reader.parseEndPosition(oldLine);
      if (endPos.isPresent() && excluded.before(endPos.get())) {
        if (!excluded.contig().equals(endPos.get().contig())) {
          throw new RuntimeException(
              "Unexpected (bug?): "
                  + "last record ends after target, but starts in a different contig?\n"
                  + oldLine
                  + "target: "
                  + excluded);
        }
        // We need to include the first half.
        parts[7] = "END=" + (excluded.pos() - 1);
        oldLine = String.join("\t", parts);
      }
      total += (oldLine.length() + 1);
      dest.write(ByteBuffer.wrap(oldLine.getBytes(UTF_8)));
      newline.position(0);
      newline.mark();
      newline.limit(1);
      dest.write(newline);
    }
    return total;
  }
}
