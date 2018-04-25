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

import com.google.common.base.Preconditions;
import com.verily.genomewarp.utils.Fasta;
import java.io.IOException;
import org.checkerframework.checker.nullness.qual.MonotonicNonNull;

/**
 * Lets you query for a reference base. It'll only load up the reference on the first query.
 */
class Reference {

  private final String fastaFile;
  @MonotonicNonNull private Fasta fasta;
  private String cachedContig = "";
  private long cachedPos;
  private String cachedResult;
  private int queryCount;

  Reference(String fastaFile) {
    this.fastaFile = fastaFile;
  }

  /**
   * Gets the reference base at the specified position. Thread-safe.
   *
   * @param contig - the contig (chromosome) of interest.
   * @param oneBasedPosition - 1-based position of interest in the contig.
   * @throws IOException
   */
  synchronized String getBase(String contig, long oneBasedPosition) throws IOException {
    Preconditions.checkArgument(oneBasedPosition > 0);
    queryCount++;
    if (null == fasta) {
      // Creating the object fills the index, so we do it here instead of in the ctor,
      // this way we don't have to pay that cost if we never need to look up the base.
      fasta = new Fasta(fastaFile);
    }
    if (oneBasedPosition == cachedPos && contig.equals(cachedContig)) {
      return cachedResult;
    }
    // If an exception hits before the end of the method we don't want the cache
    // to appear filled.
    cachedPos = -1;
    cachedResult = fasta.get(contig, oneBasedPosition - 1, oneBasedPosition);
    cachedContig = contig;
    cachedPos = oneBasedPosition;
    return cachedResult;
  }

  /**
   * Returns the number of times {@link #getBase} was called.
   */
  int getQueryCount() {
    return queryCount;
  }

}
