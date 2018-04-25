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

/**
 * A dummy reference, for testing purposes.
 */
public class TestReference extends Reference {
  int queryCount;

  public TestReference() {
    super("");
  }

  @Override
  synchronized String getBase(String contig, long pos) {
    queryCount++;
    return String.valueOf("ACGT".charAt((int) (pos % 4)));
  }

  @Override
  int getQueryCount() {
    return queryCount;
  }
}