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