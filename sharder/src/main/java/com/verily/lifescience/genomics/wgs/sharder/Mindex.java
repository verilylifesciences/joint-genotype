package com.verily.lifescience.genomics.wgs.sharder;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.LongBuffer;
import java.nio.file.Files;
import java.nio.file.Path;
import org.checkerframework.checker.nullness.qual.EnsuresNonNull;
import org.checkerframework.checker.nullness.qual.MonotonicNonNull;

/**
 * Load and query the mindex.
 *
 * <p>Additionally this offers a prefetch feature (since we usually need to find at least two
 * consecutive shard boundaries).
 */
public class Mindex {
  /** Path to the mindex file. */
  private final Path path;
  /** First shard we have data for. */
  private int shard;
  /** Number of shards we have data for. */
  private final int prefetch;
  /** Cached mindex data. */
  // This syntax indicates a Nullable array (of non-null longs).
  // Monotonic means that once it's non-null, it stays that way.
  private long @MonotonicNonNull [] cachedOffsets;

  /**
   * Creates a Mindex object.
   *
   * <p>The methods you call on it will refer to the path specified here.
   *
   * <p>Prefetch is set to the default value of 3.
   *
   * <p>(I chose 3 before it means we won't need a second access even if you access two neighboring
   * shards boundaries (say, because you're looking for the beginning and end point of your shards),
   * even if your shards are 2 wide.)
   */
  public Mindex(Path path) throws IOException {
    this(path, 3);
  }

  /**
   * Creates a Mindex object.
   *
   * <p>The methods you call on it will refer to the path specified here. It will fetch the
   * specified number of elements at a time.
   */
  public Mindex(Path path, int prefetch) throws IOException {
    this.path = path;
    this.prefetch = prefetch;
  }

  /**
   * Returns the corresponding entry from the mindex.
   *
   * <p>It will be fetched from memory if cached, from disk otherwise.
   */
  public long get(int shardToGet) throws IOException {
    int index = shardToGet - this.shard;
    if (cachedOffsets == null || index < 0 || index >= this.cachedOffsets.length) {
      loadMindex(shardToGet);
      index = 0;
    }
    return cachedOffsets[index];
  }

  @EnsuresNonNull("this.cachedOffsets")
  private void loadMindex(int shard) throws IOException {
    LongBuffer mindexBytes = ByteBuffer.wrap(Files.readAllBytes(path)).asLongBuffer();
    this.shard = shard;
    int toFetch = Math.min(prefetch, mindexBytes.limit() - shard);
    cachedOffsets = new long[toFetch];
    for (int i = 0; i < toFetch; i++) {
      cachedOffsets[i] = mindexBytes.get(shard + i);
    }
  }

}
