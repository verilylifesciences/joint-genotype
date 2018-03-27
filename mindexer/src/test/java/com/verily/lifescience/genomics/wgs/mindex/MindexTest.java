package com.verily.lifescience.genomics.wgs.mindex;

import static com.google.common.truth.Truth.assertThat;

import java.io.File;
import java.io.IOException;
import java.nio.LongBuffer;
import java.nio.file.Files;
import java.text.ParseException;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

@RunWith(JUnit4.class)
public class MindexTest {

  static final String TEST_DATA = "src/test/resources/";
  static final String POS = TEST_DATA + "positions.tsv";
  static final String IDX = TEST_DATA + "dummy.g.vcf.idx";
  static final int ROWS_IN_POS = 3;

  @Test
  public void goldenFile() throws IOException, ParseException {
    long[] expected = new long[] {0, 8225, Mindexer.PAST_EOF};
    LongBuffer lbuf = Mindexer.computeMindex(POS, IDX, ROWS_IN_POS).asLongBuffer();
    assertThat(lbuf.capacity()).isEqualTo(expected.length);
    for (int i = 0; i < expected.length; i++) {
      assertThat(expected[i]).isEqualTo(lbuf.get(i));
    }
  }

  @Test
  public void testGenMindexOutput() throws IOException, ParseException {
    File tempFile = File.createTempFile("testGenMindexOutput", ".mindex");
    innerGenMindexTest(tempFile);
  }

  @Test
  public void testGenMindexTruncates() throws IOException, ParseException {
    File tempFile = File.createTempFile("testGenMindexOutput", ".mindex");
    // larger than the target file size
    byte[] filler = new byte[8 * ROWS_IN_POS * 3];
    Files.write(tempFile.toPath(), filler);
    innerGenMindexTest(tempFile);
  }

  private void innerGenMindexTest(File tempFile) throws IOException, ParseException {
    long total = Mindexer.genMindex(POS, IDX, ROWS_IN_POS, tempFile.getAbsolutePath());
    assertThat(total).isEqualTo(8 * ROWS_IN_POS);
    assertThat(tempFile.exists()).isTrue();
    assertThat(tempFile.length()).isEqualTo(total);
  }
}
