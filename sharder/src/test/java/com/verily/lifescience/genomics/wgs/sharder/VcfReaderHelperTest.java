package com.verily.lifescience.genomics.wgs.sharder;

import static com.google.common.truth.Truth.assertThat;

import com.google.common.collect.ImmutableMap;
import com.verily.lifescience.genomics.wgs.mindex.Position;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.channels.Channels;
import java.nio.channels.WritableByteChannel;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.ParseException;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

@RunWith(JUnit4.class)
public class VcfReaderHelperTest {

  private static final String TEST_DATA =
      "src/test/resources/";

  private static final Path TEST1_VCF = Paths.get(TEST_DATA, "test1.vcf");

  private static final ImmutableMap<String, Integer> NORMAL_CONTIGS =
      ImmutableMap.of("chr1", 0, "chr2", 1, "chr3", 2, "chr4", 3, "chr5", 4);

  @Test
  public void testEndCommonCase() throws IOException, ParseException {
    VcfReader reader = new VcfReader(TEST1_VCF, NORMAL_CONTIGS, new TestReference());
    // This is a normal record. Write the normal record before it.
    Position end = Position.of("chr1", 2, NORMAL_CONTIGS);
    ByteArrayOutputStream bytes = new ByteArrayOutputStream();
    WritableByteChannel got = Channels.newChannel(bytes);
    reader.advanceTo(end);
    VcfReaderHelper.saveLastRecord(reader, end, got);
    String expected =
        "chr1\t1\tdelfirst\tGATGATGATGAT\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0"
            + ":48:1:51,51\n";
    assertThat(bytes.toString()).isEqualTo(expected);
  }

  @Test
  public void testEndInReference() throws IOException, ParseException {
    VcfReader reader = new VcfReader(TEST1_VCF, NORMAL_CONTIGS, new TestReference());
    // This is in a reference record, trim it.
    Position end = Position.of("chr1", 400, NORMAL_CONTIGS);
    ByteArrayOutputStream bytes = new ByteArrayOutputStream();
    WritableByteChannel got = Channels.newChannel(bytes);
    reader.advanceTo(end);
    VcfReaderHelper.saveLastRecord(reader, end, got);
    String expected =
        "chr1\t379\t.\tA\t<NON_REF>\t.\t.\tEND=399\tGT:DP:GQ:MIN_DP:PL\t0/0:9:24:9:0,24,360\n";
    assertThat(bytes.toString()).isEqualTo(expected);
  }

  @Test
  public void testEndAtEndOfReference() throws IOException, ParseException {
    VcfReader reader = new VcfReader(TEST1_VCF, NORMAL_CONTIGS, new TestReference());
    // This is the last base in the reference record. Keep the whole thing.
    Position end = Position.of("chr1", 999, NORMAL_CONTIGS);
    ByteArrayOutputStream bytes = new ByteArrayOutputStream();
    WritableByteChannel got = Channels.newChannel(bytes);
    reader.advanceTo(end);
    VcfReaderHelper.saveLastRecord(reader, end, got);
    String expected =
        "chr1\t379\t.\tA\t<NON_REF>\t.\t.\tEND=999\tGT:DP:GQ:MIN_DP:PL\t0/0:9:24:9:0,24,360\n";
    assertThat(bytes.toString()).isEqualTo(expected);
  }

  @Test
  public void testEndInGap() throws IOException, ParseException {
    VcfReader reader = new VcfReader(TEST1_VCF, NORMAL_CONTIGS, new TestReference());
    // This is just after the last base in the reference record. We should keep the whole record.
    Position end = Position.of("chr1", 1000, NORMAL_CONTIGS);
    ByteArrayOutputStream bytes = new ByteArrayOutputStream();
    WritableByteChannel got = Channels.newChannel(bytes);
    reader.advanceTo(end);
    VcfReaderHelper.saveLastRecord(reader, end, got);
    String expected =
        "chr1\t379\t.\tA\t<NON_REF>\t.\t.\tEND=999\tGT:DP:GQ:MIN_DP:PL\t0/0:9:24:9:0,24,360\n";
    assertThat(bytes.toString()).isEqualTo(expected);
  }
}
