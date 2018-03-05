package com.verily.lifescience.genomics.wgs.sharder;

import static com.google.common.base.Preconditions.checkArgument;
import static com.google.common.collect.ImmutableList.toImmutableList;
import static java.nio.charset.StandardCharsets.UTF_8;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.ParseException;
import com.google.common.collect.ImmutableList;
import java.io.IOException;
import java.net.URI;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.stream.Stream;

/** Writes outFolder one shardNumber shardsTotal a collection shardsTotal GVCF files */
public class Main {

  public static void main(String[] args) throws IOException, java.text.ParseException, ParseException, InterruptedException {

    Options options = new Options()
      .addOption(Option.builder().required()
        .longOpt("shards_file")
        .desc("Path to TSV file that describes each shard").build())
      .addOption(Option.builder().type(Integer.class).required()
        .longOpt("shard_number")
        .desc("Index of the shard to cut (starts at 0)").build())
      .addOption(Option.builder().type(Integer.class).required()
        .longOpt("shards_total")
        .desc("Total number of shards for the input file (even though we only write one)").build())
      .addOption(Option.builder().required()
        .longOpt("vcf_files")
        .desc("File that lists the input VCF files (one per line). "
          + "They may be on GCS.").build())
      .addOption(Option.builder().required()
        .longOpt("mindex_files")
        .desc("File that list the mindex files that correspond to the VCF files (one per line). "
          + "They may be on GCS.").build())
      .addOption(Option.builder().required()
        .longOpt("reference")
        .desc("FASTA file for the reference. Must be on the local disk.").build())
      .addOption(Option.builder().required()
        .longOpt("output_folder")
        .desc("Folder for the output files. Their names are derived from the vcf file name.").build())
      .addOption(Option.builder().type(Integer.class).required()
        .longOpt("threads")
        .desc("Number of threads for IO. Also enables parallel init if > 1.").build())
      .addOption(Option.builder().type(Integer.class).required()
        .longOpt("metrics")
        .desc("Path for the metrics file. If unset, metrics are not saved.").build());

    CommandLine cmd = new DefaultParser().parse(options, args);

    ImmutableList<Path> vcfPaths = pathsInFile(cmd.getOptionValue("vcf_files"));
    ImmutableList<Path> outPaths =
        makeOutputPaths(cmd.getOptionValue("output_folder"), vcfPaths,
          (Integer)cmd.getParsedOptionValue("shard_number"),
          (Integer)cmd.getParsedOptionValue("shards_total"));
    new Sharder(
            makePath(cmd.getOptionValue("shards_file")),
            pathsInFile(cmd.getOptionValue("mindex_files")),
            vcfPaths,
            (Integer)cmd.getParsedOptionValue("shards_total"),
            outPaths,
            makePath(cmd.getOptionValue("metrics")),
            new Reference(cmd.getOptionValue("reference")))
        .shard((Integer)cmd.getParsedOptionValue("shard_number"),
          (Integer)cmd.getParsedOptionValue("threads"));
  }

  private static ImmutableList<Path> pathsInFile(String filename) throws IOException {
    Path p = makePath(filename);
    try (Stream<String> stream = Files.lines(p, UTF_8)) {
      return stream.map(Main::makePath).collect(toImmutableList());
    }
  }

  private static ImmutableList<Path> makeOutputPaths(
      String outputFolder, ImmutableList<Path> inputs, int shardNo, int shardTotal) {
    Path out = makePath(outputFolder);
    ImmutableList.Builder<Path> builder = ImmutableList.builder();
    for (Path input : inputs) {
      String inputName = input.getFileName().toString();
      String prefix = inputName;
      String suffix = "vcf";
      String keepExtension = ".g.vcf";
      if (inputName.endsWith(keepExtension)) {
        prefix = inputName.substring(0, inputName.length() - keepExtension.length());
        suffix = inputName.substring(inputName.length() - keepExtension.length() + 1);
      }
      builder.add(
          out.resolve(
              String.format(
                  "%s.shardNumber-%05d-of-%05d.%s", prefix, shardNo, shardTotal, suffix)));
    }
    return builder.build();
  }

  private static Path makePath(String path) {
    if (path.contains("://")) {
      // Only the URI overload loads alternative filesystems
      // such as the NIO GCS provider.
      return Paths.get(URI.create(path));
    } else {
      // The URI overload fails if the string isn't a URI
      // (e.g. when it's a normal path), so use the String overload.
      return Paths.get(path);
    }
  }
}
