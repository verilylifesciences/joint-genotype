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

import static com.google.common.collect.ImmutableList.toImmutableList;
import static java.nio.charset.StandardCharsets.UTF_8;

import java.nio.file.spi.FileSystemProvider;
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
import com.google.cloud.storage.contrib.nio.CloudStorageFileSystemProvider;

/** Writes outFolder one shardNumber shardsTotal a collection shardsTotal GVCF files */
public class Main {

  public static void main(String[] args) throws IOException, java.text.ParseException, ParseException, InterruptedException {

    Options options = new Options()
      .addOption(Option.builder().required()
        .argName("shards_file")
        .longOpt("shards_file")
        .hasArg()
        .desc("Path to TSV file that describes each shard").build())
      .addOption(Option.builder().type(Integer.class).required()
        .argName("shard_number")
        .longOpt("shard_number")
        .hasArg()
        .desc("Index of the shard to cut (starts at 0)").build())
      .addOption(Option.builder().type(Integer.class).required()
        .argName("shards_total")
        .longOpt("shards_total")
        .hasArg()
        .desc("Total number of shards for the input file (even though we only write one)").build())
      .addOption(Option.builder().required()
        .argName("vcf_files")
        .longOpt("vcf_files")
        .hasArg()
        .desc("File that lists the input VCF files (one per line). "
          + "They may be on GCS.").build())
      .addOption(Option.builder().required()
        .argName("mindex_files")
        .longOpt("mindex_files")
        .hasArg()
        .desc("File that list the mindex files that correspond to the VCF files (one per line). "
          + "They may be on GCS.").build())
      .addOption(Option.builder().required()
        .argName("reference")
        .longOpt("reference")
        .hasArg()
        .desc("FASTA file for the reference. Must be on the local disk.").build())
      .addOption(Option.builder().required()
        .argName("output_folder")
        .longOpt("output_folder")
        .hasArg()
        .desc("Folder for the output files. Their names are derived from the vcf file name.").build())
      .addOption(Option.builder().type(Integer.class).required()
        .argName("threads")
        .longOpt("threads")
        .hasArg()
        .desc("Number of threads for IO. Also enables parallel init if > 1.").build())
      .addOption(Option.builder().type(Integer.class).required()
        .argName("metrics")
        .longOpt("metrics")
        .hasArg()
        .desc("Path for the metrics file.").build())
      .addOption(Option.builder().type(Integer.class)
        .argName("list_nio_providers")
        .longOpt("list_nio_providers")
        .desc("Prints the list of installed NIO providers.").build())
        .addOption(Option.builder().type(Boolean.class)
            .argName("verbose")
            .longOpt("verbose")
            // This is more of a debug option, really.
            .desc("Print a few more things along the way.").build())
      .addOption(Option.builder().type(Boolean.class)
          .argName("skip_writing")
          .longOpt("skip_writing")
          // This is more of a debug option, really.
          .desc("If specified, the cut is determined but no files are copied.").build());


    CommandLine cmd = new DefaultParser().parse(options, args);

    if (cmd.hasOption("list_nio_providers")) {
      System.out.println("Installed providers:");
      for (FileSystemProvider x : CloudStorageFileSystemProvider.installedProviders()) {
        System.out.println(x.getScheme()+"://");
      }
      System.out.println();
    }

    ImmutableList<Path> vcfPaths = pathsInFile(cmd.getOptionValue("vcf_files"));
    ImmutableList<Path> outPaths =
        makeOutputPaths(cmd.getOptionValue("output_folder"), vcfPaths,
          Integer.valueOf(cmd.getOptionValue("shard_number")),
          Integer.valueOf(cmd.getOptionValue("shards_total")));
    new Sharder(
            makePath(cmd.getOptionValue("shards_file")),
            pathsInFile(cmd.getOptionValue("mindex_files")),
            vcfPaths,
            Integer.valueOf(cmd.getOptionValue("shards_total")),
            outPaths,
            makePath(cmd.getOptionValue("metrics")),
            new Reference(cmd.getOptionValue("reference")))
        .setVerbose(cmd.hasOption("verbose"))
        .skipWriting(cmd.hasOption("skip_writing"))
        .shard(Integer.valueOf(cmd.getOptionValue("shard_number")),
          Integer.valueOf(cmd.getOptionValue("threads")));
  }

  private static ImmutableList<Path> pathsInFile(String filename) throws IOException {
    if (null == filename) {
      System.out.println("Warning: null filename to pathsInFile (shouldn't happen. Bug?)");
      return ImmutableList.<Path>of();
    }
    Path p = makePath(filename);
    try (Stream<String> stream = Files.lines(p, UTF_8)) {
      return stream
             .filter(line -> line != null)
             .map(Main::makePath)
             .collect(toImmutableList());
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
                  "%s.shard-%05d-of-%05d.%s", prefix, shardNo, shardTotal, suffix)));
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
