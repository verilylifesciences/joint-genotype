# SHARDER

The Sharder, given a shard number, will grab that shard from every input VCF.
The inputs may be on Google Cloud Storage.

So for example if you ask for shard 0, you will get the first part of each of
the inputs.

This is useful for the following reasons:

 - a shard is smaller than the full input, so GenotypeGVCF will finish sooner.
 - no need to copy potentially large inputs to more than one worker computer,
   each worker only grabs the part they need.

Sharder uses the shards files as a guide, but may choose to move the shard
boundaries a bit in order to ensure that no shard cuts a deletion variant.
Sharder guarantees that each shard starts immediately after where the previous
one left off.

The options are listed in Main.java. Here they are again:

```
    --shards_file: Path to TSV file that describes each shard
    --shard_number: Index of the shard to cut (starts at 0)
    --shards_total: Total number of shards for the input file (even though we only write one)
    --vcf_files: File that lists the input VCF files (one per line).
                 They may be on GCS.
    --mindex_files: File that list the mindex files that correspond to the VCF files (one per line).
                    They may be on GCS.
    --reference: FASTA file for the reference. Must be on the local disk.
    --output_folder: Folder for the output files. Their names are derived from the vcf file name.
    --threads: Number of threads for IO. Also enables parallel init if > 1.
    --metrics: Path for the metrics file.
```

## BUILDING

First build the non-Maven dependencies:

```
$ ./install-genomewarp.sh
$ cd ../mindexer
$ mvn install
$ cd ../sharder
```

Then build Sharder:

```
$ mvn package
```

## TESTING

To test, follow the building steps above at least once, then:

```
$ mvn test
```

## RUNNING

```
$ java -jar target/Sharder-1.0-SNAPSHOT-jar-with-dependencies.jar
```

(see src/main/java/com/verily/lifescience/genomics/wgs/sharder/Main.java for the list of parameters)

