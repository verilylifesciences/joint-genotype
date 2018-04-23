# Joint Genotyping

## ABOUT

A package to speed up GATK genotyping by sharding the inputs into tiny pieces.

## BUILDING

Follow the instructions in each folder.

## RUNNING

1. Run shards_picker to pick the tentative shard boundaries given your chosen
   number of shards.
   I suggest picking the shards such that each shard has total size on the order
   of the amount of memory available in your machines.

  The resulting output is the "shards file". Keep it locally, it's an input to
the Sharder.

2. Run mindexer on each input's index. This can be done in parallel.

  The resulting outputs are the "mindexes". Keep them on Google Cloud Storage.
Write a file that lists the mindexes in order. This is the "mindex file", it's
an input to the Sharder.

3. Run Sharder and GATK GenotypeGVCF for each shard. This can be done in parallel.

  The sharder outputs do not include headers so you'll have to put them back on
for GATK to work. Make sure that the sample name in the CHROM line is correct
for each sample.

4. Concatenate the GenotypeGVCF outputs. The resulting file is your final
   output.

  For concatenation, remove the headers from all but the first file.
