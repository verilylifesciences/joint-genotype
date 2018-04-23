# MINDEXER

Mindexer takes a shards list and a VCF index file and subsamples the latter so
it only includes the offsets that matter for the shards list.

## BUILDING

```
$ mvn package
```

## TESTING

```
$ mvn test
```

## RUNNING

```
$ java -jar target/Mindexer-1.0-SNAPSHOT-jar-with-dependencies.jar \
    $shardsfile $indexfile $shardcount $output_file
```
