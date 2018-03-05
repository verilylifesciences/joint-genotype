#!/bin/bash

# Sharder depends on GenomeWarp (https://github.com/verilylifesciences/genomewarp)
# but genomewarp is not in the Maven repository. So to be able to compile
# Sharder you need to run this script.
# Alternatively, if you prefer, you can download genomewarp from GitHub and
# run `mvn install`.

mvn install:install-file \
  -Dfile=lib/verilylifesciences-genomewarp-1.0.0.jar \
  -DgroupId=com.verily.genomewarp \
  -DartifactId=verilylifesciences-genomewarp \
  -Dversion=1.0.0 \
  -Dpackaging=jar
