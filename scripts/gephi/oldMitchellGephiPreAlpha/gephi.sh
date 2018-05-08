#!/bin/bash

#base=/net/eichler/vol21/projects/bac_assembly/nobackups/scripts/gephi
base=/net/eichler/vol2/home/mvollger/projects/abp/gephi
module load java_jdk/8u91

javac -cp .:$base:$base/gephi-toolkit-0.9.1-all.jar $base/Headless.java 

java -cp .:$base:$base/gephi-toolkit-0.9.1-all.jar Headless




