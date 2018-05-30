#!/bin/bash
base="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

javac -cp .:$base:$base/gephi-toolkit-0.9.2-all.jar $base/Headless.java 

java -Xmx2g -cp .:$base:$base/gephi-toolkit-0.9.2-all.jar Headless $1 $2




