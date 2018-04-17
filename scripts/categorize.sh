#!/bin/bash

# this adds a fake catigory on the end
cat $1 | awk '{ print $1"\t"$2"\tall"}' > $2


