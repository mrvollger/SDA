#!/bin/bash
jar tf gephi-toolkit-0.9.1-all.jar | grep "$1" | sed 's/\//./g' | sed 's/^/import /g'
