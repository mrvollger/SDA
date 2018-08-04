#!/usr/bin/env bash

cplex -c "read $1" "$2" "display solution variables *" > $1.sol
