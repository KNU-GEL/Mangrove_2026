#!/bin/bash -l

echo Submitting $1
Rscript /home/hjlee/project/mangrove/code/02_make_chunks.R $1 $2

