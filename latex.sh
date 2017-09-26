#!/bin/bash
module load latex2rtf/2.3.2a


pdflatex --draftmode $1
pdflatex $1

