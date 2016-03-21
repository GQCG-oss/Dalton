#!/bin/bash
#
RELEASE_VERSION=`cut -d ' ' -f 1 ../../VERSION`
DALTON_MANUAL=${RELEASE_VERSION}_manual.pdf
pdflatex  Master
bibtex    Master
pdflatex  Master
pdflatex  Master
makeindex Master
pdflatex  Master
mv Master.pdf $DALTON_MANUAL
echo "$DALTON_MANUAL script finished"
#
