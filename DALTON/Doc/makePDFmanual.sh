#!/bin/bash
#
RELEASE_VERSION=`cat ../../VERSION`
DALTON_MANUAL=Dalton_${RELEASE_VERSION}_manual.pdf
pdflatex  Master
bibtex    Master
pdflatex  Master
pdflatex  Master
makeindex Master
pdflatex  Master
mv Master.pdf $DALTON_MANUAL
echo "$DALTON_MANUAL script finished"
#
