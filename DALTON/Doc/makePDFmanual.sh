#!/bin/bash
#
RELEASE_VERSION=`cat ../../VERSION`
DALTON_MANUAL=Dalton_${RELEASE_VERSION}_manual.pdf
pdflatex  master
bibtex    master
pdflatex  master
pdflatex  master
makeindex master
pdflatex  master
mv master.pdf $DALTON_MANUAL
echo "$DALTON_MANUAL script finished"
#
