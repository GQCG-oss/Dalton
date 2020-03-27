#!/bin/bash
#
DALTON_VERSION=${DALTON_VERSION:-../../VERSION}
DOC_DIRECTORY=${DOC_DIRECTORY:-./}
RELEASE_VERSION=`cut -d ' ' -f 1 ${DALTON_VERSION}`
DALTON_MANUAL=${RELEASE_VERSION}_manual.pdf
if [ ! -d $DOC_DIRECTORY ]; then
   echo 'ERROR, Dalton document directory "' $DOC_DIRECTORY '" does not exist!'
   exit 1
fi

cd $DOC_DIRECTORY

if [ -z Master.tex ]; then
   echo 'ERROR, Dalton Master.tex does not exist in "' $DOC_DIRECTORY '".'
   exit 2
fi

# debug print:
#echo "current directory is `pwd` =? " $DOC_DIRECTORY
#echo "DALTON_VERSION is " $DALTON_VERSION
#echo "RELEASE_VERSION is " $RELEASE_VERSION

echo "Output from making of pdfmanual will be in Dalton_pdfmanual.log"

(pdflatex  Master; bibtex    Master; pdflatex  Master; pdflatex  Master; makeindex Master; pdflatex  Master) \
>& dalton_pdfmanual.log

mv Master.pdf $DALTON_MANUAL
echo "$DALTON_MANUAL script finished"
#
