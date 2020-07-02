#!/bin/bash
#
if [ $# -gt 0 ]; then
   echo "--> Using user specifiend PDF directory: $1"
fi

DOC_DIRECTORY=${DOC_DIRECTORY:-./}
# Declare variables
#
if [ -n "$DALTON_VERSION" ]; then
   RELEASE_VERSION=$DALTON_VERSION
   echo "--> User specified release version: $RELEASE_VERSION"
elif [ -n "$RELEASE_VERSION" ]; then
   echo "--> User specified release version: $RELEASE_VERSION"
else
   DALTON_VERSION=${DALTON_VERSION:-../../VERSION}
   RELEASE_VERSION=`cut -d ' ' -f 1 ${DALTON_VERSION}`
fi

PDFDIR=${1:-${DOC_DIRECTORY}../Dalton_PDFmanual}
if [ -d "$PDFDIR" ]; then
   echo "--> pdf output directory $PDFDIR exists"
else
   echo "--> pdf output directory $PDFDIR does not exist - creating it ..."
   mkdir $PDFDIR
fi
DALTON_MANUAL=${RELEASE_VERSION}_manual.pdf

if [ ! -d $DOC_DIRECTORY ]; then
   echo '--> ERROR, Dalton document directory "' $DOC_DIRECTORY '" does not exist!'
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

echo "--> Output from making of pdf manual will be in $PDFDIR/dalton_pdfmanual.log"

PDFOUT="-halt-on-error -output-directory $PDFDIR"

cp $DOC_DIRECTORY/*.bib $PDFDIR # needed for bibtex

(pdflatex $PDFOUT Master; (cd $PDFDIR; bibtex Master); pdflatex $PDFOUT Master; pdflatex $PDFOUT Master; makeindex $PDFOUT Master; pdflatex $PDFOUT Master) \
>& $PDFDIR/dalton_pdfmanual.log

if [ $? -ne 0 ]; then
   echo "--> latex error, see $PDFDIR/dalton_pdfmanual.log"
   exit 2
fi

mv $PDFDIR/Master.pdf $PDFDIR/$DALTON_MANUAL
echo "--> $DALTON_MANUAL script finished, manual is saved in $PDFDIR"
#
