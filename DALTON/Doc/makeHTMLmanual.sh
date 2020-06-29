#!/bin/bash
#
if [ $# -gt 0 ]; then
   echo "--> Using user specified html directory: $1"
   HTMLDIR="$1"
else
   HTMLDIR="${DOC_DIRECTORY}/Dalton_HTMLmanual"
if [ -d "$HTMLDIR" ]; then
   echo "--> html output directory $HTMLDIR exists"
else
   echo "--> html directory $HTMLDIR does not exist - creating it ..."
   mkdir $HTMLDIR
fi

fi
if [ -d "$HTMLDIR" ]; then
   echo "--> html output directory $HTMLDIR exists"
else
   echo "--> html directory $HTMLDIR does not exist - creating it ..."
   mkdir $HTMLDIR
fi

DOC_DIRECTORY=${DOC_DIRECTORY:-.}
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

if [ ! -d $DOC_DIRECTORY ]; then
   echo '--> ERROR, Dalton document directory "' $DOC_DIRECTORY '" does not exist!'
   exit 1
fi

cp -t $HTMLDIR $DOC_DIRECTORY/*.tex
cp -t $HTMLDIR $DOC_DIRECTORY/*.bib

cd $HTMLDIR
HTMLDIR=`pwd`

if [ -z Master.tex ]; then
   echo 'ERROR, Dalton Master.tex does not exist in "' $HTMLDIR '".'
   exit 2
fi

# debug print:
#echo "current directory is `pwd` =? " $DOC_DIRECTORY
#echo "DALTON_VERSION is " $DALTON_VERSION
#echo "RELEASE_VERSION is " $RELEASE_VERSION
#echo "HTMLDIR is " $HTMLDIR

echo "--> Generation of .aux files with latex is needed for latex2html"

#
# Run latex2html
#
(latex Master; bibtex Master; latex Master; latex Master; makeindex Master; latex Master) \
>& dalton_latex.log

echo "--> Output from latex2html will be in dalton_latex2html.log"

latex2html -transparent -font_size 11pt -no_antialias -local_icons \
    -address "Dalton Manual - Release ${RELEASE_VERSION}" -dir $HTMLDIR Master.tex >& dalton_latex2html.log
mv Master.html Dalton_Master.html

# clean-up
rm *.aux *.ind *.blg *.bbl *.toc *.tex *.bib

echo "--> Dalton html manual is now created in $HTMLDIR; html root file is Dalton_Master.html"
#
# Post-processing from release 1.2.1
#
# sed 's/<PRE>/<\!--/' $HTMLDIR/footnode.html > $HTMLDIR/footnode.tmp
# sed 's/<\/PRE>/-->/' $HTMLDIR/footnode.tmp > $HTMLDIR/footnode.tmp2
# mv -f $HTMLDIR/footnode.tmp2 $HTMLDIR/footnode.html
# rm -f $HTMLDIR/footnode.tm*
# sed 's/http\:\/\/www.emsl.pnl.gov:2080\/forms\/basisform.html/<a target\=\"ext\" href=\"http\:\/\/www.emsl.pnl.gov:2080\/forms\/basisform.html\">http\:\/\/www.emsl.pnl.gov:2080\/forms\/basisform.html<\/a>/' $HTMLDIR/node118.html > $HTMLDIR/node118.tmp
# mv -f $HTMLDIR/node118.tmp $HTMLDIR/node118.html
# sed 's/http\:\/\/garm.teokem.lu.se\/MOLCAS\//<a target\=\"ext\" href=\"http\:\/\/garm.teokem.lu.se\/MOLCAS\/\">http\:\/\/garm.teokem.lu.se\/MOLCAS\/<\/a>/' $HTMLDIR/node118.html > $HTMLDIR/node118.tmp
# mv -f $HTMLDIR/node118.tmp $HTMLDIR/node118.html
# echo "BODY {background-color: #ffffff; color: #000000}" >> $HTMLDIR/Master.css
# echo "A:link {color: #6060bb}" >> $HTMLDIR/Master.css
# echo "A:visited {color: #999999}" >> $HTMLDIR/Master.css
