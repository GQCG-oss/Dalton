#!/bin/bash
#
# Declare variables
#
HTMLDIR="HTMLmanual"
export HTMLDIR
if [ -d "$HTMLDIR" ]; then
   echo "directory $HTMLDIR exists"
else
   echo "directory $HTMLDIR does not exist - creating it ..."
   mkdir $HTMLDIR
fi
#
# Run latex2html
#
latex2html -transparent -font_size 11pt -no_antialias -local_icons \
    -address "Dalton Manual - Release 1.2.1" -dir $HTMLDIR Master.tex
#
# Post-processing
#
sed 's/<PRE>/<\!--/' $HTMLDIR/footnode.html > $HTMLDIR/footnode.tmp
sed 's/<\/PRE>/-->/' $HTMLDIR/footnode.tmp > $HTMLDIR/footnode.tmp2
mv -f $HTMLDIR/footnode.tmp2 $HTMLDIR/footnode.html
rm -f $HTMLDIR/footnode.tm*
sed 's/http\:\/\/www.emsl.pnl.gov:2080\/forms\/basisform.html/<a target\=\"ext\" href=\"http\:\/\/www.emsl.pnl.gov:2080\/forms\/basisform.html\">http\:\/\/www.emsl.pnl.gov:2080\/forms\/basisform.html<\/a>/' $HTMLDIR/node118.html > $HTMLDIR/node118.tmp
mv -f $HTMLDIR/node118.tmp $HTMLDIR/node118.html
sed 's/http\:\/\/garm.teokem.lu.se\/MOLCAS\//<a target\=\"ext\" href=\"http\:\/\/garm.teokem.lu.se\/MOLCAS\/\">http\:\/\/garm.teokem.lu.se\/MOLCAS\/<\/a>/' $HTMLDIR/node118.html > $HTMLDIR/node118.tmp
mv -f $HTMLDIR/node118.tmp $HTMLDIR/node118.html
echo "BODY {background-color: #ffffff; color: #000000}" >> $HTMLDIR/Master.css
echo "A:link {color: #6060bb}" >> $HTMLDIR/Master.css
echo "A:visited {color: #999999}" >> $HTMLDIR/Master.css
