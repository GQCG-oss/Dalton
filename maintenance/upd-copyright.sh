#!/bin/tcsh
#
# hjaaj Aug. 2014
#
# Utility used to update all DALTON copyright statements
# from DALTON2013 to DALTON2014 as follows:
#
#   cd DALTON
#   ../maintenance/upd-copyright.sh '*/*.F'
#   ../maintenance/upd-copyright.sh '*/*.F90'
#   ../maintenance/upd-copyright.sh '*/*.c'
#
# OBS! remember to update the "/2013/2014/" next year!
#
foreach ff ( $1 )
   echo Updating $ff
   sed -i '' -e '1,20s/2013/2014/g' $ff
end
