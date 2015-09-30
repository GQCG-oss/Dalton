#Author:  Patrick Ettenhuber, pettenhuber@gmai.com, 2015
#Purpose: Run a set of build instrcutctions stored in bash scripts in the folder hierarchy
# ./machinename/<nightly|continuous>/*bash
# which take as input the temporary folder <tmpd> where the script is executed. This script assumes
# that the build instruction will make a folder in <tmpd> which contains the name of the 
# buildinstrcution. If the test is set to "nightly" the build folder is removed immediately
#
#
#INPUT
#host name first
host=$1
#this has to be nightly or continuous, maybe an experimental can be added
buil=$2
#full path to the dir where the repos shall be stored
tmpd=$3
##################################

if [[ $buil != "continuous" && $buil != "nightly" ]]
then
   echo "Only continuous and nightly are supported"
   exit
fi

if [[ -z $tmpd ]]
then
   echo "Tmdpdir not given"
   exit
fi
if [[ ! -d $tmpd ]]
then
   echo "Tmpdir does not exist!!"
fi
#
sname=$(basename "$0")
pname="$0"
aname=$(pwd)
#
if [[ -z ${pname/"$sname"/} ]]
then
   p="$aname"/$host/$buil
else
   p=${pname/"$sname"/}$host/$buil
fi
#check if info for host exists
if [[ ! -d $p ]]
then
   echo "Information for host \"$host\" and build \"$buil\" not found"
   exit
fi
#
for name in $p/*bash
do
   echo $name $tmpd
   bash $name $tmpd
   #if nighly remove the repo from the tmp dir  
   if [[ "$buil" == "nightly" ]]
   then
      echo "rm -rf $tmpd/*"$(basename ${name/".bash"/})"*"
      rm -rf $tmpd/*"$(basename ${name/".bash"/})"*
   fi
done
