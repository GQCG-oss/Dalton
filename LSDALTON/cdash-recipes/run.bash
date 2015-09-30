#INPUT
#host name first
host=$1
#this has to be nightly or continuous, maybe an experimental can be added
buil=$2
#full path to the dir where the repos shall be stored
tmpd=$3
##################################

sname=`basename "$0"`
pname="$0"
aname=$(pwd)

if [ -z ${pname/"$sname"/} ]
then
   p="$aname"/$host/$buil
else
   p=${pname/"$sname"/}$host/$buil
fi

for name in $p/*bash
do
   bash $name $tmpd
   #if nighly remove the repo from the tmp dir  
   if [ "$buil" == "nightly" ]
   then
      rm -rf "$tmpd"/*$name*
   fi
done
