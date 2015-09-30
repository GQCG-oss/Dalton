sname=`basename "$0"`
pname="$0"
aname=$(pwd)

if [ -z ${pname/"$sname"/} ]
then
   p="$aname"/$1/$2
else
   p=${pname/"$sname"/}$1/$2
fi

##################################

for name in $p/*bash
do
   bash $name $3
done
