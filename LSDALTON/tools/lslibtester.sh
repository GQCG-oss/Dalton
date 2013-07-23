log=$1

if [ `uname` = Linux ]; then
   GREP="egrep -a"
else
   GREP="egrep"
fi

typeset -a STRINGS
STRINGS[ 1]="                       One-electron AO-matrix matrix: RMS and index-weighted sum"
STRINGS[ 2]="                       One-electron AO-matrix matrix: RMS and index-weighted sum"
STRINGS[ 3]="                   Coulomb AO-matrix matrix number 1: RMS and index-weighted sum"
STRINGS[ 4]="                   Coulomb AO-matrix matrix number 2: RMS and index-weighted sum"
STRINGS[ 5]="                  Exchange AO-matrix matrix number 1: RMS and index-weighted sum"
STRINGS[ 6]="                  Exchange AO-matrix matrix number 2: RMS and index-weighted sum"
STRINGS[ 7]="                        XC AO-matrix matrix number 1: RMS and index-weighted sum"
STRINGS[ 8]="                        XC AO-matrix matrix number 2: RMS and index-weighted sum"
STRINGS[ 9]="                   Fock/KS AO-matrix matrix number 1: RMS and index-weighted sum"
STRINGS[10]="                   Fock/KS AO-matrix matrix number 2: RMS and index-weighted sum"
STRINGS[11]="                                         nn gradient: RMS and index-weighted sum"
STRINGS[12]="                                    kinetic gradient: RMS and index-weighted sum"
STRINGS[13]="                                         ne gradient: RMS and index-weighted sum"
STRINGS[14]="                                 1-electron gradient: RMS and index-weighted sum"
STRINGS[15]="                                    Coulomb gradient: RMS and index-weighted sum"
STRINGS[16]="                                   Exchange gradient: RMS and index-weighted sum"
STRINGS[17]="                                         XC gradient: RMS and index-weighted sum"
STRINGS[18]="                                 2-electron gradient: RMS and index-weighted sum"
STRINGS[19]="                       Reorthonormalization gradient: RMS and index-weighted sum"
STRINGS[20]="        Coulomb AO-matrix matrix from eri .Mulliken.: RMS and index-weighted sum"
STRINGS[21]="       Exchange AO-matrix matrix from eri .Mulliken.: RMS and index-weighted sum"
STRINGS[22]="           Coulomb AO-matrix matrix from eri .Dirac.: RMS and index-weighted sum"
STRINGS[23]="          Exchange AO-matrix matrix from eri .Dirac.: RMS and index-weighted sum"
STRINGS[24]="           Coulomb gradient from diff-eri .Mulliken.: RMS and index-weighted sum"


typeset -a TOLERANCES
TOLERANCES[ 1]=1e-9
TOLERANCES[ 2]=1e-9
TOLERANCES[ 3]=1e-9
TOLERANCES[ 4]=1e-9
TOLERANCES[ 5]=1e-9
TOLERANCES[ 6]=1e-9
TOLERANCES[ 7]=1e-9
TOLERANCES[ 8]=1e-9
TOLERANCES[ 9]=1e-9
TOLERANCES[10]=1e-9
TOLERANCES[11]=1e-9
TOLERANCES[12]=1e-9
TOLERANCES[13]=1e-9
TOLERANCES[14]=1e-9
TOLERANCES[15]=1e-9
TOLERANCES[16]=1e-9
TOLERANCES[17]=1e-9
TOLERANCES[18]=2e-9
TOLERANCES[19]=1e-9
TOLERANCES[20]=1e-9
TOLERANCES[21]=1e-9
TOLERANCES[22]=1e-9
TOLERANCES[23]=1e-9
TOLERANCES[24]=1e-9

passed=1
for i in {1..24}
do 
   STRING=${STRINGS[$i]}
   TOLERANCE=${TOLERANCES[$i]}
   LOGSTR=`$GREP "${STRING}" $log | cut -c 81-`
   base=`echo ${log%.*}`
   REFSTR=`$GREP "${STRING}" ${base}.ref | cut -c 81-`
   LOGRMS=`echo $LOGSTR | awk '{printf("%25.15f",$1)}'`
   LOGIWS=`echo $LOGSTR | awk '{printf("%25.15f",$2)}'`
   REFRMS=`echo $REFSTR | awk '{printf("%25.15f",$1)}'`
   REFIWS=`echo $REFSTR | awk '{printf("%25.15f",$2)}'`

   #Test for NaN
   isnan=`echo $LOGRMS | awk '{print ( $1 == "nan" ? 1 : 0 )}'`
   if [ $isnan -eq 1 ]
   then
     echo "RSM error in, NaN in" $STRING
     passed=0
   fi
   isnan=`echo $LOGIWS | awk '{print ( $1 == "nan" ? 1 : 0 )}'`
   if [ $isnan -eq 1 ]
   then
     echo "IWS error in, NaN in" $STRING
     passed=0
   fi

   ABSDIFFRMS=`echo $LOGRMS $REFRMS | awk '{printf("%25.15f",( $1-$2 >=0 ? $1-$2 : $2-$1 ))}'`
   ABSDIFFIWS=`echo $LOGIWS $REFIWS | awk '{printf("%25.15f",( $1-$2 >=0 ? $1-$2 : $2-$1 ))}'`

   passrms=`echo $TOLERANCE $ABSDIFFRMS | awk '{print ( $1-$2 >=0 ? 1 : 0 )}'`
   passiws=`echo $TOLERANCE $ABSDIFFIWS | awk '{print ( $1-$2 >=0 ? 1 : 0 )}'`

   if [ $passrms -eq 0 ]
   then
     echo "RSM error in" $STRING $LOGRMS $REFRMS $TOLERANCE
     passed=0
   fi
   if [ $passiws -eq 0 ]
   then
     echo "IWS error in" $STRING $LOGIWS $REFIWS $TOLERANCE
     passed=0
   fi
done

completed=`$GREP "\*\*\* LSlib tester completed \*\*\*" $log | wc -l`
if [ $completed -ne 1 ]
then
   echo "LSlib tester did not complete properly"
   passed=0
fi

if [ $passed -eq 1 ]
then
  echo "TEST ENDED PROPERLY"
  exit 0
else
  echo "THERE IS A PROBLEM in testcase " $log
  exit 1
fi
