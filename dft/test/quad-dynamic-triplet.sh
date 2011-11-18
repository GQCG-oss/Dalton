#!/bin/sh
# tests dynamics polariabilities.
# uses H=H_0 + e C -> d/de <<A;B>>(w) = <<A;B,C>>(w,0)
# for testing. 
# Pawel Salek, pawsa@theochem.kth.se, 2001.02.02

fld=2e-5
freq=${1:-0.5}
c_op=XDIPLEN
run_linear=1
grid=1e-6
direct='!.DIRECT'

# OTHER CONFIGURATION VARIABLES ------------------------------------
dalton=`pwd`/../../dalton.x
#dalton="/home/pawsa/dalton/head/dalton.x"
tmp=${TMPDIR:-/tmp}/${USER}

# MOLECULE GENERATION ----------------------------------------------

gen_LiH() {
cat > MOLECULE.INP <<EOF
BASIS
STO-2G
LiH STO-2G or AhlrichsVDZ
--------
    2    0
        2.    1
He     .0000     .0000     .0000
        1.    1
H    1.41421356237309504880 1.41421356237309504880    .0000
EOF
}

gen_CO() {
cat > MOLECULE.INP <<EOF
BASIS
cc-pVDZ
CO with STO-2G basis set
--------
    2    0  Z
        6.    1
C     .0000     .0000     .0000
        8.    1
O    1.41421356 1.41421356  .0000
EOF
}

gen_CO2() {
cat > MOLECULE.INP <<EOF
BASIS
cc-pVTZ
CO2 STO-2G or AhlrichsVDZ
--------
    2    0
        6.    1
C     .0000     .0000     .0000
        8.    2
O1    2.3 0.0 0.0
O2   -2.3 0.1 0.0
EOF
}

gen_C3H4() {
cat > MOLECULE.INP <<EOF
BASIS
DunningDZ
Testjob for ROA calculation: Allene (H2C=C=CH2) at optim. geometry SCF-DZ
DZ(Dunning)  Symmetrie C1
    2  0 0      
        6.    3
C1    0.0            0.0           0.0
C2    0.0	     2.4743846667  0.0   
C3    0.0           -2.4743846667  0.0   
        1.    4
H4    1.7367853003   3.5259125935  0.0     
H5   -1.7367853003   3.5259125935  0.0     
H6    0.0           -3.5259125935  1.7367853003
H7    0.0           -3.5259125935 -1.7367853003
EOF
}


# LINEAR FINITE -----------------------------------------------------

finite_run() {
dftwg="$1"
field="$2"

freq_cnt=`echo $freq | wc -w`
if [ -f SIRIUS.RST ]; then mostart=NEWORB; else mostart=HUCKEL; fi

[ `pwd` = "$tmp" ] && rm -f A* D* U* R*
cat > DALTON.INP <<EOF
**DALTON INPUT
.RUN RESPONSE
$direct
**INTEGRALS
.DIPLEN
.NOSUP
**WAVE FUNCTIONS
.HF
*HF INPUT
.DFTFUNCTIONAL
$dftwg
.RADINT
$grid
.HF OCCU
$hfocc
.THRESHOLD
1e-9
*HAMILTONIAN
.FIELD
$field
 $c_op
*ORBITAL INPUT
.MOSTART
$mostart
**RESPONSE
.TRPFLG
*LINEAR
.PROPRT
YDIPLEN
.THCLR 
1e-8
.FREQUENCIES
$freq_cnt
$freq
*END OF INPUT
EOF

$mol
$dalton
}

# QUADRATIC --------------------------------------------------------
quadratic_run() {
dftw="$1"

freq_cnt=`echo $freq | wc -w`
if [ -f SIRIUS.RST ]; then mostart=NEWORB; else mostart=HUCKEL; fi

[ `pwd` = "$tmp" ] && rm -f A* D* U* R*
cat > DALTON.INP <<EOF
**DALTON INPUT
.RUN RESPONSE
$direct
**INTEGRALS
.DIPLEN
.NOSUP
**WAVE FUNCTIONS
.HF
*HF INPUT
.DFTFUNCTIONAL
$dftw
.RADINT
$grid
.HF OCCU
$hfocc
*ORBITAL INPUT
.MOSTART
$mostart
**RESPONSE
.TRPFLG
*QUADRATIC
.E3TEST
.PRINT
 11
.THCLR 
1e-7
.BFREQ
$freq_cnt
$freq
!.PRINT LEV
!301
.ISPABC
1 1 0
.APROP
YDIPLEN
.BPROP
YDIPLEN
.CPROP
$c_op
*END OF INPUT
EOF

$mol
$dalton
}

extr_lin() {
awk '/^@YDIPLEN/ {gsub("[dD]","e",$3);print $3}' $1
}
checkfound() {
#echo "check found: $1"
if [ -z "$1" ] ; then
    echo $2
    exit 1
fi
}

run_test() {
    mol="$1"
    hfocc="$2"
    dft_weight="$3"
    echo "===================================================="
    echo "Running calculations with DFT weights=$dft_weight molecule: $mol"
    echo "for frequencies $freq"
    #rm -f SIRIUS.RST
    if [ "$run_linear" != 0 ]; then 
	finite_run "$dft_weight" -$fld
	min=`extr_lin DALTON.OUT`
	checkfound "$min" "Finite field (-field) run returned nothing. I quit."
	finite_run "$dft_weight" +$fld
	plu=`extr_lin DALTON.OUT`
	checkfound "$plu" "Finite field run returned nothing. I quit."
	
	echo "min=$min plu=$plu 2fld=$fld"
	echo "The expected value for the polarizability:" \
	    `awk "BEGIN{print -(($min)-($plu))/(2*($fld))}" </dev/null`
    else
	echo "Running directly the quadratic response"
    fi
    quadratic_run "$dft_weight"
    awk '/^@.*beta/{print "Quadratic response returns               : "$10}'\
	DALTON.OUT
}

# main -----------------------------------------------------------
mkdir -p $tmp || exit 1

cd $tmp || exit 1

rm SIRIUS.RST
run_test gen_LiH  1  "GGAKey hf=.7"
run_test gen_LiH  1  Dirac
#run_test gen_LiH  1  Becke
run_test gen_LiH  1  "GGAKey becke=1 dirac=1 vwn=1"
run_test gen_LiH  1  "LDA"
rm SIRIUS.RST
#run_test gen_CO   7  LDA
#run_test gen_CO   7  "GGAKey becke=1 dirac=1 vwn=1"
run_test gen_CO   7  BLYP
exit 0
run_test gen_CO   7  "GGAKey hf=0.2 dirac=0.8 vwn=0.19 becke=0.72 lyp=0.81"
run_test gen_CO   7  B3LY
run_test gen_CO   7  Example2
exit 0
run_test gen_CO   7  Example3
run_test gen_CO   7  Example4
run_test gen_CO   7  Becke
run_test gen_CO2  11 GGAKey
run_test gen_CO2  11 Example4
run_test gen_CO2  11 Becke
run_test gen_LiH  1  "GGA 0 1 0 0 0"
run_test gen_LiH  1  "0 0 1 0 0"
run_test gen_CO2  11 "1 0 0 0 0"
run_test gen_CO2  11 "0 1 1 0 0"
run_test gen_C3H4 11 "1 0 0 0 0"
run_test gen_C3H4 11 "0 1 1 0 0"
