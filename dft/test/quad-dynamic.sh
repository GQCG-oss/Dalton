#!/bin/sh
# tests dynamics polariabilities.
# uses H=H_0 + e C -> d/de <<A;B>>(w) = <<A;B,C>>(w,0)
# for testing. 
# Pawel Salek, pawsa@theochem.kth.se, 2001.02.02

fld=2e-5
freq=${1:-0.5}
c_op=ZDIPLEN
run_linear=1
grid=1e-9
direct='.DIRECT'
thresholds="1e-16 1e-16 1e-17"
thresholds="1e-10 1e-12 1e-14"

# OTHER CONFIGURATION VARIABLES ------------------------------------
dalton=`pwd`/../../dalton.x
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
H    1.41421356237309504880 0.0 1.41421356237309504880
EOF
}

gen_CO() {
cat > MOLECULE.INP <<EOF
BASIS
AhlrichsVDZ
CO molecule (geometry from O. Christiansen et al. /CPL 305 (1999) 147-155.)
------------
    2    2  X  Y
        6.    1 
C    0.000  0.000   0.000
        8.    1 
O    0.000  0.000   2.132
EOF
}

gen_H2O() {
cat > MOLECULE.INP <<EOF
BASIS
daug-cc-pVTZ
H2O molecule (geometry from O. Christiansen et al. /CPL 305 (1999) 147-155.)
------------
    2    2  X  Y
        8.    1 
O    0.0000  0.0000   0.0000
        1.    1 
H    0.0000  1.42985  1.1071 
EOF
}

gen_HF() {
cat > MOLECULE.INP <<EOF
BASIS
taug-cc-pVTZ
HF molecule (geometry from  J Gauss et al. /CPL 296 (1998) 117-124.)
------------
    2    2  X  Y
        9.    1 
F    0.0000  0.0000   0.0000
        1.    1 
H    0.0000  0.0000   1.7328
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
.NOSUP
.DIPLEN
**WAVE FUNCTION
.DFT
$dftwg
*SCF INPUT
$hfocccommand
$hfocc
.THRESHOLD
1e-9
*DFT INPUT
.DFTTHR
$thresholds
.RADINT
$grid
*HAMILTONIAN
.FIELD
$field
 $c_op
*ORBITAL INPUT
.MOSTART
$mostart
**RESPONSE
*LINEAR
.DIPLNZ
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
.NOSUP
*READIN
**WAVE FUNCTIONS
.DFT
$dftw
*SCF INPUT
$hfocccommand
$hfocc
.THRESHOLD
1e-9
*DFT INPUT
.RADINT
$grid
.DFTTHR
$thresholds
*ORBITAL INPUT
.MOSTART
$mostart
**RESPONSE
*QUADRATIC
.DIPLNZ
.THCLR 
1e-8
.BFREQ
$freq_cnt
$freq
!.PRINT LEV
!301
*END OF INPUT
EOF

$mol
$dalton
}

extr_lin() {
awk '/^@.*ZDIPLEN/ {gsub("[dD]","e",$8);print $8}' $1
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
    if [ "$hfocc" = "" ]; then
        hfocccommand='!no enforced occupation'
        hfocc='!.'
    else
        hfocccommand='.HF OCCUPATION'
    fi
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
#run_test gen_LiH  1  "GGAKey"
run_test gen_LiH  1  "GGAKey PBEx=1 vwn=1"
rm SIRIUS.RST
run_test gen_CO   ""  "PBE"
rm SIRIUS.RST
run_test gen_HF    ""  "PBE"
run_test gen_HF    ""  "B3LYP"
rm SIRIUS.RST
run_test gen_CO2  11 B3LYP
run_test gen_CO2  11 CAMB3LYP
