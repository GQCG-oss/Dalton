#!/bin/sh
dalton=`pwd`/../../dalton.x
tmp=/tmp/${USER}
fld=2e-4
runenergy=1
grid=1e-5

gen_LiH() {
cat > MOLECULE.INP <<EOF
BASIS
STO-2G
LiH STO-2G
--------
    2    0
        3.    1
Li     .0000     .0000     .0000
        1.    1
H1    1.41421356237309504880 1.41421356237309504880    .0000
EOF
}

gen_O2() {
cat > MOLECULE.INP <<EOF
BASIS
cc-pVDZ
Test1 i test-suite: Geometrioptimering med symmetri, med beregning av
egenskaper i foerste og siste punkt. Generelt kontraherte basissett
    1    0   
        8.    2
O     0.00000000     0.00000000     1.1405
O     0.00000000     0.00000000    -1.1405
EOF
}



gen_BeH() {
cat > MOLECULE.INP <<EOF
BASIS
AhlrichsVTZ
BeH AhlrichsVDZ
--------
    2    0
        4.    1
Be   .0000     .0000     .0000
        1.    1
H1    1.41421356237309504880 1.41421356237309504880    .0000
EOF
}

gen_H2Oion() {
cat > MOLECULE.INP <<EOF
BASIS
aug-cc-pVDZ
Water                                                                   
cc-pVDZ basis
    2  1 2  X  Z
        8.    1
O      .0000   -.2249  0.0000000000
        1.    1                                             
H     1.4523    .8996   .0000000000
FINISH
EOF
}

gen_CO() {
cat > MOLECULE.INP <<EOF
BASIS
cc-pVDZ
CO with STO-2G basis set
--------
    2    0
        6.    1
C     .0000     .0000     .0000
        8.    1
O    1.41421356 1.41421356  .0000
EOF
}

gen_COs() {
cat > MOLECULE.INP <<EOF
BASIS
aug-cc-pVDZ
CO with STO-2G basis set
--------
    2    1  Z
        6.    1
C     .0000     .0000     .0000
        8.    1
O    1.41421356 1.41421356  .0000
EOF
}

gen_HF() {
cat > MOLECULE.INP <<EOF
BASIS
taug-cc-pVTZ
HF molecule (geometry from  J Gauss et al. /CPL 296 (1998) 117-124.)
------------
    2    2  X  Z
        9.    1 
F    0.0000  0.0000   0.0000
        1.    1 
H    0.0000  1.7328   0.0000   
EOF
}

# FINITE -----------------------------------------------------------

finite_run() {
dftwg="$1"
field=$2

[ `pwd` = "$tmp" ] && rm -f AO* D* R* U*
cat > DALTON.INP <<EOF
**DALTON INPUT
.RUN PROPERTIES
.DIRECT
**INTEGRALS
.NOSUP
.DIPLEN
**WAVE FUNCTIONS
.HF
*HF INPUT
.DFTFUNCTIONAL
$dftwg
.THRESHOLD
1e-9
.HF OCC
 2 
.HSROHF
 2 
.RADINT
$grid
.DFTTHRESHOLD
0 0 0
*HAMILTONIAN
.FIELD
$field
 YDIPLEN
*ORBITAL INPUT
.MOSTART
$mostart
*END OF INPUT
EOF

$mol
$dalton
}

# LINEAR -----------------------------------------------------------
linear_run() {
dftw="$1"

[ `pwd` = "$tmp" ] && rm -f AO* D* R* U*
cat > DALTON.INP <<EOF
**DALTON INPUT
.RUN RESPONSE
.DIRECT
**INTEGRALS
.NOSUP
.DIPLEN
**WAVE FUNCTIONS
.HF
*HF INPUT
.DFTFUNCTIONAL
$dftw
.HF OCC
 2                 
.HSROHF
 2               
.RADINT
$grid
.THRESHOLD
1e-9
.DFTTHRESHOLD
0 0 0
*ORBITAL INPUT
.MOSTART
$mostart
**RESPONSE
*LINEAR
.DIPLNY
.THCLR 
1e-7
*END OF INPUT
EOF

$mol
$dalton
}

extr_dip() {
perl -e 'undef $/;$_=<>;print(($_ =~/Dipole moment components.*?y\s+(-?[\d.]+)/gs),"\n");'  $1
}
# main -----------------------------------------------------------
mkdir -p $tmp || exit 1

cd $tmp || exit 1

run_test() {
    field="$1"
    mol="$2"
    hfocc="$3"
    weights="$4"
    if [ "$hfocc" = "" ]; then
        hfocccommand='!no enforced occupation'
        hfocc='!.'
    else
        hfocccommand='.HF OCCUPATION'
    fi
    echo ""
    echo "*** RUNNING $mol with functional '$weights'"
    rm -f SIRIUS.RST; mostart=HUCKEL
    if [ "$runenergy" = 1 ]; then 
	finite_run "$weights" -$fld
	mostart=NEWORB
	min=`extr_dip DALTON.OUT`
	if [ -z "$min" ]; then echo "failed"; exit 1; fi
	finite_run "$weights" +$fld
	plu=`extr_dip DALTON.OUT`
	if [ -z "$plu" ]; then echo "failed"; exit 1; fi
	str="-(($min)-($plu))/(2*($fld))"
	echo "(min=$min plu=$plu)"
	echo "Finite field returns   :" `awk "BEGIN{print $str};"</dev/null`
    fi
    linear_run "$weights"
    awk '/@YDIPLEN/ {gsub("[dD]","e",$3);print "Linear response returns: " $3+0}' DALTON.OUT
}

# -------------------------------------------------------------------
run_test 1e-4 gen_BeH "" "GGAKey hf=1"  || exit 1
run_test 1e-4 gen_BeH "" "LDA"           || exit 1
run_test 1e-4 gen_BeH "" "BLYP"          || exit 1
run_test 1e-4 gen_BeH "" "B3LYP"         || exit 1
exit 0
