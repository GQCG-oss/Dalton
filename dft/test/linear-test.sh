#!/bin/sh
#dalton="mpirun -np 4 `pwd`/../../dalpar.x"; par='.PARALLEL'
dalton="`pwd`/../../bin/dalton.x"; par='!.PARALLEL'
tmp=/scratch/${USER}
runenergy=1
grid=1e-7

gen_LiH() {
cat > MOLECULE.INP <<EOF
BASIS
STO-3G
LiH STO-2G
--------
    2    0
        3.    1
Li     .0000     .0000     .0000
        1.    1
H1    1.41421356237309504880 1.41421356237309504880    .0000
EOF
}

gen_BeH() {
cat > MOLECULE.INP <<EOF
BASIS
Ahlrichs-VTZ
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
taug-cc-pVTZ
CO with STO-2G basis set
--------
    2    0
        6.    1
C     .0000     .0000     .0000
        8.    1
O    1.01421356 1.41421356  .0000
EOF
}

gen_CO2() {
cat > MOLECULE.INP <<EOF
BASIS
taug-cc-pVTZ
CO2 STO-2G or Ahlrichs-VDZ
--------
    2    0
        6.    1
C     .0000     .0000     .0000
        8.    2
O1    2.0 0.0 0.0
O2   -2.0 0.1 0.0
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
finite_run_field=$2

[ `pwd` = "$tmp" ] && rm -f AO* D* R* U*
cat > DALTON.INP <<EOF
**DALTON INPUT
.RUN PROPERTIES
.DIRECT
$par
**INTEGRALS
.NOSUP
.DIPLEN
**WAVE FUNCTIONS
.DFT
$dftwg
*SCF INPUT
.THRESHOLD
1e-8
$hfocccommand
$hfocc
*DFT INPUT
.RADINT
$grid
.DFTTHRESHOLD
0 0 0
.DFTELS
1
*HAMILTONIAN
.FIELD
$finite_run_field
YDIPLEN
*ORBITAL INPUT
.MOSTART
$mostart
.AO DELETE
1e-3
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
$par
**INTEGRALS
.NOSUP
.DIPLEN
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
.DFTTHRESHOLD
0 0 0
.DFTELS
1
*ORBITAL INPUT
.MOSTART
$mostart
.AO DELETE
1e-4
**RESPONSE
*LINEAR
.DIPLNY
.THCLR 
1e-3
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
        hfocccommand='.DOUBLY OCCUPIED'
    fi
    echo ""
    echo "*** RUNNING $mol with functional '$weights'"
    rm -f SIRIUS.RST; mostart=HUCKEL
    if [ "$runenergy" = 1 ]; then 
	finite_run "$weights" -$field
	#mostart=NEWORB
	min=`extr_dip DALTON.OUT`
	if [ -z "$min" ]; then echo "failed"; exit 1; fi
	finite_run "$weights" +$field
	plu=`extr_dip DALTON.OUT`
	if [ -z "$plu" ]; then echo "failed"; exit 1; fi
	str="-(($min)-($plu))/(2*($field))"
	echo "(min=$min plu=$plu)"
	echo "Finite field returns   :" `awk "BEGIN{print $str};"</dev/null`
    fi
    linear_run "$weights"
    awk '/@.*YDIPLEN/ {gsub("[dD]","e",$8);print "Linear response returns: " $8+0}' DALTON.OUT
}

run_test 1e-4 gen_CO2    ""  "B3LYP"  || exit 1
#run_test 2e-4 gen_BeH    "" "BP86"           || exit 1
run_test 1e-5 gen_CO    "" "CAMb3LYP"             || exit 1
#run_test 5e-5 gen_H2Oion "" "BLYP"            || exit 1
#run_test 2e-4 gen_H2Oion "" "B3LYP"           || exit 1
run_test 1e-5 gen_CO  "" "PBE"     || exit 1
exit 0
run_test 1e-4 gen_CO  "" B3LYP                || exit 1
run_test 1e-4 gen_LiH "" "LDA"                || exit 1
exit 0
#run_test 3e-4 gen_LiH    1  "GGAKey p86c=1"  || exit 1
run_test 2e-4 gen_BeH   ""  "GGAKey slater=1 pw91cl=1"  || exit 1
exit 0
run_test 1e-4 gen_BeH    "" "LDA"             || exit 1
run_test 5e-5 gen_H2Oion "" "BLYP"            || exit 1
run_test 2e-4 gen_H2Oion "" "B3LYP"           || exit 1
run_test 2e-4 gen_H2Oion "" "BP86"            || exit 1
run_test 2e-4 gen_H2Oion "" "B3P86"           || exit 1
exit 0
run_test 2e-4 gen_CO  "" LDA                  || exit 1
run_test 2e-4 gen_COs "" B3LYP                || exit 1
run_test 2e-4 gen_LiH "" "LDA"                || exit 1
run_test 2e-4 gen_LiH 1 "LDA"       || exit 1
run_test 2e-4 gen_LiH 1 "Example4"  || exit 1
run_test 2e-4 gen_LiH 2 "LDA"       || exit 1
run_test 2e-4 gen_LiH 2 "Example"   || exit 1
run_test 2e-4 gen_LiH 2 "Example2"  || exit 1
run_test 2e-4 gen_LiH 2 "Example3"  || exit 1
run_test 2e-4 gen_LiH 2 "Example4"  || exit 1
run_test 2e-4 gen_LiH 2 "B3LYP"     || exit 1
exit 0
