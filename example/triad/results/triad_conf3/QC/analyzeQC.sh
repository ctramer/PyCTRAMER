echo Analysis of CT states using TheoDORE software

export ORBKITPATH=/gpfsnyu/xspace/sungroup/software/orbkit
export PYTHONPATH=$PYHONPATH:$ORBKITPATH

WORKDIR=/gpfsnyu/xspace/projects/PCT-example/triad/results/triad_conf3/QC/
SRCDIR=/gpfsnyu/xspace/projects/PCT-example/triad/template
MYQC=/gpfsnyu/packages/qchem/4.4
jobname=TRIAD
system=triad_    
START=
END= 
email=netid@nyu.edu
module_loc=/gpfsnyu/xspace/sungroup/modules
functional=BNL
method=BNL
basis_set=SV
STATE_NUM=25
GIVEN_STRUC=conf3
LE_fragment=2
omega=157
TD_approx=True

theodir=/gpfsnyu/xspace/sungroup/software/TheoDORE_2.4
fraginfo='[[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60],[61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135],[136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207]]'  #'[[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60], [61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135], [136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207]]'
# fraginfo is generated by python module , integrate and automate in the future.

cd $WORKDIR

for structure  in $GIVEN_STRUC
do 

for n_states in  $STATE_NUM
do 

for w in $omega 
do

#rm -r cc${structure}_BNL_${n_states}states_w${w}

#cp -r cc${structure}_BNL cc${structure}_BNL_${n_states}states_w${w}

RUNDIR=/gpfsnyu/xspace/projects/PCT-example/triad/results/triad_conf3/QC/
#mkdir -p ${RUNDIR}
cd ${theodir} 
source  ${theodir}/setpaths.bash
cd $RUNDIR


# write input file for tddft analysis using theodore
cat > dens_ana.in <<EOF
rtype='qctddft'
rfile="${system}${structure}.out"
read_libwfa=True
TDA=$TD_approx
at_lists=${fraginfo}
eh_pop=1
coor_file='${system}${structure}.out'
coor_format='out'
prop_list=['Om', 'POS', 'PR', 'CT', 'COH', 'CTnt', 'RMSeh']
EOF

analyze_tden.py

cp tden_summ.txt tden_summ_tddft.txt

cp dens_ana.in dens_ana_tddft.in

cp ${system}${structure}.FChk ${system}${structure}.fchk

# write input file for fchk analysis using theodore
cat > dens_ana.in <<EOF
rtype='fchk'
rfile='${system}${structure}.fchk'
at_lists=${fraginfo}
Om_formula=1
eh_pop=1
comp_ntos=True
comp_dntos=False
jmol_orbitals=True
molden_orbitals=False
coor_file='${system}${structure}.out'
coor_format='out'
prop_list=['Om', 'POS', 'PR', 'CT', 'COH', 'CTnt', 'PRNTO', 'Z_HE', 'RMSeh']
EOF

analyze_tden.py

cp tden_summ.txt tden_summ_fchk.txt

mv dens_ana.in dens_ana_fchk.in

cd $WORKDIR
done
done
done
