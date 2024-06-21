echo Preparing files for QC calculation

WORKDIR=/gpfsnyu/xspace/projects/PCT-example/triad/results/triad_conf3/QC/ 
SRCDIR=/gpfsnyu/xspace/projects/PCT-example/triad/template
MYQC=/gpfsnyu/packages/qchem/4.4
partition=argon  # parallel
core=96
jobname=PCT-triad
system=triad_
email=netid@nyu.edu
module_loc=/gpfsnyu/xspace/sungroup/modules
functional=
basis_set=SV
START=
END=
STATE_NUM=25
GIVEN_STRUC=conf3
D_range=61-207
B_range=61-135
A_range=1-60
method=BNL
omega=157
natoms=
seg0=
seg1=
seg2= 
frag=
omega=157
TD_approx=false

cd $WORKDIR
for structure in $GIVEN_STRUC 
do
#for structure in {0..2..1}    #Bent Linear ##{4..4..1} 
for n_states in $STATE_NUM 
do 

for w in $omega  
do

#rm -r cc${structure}_BNL_${n_states}states_w${w}

#cp -r cc${structure}_BNL cc${structure}_BNL_${n_states}states_w${w}

#RUNDIR=/gpfsnyu/xspace/projects/PCT-example/triad/results/triad_conf3/QC/${WORKDIR}/${system}${structure}_${functional}_${n_states}states_w${w}
RUNDIR=/gpfsnyu/xspace/projects/PCT-example/triad/results/triad_conf3/QC/
mkdir -p ${RUNDIR}
mkdir -p ${RUNDIR}/mulliken/
#cp *${system}${structure}.* ${RUNDIR}      
#cp *${system}${structure}.* ${RUNDIR}/mulliken

cd ${SRCDIR} 
cp charge0.txt  ${RUNDIR}
cp charge0.txt  ${RUNDIR}/mulliken

cd $RUNDIR
#rm *.e *.o *.dat* *.fchk *.out *.inp 
#rm -r tmp/

# write slurm submission script
cat > job.slurm <<EOF
#!/bin/bash

#SBATCH --time=10-00:00:00
#SBATCH --mem=100gb
#SBATCH --job-name="${jobname}"
#SBATCH --mail-type=END
#SBATCH --mail-user=${email}
#SBATCH --output=%j.o
#SBATCH --error=%j.e
#SBATCH --partition=${partition}
#SBATCH -c ${core}
##SBATCH --qos=argon
##SBATCH --constraint=g6132

TEMPDIR=/tmp/\${SLURM_JOB_USER}_\${SLURM_JOB_ID}
mkdir -p \$TEMPDIR
case=${structure}
module use ${module_loc}
module load qchem
module load openmpi/gnu/1.10.7 gcc/7.3
source /gpfsnyu/packages/qchem/4.4/qcenv.sh
export QC=$MYQC
export QCSCRATCH=\$TEMPDIR
export QCAUX=$MYQC/qcaux

qchem -nt ${core} ${system}${structure}.inp > ${system}${structure}.out


#remove TEMPDIR
rm -rf \$TEMPDIR

cd ${RUNDIR}
EOF


# ---------- write QChem input file ----------- #
cat > ${system}${structure}.inp <<EOF
\$molecule
read                  ${system}${structure}.txt
\$end

\$rem
JOBTYPE               sp        ! single point charge calculation
EXCHANGE              general   !
CORRELATION           none      !
BASIS                 ${basis_set}       !
SYM_IGNORE            TRUE      !
purecart              2111      !
separate_jk           TRUE      !
derscreen             FALSE     !
ideriv                0         !
OMEGA                 ${w}       ! geometric parameter tuned for geometry
MAX_SCF_CYCLES        5000      !
scf_convergence       8         ! for SP, default is 5 ; for opt/freq, default is 8
CIS_DYNAMIC_MEM       True      ! using dynamic memory for large TDDFT calculation
MEM_STATIC            16000     ! static memory for storing data.
MEM_TOTAL             180000    !
max_cis_cycles        500000    !
CIS_N_ROOTS           ${n_states}        ! number of excited states
CIS_SINGLETS          true      ! compute singlet excited states, default is true
CIS_TRIPLETS          false     ! compute triplet excited states, default is true
RPA                   false     ! false: TDA-DFT, true: TDA- and TDDFT , 2: TDDFT
sts_fcd               true      ! compute electronic couplings with FCD method
sts_donor             ${D_range}    ! donor atom index
sts_acceptor          ${A_range}      ! acceptor atom index
POP_MULLIKEN          -1        ! get Mulliken Charge for ground and excited states
GUI                   2         ! for producing wavefunction .fchk file
state_analysis        true      ! state analysis module
MOLDEN_FORMAT         true      ! generate molden formatted NTO
\$end

\$xc_functional
X HF 1.0
X BNL 1.0
C LYP 1.0
\$end

EOF

sbatch job.slurm
echo submitted job.slurm for TheoDORE analysis
#write another qchem input for mulliken

cp job.slurm $RUNDIR/mulliken
cd $RUNDIR/mulliken
cp ../${system}${structure}.txt .

cat > ${system}${structure}.inp <<EOF
\$molecule
read                  ${system}${structure}.txt
\$end

\$rem
JOBTYPE               sp        ! single point charge calculation
EXCHANGE              general   !
CORRELATION           none      !
BASIS                 ${basis_set}        !
SYM_IGNORE            TRUE      !
purecart              2111      !
separate_jk           TRUE      !
derscreen             FALSE     !
ideriv                0         !
OMEGA                 ${w}       ! geometric parameter tuned for geometry
MAX_SCF_CYCLES        5000      !
scf_convergence       8         ! for SP, default is 5 ; for opt/freq, default is 8
CIS_DYNAMIC_MEM       True      ! using dynamic memory for large TDDFT calculation
MEM_STATIC            16000     ! static memory for storing data.
MEM_TOTAL             180000    !
max_cis_cycles        500000    !
CIS_N_ROOTS           ${n_states}        ! number of excited states
CIS_SINGLETS          true      ! compute singlet excited states, default is true
CIS_TRIPLETS          false     ! compute triplet excited states, default is true
RPA                   false     ! false: TDA-DFT, true: TDA- and TDDFT , 2: TDDFT
sts_fcd               true      ! compute electronic couplings with FCD method
sts_donor             ${D_range}    ! donor atom index
sts_acceptor          ${A_range}      ! acceptor atom index
POP_MULLIKEN          -1        ! get Mulliken Charge for ground and excited states
\$end

\$xc_functional
X HF 1.0
X BNL 1.0
C LYP 1.0
\$end

EOF



cat > input <<EOF
NUM_TOT                ${natoms}             # Number of total atoms in system
NUM_of_ACCEPTOR        ${ndonor}            # Total number of acceptor fragments
NUM_of_DONOR           ${nacceptor}               # Total number of donor fragments
ATOMS_of_FRAGMENT      ${frag}        # Total number of atoms for each acceptor and donor fragment,(such as, Acc,Don,Don)
DFT_LEVEL              1               # BNL:1, B3LYP:2, wPBEh:3, SRSH:4
PHASE                  0               # GAS:0, PCM:1
NUMEXS                 ${n_states}              # Desired number of excited states
BASISSET               SV              # SV, 6-31G*... consult to the literature
DIEL                   7.6             # Dielectric value in PCM section; THF:7.6
OPTDIEL                1.78            # Optical dielectric value in PCM section
NOMEGA                 ${w}             # If known! 157 for linear and 154 for bent. For now use 157 for all
TUNOMEGA               0               # Enter 1 for tuning the omega for each geometry. Not implemented yet!!!
EOF


echo created input file and job.slurm for case ${system}${structure} with ${n_states} states and w = ${w}

sbatch job.slurm
#sbatch  --parsable job.slurm 
#@echo $jobID
echo job.slurm submitted for mulliken charge calculation
cd $WORKDIR

done
done
done
