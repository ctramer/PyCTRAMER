#!/bin/sh
#SBATCH --job-name=PCT-GPU
#SBATCH --output=%j.o
#SBATCH --error=%j.e
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --mem=10GB
#SBATCH -p argon 
#SBATCH -c 1
#SBATCH --gres=gpu:1
#SBATCH --dependency=afterok:5776454
#SBATCH --constraint=g6230|2630v4

#load amber20z module automatically according to gpu model
module purge
export MODULEPATH=$MODULEPATH:/xspace/sungroup/modules
if [ "$SLURM_JOB_NODELIST" == "gpu7" ];then
# for 2080Ti, load avx512_cuda11.02 version
module load amber/20z
elif [ "$SLURM_JOB_NODELIST" == "gpu6" ];then
# for P100, load avx2_cuda11.02 version
module load amber/20z_avx2
elif [ "$SLURM_JOB_NODELIST" == "agpu8" ];then
# for 2080Ti, load avx512_cuda11.02 version
module load amber/20z
else
# for K10, load avx2_cuda10.2 version
module load amber/20z_avx2_cuda10.2
fi
module load amber/22gpu

echo The simulation will run in GPU Platform using PMEMD.cuda
#PLATFORM=GPU


#HQDIR is where all the .in file will be placed
# This section is for input file will be included from external control 
THE_DIR=/gpfsnyu/xspace/projects/PCT-example/DBP-C70/results/DBP-C70_conf11/MD/
HQDIR=/gpfsnyu/xspace/projects/PCT-example/DBP-C70/results/DBP-C70_conf11/MD/
SYSTEM=DBP-C70_
ID_NO=conf11
SOLVENT=DBPC70  #DDDAAA
CONFIG=${ID_NO}_${SOLVENT}
CONFIG_GAS=${ID_NO}_no${SOLVENT}
STATE=LE
DONOR=LE
d_resname=DBP
a_resname=C70
JOBDIR=${THE_DIR} # /${SYSTEM}${ID_NO}
TRAJDIR=${JOBDIR}/Traj_${STATE}
CASE=0
NTRAJ=1
CONTINUE_FROM=0
min_inputfile=min.in
heat_inputfile=heat.in
npt_eq_inputfile=equil_NPT.in #npt_eq_fixed_solute.in
nvt_eq_inputfile=equil_NVT.in #nvt_eq_fixed_solute.in
nvt_sampling_inputfile=nvt_sampling_fixed_solute.in
nve_relax_inputfile=nve_relax_fixed_solute.in
nve_prod_inputfile=nve_prod_fixed_solute.in
state_list='GR LE CT1 CT2 CT3 CT4 CT5 CT6 '
for inp in $min_inputfile $heat_inputfile $npt_eq_inputfile $nvt_eq_inputfile $nvt_sampling_inputfile $nve_relax_inputfile $nve_prod_inputfile
do
cp $HQDIR/$inp $JOBDIR/.
done
mkdir -p $TRAJDIR
TOP=${SYSTEM}${CONFIG}_${STATE}.prmtop 
INIT_CONFIG=${SYSTEM}${CONFIG}.inpcrd
NT=48 #48 96
md_control=pmemd # only pmemd is supported for GPU run
run_on_tmp=yes
min=no # GPU Run may not produce appropriate minimization due to single precision limitation, use CPU instead
heat=no #yes
eq_npt=no #yes
eq_nvt=yes
sample_nvt=yes
relax_nve=yes
prod_nve=yes
get_correction=yes
cd $TRAJDIR
# 6 is for starting from previous run
#Control for loop job
START=$((CASE+ $CONTINUE_FROM + 1))
STOP=$((CASE+ $CONTINUE_FROM + $NTRAJ))

# before we start with everything else, create a directory 
mkdir -p $TRAJDIR
cp $THE_DIR/*inpcrd  $TRAJDIR/.
cp $THE_DIR/$TOP    $JOBDIR/.
if [ "$run_on_tmp" == "yes" ];then 
	echo Calculation were performed on tmp directory on computation node
	echo NOTE: if job failure happen make sure you clean up any cookie crumble of your calculation may left on the computation node. 
	TEMPDIR=/tmp/${SLURM_JOB_USER}_${SLURM_JOB_ID}
	mkdir $TEMPDIR
	cd $TEMPDIR
else 
	echo Calculation were performed on current node
	echo NOTE: Make sure that this calculation is not heavy one that read and write frequently otherwise this job might influence cluster performance.
fi

#############################################################################
if [ "$min" == "yes" ];then 
echo Starting Energy Minimization... `date` 
###[1] minimization run
	pmemd.cuda -O -i ${HQDIR}/${min_inputfile} \
           -o min.out \
           -p ${JOBDIR}/${TOP} \
           -c ${JOBDIR}/${INIT_CONFIG} -ref ${JOBDIR}/${INIT_CONFIG} \
           -r min.rst \
           -inf min.mdinfo 
echo Job minimization finished on `date` 
else
	echo WARNING: No minimization were performed. 
	echo min.rst is to be fetch from ${TRAJDIR}
fi
###########################################################################

if [ "$heat" == "yes" ];then
echo Start Heating ... `date` 
## [2] heat
	pmemd.cuda -O -i ${HQDIR}/${heat_inputfile} \
           -o heat.out \
           -p ${JOBDIR}/${TOP} \
           -c ${TRAJDIR}/min.rst -ref ${TRAJDIR}/min.rst \
           -r heat.rst \
           -x heat.mdcrd.nc \
           -inf heat.mdinfo

echo Job heat finished on `date`
else 
	echo WARNING: Proper heating to target temperature were not performed.
fi
############################################################################

if [ "$eq_npt" == "yes" ];then
echo Starting NPT Equilibration ... `date` 
### [3] 1 ns NPT equilibriation
	pmemd.cuda -O -i ${HQDIR}/${npt_eq_inputfile} \
           -o equil_NPT.out \
           -p ${JOBDIR}/${TOP} \
           -c heat.rst -ref heat.rst \
           -r equil_NPT.rst \
           -x equil_NPT.mdcrd.nc \
           -inf equil_NPT.mdinfo
   
echo Job NPT_eq finished on `date`
else 
	echo WARNING: No NPT equilibration were performed.
fi



############################################################################

if [ "$eq_nvt" == "yes" ];then
echo Starting NVT Equilibration ... `date`
cp $TRAJDIR/equil_NPT.rst  $TEMPDIR
echo `pwd; ls`
## [4] 1 ns NVT equilibriation
	pmemd.cuda -O -i ${HQDIR}/${nvt_eq_inputfile} \
           -o ${SYSTEM}${ID_NO}_NVT_${CONTINUE_FROM}.out \
           -p ${JOBDIR}/${TOP} \
           -c equil_NPT.rst -ref equil_NPT.rst \
           -r ${SYSTEM}${ID_NO}_NVT_${CONTINUE_FROM}.rst \
           -x ${SYSTEM}${ID_NO}_NVT_${CONTINUE_FROM}.mdcrd.nc \
           -inf ${SYSTEM}${ID_NO}_NVT_${CONTINUE_FROM}.mdinfo
## NOTE: _0 is indexing to represent this is initial configuration for further MD Run.
echo Finished NVT equilibration run
# to prevent memory overload move everything to job dir every timestep, leave only rst behind
cp *.out* *.mdcrd* *.mdinfo  $TRAJDIR
# perform backup 
cp *.rst $TRAJDIR 
else 
	echo Warning: No NVT equilibration were performed.
fi 

###########################################################################
if [ "$get_correction" == "yes" ];then
echo prepare equilibrated NVT .rst without THF
cat > stripTHFrst.cpptraj <<EOF
parm ${JOBDIR}/${TOP}
trajin ${TRAJDIR}/${SYSTEM}${ID_NO}_NVT_0.rst
autoimage
strip :${d_resname} #${SOLVENT}
strip :${a_resname}
trajout ${SYSTEM}${ID_NO}_no${SOLVENT}_NVT_0.rst
EOF
cpptraj -i stripTHFrst.cpptraj
echo created a NVT equilibrated .rst file for the solute only on `date`
echo performing a single frame MD energy calculation for all states... 

# perform calculation for gas phase MD to obtained solute only energy needed for energy correction
for CHARGE_STATE in ${state_list} #GR PI CT #1 CT2
do

cat > esander_gas_${CHARGE_STATE}.in <<EOF
parm ${JOBDIR}/${SYSTEM}${CONFIG_GAS}_${CHARGE_STATE}.prmtop  
trajin ${SYSTEM}${ID_NO}_no${SOLVENT}_NVT_0.rst
esander ${CHARGE_STATE} out ${SYSTEM}${ID_NO}_GAS_${CHARGE_STATE}_${DONOR}_NVT.dat ntb 1 cut 9 ntf 2 ntc 2
EOF
mpirun -np 1 cpptraj.MPI -i esander_gas_${CHARGE_STATE}.in > esander_gas_${CHARGE_STATE}.out
echo Finished gas phase energy calculation on state ${CHARGE_STATE} on `date`

###########################################################################
# We retain only at max 1 big .nc trajectories at a given time 
done
mv *.in *.out *.dat $TRAJDIR
echo Finished all gas phase nergy calculation on `date`

fi


if [ "$sample_nvt" == "yes" ];then

echo Starting NVT sampling from equilibrated NVT ... `date` 
cp  ${SYSTEM}${ID_NO}_NVT_${CONTINUE_FROM}.rst ${TRAJDIR}
# backup equilibrated .rst just in case



for  (( j = ${START}; j <= ${STOP}; j++))
do
beforej=$((j-1))
# [5]  Run NVT sampling from equilibrated NVT
     pmemd.cuda -O -i ${HQDIR}/${nvt_sampling_inputfile} \
           -o ${SYSTEM}${ID_NO}_NVT_${j}.out \
           -p ${JOBDIR}/${TOP} \
           -c ${SYSTEM}${ID_NO}_NVT_${beforej}.rst -ref ${SYSTEM}${ID_NO}_NVT_${beforej}.rst \
           -r ${SYSTEM}${ID_NO}_NVT_${j}.rst \
           -x ${SYSTEM}${ID_NO}_NVT_${j}.mdcrd.nc \
           -inf ${SYSTEM}${ID_NO}_NVT_${j}.mdinfo
done
# to prevent memory overload move everything to job dir every timestep, leave only rst behind
cp *.out* *.mdcrd* *.mdinfo  $TRAJDIR
echo Finished all NVT sampling on  `date` 
fi


##############################################################################################33
if [ "$relax_nve" == "yes" ];then
echo Starting NVE relaxation ... `date` 
for  (( j = ${START}; j <= ${STOP}; j++)) 
do
# [6]  Run NVE relaxation run
     pmemd.cuda -O -i ${HQDIR}/${nve_relax_inputfile} \
           -o ${SYSTEM}${ID_NO}_relax_NVE_${j}.out \
	   -p ${JOBDIR}/${TOP} \
           -c ${SYSTEM}${ID_NO}_NVT_${j}.rst -ref ${SYSTEM}${ID_NO}_NVT_${j}.rst \
           -r ${SYSTEM}${ID_NO}_relax_NVE_${j}.rst \
           -x ${SYSTEM}${ID_NO}_relax_NVE_${j}.mdcrd.nc \
           -inf ${SYSTEM}${ID_NO}_relax_NVE_${j}.mdinfo
# to prevent memory overload move everything to job dir every timestep, leave only rst behind
cp *.out* *.mdcrd* *.mdinfo  $TRAJDIR
done
echo Finished all NVE relaxation run on `date` 
fi 
######################################################################################################

if [ "$prod_nve" == "yes" ];then
echo Starting NVE production ... `date` 
for  (( j = ${START}; j <= ${STOP}; j++))
do
# [7] Run NVE production run
     pmemd.cuda -O -i ${HQDIR}/${nve_prod_inputfile} \
           -o ${SYSTEM}${ID_NO}_prod_NVE_${j}.out \
           -p ${JOBDIR}/${TOP} \
           -c ${SYSTEM}${ID_NO}_relax_NVE_${j}.rst -ref ${SYSTEM}${ID_NO}_relax_NVE_${j}.rst \
           -r ${SYSTEM}${ID_NO}_prod_NVE_${j}.rst \
           -x ${SYSTEM}${ID_NO}_prod_NVE_${j}.mdcrd.nc \
           -inf ${SYSTEM}${ID_NO}_prod_NVE_${j}.mdinfo
# to prevent memory overload move everything to job dir every timestep, leave only rst behind
cp *.out* *.mdcrd* *.mdinfo  $TRAJDIR
done
echo Finished all NVE production run on `date` 
fi
###############################################################################################3
if [ "$run_on_tmp" == "yes" ];then 
	echo cleaning up ... 
	echo NOTE: double check if necessary, make sure no one got left behind.
	# This time around we take everything with us no one got left behind.
	cp *.out* *.mdcrd* *.rst *.mdinfo  $TRAJDIR
	rm -rf $TEMPDIR
else 
	echo Check your result and have fun ! `date` 
fi 

echo Check your result and have fun ! 
echo All job is finished on `date` 
















