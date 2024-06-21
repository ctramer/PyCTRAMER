#!/bin/sh
#SBATCH --job-name=PCT
#SBATCH --output=%j.o
#SBATCH --error=%j.e
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --mem=100GB
#SBATCH -p parallel #argon 
##SBATCH -w argon9
#SBATCH -c 40 #96
##SBATCH --qos=argon
##SBATCH  --exclude=compute[31,41,51,52]
##SBATCH --constraint=g6132
##SBATCH --constraint=g6146|9242
##SBATCH --dependency=afterok:5777088
echo The simulation will run in CPU Platform using PMEMD.MPI
# load amber20z module automatically according to cpu model
module purge
export MODULEPATH=$MODULEPATH:/gpfsnyu/xspace/sungroup/modules
# for e5 series cpu node (compute1-27), load avx2 version
if [ -n ${SLURM_JOB_NODELIST:7} ] && [ ${SLURM_JOB_NODELIST:7} -lt 28 ];then
module load amber/20z_avx2
else
# for gold/platinum seris cpu node, load avx512 version
module load amber/20z
fi
#@ module load amber
module load amber/22z
module list 

#HQDIR is where all the .in file will be placed
MOLECULE=
THE_DIR=
HQDIR=
SYSTEM=
ID_NO=
SOLVENT=
d_resname=
a_resname=
CONFIG=${ID_NO}_${SOLVENT}
CONFIG_GAS=${ID_NO}_no${SOLVENT}
STATE=PI
DONOR=PI
TRAJ_STATE=NVE
JOBDIR=${THE_DIR} #/${SYSTEM}${ID_NO}
TRAJDIR=${JOBDIR}/Traj_${STATE}
NTRAJ=
CONTINUE_FROM=  # if CONTINUE_FROM = 0 ; starts from beginning i.e. NPT equilibrated and start sampling.
min_inputfile=min.in
heat_inputfile=heat.in
npt_eq_inputfile=equil_NPT.in #npt_eq_fixed_solute.in
nvt_eq_inputfile=equil_NVT.in #nvt_eq_fixed_solute.in
nvt_sampling_inputfile=nvt_sampling_fixed_solute.in
nve_relax_inputfile=nve_relax_fixed_solute.in
nve_prod_inputfile=nve_prod_fixed_solute.in
state_list=
mkdir -p $TRAJDIR

for inp in $min_inputfile $heat_inputfile $npt_eq_inputfile $nvt_eq_inputfile $nvt_sampling_inputfile $nve_relax_inputfile $nve_prod_inputfile
do
cp $HQDIR/$inp $JOBDIR/.
done

TOP=${SYSTEM}${CONFIG}_${STATE}.prmtop  
INIT_CONFIG=${SYSTEM}${CONFIG}.inpcrd
NT=96
md_control=pmemd
run_on_tmp=yes
min=no
heat=no
eq_npt=no
eq_nvt=no
sample_nvt=no
relax_nve=no
prod_nve=no
recalculate=no
energy_analysis=no
get_correction=no

#Control for loop job
START=$(($CONTINUE_FROM + 1))
STOP=$(($CONTINUE_FROM +  $NTRAJ))

# before we start with everything else, create a directory 
mkdir -p $TRAJDIR
cp $THE_DIR/*inpcrd  $JOBDIR/.
cp $THE_DIR/$TOP    $JOBDIR/.
cp $THE_DIR/*prmtop $JOBDIR/.
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
# cd $TRAJDIR
#############################################################################
if [ "$min" == "yes" ];then
echo Starting Energy Minimization... `date` 
echo `pwd` 
###[1] minimization run
       mpirun -np ${NT} ${md_control}.MPI -O -i ${HQDIR}/${min_inputfile} \
           -o min.out \
           -p ${JOBDIR}/${TOP} \
           -c ${JOBDIR}/${INIT_CONFIG} -ref ${JOBDIR}/${INIT_CONFIG} \
           -r min.rst \
           -inf min.mdinfo 
echo Job minimization finished on `date` 
else
	echo WARNING: No minimization were performed. 
fi
###########################################################################

if [ "$heat" == "yes" ];then
echo Start Heating ... `date` 
## [2] heat
       mpirun -np ${NT} ${md_control}.MPI -O -i ${HQDIR}/${heat_inputfile} \
           -o heat.out \
           -p ${JOBDIR}/${TOP} \
           -c min.rst -ref min.rst \
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
       mpirun -np ${NT} ${md_control}.MPI -O -i ${HQDIR}/${npt_eq_inputfile} \
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


echo Starting NVT Equilibration ... `date` 

############################################################################

if [ "$eq_nvt" == "yes" ];then
echo Starting NVT Equilibration ...
       mpirun -np ${NT} ${md_control}.MPI -O -i ${HQDIR}/${nvt_eq_inputfile} \
           -o ${SYSTEM}${ID_NO}_NVT_0.out \
           -p ${JOBDIR}/${TOP} \
           -c equil_NPT.rst -ref equil_NPT.rst \
           -r ${SYSTEM}${ID_NO}_NVT_0.rst \
           -x ${SYSTEM}${ID_NO}_NVT_0.mdcrd.nc \
           -inf ${SYSTEM}${ID_NO}_NVT_0.mdinfo
## NOTE: _0 is indexing to represent this is initial configuration for further MD Run.
echo Finished NVT equilibration run on `date` 
# to prevent memory overload move everything to job dir every timestep, leave only rst behind
mv *.out* *.mdcrd* *.mdinfo  $TRAJDIR
cp *.rst $TRAJDIR 
	echo Warning: No NVT equilibration were performed.
fi 

###########################################################################


if [ "$sample_nvt" == "yes" ];then

echo Starting NVT sampling from equilibrated NVT ... `date` 

TRAJDIR=${JOBDIR}/Traj_${STATE}
mkdir -p $TRAJDIR
cp  ${SYSTEM}${ID_NO}_NVT_0.rst ${TRAJDIR}
# backup equilibrated .rst just in case



for  (( j = ${START}; j <= ${STOP}; j++))
do
beforej=$((j-1))
# [5]  Run NVT sampling from equilibrated NVT
    mpirun -np ${NT} ${md_control}.MPI -O -i ${HQDIR}/${nvt_sampling_inputfile} \
           -o ${SYSTEM}${ID_NO}_NVT_${j}.out \
           -p ${JOBDIR}/${TOP} \
           -c ${SYSTEM}${ID_NO}_NVT_${beforej}.rst -ref ${SYSTEM}${ID_NO}_NVT_${beforej}.rst \
           -r ${SYSTEM}${ID_NO}_NVT_${j}.rst \
           -x ${SYSTEM}${ID_NO}_NVT_${j}.mdcrd.nc \
           -inf ${SYSTEM}${ID_NO}_NVT_${j}.mdinfo
done
# to prevent memory overload move everything to job dir every timestep, leave only rst behind
mv *.out* *.mdcrd* *.mdinfo  $TRAJDIR
echo Finished performing NVT sampling on `date` 
fi


##############################################################################################33
if [ "$relax_nve" == "yes" ];then
echo Starting NVE relaxation ... `date` 
for  (( j = ${START}; j <= ${STOP}; j++)) 
do
# [6]  Run NVE relaxation run
    mpirun -np ${NT} ${md_control}.MPI  -O -i ${HQDIR}/${nve_relax_inputfile} \
           -o ${SYSTEM}${ID_NO}_relax_NVE_${j}.out \
	   -p ${JOBDIR}/${TOP} \
           -c ${SYSTEM}${ID_NO}_NVT_${j}.rst -ref ${SYSTEM}${ID_NO}_NVT_${j}.rst \
           -r ${SYSTEM}${ID_NO}_relax_NVE_${j}.rst \
           -x ${SYSTEM}${ID_NO}_relax_NVE_${j}.mdcrd.nc \
           -inf ${SYSTEM}${ID_NO}_relax_NVE_${j}.mdinfo
# to prevent memory overload move everything to job dir every timestep, leave only rst behind
mv *.out* *.mdcrd* *.mdinfo  $TRAJDIR
done
echo Finished all NVE relaxation run on `date` 
fi 
######################################################################################################

if [ "$prod_nve" == "yes" ];then
echo Starting NVE production ... `date` 
for  (( j = ${START}; j <= ${STOP}; j++))
do
# [7] Run NVE production run
     mpirun -np ${NT} ${md_control}.MPI  -O -i ${HQDIR}/${nve_prod_inputfile} \
           -o ${SYSTEM}${ID_NO}_prod_NVE_${j}.out \
           -p ${JOBDIR}/${TOP} \
           -c ${SYSTEM}${ID_NO}_relax_NVE_${j}.rst -ref ${SYSTEM}${ID_NO}_relax_NVE_${j}.rst \
           -r ${SYSTEM}${ID_NO}_prod_NVE_${j}.rst \
           -x ${SYSTEM}${ID_NO}_prod_NVE_${j}.mdcrd.nc \
           -inf ${SYSTEM}${ID_NO}_prod_NVE_${j}.mdinfo
# to prevent memory overload move everything to job dir every timestep, leave only rst behind
mv *.out* *.mdcrd* *.mdinfo  $TRAJDIR
done
echo Finished all NVE production run on `date` 
fi
###############################################################################################3

# [8] Energy Recalculation 
#esander Recalculation running in TEMPDIR
if [ "$recalculate" == "yes" ];then
for (( j = ${START}; j <= ${STOP}; j++))
do
echo Starting Energy Recalculation ... `date` 
echo Transporting .mdcrd.nc to a save tempdir ... 
cp $TRAJDIR/${SYSTEM}${ID_NO}_prod_NVE_${j}.mdcrd.nc $TEMPDIR
cd $TEMPDIR
echo `pwd; ls`


for CHARGE_STATE in ${state_list} #GR PI CT #1 CT2 
do
cat > esander_${CHARGE_STATE}.in <<EOF
parm ${JOBDIR}/${SYSTEM}${CONFIG}_${CHARGE_STATE}.prmtop  
trajin ${SYSTEM}${ID_NO}_prod_NVE_${j}.mdcrd.nc
esander ${CHARGE_STATE} out ${SYSTEM}${ID_NO}_REC_${CHARGE_STATE}_${DONOR}_${TRAJ_STATE}_${j}.dat ntb 1 cut 9 ntf 2 ntc 2
EOF
mpirun -np ${NT} cpptraj.MPI -i esander_${CHARGE_STATE}.in > esander_${CHARGE_STATE}.out
echo Finished Energy Recalculation on state ${CHARGE_STATE} on `date` 
# We retain only at max 1 big .nc trajectories at a given time 
done
mv *.in *.out *.dat $TRAJDIR

# rm *.nc 
done
echo Finished all Energy Recalculation on `date` 
fi 


###########################################################################
if [ "$get_correction" == "yes" ];then
# [9] Get Energy Correction for gas phase MD 
cd $TRAJDIR

echo prepare equilibrated NVT .rst without THF
cat > stripTHFrst.cpptraj <<EOF
parm ${JOBDIR}/${TOP}
trajin ${TRAJDIR}/${SYSTEM}${ID_NO}_NVT_0.rst
autoimage
strip :${d_resname}
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
echo Finished all gas phase energy calculation on `date`

fi


############################################################################

# [10] Energy Analysis 
# extract only total energy from all the esander calculated energy

if [ "$energy_analysis" == "yes" ];then
echo Starting Energy analysis ... `date`

echo Start with energy analysis for gas phase ...

for CHARGE_STATE in  $state_list  #GR PI CT # 1 CT2
do

# this is always on the NVT traj_state
awk '{print $2}' ${TRAJDIR}/${SYSTEM}${ID_NO}_GAS_${CHARGE_STATE}_${DONOR}_NVT.dat | sed 1d > ${TRAJDIR}/energy_gas_${CHARGE_STATE}_TRAJ_NVT.dat


for (( j = ${START}; j <= ${STOP}; j++))
do
echo Continue with energy analysis of the solvated system

awk '{print $2}' ${TRAJDIR}/${SYSTEM}${ID_NO}_REC_${CHARGE_STATE}_${DONOR}_${TRAJ_STATE}_${j}.dat | sed 1d > ${TRAJDIR}/energy_${CHARGE_STATE}_TRAJ_${TRAJ_STATE}_${j}.dat

done

done
# stitch everything into 1
for REC_STATE in ${state_list} # GR PI CT # 1 CT2
do
rm ${MOLECULE}_E_${REC_STATE}_TRAJ_${STATE}_${TRAJ_STATE}.dat
done

for  (( j = ${START}; j <= ${STOP}; j++))
do
        for REC_STATE in ${state_list} # GR PI CT #1 CT2
        do
        cat energy_${REC_STATE}_TRAJ_${TRAJ_STATE}_${j}.dat >> ${MOLECULE}_E_${REC_STATE}_TRAJ_${STATE}_${TRAJ_STATE}.dat
        done
done


echo Finished all energy analysis on `date`

fi


if [ "$run_on_tmp" == "yes" ];then 
	echo cleaning up ... 
	echo NOTE: double check if necessary, make sure no one got left behind.
	# This time around we take everything with us no one got left behind.
	mv *.out* *.mdcrd* *.rst *.mdinfo  $TRAJDIR
	rm -rf $TEMPDIR
else 
	echo Check your result and have fun !
fi 
echo Check your result and have fun ! 
echo All job is finished on `date` 
