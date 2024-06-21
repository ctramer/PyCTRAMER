#!/bin/bash 

Template_Dir=/gpfsnyu/xspace/projects/PCT-example/DBP-C70/template
Target_Dir=/gpfsnyu/xspace/projects/PCT-example/DBP-C70/results/DBP-C70_conf11/MD/

# Warning! Warning! Warning!
# one need to use different template dirs with different simulation parameters in order to avoid misuse of MD input files 

MD_steps=100000      
MD_dt=0.001          
sample_steps=5 
MD_temperature=300.0

# 

cd ${Template_Dir}


for simu in  heat  equil_NPT  equil_NVT  nvt_eq_fixed_solute  nvt_sampling_fixed_solute  nve_prod_fixed_solute  nve_relax_fixed_solute
do

ori_file=${Template_Dir}/${simu}.temp 
tar_file=${simu}.in 

sed     "s/MYDT/${MD_dt}/"  $ori_file >    $tar_file
sed -i  "s/TARTEMP/${MD_temperature}/g"     $tar_file
sed -i  "s/MYSAMPLESTEP/${sample_steps}/g"  $tar_file
sed -i  "s/MYMDSTEP/${MD_steps}/"          $tar_file

cp $tar_file $Target_Dir

done  

cp ${Template_Dir}/min.in  $Target_Dir
