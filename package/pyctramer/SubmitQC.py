import numpy as np 
import pandas as pd
import sys
import os
from .Utilities import * 



def run_step1(i, job_control, dict_of_simulation):
    #print(arrow_str,"Started Step 1: SubmitQC") 
    project       = dict_of_simulation['project']
    work_dir      = dict_of_simulation['work_dir']
    structure_dir = dict_of_simulation['structure_dir']
    print(arrow_str, 'caseid: ' ,job_control[i]['caseid'] )
    caseid =  job_control[i]['caseid'] #str(i)
    case_dir = work_dir + '/' + project +"_"+ caseid
    # Quantum chemistry set up
    qc_dir,md_dir = prepare_folder_structure(case_dir)
    ## charge of the molecule of interest
    charge = dict_of_simulation['charge']
    ## mulitiplicity
    multiplicity = dict_of_simulation['multiplicity']
    # write the input structure location
    input_structure = structure_dir + '/' + project + "_" + caseid + ".xyz"
    ## write the structure file into the target QC folder
    write_txt_for_qchem(structure_dir,input_structure, charge, multiplicity)
    input_structure_qchem = structure_dir + '/' + project + "_" + caseid + ".txt"
    # prepare pdb in target folder
    input_structure_pdb = structure_dir + '/' + project + "_" + caseid + ".pdb"
    os.popen("mkdir -p " +qc_dir+"/mulliken").read()
    prepare_QC(input_structure_qchem, qc_dir)
    prepare_QC(input_structure_qchem, qc_dir+"/mulliken/.")
    prepare_QC(input_structure, qc_dir)
    prepare_QC(input_structure_pdb, qc_dir) 

    # print(dict_of_simulation) 
    for dosname in dict_of_simulation:
         print(dosname,': ', dict_of_simulation[dosname])
    template_dir =  dict_of_simulation['template_dir']

    # loading quantum chemistry parameters
    N_atom = dict_of_simulation['N_atom']
    N_cis  = dict_of_simulation['N_cis']

    QC_method   = dict_of_simulation['QC_method']
    basis_set   = dict_of_simulation['basis_set']
    RSH_omega = dict_of_simulation['RSH_omega']
    TD_approx = dict_of_simulation['TD_approx'] 
    charge_type = dict_of_simulation['charge_type']

    DA_type  = dict_of_simulation['DA_type']
    D_range  =  dict_of_simulation['D_range']
    B_range  =  dict_of_simulation['B_range']
    A_range  =  dict_of_simulation['A_range']
    A1_range = dict_of_simulation['A1_range']
    A2_range = dict_of_simulation['A2_range']
    CT_donor_fragment = dict_of_simulation['CT_donor_fragments']
    CT_acceptor_fragment = dict_of_simulation['CT_acceptor_fragments'] 

    

    jobs = []
    returning_info = submitQC(template_dir, qc_dir, QC_method, basis_set, caseid, D_range, B_range, A_range, A1_range, A2_range,CT_donor_fragment,CT_acceptor_fragment,N_cis, project, qc_dir,TD_approx,RSH_omega,DA_type,charge_type,N_atom)

    find_job_id(returning_info, jobs)
    dependency_at_this_step = list_to_dep_line(jobs)
    job_control[i]['QC_jobids'] = dependency_at_this_step

    
    stroutput = "\n" + arrow_str + "STEP 1: SubmitQC\n" 
    stroutput = stroutput + "Input: \nsolute structure:  "+ input_structure + "\n"
    stroutput = stroutput + "Job submit script:\nQuantum chemistry: " + qc_dir + "\n" 
    stroutput = stroutput + "SLURM script:\n" + qc_dir + "submitQC_flexible*.slurm \n"
    stroutput = stroutput + "QChem input file:\n " + qc_dir +  project +"_"+ caseid + ".inp \n"
    stroutput = stroutput + qc_dir + 'mulliken/' +  project +"_"+ caseid + ".inp \n"
    stroutput = stroutput + "Output: \n    QChem output file: \n" + qc_dir +  project +"_"+ caseid + ".out \n"
    stroutput = stroutput + qc_dir +  project +"_"+ caseid + ".fchk \n"
    stroutput = stroutput + qc_dir + "mulliken/" +  project +"_"+ caseid + ".out \n"
    stroutput = stroutput + qc_dir + "mulliken/" +  project +"_"+ caseid + ".fchk \n"
    stroutput = stroutput + "#######################################################\n"
    write_to_output(i, job_control, dict_of_simulation, stroutput) 
    return job_control

   
def submitQC(template_dir,running_dir,QCmethod,basis_set,testing_str,D_range,B_range,A_range,A1_range,A2_range,CT_donor_fragment,CT_acceptor_fragment,N_cis, project,basedir,TD_approx,RSH_omega,DA_type,charge_type,N_atom):
    """
    template_dir: slrum script template for the project 
    running_dir:  where to run this quantum chemistry jobs
    QCmethod: quantum chemistry calculation method like HF or DFT method like b3lyp, wb97x-b, should be a string 
    basis_set: basis set for quantum chemistry electronic structure calculation 
    testing_str: a string containing several space separated values utilized as part of the name of quantum chemistry jobs
    Drange: string like '1-60'
    Arange: string like '61-207'
    N_cis: number of singlets are calculated, int 
    case: series name for the calculation 
    
    This function prepare and submit QC job based on a submitQC_template.sh 
    
    check the template for platform and cluster specific parameters
    in the next release the cluster control can be performed from external function ctr.job_submit_setup()
    """

    os.chdir(running_dir)

    if QCmethod == "BNL":
        submitQC_file = shell_copy_template(template_dir,running_dir, filename='submitQC_general')
    else :
        submitQC_file = shell_copy_template(template_dir,running_dir, filename='submitQC_preset')
    
    if TD_approx == "TDA": 
        TD_token = "false"
    elif TD_approx == "RPA":
        TD_token = "true"
    else:
        TD_token = "2"
        print("Full TDDFT would be running ")
    if charge_type == "Mulliken":
        charge_token = "Mulliken"
    elif charge_type == "RESP":
        raise ValueError("charge_type RESP is not implemented yet")
    else:
        raise ValueError("charge_type should be Mulliken or RESP")
    
    if RSH_omega == "":
        sed_inplace(submitQC_file,"omega=","omega=114514")
        sed_inplace(submitQC_file,"OMEGA", "! OMEGA")
    else:
        sed_inplace(submitQC_file,"omega=","omega="+RSH_omega)

    ordersD  = np.fromstring(D_range, dtype=int, sep='-')
    if B_range != "":
        ordersB  = np.fromstring(B_range, dtype=int, sep='-') 
    ordersA  = np.fromstring(A_range, dtype=int, sep='-')
    if A1_range != "":
        ordersA1  = np.fromstring(A1_range, dtype=int, sep='-') 
    if A2_range != "":
        ordersA2  = np.fromstring(A2_range, dtype=int, sep='-')
    
    dsegmentnum = len(CT_donor_fragment)
    asegmentnum = len(CT_acceptor_fragment) 
    for i in [dsegmentnum,asegmentnum]:
        print(i)
        if i > 1:
            i = i - 1
        else:
            pass # doing nothing
    
    donor_fragment = np.ones(dsegmentnum)
    acceptor_fragment = np.ones(asegmentnum)

    if len(CT_donor_fragment) != 1:
        # print(CT_donor_fragment, type(CT_donor_fragment))
        donor_fragment = np.fromstring(CT_donor_fragment, sep=',',dtype=int)
    if len(CT_acceptor_fragment) != 1:
        acceptor_fragment = np.fromstring(CT_acceptor_fragment, sep=',',dtype=int)

    if DA_type == "DBA":
        orders = np.array([ordersD[0],ordersD[1],ordersB[0],ordersB[1],ordersA[0],ordersA[1]])
        if len(donor_fragment) == 2:
            D_range = str( min([min(ordersD),min(ordersB)]) ) + "-" + str( max([max(ordersB),max(ordersD)] ))

        if len(acceptor_fragment) == 2:
            A_range = str( min([min(ordersA),min(ordersB)]) ) + "-" + str( max([max(ordersB),max(ordersA)] ))


    os.chdir(running_dir)
    QC_DIR = f'{running_dir}'
    #QC_DIR=f'{project_dir}/'+case+'/QC/'
    
    sed_inplace(submitQC_file,'WORKDIR=','WORKDIR='+QC_DIR)
    sed_inplace(submitQC_file,'SRCDIR=','SRCDIR='+template_dir)
    sed_inplace(submitQC_file,'RUNDIR=','RUNDIR='+QC_DIR)
    sed_inplace(submitQC_file,'system=','system='+project +'_')
    sed_inplace(submitQC_file,'basis_set=','basis_set='+basis_set)
    sed_inplace(submitQC_file,'method=','method='+QCmethod)


    sed_inplace(submitQC_file,'D_range=','D_range='+D_range)
    sed_inplace(submitQC_file,'B_range=','B_range='+B_range)
    sed_inplace(submitQC_file,'A_range=','A_range='+A_range)
    sed_inplace(submitQC_file,'STATE_NUM=','STATE_NUM='+str(N_cis))
    sed_inplace(submitQC_file,'GIVEN_STRUC=','GIVEN_STRUC='+ testing_str)
    sed_inplace(submitQC_file, "TD_approx=","TD_approx="+TD_token)
    # if RSH_omega == "":
    #     sed_inplace(submitQC_file,"omega=","omega=114514")
    #     sed_inplace(submitQC_file,"OMEGA", "!OMEGA")
    # else:
    #     sed_inplace(submitQC_file,"omega=","omega="+RSH_omega)

    os.chdir(running_dir)
    if QCmethod == "BNL":
        jobid = os.popen('sh submitQC_general.sh').read()
    else: 
        jobid = os.popen('sh submitQC_preset.sh').read()
    
    return jobid

def read_xyz(filename, debug=False):
    # read xyz file fr
    with open(filename) as fi:
        lines = fi.readlines()

    if debug == True:
        print('DEBUG MODE ON: def read_xyz')
        print(lines)
    return lines

def write_txt_for_qchem(structure_dir,filename,charge,multiplicity,debug=False):
    # read files of a given name and write the txt files for Q-Chem calculation in txt format
    # return filename of the xyz and txt in string
    lines_in_txt_files = read_xyz(filename,debug=debug)
    Output_File_name = filename[0:-3]+"txt"
    # edit Output File
    lines_in_txt_files[0] = '$molecule \n'
    lines_in_txt_files[1] = str(charge) + ' ' + str(multiplicity) + ' \n'
    lines_in_txt_files[-1]=('$end')
    os.chdir(structure_dir)
    os.popen("cp "+filename + "  " +Output_File_name ).read()
    os.popen("sed -i '1c " +  lines_in_txt_files[0] + "' " + Output_File_name).read()
    os.popen("sed -i '2c " +  lines_in_txt_files[1] + "' " + Output_File_name).read()
    os.popen("sed -i '$a " +  lines_in_txt_files[-1] + "' " + Output_File_name).read()
    
    return filename.split('/')[-1][0:-4]

def prepare_QC(input_structure, QC_path):
    os.popen('cp '+input_structure+'  '+QC_path+' ').read()
    return 


def main_SubmitQC():
    dict_of_simulation =  parse_input()
    # prepare for quantum chemistry simulation 
    work_dir = dict_of_simulation['work_dir']
    # where we have the structure directory
    structure_dir = dict_of_simulation['structure_dir']
    # load the project folder 
    project = dict_of_simulation['project']
    # get case id for the simulation
    caseid = dict_of_simulation['caseid']
    # get job list
    case_id_list = separate_idlist(caseid)
    # get list of job control token 
    job_control  = init_job_control(case_id_list)

    startime = time.time 
    print("SubmitQC begins at "+str(startime))
    # print("Input parameters ")
    # print(simulation_parameter)
    for ind, ctr_ele in enumerate(job_control): 
        print(ind, job_control)
        run_step1(ind, job_control, dict_of_simulation)

    # run_submitQC(simulation_parameter)
    end_time = time.time  
    print("SubmitQC ends at "+str(end_time))

if __name__ == "__main__":
    main_SubmitQC()
