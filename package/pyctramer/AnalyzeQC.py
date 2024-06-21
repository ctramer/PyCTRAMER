from .Utilities import *  
import os
import numpy as np 

def run_step2(i, job_control, dict_of_simulation):
    #print(arrow_str,"Started Step 2: AnalyzeQC")
    project  = dict_of_simulation['project']
    work_dir = dict_of_simulation['work_dir']
    caseid =  job_control[i]['caseid'] 
    case_dir = work_dir + '/' + project +"_"+ caseid
    qc_dir,md_dir =  get_folder_path(case_dir)
    
    template_dir = dict_of_simulation['template_dir']

    N_atom = dict_of_simulation['N_atom']
    N_cis  = dict_of_simulation['N_cis']

    QC_method   = dict_of_simulation['QC_method']
    basis_set   = dict_of_simulation['basis_set']
    RSH_omega = dict_of_simulation['RSH_omega']
    TD_approx = dict_of_simulation['TD_approx'] 
    charge_type = dict_of_simulation['charge_type']

    DA_type     = dict_of_simulation['DA_type']
    D_range  =  dict_of_simulation['D_range']
    B_range  =  dict_of_simulation['B_range']
    A_range  =  dict_of_simulation['A_range']
    A1_range = dict_of_simulation['A1_range']
    A2_range = dict_of_simulation['A2_range']
    CT_donor_fragment = dict_of_simulation['CT_donor_fragments']
    CT_acceptor_fragment = dict_of_simulation['CT_acceptor_fragments'] 

    LE_fragment = dict_of_simulation['LE_fragment']
    CT_amount_threshold = dict_of_simulation['CT_amount_threshold']
    OS_threshold = dict_of_simulation['OS_threshold']
    LE_position_threshold = dict_of_simulation['LE_position_threshold']
    TD_approx = dict_of_simulation["TD_approx"] 
    analyzeQC(template_dir, qc_dir, QC_method, basis_set, caseid, D_range,B_range,A_range,A1_range,A2_range,LE_fragment, N_cis, project, qc_dir, DA_type, RSH_omega,TD_approx)

    stroutput = "\n" + arrow_str
    stroutput = stroutput + "STEP 2: AnalyzeQC \n"
    stroutput = stroutput + "Input: \n" + qc_dir + project +"_"+ caseid + ".out \n" 
    stroutput = stroutput + qc_dir + project +"_"+ caseid + ".fchk \n" 
    stroutput = stroutput + qc_dir + "dens_ana.in \n"
    stroutput = stroutput + qc_dir + "dens_ana_fchk.in \n"
    stroutput = stroutput + qc_dir + "dens_ana_tddft.in \n\n"
    stroutput = stroutput + "Job submit script: \n"+ qc_dir +"analyzeQC_flexible1.sh \n\n" 
    stroutput = stroutput + "Output: \n" + qc_dir + "tden_summ.txt \n"
    stroutput = stroutput + qc_dir + "tden_summ_fchk.txt \n"
    stroutput = stroutput + qc_dir + "tden_summ_tddft.txt \n"
    stroutput = stroutput + "#######################################################\n"
    write_to_output(i, job_control, dict_of_simulation, stroutput) 
    return job_control

def analyzeQC(template_dir,running_dir,QC_method,basis_set,testing_str,D_range,B_range,A_range,A1_range,A2_range,LE_fragment,N_cis,project,basedir,DA_type,RSH_omega,TD_approx,dependency=[]):
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
    #os.system('mkdir -p '+running_dir)

    os.chdir(running_dir)
    submitQC_file = shell_copy_template(template_dir,running_dir, filename='analyzeQC')
    
    os.chdir(running_dir)
    QC_DIR = f'{running_dir}'
    #QC_DIR=f'{project_dir}/'+case+'/QC/'
    
    ordersD  = np.fromstring(D_range, dtype=int, sep='-') 
    ordersA  = np.fromstring(A_range, dtype=int, sep='-')
    if B_range != "":
        ordersB  = np.fromstring(B_range, dtype=int, sep='-') 
    if A1_range != "":
        ordersA1  = np.fromstring(A1_range, dtype=int, sep='-') 
    if A2_range != "":
        ordersA2  = np.fromstring(A2_range, dtype=int, sep='-')
    #ordersB  = np.fromstring(B_range, dtype=int, sep='-')
    #ordersA1  = np.fromstring(A1_range, dtype=int, sep='-') 
    #ordersA2  = np.fromstring(A2_range, dtype=int, sep='-')

    if RSH_omega == "":
        sed_inplace(submitQC_file,"omega=","omega=114514")
        # sed_inplace(submitQC_file,"OMEGA", "!OMEGA")
    else:
        sed_inplace(submitQC_file,"omega=","omega="+RSH_omega)
    
    if TD_approx == "TDA":
        TD_approx = "True"
    else:
        TD_approx = "False"

    sed_inplace(submitQC_file,'WORKDIR=','WORKDIR='+QC_DIR)
    sed_inplace(submitQC_file,'SRCDIR=','SRCDIR='+template_dir)
    sed_inplace(submitQC_file,'RUNDIR=','RUNDIR='+QC_DIR)
    sed_inplace(submitQC_file,'system=','system='+project+'_')
    sed_inplace(submitQC_file,'basis_set=','basis_set='+basis_set)
    sed_inplace(submitQC_file,'method=','method='+QC_method)
    sed_inplace(submitQC_file,'D_range=','D_range='+D_range)
    sed_inplace(submitQC_file,'B_range=','B_range='+B_range)
    sed_inplace(submitQC_file,'A_range=','A_range='+A_range)
    sed_inplace(submitQC_file,'STATE_NUM=','STATE_NUM='+str(N_cis))
    sed_inplace(submitQC_file,'GIVEN_STRUC=','GIVEN_STRUC='+ testing_str)
    sed_inplace(submitQC_file,'LE_fragment=','LE_fragment='+LE_fragment)
    sed_inplace(submitQC_file,'TD_approx=','TD_approx='+TD_approx)
    
    frag_str = "'[["
    if DA_type == "DA" and max(ordersD) < max(ordersA):
        int1 = max(ordersD)-min(ordersD)+1
        for i in range(int1):
            frag_str = frag_str + str(i+1) 
            if i != (int1-1):
                frag_str = frag_str + ","
            else :  
                frag_str = frag_str + "],["
        # frag_str = frag_str + "],["
        
        int1 = max(ordersA)-min(ordersA)+1
        for i in range(int1): 
            frag_str = frag_str + str(max(ordersD)-min(ordersD) +2+i) 
            if i != (max(ordersA)-1 ):
                frag_str = frag_str + ","
        frag_str = frag_str +"]]'  #"
        # print(int1,max(ordersD)-min(ordersD),max(ordersA)-min(ordersA),'max(ordersA)-min(ordersA)')
        
    elif DA_type == "DA" and max(ordersA) < max(ordersD):
        int1 = max(ordersA)-min(ordersA)+1
        for i in range(int1): 
            frag_str = frag_str + str(1+i) 
            if i != (int1-1):
                frag_str = frag_str + ","
            else: 
                frag_str = frag_str + "],["
                
        int1 = max(ordersD)-min(ordersD)+1
        for i in range(int1):
            frag_str = frag_str + str(i+1+max(ordersA)-min(ordersA) +1) 
            if i != (int1-1):
                frag_str = frag_str + ","    
        frag_str = frag_str +"]]'  #"
        # print(int1,max(ordersD)-min(ordersD),max(ordersA)-min(ordersA),'max(ordersA)-min(ordersA)')
    
    
    elif DA_type == "DBA":
        # seglist should be DBA order
        # for example triad: 60+75+72 = 207, its corresponding segorder should be 60,75,72

        daseg = np.array([ordersA[1]-ordersA[0]+1,ordersB[1]-ordersB[0]+1,ordersD[1]-ordersD[0]+1  ] )
        
        accseg = 0
        for ind,seg in enumerate(daseg):

            for i in range(seg):
                if ind == 0:
                    frag_str = frag_str+str(i+1)
                else:
                    frag_str = frag_str+str(i+accseg+1)
                if i!= seg-1:
                    frag_str = frag_str + ","
                else:
                    if ind == 2:
                        frag_str = frag_str + "]]'  #"
                    else:
                        frag_str = frag_str + "],["
            accseg += seg 
    # elif DA_type == "DBA":
    #     # seglist should be DBA order
    #     # for example triad: 60+75+72 = 207, its corresponding segorder should be 60,75,72
    #     daseg = np.array([ordersD[1]-ordersD[0]+1,ordersB[1]-ordersB[0]+1,ordersA[1]-ordersA[0]+1  ] )
    #     
    #     for ind,seg in enumerate(daseg):
    #         for i in range(seg):
    #             if ind == 0:
    #                 frag_str = frag_str+str(i+1)
    #             else:
    #                 frag_str = frag_str+str(i+1+daseg[ind-1])
    #             if i!= seg-1:
    #                 frag_str = frag_str + ","
    #             else:
    #                 if ind == 2:
    #                     frag_str = frag_str + "]]'"
    #                 else:
    #                     frag_str = frag_str + "],["
            
        
    elif DA_type == "ADA":
        # seglist should be D-A1-A2 order
        # daseg  = np.fromstring(daseg, dtype=int, sep=',')
        daseg = np.array([ordersD[1]-ordersD[0]+1,ordersA1[1]-ordersA1[0]+1,ordersA2[1]-ordersA2[0]+1  ] )
         
        for ind,seg in enumerate(daseg):
            for i in range(seg):
                if ind == 0:
                    frag_str = frag_str+str(i+1)
                else:
                    frag_str = frag_str+str(i+1+daseg[ind-1])
                if i!= seg-1:
                    frag_str = frag_str + ","
                else:
                    if ind == 2:
                        frag_str = frag_str + "]]'"
                    else:
                        frag_str = frag_str + "],["
    
    sed_inplace(submitQC_file,"fraginfo=","fraginfo="+frag_str)
    os.chdir(running_dir)
     
    jobinfo =  os.popen('sh analyzeQC.sh').read()
    
    return jobinfo

def main_AnalyzeQC():
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
    print("AnalyzeQC begins at "+str(startime))
    # print("Input parameters ")
    # print(simulation_parameter)
    for ind, ctr_ele in enumerate(job_control): 
        print(ind, job_control)
        run_step2(ind, job_control, dict_of_simulation)

    # run_submitQC(simulation_parameter)
    end_time = time.time  
    print("AnalyzeQC ends at "+str(end_time))

if __name__ == "__main__":
    main_AnalyzeQC()
