from .SubmitQC import run_step1
from .AnalyzeQC import run_step2
from .ConstructFF import run_step3
from .RunMD import run_step4
from .MarcusAA import run_step5
from .SpinBoson import run_step6
from .FGR import run_step7
from .Utilities import * 

def run_Marcus(dict_of_simulation):
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
    # print(job_control)
    # job infomation
    job_info = init_job_info(case_id_list)
    ctr_str = "Marcus"
    count = 0
    while check_job_unfinished(job_control,ctr_str):
        count +=1 
        job_control = check_job_status(job_control, dict_of_simulation, job_info)
        if count % 1000 ==0:
            print(count,job_control)
        time.sleep(1)
    return 
        
def run_FGR(dict_of_simulation):
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
    # print(job_control)
    ctr_str = "FGR"
    # job infomation
    job_info = init_job_info(case_id_list)
    count = 5
    for i in range(len(job_control)):
        job_control[i]['status'][0] = 5
        job_control[i]['status'][5] = 'finished'
    while check_job_unfinished(job_control,ctr_str):
        count +=1 
        job_control = check_job_status(job_control, dict_of_simulation, job_info)
def check_job_status(job_control, dict_of_simulation, job_info):
    for i in range(len(job_control)):
        curr_step = job_control[i]['status'][0]
        if curr_step == 0 and job_control[i]['status'][curr_step+1] == 'unstarted':
            if dict_of_simulation['skip_QC'] == 'True' or dict_of_simulation['skip_QC'] == 'true':
                curr_step = 1
                job_control[i]['status'][curr_step] = 'finished'
            else:
                job_control = run_step1(i,job_control,dict_of_simulation)
                curr_step = 1
                job_control[i]['status'][curr_step] = 'running'
                print("==> Started current step # ", curr_step, ": SubmitQC")
            
        if curr_step == 1 and job_control[i]['status'][curr_step] == 'running':
            if check_slurm_status(job_control[i]['QC_jobids']):
                job_control[i]['status'][curr_step] = 'finished'
                job_control[i]['status'][0] = curr_step
            else: 
                pass
            if dict_of_simulation['skip_QC'] == 'True' or dict_of_simulation['skip_QC'] == 'true' :
                job_control[i]['status'][curr_step] = 'finished'
                job_control[i]['status'][0] = curr_step
            else:
                pass
              
            
        if curr_step == 1 and job_control[i]['status'][curr_step] == 'finished':
            if dict_of_simulation['skip_FF'] == 'True' or dict_of_simulation['skip_FF'] == 'true':
                curr_step = 3
                job_control[i]['status'][curr_step-1] = 'finished'
                job_control[i]['status'][curr_step] = 'finished'
                job_control[i]['status'][0] = curr_step
            else: 
                print("==> Started current step # ", curr_step+1, ": AnalyzeQC")
                job_control = run_step2(i,job_control,dict_of_simulation)
                curr_step = 2
                job_control[i]['status'][curr_step] = 'finished'
                print("==> Finished current step # ", curr_step)
                job_control[i]['status'][0] = curr_step
 
        if curr_step == 2 and job_control[i]['status'][curr_step] == 'finished':
            print("==> Started current step # ", curr_step+1, ": ConstructFF")
            job_control = run_step3(i,job_control,dict_of_simulation, job_info)
            curr_step = 3
            job_control[i]['status'][curr_step] = 'finished'
            print("==> Finished current step # ", curr_step)
            job_control[i]['status'][0] = curr_step

        if curr_step == 3 and job_control[i]['status'][curr_step] == 'finished':
            if dict_of_simulation['skip_MD'] == 'True' or dict_of_simulation['skip_MD'] == 'true':
                curr_step = 4
                job_control[i]['status'][curr_step] = 'finished'
                job_control[i]['status'][0] = curr_step
            else: 
                print("==> Started current step # ", curr_step+1, ": RunMD")
                job_control = run_step4(i,job_control,dict_of_simulation, job_info)
                curr_step = 4
                job_control[i]['status'][curr_step] = 'running' 
 
        
        if curr_step == 4 and job_control[i]['status'][curr_step] == 'running':
            if check_slurm_status(job_control[i]['MD_jobids']):
                print("job_control[i]['MD_jobids']",job_control[i]['MD_jobids'])
                job_control[i]['status'][curr_step] = 'finished'
                job_control[i]['status'][0] = curr_step
            else: 
                # print("MD_jobids",job_control[i]['MD_jobids'] )
                pass
                    
        if curr_step == 4 and job_control[i]['status'][curr_step] == 'finished':
            print("==> Finished current step # ", curr_step)
            print("==> Started current step # ", curr_step+1, ": MarcusAA")
            if dict_of_simulation['workflow'] == "Marcus":
                print("==> Finished current step # ", curr_step)
                job_control[i]['status'][0] = curr_step
                break
            job_control = run_step5(i,job_control,dict_of_simulation,job_info)
            curr_step = 5
            job_control[i]['status'][curr_step] = 'finished'
            print("==> Finished current step # ", curr_step)
            job_control[i]['status'][0] = curr_step
        # print(job_control) 
        
        if curr_step == 5 and job_control[i]['status'][curr_step] == 'finished':
            print("==> Started current step # ", curr_step+1, ": SpinBoson")
            job_control = run_step6(i,job_control,dict_of_simulation,job_info)
            curr_step = 6
            # print(job_control[i])
            job_control[i]['status'][curr_step] = 'finished'
            print("==> Finished current step # ", curr_step)
            job_control[i]['status'][0] = curr_step

        if curr_step == 6 and job_control[i]['status'][curr_step] == 'finished':
            print("==> Started current step # ", curr_step+1, ": FGR")
            job_control = run_step7(i,job_control,dict_of_simulation,job_info)
            curr_step = 7
            job_control[i]['status'][curr_step] = 'finished'
            print("==> Finished current step # ", curr_step)
            job_control[i]['status'][0] = curr_step

        job_control[i]['status'][0] = curr_step
        
    return job_control
