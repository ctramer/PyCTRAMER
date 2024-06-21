import shutil
from shutil import copyfile
import re, os, sys
import tempfile
import argparse
import fileinput
import time

def shell_copy_template(template_dir, running_dir, filename):
    """
    check into srcdir,
    then copy the file called filename_template.sh to filename.sh
    return filename.sh as a python string
    """
    the_template =f'{template_dir}/{filename}_template.sh'
    the_file =f'{running_dir}/{filename}.sh'
    copyfile(the_template,the_file)
    return the_file

def parse_input_for_jupyter(file_path):
    Simulation_Param = {'PyCTRAMER': '0.1'}
    f = open(file_path)
    inputfile = f.readlines() #fileinput.input()
    for line in inputfile:
        # print("line ",line,len(line))
        for t_len in range(len(line)):
            #print('tlen',t_len, line[t_len+1])
            if line[t_len] == "#":
                line = line[:t_len]
                break
        a = line.split()
        if len(a) >= 2 :
            if a[0][0] == "#":
                pass
            else:
                Simulation_Param.update({a[0]:a[1]})
        elif len(a) >= 1 :
           Simulation_Param.update({a[0]:""})
    return Simulation_Param


def separate_idlist(string):
    id_list = []
    string_list = string.split(',')
    for num in string_list:
        elements = num.split('-')
        if len(elements) == 2:
            for i in range(int(elements[0]),int(elements[1])+1):
                id_list.append(str(i))
        elif len(elements) == 1:
            id_list.append(elements[0])
        else :
            return "error"
    return id_list


def get_output_perfix(string1, string2, argu):
    return string1 + "_" + string2

def sed_inplace(filename, pattern, repl):
    '''
    repl: string to be replaced
    pattern: pattern to be found and matched
    Perform the pure-Python equivalent of in-place `sed` substitution: 
    e.g.,
    `sed -i -e 's/'${pattern}'/'${repl}' "${filename}"`.
    '''
    # For efficiency, precompile the passed regular expression.
    pattern_compiled = re.compile(pattern)

    # For portability, NamedTemporaryFile() defaults to mode "w+b" (i.e., binary
    # writing with updating). This is usually a good thing. In this case,
    # however, binary writing imposes non-trivial encoding constraints trivially
    # resolved by switching to text writing. Let's do that.
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp_file:
        with open(filename) as src_file:
            for line in src_file:
                tmp_file.write(pattern_compiled.sub(repl, line))

    # Overwrite the original file with the munged temporary file in a
    # manner preserving file attributes (e.g., permissions).
    shutil.copystat(filename, tmp_file.name)
    shutil.move(tmp_file.name, filename)

def find_job_id(token, jobs_list):
    token_lst = token.split()
    print(token)
    ind = 0
    for lines in token_lst:
        ind+=1
        if lines == 'Submitted':
            jobs_list.append(token_lst[ind+2])
    return jobs_list

def prepare_folder_structure(datapath):
    QC_dir = datapath+'/QC/'
    MD_dir = datapath+"/MD/"
    
    os.popen('mkdir -p '+QC_dir).read()
    os.popen('mkdir -p '+MD_dir).read()
    return QC_dir,MD_dir
    
def get_folder_path(datapath):
    QC_dir = datapath+'/QC/'
    MD_dir = datapath+"/MD/"
    return QC_dir,MD_dir

def list_to_dep_line(jobs_list):
    # print(len(jobs_list))
    inner_str = jobs_list[0]
    for j in range(len(jobs_list)-1):
        inner_str = inner_str + ','+ jobs_list[j+1]
    return inner_str

def parse_input():
    inputfile = fileinput.input()
    Simulation_Param = {'PyCTRAMER': '0.1'}
    #f = open(file_path)
    #inputfile = f.readlines() #fileinput.input()
    for line in inputfile:
        # print("line ",line,len(line))
        for t_len in range(len(line)):
            #print('tlen',t_len, line[t_len+1])
            if line[t_len] == "#":
                line = line[:t_len]
                break
        a = line.split()
        if len(a) >= 2 :
            if a[0][0] == "#":
                pass
            else:
                Simulation_Param.update({a[0]:a[1]})
        elif len(a) >= 1 :
           Simulation_Param.update({a[0]:""})
    return Simulation_Param



#def parse_input():
#    inputfile = fileinput.input()
#    Simulation_Param= {}
#
#    for line in inputfile: 
#        a = line.split()
#        if len(a) >= 2 :
#            if a[0][0] == "#":
#                pass
#            else:
#                Simulation_Param.update({a[0]:a[1]})
#        if len(a) == 1:
#            Simulation_Param.update({a[0]:""})
#    return Simulation_Param

def sort_and_index(df,column_name):
    """
    df : dataframe
    column_name : string , header of column to be sorted
    sort data
    """
    final_list=[]
    mylist= list(df[column_name])
    sorted_list = sorted(mylist,reverse=True)
    for i in sorted_list:
        a = mylist.index(i) + 1 # 1 because python index start from 0 and excited states = from 1 
        final_list.append(a)
    return final_list, sorted_list

def write_output(filename, new_line_of_string):
    f = open(filename,'w')
    f.write(new_line_of_string + '\n')
    f.close()    
    return


def init_job_control(case_id_list):
    job_control = []
    for caseid in case_id_list:
        current_job_control = {}
        current_job_control['caseid']    = caseid
        current_job_control['QC_jobids'] = []
        current_job_control['MD_jobids'] = []
        current_job_control['status'] = [0,'unstarted','unstarted','unstarted','unstarted','unstarted','unstarted','unstarted']
        job_control.append(current_job_control)
    
    return job_control 

def init_job_info(case_id_list):
    job_info = []
    for caseid in case_id_list:
        current_job_info = {"namemax":""}
        job_info.append(current_job_info)
    
    return job_info

def check_slurm_status(job_ids_list):
    """ if all jobs are finished, return True; else, return false  """
    list_of_str_job_ids = job_ids_list.split(",")
    for jobid in list_of_str_job_ids:
        one_piece = os.popen('sacct -j ' + jobid).read()
        num_Complete = one_piece.count("COMPLETED")
        num_Line     = one_piece.count("\n")
        if num_Complete != num_Line-2 or num_Line == 2:
            return False
    return True

def check_job_unfinished(job_control,ctr_str):
    if ctr_str == "Marcus":
        final_step = 5
    elif ctr_str == "FGR":
        final_step = 7
    for i in range(len(job_control)):
        if job_control[i]['status'][final_step] != 'finished':
            return True
    return False


def check_job_status(job_control, dict_of_simulation, job_info):
    for i in range(len(job_control)):
        curr_step = job_control[i]['status'][0]
        if curr_step == 0 and job_control[i]['status'][curr_step+1] == 'unstarted':
            if dict_of_simulation['skip_QC'] == 'True' or dict_of_simulation['skip_QC'] == 'true':
                curr_step = 1
                job_control[i]['status'][curr_step] = 'finished'
            else: 
                print("==> Start current step # ", curr_step+1, ": SubmitQC")
                job_control = run_step1(i,job_control,dict_of_simulation)
                curr_step = 1
                job_control[i]['status'][curr_step] = 'running'
            
        if curr_step == 1 and job_control[i]['status'][curr_step] == 'running':
            if check_slurm_status(job_control[i]['QC_jobids']):
                job_control[i]['status'][curr_step] = 'finished'
            else: 
                pass
            if dict_of_simulation['skip_QC']=='True' or dict_of_simulation['skip_QC'] == 'true':
                job_control[i]['status'][curr_step] = 'finished'
            else:
                pass
              
            
        if curr_step == 1 and job_control[i]['status'][curr_step] == 'finished':
            print("==> Finished current step # ", curr_step, '\n')

            print("\n==> Start current step # ", curr_step+1,": AnalyzeQC")
            job_control = run_step2(i,job_control,dict_of_simulation)
            curr_step = 2
            job_control[i]['status'][curr_step] = 'finished'
        
        if curr_step == 2 and job_control[i]['status'][curr_step] == 'finished':
            print("\n==> Finished current step # ", curr_step, '\n')
            print(' ') 
            print("==> Start current step # ", curr_step+1, ": ConstructFF")
            job_control = run_step3(i,job_control,dict_of_simulation, job_info)
            curr_step = 3
            job_control[i]['status'][curr_step] = 'finished'
        
        if curr_step == 3 and job_control[i]['status'][curr_step] == 'finished':    
            print("\n==> Finished current step # ", curr_step, '\n')  
            print(' ')           
            print("==> Start current step # ", curr_step+1, ": RunMD")
            job_control = run_step4(i,job_control,dict_of_simulation)
            curr_step = 4
            job_control[i]['status'][curr_step] = 'running'  
            print("MD_jobids",job_control[i]['MD_jobids'] )
        
        if curr_step == 4 and job_control[i]['status'][curr_step] == 'running':
            if check_slurm_status(job_control[i]['MD_jobids']):
                print("job_control[i]['MD_jobids']",job_control[i]['MD_jobids'])
                job_control[i]['status'][curr_step] = 'finished'
                print("\n==> Finished current step # ", curr_step, '\n')
                 
            else: 
                # print("MD_jobids",job_control[i]['MD_jobids'] )
                pass
            
        if curr_step == 4 and job_control[i]['status'][curr_step] == 'finished':
            print("==> Start current step # ", curr_step+1, ": MarcusAA")
            job_control = run_step5(i,job_control,dict_of_simulation,job_info)
            curr_step = 5
            job_control[i]['status'][curr_step] = 'finished'
            print("\n==> Finished current step # ", curr_step, '\n')
            

        if curr_step == 5 and job_control[i]['status'][curr_step] == 'finished':
            print("==> Start current step # ", curr_step+1, ": SpinBoson")
            job_control = run_step6(i,job_control,dict_of_simulation,job_info)
            curr_step = 6 
            job_control[i]['status'][curr_step] = 'finished'
            print("\n==> Finished current step # ", curr_step, '\n')
             

        if curr_step == 6 and job_control[i]['status'][curr_step] == 'finished':
            print("==> Start current step # ", curr_step+1, ": FGR")
            job_control = run_step7(i,job_control,dict_of_simulation,job_info)
            curr_step = 7 
            job_control[i]['status'][curr_step] = 'finished'       
            print("\n==> Finished current step # ", curr_step, '\n')
             
        job_control[i]['status'][0] = curr_step
        
    return job_control

def write_to_output(i, job_control, dict_of_simulation, strr) :
    outdir = dict_of_simulation['work_dir'] +'/'+ dict_of_simulation['project']+'_'+dict_of_simulation['caseid'] 
    
    filename = dict_of_simulation['output_filename']
    if len(filename) == 0:
        filename = dict_of_simulation['project']+'_'+dict_of_simulation['caseid']+ '_PyCTRAMER.out'
    
    outfile = outdir + '/' + filename
    #print(len(job_control),job_control, i)
    #print("dict_of_simulation['workflow']" )
    if job_control[i]['status'][0] == 0 and (dict_of_simulation['workflow'] == "alltheway" or dict_of_simulation['workflow'] == "Marcus") :
        fo = open(outfile, "w+")
        fo.write(header)
        fo.write(strr)
        fo.close()
    elif job_control[i]['status'][0] == 5 and dict_of_simulation['workflow'] == "FGR" :
        # print(dict_of_simulation['workflow'], "FGR", job_control[i]['status'] )
        with open(outfile, "w+") as fo:
            fo.write(header)
            fo.write(strr)
            fo.close()
    else: 
        # print("job_control[i]['status']",job_control[i]['status'])
        with open(outfile,'a') as fo:
            fo.write(strr)
            fo.close()
    # fo.write(strr)
    #fo.close()
    return



cm2au =  4.55633525e-6
au2cm = 219474.6313632
ps2au = 41341.373335
au2ps = 2.4188843266e-5
fs2au = 41.341373335
au2fs = 2.4188843266e-3
au2kT = 315775.024804
kT2au =  3.16681156e-6
au2eV = 27.211386246
eV2au = 0.0367493221
s2au  =4.1341374575751E+16
space_str = "    "
arrow_str = "==> "  
header = "        #######################################\n        #              PyCTRAMER              #\n        #######################################\n"


