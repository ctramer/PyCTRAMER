# FGR
import numpy as np
import pandas as pd
from .Utilities import *

def coth(x):
    return 1. / np.tanh(x)

sin = np.sin 
cos = np.cos 
exp = np.exp 
I = 1j 

def run_step7(i, job_control, dict_of_simulation, job_info): 
    # print('')
    #print(arrow_str,"Started Step 7: FGR")
    project  = dict_of_simulation['project']
    work_dir = dict_of_simulation['work_dir']
    caseid =  job_control[i]['caseid'] # caseid = str(i)
    case_dir = work_dir + '/' + project +"_"+ caseid
    qc_dir,md_dir =  get_folder_path(case_dir)
    template_dir = dict_of_simulation['template_dir']
    temperature  = float(dict_of_simulation['MD_temperature'])
    beta = 1. / (temperature * kT2au)
    MD_dyn_state = dict_of_simulation['MD_traj_state']

    sbmsym = np.loadtxt(case_dir+'/H_sys.dat')
    sbmbath = np.loadtxt(work_dir+'/SBM'+MD_dyn_state+'CT.dat')
    omega = sbmbath[:,1]
    req = sbmbath[:,3]

    inttt = 9
    points = 2**inttt 
    point2 = 2**(inttt-1)

    dt_tcf_sbm = 30 #$1 #0.5 # in au 
    dw_QM = 2. * np.pi / points / dt_tcf_sbm # in au 
    w_QM = np.fft.fftfreq(points)*points*dw_QM
    t = np.arange(points) * dt_tcf_sbm
    t = t - dt_tcf_sbm * point2
    dt_au = t[1] - t[0]
    Gamma_DA = sbmsym[0,1]
    hw_DA = sbmsym[0,0] - sbmsym[1,1]


    #print('domega (eV) = ', dw_QM * au2eV )
    #print("omega span (eV) = ", dw_QM * au2eV * point2)
    print('-DeltaE = ', hw_DA * au2eV, ' eV\nGamma_DA = ', Gamma_DA*au2eV, ' eV\nEr = ', np.sum(omega**2*req**2)/2.*au2eV, ' eV')
    print("The following 2 rate constants are from direct integral and FFT, respectively. ")
    print('QM  ', ApproxRates(Ct_exact,omega,req,beta,t,hw_DA,Gamma_DA)*ps2au*dt_au*1e12, 'Hz' )
    print('W0  ', ApproxRates(Ct_W0,omega,req,beta,t,hw_DA,Gamma_DA)*ps2au*dt_au*1e12, 'Hz' )
    #print('Marcus-Levitch ', MarcusLevichRate(omega,req,beta,t,hw_DA,Gamma_DA)*ps2au*dt_au,'THz')
    print('CAV ', ApproxRates(Ct_CAV,omega,req,beta,t,hw_DA,Gamma_DA)*ps2au*dt_au*1e12, 'Hz' )
    print('CD  ', ApproxRates(Ct_CD,omega,req,beta,t,hw_DA,Gamma_DA)*ps2au*dt_au*1e12, 'Hz' )
    print('C0  ', ApproxRates(Ct_C0,omega,req,beta,t,hw_DA,Gamma_DA)*ps2au*dt_au*1e12, 'Hz' )
    print('Marcus-Levitch ', MarcusLevichRate(omega,req,beta,t,hw_DA,Gamma_DA)*ps2au*1e12,'Hz' )
    print('Marcus Theory  ', MarcusRate(omega,req,beta,t,hw_DA,Gamma_DA)*ps2au*1e12,'Hz') 
    

    stroutput = ""
    stroutput = stroutput + "\n" + arrow_str + "STEP 7: FGR \n"
    stroutput = stroutput + "Spin Boson model for Linearized Semiclassical FGR Rates \n"
    stroutput = stroutput + "Note: The 2 rate constants are from direct integral and FFT by order.\n"
    stroutput = stroutput + 'QM  ' +str( ApproxRates(Ct_exact,omega,req,beta,t,hw_DA,Gamma_DA)*ps2au*dt_au*1e12 )  + " Hz \n"
    stroutput = stroutput + 'W0  ' +str( ApproxRates(Ct_W0,omega,req,beta,t,hw_DA,Gamma_DA)*ps2au*dt_au*1e12 )  + " Hz \n"
    # stroutput = stroutput + 'Marcus-Levitch rates ' +str( MarcusLevichRate(omega,req,beta,t,hw_DA,Gamma_DA)*ps2au*dt_au*1e12)  + " Hz \n"
    stroutput = stroutput + 'CAV ' +str( ApproxRates(Ct_CAV,omega,req,beta,t,hw_DA,Gamma_DA)*ps2au*dt_au*1e12 )  + " Hz \n"
    stroutput = stroutput + 'CD  ' +str( ApproxRates(Ct_CD,omega,req,beta,t,hw_DA,Gamma_DA)*ps2au*dt_au*1e12 )  + " Hz \n"
    stroutput = stroutput + 'C0  ' +str( ApproxRates(Ct_C0,omega,req,beta,t,hw_DA,Gamma_DA)*ps2au*dt_au*1e12 )  + " Hz \n"
    stroutput = stroutput + 'Marcus-Levitch ' +str( MarcusLevichRate(omega,req,beta,t,hw_DA,Gamma_DA)*ps2au*1e12)  + " Hz \n"
    stroutput = stroutput + 'Marcus Theory  ' +str( MarcusRate(omega,req,beta,t,hw_DA,Gamma_DA)*ps2au*1e12)  + " Hz \n"
    # stroutput =  stroutput + '\n\n'
    stroutput = stroutput + "#######################################################\n"
    stroutput = stroutput + "\nJob Ends at " + time.ctime(time.time()) + "\n"
    write_to_output(i, job_control, dict_of_simulation, stroutput) 

    return job_control 


def ApproxRates(Ct,omega,req,beta,t,hw_DA,Gamma_DA):
    CqmDA = Ct(omega,req,beta,t,hw_DA,Gamma_DA)
    Cqm0 = Ct(omega,req,beta,t,0,Gamma_DA)
    dt_tcf_sbm = t[1] - t[0]
    # rateQM = (np.sum(CqmDA) - CqmDA[0]/2. - CqmDA[-1]/2. ) * dt_tcf_sbm
    rateQM_approach1 = (np.sum(CqmDA) - 0*CqmDA[0]  - 0*CqmDA[-1]  ) # * dt_tcf_sbm  
    rateQM_approach2 = np.fft.ifft(CqmDA)[0] * len(t) 
    return np.array([rateQM_approach1.real,rateQM_approach2.real])

#############################################################
#
# Ref: J. Chem. Phys. 153, 044105 (2020) 
# Eq (26) to Eq (30)
#
#############################################################
 
def Ct_exact(omega,Req,beta,t,omega_DA,Gamma_DA):
    # all in a.u. 
    hbar = 1. 
    shift = -omega * Req * Req / 2. / hbar 
    bhw2 = 0.5 * beta * hbar * omega 
    size = len(t) 
    Ct = np.zeros(size, dtype=complex) 
    for ind, w in enumerate(omega):
        wt = w * t 
        Ct += shift[ind] * (coth(bhw2[ind])*(1 - cos(wt)) + I*sin(wt))
    
    # print(np.shape(Ct),np.shape(omega_DA),np.shape(I))
    Ct += I * omega_DA  * t 
    Ct = Gamma_DA * Gamma_DA * np.exp(Ct)
    return Ct

def Ct_W0(omega,Req,beta,t,omega_DA,Gamma_DA):
    # all in a.u. 
    hbar = 1. 
    shift = -omega * Req * Req / 2. / hbar 
    bhw2 = 0.5 * beta * hbar * omega 
    size = len(t) 
    Ct = np.zeros(size, dtype=complex) 
    for ind, w in enumerate(omega):
        wt = w * t 
        Ct += shift[ind] * (coth(bhw2[ind])*wt*wt/2. + I*wt)
    
    # print(np.shape(Ct),np.shape(omega_DA),np.shape(I))
    Ct += I * omega_DA  * t 
    Ct = Gamma_DA * Gamma_DA * np.exp(Ct)
    return Ct

def MarcusLevichRate(omega,Req,beta,t,omega_DA,Gamma_DA):
    # all in a.u. 
    hbar = 1. 
    shift = -omega * Req * Req / 2. / hbar 
    bhw2 = 0.5 * beta * hbar * omega 
    size = len(t) 
    Ct = np.zeros(size, dtype=complex) 
    a = 0.25 * np.sum(omega * omega * omega * Req * Req / np.tanh(0.5 * beta * omega * hbar) )
    # print('a_parameter = ', a)
    pref = np.sqrt(np.pi/a) 

    Er = 0.5 * np.sum(omega**2*Req**2) 
    DeltaE = - hbar * omega_DA
    # print(Er, DeltaE, Gamma_DA)
    return Gamma_DA**2 * pref * np.exp( - (Er+DeltaE)**2/4. / a )
    #for ind, w in enumerate(omega):
    #    wt = w * t 
    #    Ct += shift[ind] * (coth(bhw2[ind])*wt*wt*2 + I*sin(wt))
    
    # print(np.shape(Ct),np.shape(omega_DA),np.shape(I))
    #Ct += I * omega_DA  * t 
    # Ct = Gamma_DA * Gamma_DA * np.exp(Ct)
    # return Ct
def MarcusRate(omega,Req,beta,t,omega_DA,Gamma_DA): 
    hbar = 1. 
    Er = np.sum(omega**2*Req**2*0.5)
    U = Er - omega_DA * hbar 
    # print(Er*au2eV,Gamma_DA*au2eV, U*au2eV, beta)
    return Gamma_DA ** 2 * np.sqrt(np.pi*beta / Er) * np.exp( - U**2 * beta * 0.25 / Er)

def Ct_CAV(omega,Req,beta,t,omega_DA,Gamma_DA):
    # all in a.u. 
    hbar = 1. 
    shift = -omega * Req * Req / 2. / hbar 
    bhw2 = 0.5 * beta * hbar * omega 
    size = len(t) 
    Ct = np.zeros(size, dtype=complex) 
    for ind, w in enumerate(omega):
        wt = w * t 
        Ct += shift[ind] * (( 1 / bhw2[ind] )*(1 - cos(wt)) + I*sin(wt))
    Ct += I * t * omega_DA 
    Ct = Gamma_DA * Gamma_DA * np.exp(Ct)
    return Ct

def Ct_CD(omega,Req,beta,t,omega_DA,Gamma_DA):
    # all in a.u. 
    hbar = 1. 
    shift = -omega * Req * Req / 2. / hbar 
    bhw2 = 0.5 * beta * hbar * omega 
    size = len(t) 
    Ct = np.zeros(size, dtype=complex) 
    for ind, w in enumerate(omega):
        wt = w * t 
        Ct += shift[ind] * (( 1 / bhw2[ind] )*(1 - cos(wt)) + I*wt)
    Ct += I * t * omega_DA 
    Ct = Gamma_DA * Gamma_DA * np.exp(Ct)
    return Ct

def Ct_C0(omega,Req,beta,t,omega_DA,Gamma_DA):
    # all in a.u. 
    hbar = 1. 
    shift = -omega * Req * Req / 2. / hbar 
    
    size = len(t) 
    Ct = np.zeros(size, dtype=complex) 
    for ind, w in enumerate(omega):
        wt = w * t 
        wttbetahbar = w * t * t / beta / hbar 
        Ct += shift[ind] * ( wttbetahbar + I * wt)
    Ct += I * t * omega_DA 
    Ct = Gamma_DA * Gamma_DA * np.exp(Ct)
    return Ct



import pickle 

def main_FGR():
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

    # job_info = init_job_info(case_id_list)
    startime = time.time ()

    # Read list to memory
    def read_list():
        # for reading also binary mode is important
        with open(work_dir+'/listfile', 'rb') as fp:
            n_list = pickle.load(fp)
            return n_list 

    def write_list(a_list):
        # store list in binary file so 'wb' mode
        with open(work_dir+'/listfile', 'wb') as fp:
            pickle.dump(a_list, fp)
            print('Done writing list into a binary file')
        
    job_info = read_list()
    print("FGR begins at "+str(startime))
    # print("Input parameters ")
    # print(simulation_parameter)
    for ind, ctr_ele in enumerate(job_control): 
        print(ind, job_control)
        run_step7(ind, job_control, dict_of_simulation,job_info)

    write_list(job_info)
    end_time = time.time  ()
    print("FGR ends at "+str(end_time))

if __name__ == "__main__":
    main_FGR()
