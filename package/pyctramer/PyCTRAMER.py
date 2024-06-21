# PyCTRAMER
#standard library import
import numpy as np 
import pandas as pd
import sys

import tempfile
import parmed
from mdtraj.reporters import HDF5Reporter

import math
import os
import re
import shutil
from shutil import copyfile
import tempfile
import json
import pickle as pkl
import matplotlib
import matplotlib.pyplot as plt
 
import mdtraj as md
import argparse
import fileinput
import time

# SubmitQC
from .JobControl import *
from .Utilities import * 


def simulate(input) :
    simulation_parameter =  parse_input_for_jupyter(input)
    startime = time.time()
    print("simulation begins at "+ time.ctime(startime) )
    print("Input parameters ")
    print(simulation_parameter)
    run_simulation(simulation_parameter)
    end_time = time.time()
    print("simulation ends at "+ time.ctime(end_time) )
    timecosu = end_time - startime 
    hour = timecosu // 3600 
    minute = (timecosu - hour * 3600) // 60 
    second = timecosu - 3600 * hour - 60 * minute 
    print("Computational time "+ str(hour) + " hour " + str(minute) + " minutes " + second + " second "  )
    return 0

def main() :
    simulation_parameter =  parse_input()
    startime = time.time 
    print("simulation begins at "+str(startime ))
    print("Input parameters ")
    print(simulation_parameter)
    run_simulation(simulation_parameter)
    end_time = time.time  
    print("simulation ends at "+str(end_time ))
    return 0
    
def run_simulation(dict_of_simulation):
    if dict_of_simulation['workflow'] == 'Marcus':
        run_Marcus(dict_of_simulation)
    elif dict_of_simulation['workflow'] == 'FGR':
        run_FGR(dict_of_simulation)
    elif dict_of_simulation['workflow'] == 'alltheway':
        run_Marcus(dict_of_simulation)
        run_FGR(dict_of_simulation)
        # OR define a run_all_the_way(dict_of_simulation)
    else :
        return 'Not defined yet'   
#if __name__ == "__main__":
#    main(); 


