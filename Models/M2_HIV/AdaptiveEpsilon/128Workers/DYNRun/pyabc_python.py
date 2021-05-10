# -*- coding: utf-8 -*-
import pyabc
import numpy as np
from pyabc.external import R
from pyabc import Distribution, RV, ABCSMC
from pyabc import populationstrategy
from pyabc.distance.scale import median
import fitmulticell as fmc
import os
from pyabc.sampler import RedisEvalParallelSampler
import argparse
import glob
import pandas as pd

parser = argparse.ArgumentParser(description='Parse necessary arguments')
parser.add_argument('-pt', '--port', type=str, default="50004",
                    help='Which port should be use? Nils has 50004-50007')
parser.add_argument('-ip', '--ip', type=str,
                    help='BQ: Dynamically passed - BW: Login Node 3')
args = parser.parse_args()

# Password for the redis_server
port_to_pw = pd.read_csv('/home/hd/hd_hd/hd_ty170/port_pw.csv',
                         squeeze=True,
                         index_col=0,
                         header=None).to_dict()
password=port_to_pw[int(args.port)]

# Morpheus Model
cur_path = os.getcwd()
filename = glob.glob('*.xml')[0]
file_ = os.path.join(cur_path, filename)
par_map = {'PS_tar': './Global/Constant[@symbol="PS_tar"]',
           'DT_tar': './Global/Constant[@symbol="DT_tar"]',
           'PS_inf': './Global/Constant[@symbol="PS_inf"]',
           'DT_inf': './Global/Constant[@symbol="DT_inf"]',
           'JTT': './Global/Constant[@symbol="JTT"]',
           'JTI': './Global/Constant[@symbol="JTI"]',
           'JIM': './Global/Constant[@symbol="JIM"]',
           'JTM': './Global/Constant[@symbol="JTM"]',
           'JTC': './Global/Constant[@symbol="JTC"]',
           'JIC': './Global/Constant[@symbol="JIC"]',
           'JII': './Global/Constant[@symbol="JII"]',
           'C_scs': './Global/Constant[@symbol="C_scs"]',
           'C_vcs': './Global/Constant[@symbol="C_vcs"]'}

model = fmc.model.MorpheusModel(model_file =  file_, 
                                sumstat_funs = [],
                                par_map=par_map,
                                show_stdout=False,
                                show_stderr=False,
                                raise_on_error=False)

# R - Distances and sumstats
r = R(os.path.join(cur_path, 'abc_r_blueprint.R'))

# Target distances
dist_speed_tar = r.distance("distance_speed_tar")
dist_trn_tar  = r.distance("distance_trn_tar")
dist_arr_tar   = r.distance("distance_arr_tar")
dist_str_tar   = r.distance("distance_str_tar")
dist_msd_tar   = r.distance("distance_msd_tar")

# Infected distances
dist_speed_inf = r.distance("distance_speed_inf")
dist_trn_inf  = r.distance("distance_trn_inf")
dist_arr_inf   = r.distance("distance_arr_inf")
dist_str_inf   = r.distance("distance_str_inf")
dist_msd_inf   = r.distance("distance_msd_inf")

# Aggregate distances
distance = pyabc.AggregatedDistance([dist_speed_tar,
                                             dist_trn_tar,
                                             dist_arr_tar,
                                             dist_str_tar,
                                             dist_msd_tar,
                                             dist_speed_inf,
                                             dist_trn_inf,
                                             dist_arr_inf,
                                             dist_str_inf,
                                             dist_msd_inf])

# Sumstats
sum_stat = r.summary_statistics("mySummaryStatistics", is_py_model=True) # If pyabc ver 0.9.22 add is_py_model=True

dt_domain = np.arange(60)

prior = Distribution(PS_tar=RV("uniform", 0, 100),  # low,low+interval
                     DT_tar=RV("rv_discrete", values = (dt_domain, [1/60]*60)),
                     PS_inf=RV("uniform", 0, 100),
                     DT_inf=RV("rv_discrete", values = (dt_domain, [1/60]*60)),
                     JTT=RV("uniform", 0, 500),  
                     JTI=RV("uniform", 0, 500),  
                     JIM=RV("uniform", 0, 500),
                     JTM=RV("uniform", 0, 500),
                     JTC=RV("uniform", 0, 500), 
                     JIC=RV("uniform", 0, 500),
                     JII=RV("uniform", 0, 500),
                     C_scs=RV("uniform", 0, 100),
                     C_vcs=RV("uniform", 0, 100)) 


transition = pyabc.transition.AggregatedTransition(mapping={
    'PS_tar': pyabc.MultivariateNormalTransition(),
    'DT_tar': pyabc.DiscreteJumpTransition(domain=dt_domain, p_stay=0.7),
    'PS_inf': pyabc.MultivariateNormalTransition(),
    'DT_inf': pyabc.DiscreteJumpTransition(domain=dt_domain, p_stay=0.7),
    'JTT': pyabc.MultivariateNormalTransition(),
    'JTI': pyabc.MultivariateNormalTransition(),
    'JIM': pyabc.MultivariateNormalTransition(),
    'JTM': pyabc.MultivariateNormalTransition(),
    'JTC': pyabc.MultivariateNormalTransition(),
    'JIC': pyabc.MultivariateNormalTransition(),
    'JII': pyabc.MultivariateNormalTransition(),
    'C_scs': pyabc.MultivariateNormalTransition(),
    'C_vcs': pyabc.MultivariateNormalTransition()})


redis_sampler = RedisEvalParallelSampler(host=args.ip, port=args.port, password=password, look_ahead=False, look_ahead_delay_evaluation=False, log_file='look_ahead.csv')

import time
time.sleep(10)
redis_sampler.stop()

# ABCSMC
pop_strategy = pyabc.ConstantPopulationSize(nr_particles=256) 
eps = pyabc.QuantileEpsilon(65, alpha=0.7)  # Defines which distances is needed in next populatio

abc = ABCSMC(model, prior, distance, summary_statistics=sum_stat, transitions=transition, 
             sampler=redis_sampler, population_size=pop_strategy,
             eps=eps)

# Database
db_name = str.split(filename, sep = '.')[0]
db = "sqlite:///" + os.path.join(cur_path, db_name + '.db') # Define name of the database

abc.new(db, r.observation("Exp"))
#abc.load(db, 1, r.observation("Exp"))

# Abbruchkriterium
history = abc.run(minimum_epsilon=5)
