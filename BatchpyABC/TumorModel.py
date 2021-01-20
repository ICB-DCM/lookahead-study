#!/usr/bin/env python
# coding: utf-8

# Multi-scale model: Tumor spheroid growth
# ======================

# In[1]:


from time import time

import matplotlib.pyplot as plt
from string import capwords
import os
import tempfile

import pyabc
import tumor2d

import argparse

parser = argparse.ArgumentParser(description='Parse necessary arguments')
parser.add_argument('-pt', '--port', type=str, default="50004",
                    help='Which port should be used?')
parser.add_argument('-ip', '--ip', type=str,
                    help='Dynamically passed - BW: Login Node 3')
args = parser.parse_args()

port = args.port
host = args.ip

pop_size = 1000
min_eps = 700
min_eps_ori = min_eps
max_nr_pop = 100
logfilepath = "/p/home/jusers/reck1/juwels/scripts/programs/results/TumorStats.csv"
resultfile = open("/p/home/jusers/reck1/juwels/scripts/programs/results/TumorRuntimes.txt", "a")

db_path = "sqlite:///" + os.path.join("/p/home/jusers/reck1/juwels/scripts/programs/results", "TumorRes.db")

db_path_ori = "sqlite:///" + os.path.join("/p/home/jusers/reck1/juwels/scripts/programs/results", "TumorRes_ori.db")

# In[2]:


start_time = time()
observation = tumor2d.simulate(division_rate=4.17e-2,
                       initial_spheroid_radius=1.2e1,
                       initial_quiescent_cell_fraction=7.5e-1,
                       division_depth=100,
                       ecm_production_rate=5e-3,
                       ecm_degradation_rate=8e-4,
                       ecm_division_threshold=1e-2)
print(f"Simulation took {time() - start_time:.2f}s")

# In[3]:


fig, axes = plt.subplots(ncols=3)
fig.set_size_inches((16, 5))

color = {"growth_curve":
             "k",
         "extra_cellular_matrix_profile":
             "green",
         "proliferation_profile":
             "orange"}

x_label = {"growth_curve":
               "Time (d)",
           "extra_cellular_matrix_profile":
               "Distance to rim (μm)",
           "proliferation_profile":
               "Distance to rim (μm)"}

y_label = {"growth_curve":
               "Radius (μm)",
           "extra_cellular_matrix_profile":
               "Extracellular matrix intensity",
           "proliferation_profile":
               "Fraction proliferating cells"}

for ax, (key, val) in zip(axes, observation.items()):
    ax.plot(val, color=color[key])
    ax.set_title(capwords(key.replace("_", " ")))
    ax.set_xlabel(x_label[key])
    ax.set_ylabel(y_label[key])
    if key.endswith("profile"):
        ax.set_xlim(0, 600)

# In[4]:


limits = dict(log_division_rate=(-3, -1),
              log_division_depth=(1, 3),
              log_initial_spheroid_radius=(0, 1.2),
              log_initial_quiescent_cell_fraction=(-5, 0),
              log_ecm_production_rate=(-5, 0),
              log_ecm_degradation_rate=(-5, 0),
              log_ecm_division_threshold=(-5, 0))

prior = pyabc.Distribution(**{key: pyabc.RV("uniform", a, b - a)
                        for key, (a,b) in limits.items()})

data_mean = tumor2d.load_default()[1]  # (raw, mean, var)

# In[6]:


redis_sampler = pyabc.sampler.RedisEvalParallelSampler(host=host, port=port, look_ahead = False)

starttime=time()
abc = pyabc.ABCSMC(models=tumor2d.log_model, 
                   parameter_priors=prior,
                   distance_function=tumor2d.distance, 
                   population_size=pop_size, 
                   sampler=redis_sampler)

abc.new(db_path_ori, data_mean)
history_f = abc.run(max_nr_populations=max_nr_pop, minimum_epsilon=min_eps_ori)
endtime=time()

resultfile.write("Ori, " + str(endtime-starttime)+", " + str(pop_size) + ", " + str(min_eps))

# In[7]:


df, w = history_f.get_distribution(m=0,t=history_f.max_t)
pyabc.visualization.plot_kde_matrix(df, w, limits=limits);
plt.savefig("/p/home/jusers/reck1/juwels/scripts/programs/img/TumorResOri.jpg")


# In[8]:


redis_sampler = pyabc.sampler.RedisEvalParallelSampler(host=host, port=port, look_ahead = True, look_ahead_delay_evaluation=True, log_file=logfilepath)

starttime=time()
abc = pyabc.ABCSMC(tumor2d.log_model, 
                   prior, 
                   tumor2d.distance, 
                   population_size=pop_size, 
                   sampler=redis_sampler)

abc.new(db_path, data_mean)
history = abc.run(max_nr_populations=max_nr_pop, minimum_epsilon=min_eps)
endtime=time()

resultfile.write("LA, " + str(endtime-starttime)+", " + str(pop_size) + ", " + str(min_eps))

# In[10]:


df, w = history.get_distribution(m=0, t=history.max_t)
pyabc.visualization.plot_kde_matrix(df, w, limits=limits);
plt.savefig("/p/home/jusers/reck1/juwels/scripts/programs/img/TumorResPPP.jpg")


