"""Same as the .ipynb file; converted to a python file to run it a cluster"""

#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pyabc
import time
import os
import tempfile

db_path = ("sqlite://")


# In[2]:


redis_sampler = pyabc.sampler.RedisEvalParallelSampler(host="131.220.224.226", port=6358, look_ahead = False)


eps_list=[5, 3, 2, 1, 0.75, 0.5, 0.33, 0.25, 0.2, 0.15]
#eps = pyabc.ListEpsilon(eps_list)
max_nr_pop=5
eps=pyabc.MedianEpsilon(500, median_multiplier=0.8)
min_eps=0.01
pop_size = 500
noise_factor=0.025
iters=100
sleeptime=0.0
noisefactor = 0.4


# In[3]:


def model(pars):
    theta = pars['theta1']
    sol = theta**2 + 0.4 * np.random.randn()
    
    # Mimic a model with long runtimes for some parameters
    if theta < 0:
        time.sleep(sleeptime)
            
    return {"X_2": sol}

parameter_prior = pyabc.Distribution(theta1=pyabc.RV("uniform", -2, 4))

def distance(simulation, data):
    return abs(data["X_2"] - simulation["X_2"])


# In[4]:


abc = pyabc.ABCSMC(models=model,
        parameter_priors=parameter_prior,
        distance_function=distance,
        population_size=pop_size,
        sampler=redis_sampler,
        eps=eps)

log_file="/home/freck/sampling/log/PreliminaryAcceptances.csv"

histories=[]
runtimes=np.zeros(iters)
for i in range(iters):
    starttime=time.time()
    abc.new(db_path, {"X_2": 1});
    endtime=time.time()
    h = abc.run(minimum_epsilon=min_eps, max_nr_populations=max_nr_pop, log_file	=log_file)
    histories.append(h)
    runtimes[i]=endtime-starttime


# In[5]:


tmax=histories[0].max_t
fig,ax = plt.subplots()

kernel = pyabc.transition.GridSearchCV()

for i in range(iters):
    df,w = histories[i].get_distribution(m=0,t=histories[i].max_t)
    pyabc.visualization.plot_kde_1d(df, w, x='theta1',ax=ax, color='grey', kde = kernel, alpha = 0.2)

plt.savefig("/home/freck/sampling/log/img/UnbalancedDistributions.jpg")


# In[6]:


t1_quantiles_10 = np.zeros(iters)
t1_medians = np.zeros(iters)
t1_quantiles_90 = np.zeros(iters)
t1_means = np.zeros(iters)
t1_stds = np.zeros(iters)


for i in range(iters):
    df,w= histories[i].get_distribution(m=0,t=histories[i].max_t)
    points = df['theta1'].values
    t1_quantiles_10[i]=pyabc.weighted_statistics.weighted_quantile(points, w, alpha=0.1)
    t1_medians[i]=pyabc.weighted_statistics.weighted_quantile(points, w, alpha=0.5)
    t1_quantiles_90[i]=pyabc.weighted_statistics.weighted_quantile(points, w, alpha=0.9)
    t1_means[i]=pyabc.weighted_statistics.weighted_mean(points,w)
    t1_stds[i]=pyabc.weighted_statistics.weighted_std(points,w)
    
    

allnames=[]
allnames.append("t1_means")
allnames.append("t1_stds")
allnames.append("t1_quantiles_10")
allnames.append("t1_medians")
allnames.append("t1_quantiles_90")

alldata = []
for i in range(len(allnames)):
    alldata.append(eval(allnames[i]))


# In[7]:


fig = plt.figure(figsize=(10,10))
nx, ny=5, 1

for i in range(nx*ny):
    ax = fig.add_subplot(nx,ny,i+1)
    ax.hist(alldata[i], bins = 10, range=(alldata[i].min(), alldata[i].max()))
    ax.set_title(allnames[i])
    
fig.tight_layout()
plt.savefig("/home/freck/sampling/log/img/MeansIEUnbalanced.jpg")
