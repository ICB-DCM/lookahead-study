#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pyabc
import matplotlib.pyplot as plt
import os
import tempfile
import numpy as np
import scipy as sp
import time


db_path = ("sqlite://")
logfile_path="/home/felipe/testresults/MJPAcceptanceRates.csv"


# In[ ]:


port=6358

eps_list=[10, 5, 3, 2, 1, 0.8, 0.66]
eps=pyabc.ListEpsilon(eps_list)
#eps=pyabc.MedianEpsilon(500, median_multiplier=0.8)

max_nr_pop=len(eps_list)
min_eps=0.66
pop_sizes = [32, 64, 128, 256, 512, 1024, 2048]
iters_PPP = 15
iters_ori = 10
resultfilepath = "/home/freck/sampling/log/MJPruntimeresults.txt"
resultfile = open(resultfilepath, "w")
resultfile.write("Pop size, Look_ahead, Repetitions, Runtime Expectation, Runtime Variance, total Walltime\n")
resultfile.close()


# In[ ]:


def h(x, pre, c):
    return (x**pre).prod(1) * c

def gillespie(x, c, pre, post, max_t):
    """
    Gillespie simulation
    
    Parameters
    ----------
    
    x: 1D array of size n_species
        The initial numbers.
    
    c: 1D array of size n_reactions
        The reaction rates.
    
    pre: array of size n_reactions x n_species
        What is to be consumed.
    
    post: array of size n_reactions x n_species
        What is to be produced
    
    max_t: int
        Timulate up to time max_t
        
    Returns
    -------
    t, X: 1d array, 2d array
        t: The time points.
        X: The history of the species.
           ``X.shape == (t.size, x.size)``
    
    """
    t = 0
    t_store = [t]
    x_store = [x.copy()]
    S = post - pre

    while t < max_t:
        h_vec = h(x, pre, c)
        h0 = h_vec.sum()
        if h0 == 0:
            break
        delta_t = np.random.exponential(1 / h0)
        # no reaction can occur any more
        if not np.isfinite(delta_t):
            t_store.append(max_t)
            x_store.append(x)
            break
        reaction = np.random.choice(c.size, p=h_vec/h0)
        t = t + delta_t
        x = x + S[reaction]
        
        t_store.append(t)
        x_store.append(x)

    return np.array(t_store), np.array(x_store)


# In[ ]:


MAX_T = 0.1

class Model1:
    __name__ = "Model 1"
    x0 = np.array([40, 3])   # Initial molecule numbers
    pre = np.array([[1, 1]], dtype=int)
    post = np.array([[0, 2]])
    
    
    def __call__(self, par):
        t, X = gillespie(self.x0,
                         np.array([float(par["rate"])]),
                         self.pre, self.post,
                         MAX_T)
        return {"t": t, "X" : X}
    
    
true_rate = 2.3
observations = Model1()({"rate": true_rate})

N_TEST_TIMES = 20

t_test_times = np.linspace(0, MAX_T, N_TEST_TIMES)

def distance(x, y):
    xt_ind = np.searchsorted(x["t"], t_test_times) - 1
    yt_ind = np.searchsorted(y["t"], t_test_times) - 1
    error = (np.absolute(x["X"][:,1][xt_ind]
                        - y["X"][:,1][yt_ind]).sum()
             / t_test_times.size)
    return error

prior = pyabc.Distribution(rate=pyabc.RV("uniform", 0, 100))


# In[ ]:


redis_sampler = pyabc.sampler.RedisEvalParallelSampler(host="131.220.224.226",
                                                       port=port,
                                                       look_ahead = True, 
						       look_ahead_delay_evaluation = True)

histories=[]
runtimes=np.zeros(iters_PPP)

for pop_size in pop_sizes:
    totalstarttime =time.time()
    for i in range(iters_PPP):
        starttime=time.time()

        abc = pyabc.ABCSMC(models = Model1(),
                   parameter_priors = prior,
                   distance_function = distance,
                   population_size = pop_size,
                   sampler = redis_sampler,
                   eps = eps)

        abc.new(db_path, observations)
        history = abc.run(minimum_epsilon=min_eps, max_nr_populations=max_nr_pop)

        endtime=time.time()

        #histories.append(history)
        runtimes[i]=endtime-starttime
        
    totalendtime = time.time()
    walltime = totalendtime-totalstarttime

    resultfile = open(resultfilepath, "a")
    resultfile.write(str(pop_size)+", ")
    resultfile.write("True, ")
    resultfile.write(str(iters_PPP)+", ")
    resultfile.write(str(runtimes.mean())+", ")
    resultfile.write(str(runtimes.var())+", ")
    resultfile.write(str(walltime)+"\n")
    resultfile.close()


# In[ ]:

redis_sampler = pyabc.sampler.RedisEvalParallelSampler(host="131.220.224.226",
                                                       port=port,
                                                       look_ahead = False)

histories_ori=[]
runtimes_ori=np.zeros(iters_ori)

for pop_size in pop_sizes:
    totalstarttime =time.time()
    for i in range(iters_ori):
        starttime=time.time()

        abc = pyabc.ABCSMC(models = Model1(),
                   parameter_priors = prior,
                   distance_function = distance,
                   population_size = pop_size,
                   sampler = redis_sampler,
                   eps = eps)

        abc.new(db_path, observations)
        history = abc.run(minimum_epsilon=min_eps, max_nr_populations=max_nr_pop)

        endtime=time.time()

        #histories_ori.append(history)
        runtimes_ori[i]=endtime-starttime
    
    totalendtime = time.time()
    walltime = totalendtime-totalstarttime

    resultfile = open(resultfilepath, "a")
    resultfile.write(str(pop_size)+", ")
    resultfile.write("False, ")
    resultfile.write(str(iters_ori)+", ")
    resultfile.write(str(runtimes_ori.mean())+", ")
    resultfile.write(str(runtimes_ori.var())+", ")
    resultfile.write(str(walltime)+"\n")
    resultfile.close()

