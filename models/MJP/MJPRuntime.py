import pyabc
import matplotlib.pyplot as plt
import os
import tempfile
import numpy as np
import scipy as sp
import time

path="/p/project/fitmulticell/felipe/scripts/Batch_pyABC/programs/MJP"


import argparse

parser = argparse.ArgumentParser(description='Parse necessary arguments')
parser.add_argument('-pt', '--port', type=str, default="50004",
                    help='Which port should be used?')
parser.add_argument('-ip', '--ip', type=str,
                    help='Dynamically passed - BW: Login Node 3')
parser.add_argument('-nd','--nodes', type=int, default=8, help='How many nodes are used')

args = parser.parse_args()

port = args.port
host = args.ip
nodes=args.nodes

# Set constants

eps_list=[10, 5, 3, 2, 1, 0.85, 0.75, 0.66]
eps=pyabc.ListEpsilon(eps_list)
max_nr_pop=len(eps_list)
min_eps=0.66
pop_sizes = [64, 256, 1024, 4096]
iters_PPP = 25
iters_ori = 25
resultfilepath = os.path.join(path, "results/MJPruntimeresultsN"+str(nodes))
resultfile = open(resultfilepath, "w")
resultfile.write("Pop size, Look_ahead, Repetitions, Runtime Expectation, Runtime Variance, total Walltime\n")
resultfile.close()


# Define the model


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


# Run inference in Look-Ahead mode


redis_sampler = pyabc.sampler.RedisEvalParallelSampler(host=host,
                                                       port=port,
                                                       look_ahead = True,
                                                       look_ahead_delay_evaluations=False)

histories=[]
runtimes=np.zeros(iters_PPP)

for pop_size in pop_sizes:
    totalstarttime =time.time()
    for i in range(iters_PPP):
        starttime=time.time()
        
        db_path = "sqlite:///" + os.path.join(path,
                                              "results/database",
                                              "MJPLAdatabase"+str(i))
        
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


# Run in Standard mode

redis_sampler = pyabc.sampler.RedisEvalParallelSampler(host=host,
                                                       port=port,
                                                       look_ahead = False)

histories_ori=[]
runtimes_ori=np.zeros(iters_ori)

for pop_size in pop_sizes:
    totalstarttime =time.time()
    for i in range(iters_ori):
        starttime=time.time()
        
        db_path_ori = "sqlite:///" + os.path.join(path,
                                              "results/database",
                                              "MJPORIdatabase"+str(i))
        
        abc = pyabc.ABCSMC(models = Model1(),
                   parameter_priors = prior,
                   distance_function = distance,
                   population_size = pop_size,
                   sampler = redis_sampler,
                   eps = eps)

        abc.new(db_path_ori, observations)
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

