import pyabc
import matplotlib.pyplot as plt
import os
import tempfile
import numpy as np
import scipy as sp
import time

import argparse

parser = argparse.ArgumentParser(description='Parse necessary arguments')
parser.add_argument('-pt', '--port', type=str, default="50004",
                    help='Which port should be use?')
parser.add_argument('-ip', '--ip', type=str,
                    help='Dynamically passed - BW: Login Node 3')
parser.add_argument('-nd','--nodes', type=int, default=8, help='How many nodes are used')
parser.add_argument('-sv', '--var', type=int, default=1, help='Lognormal Variance value')
args = parser.parse_args()

port=args.port
host=args.ip
nodes=args.nodes
sleepvar=args.var


# Constants

noisevar = 0.03

basepath="/p/home/jusers/reck1/juwels/scripts/Batch_pyABC/programs/ODE"

pop_sizes = [64, 256, 1024, 4096]
eps_list=[8, 4, 2, 1, 0.7,  0.5, 0.33, 0.25]
eps = pyabc.ListEpsilon(eps_list)

iters_PPP=min(5*nodes+20, 50)
iters_ori=iters_PPP
iters_stat=min(2+nodes+10, 25)

resultfilepath=os.path.join(basepath, "results/Var"+str(sleepvar), "sleeptimeresults"+str(nodes)+".txt")
resultfile = open(resultfilepath, "w")
resultfile.write("Pop size, Look_ahead, Repetitions, Runtime Expectation, Runtime Variance, total Walltime\n")
resultfile.close()


# Model

theta1_true, theta2_true = np.exp([-2.5, -2])

measurement_data = np.array([0.0244, 0.0842, 0.1208,
                             0.1724, 0.2315, 0.2634,
                             0.2831, 0.3084, 0.3079,
                             0.3097, 0.3324])

measurement_times = np.arange(len(measurement_data))

init = np.array([1, 0])

def f(y, t0, theta1, theta2):
    x1, x2 = y
    dx1 = - theta1 * x1 + theta2 * x2
    dx2 =   theta1 * x1 - theta2 * x2
    return dx1, dx2


def model(pars):
    sol = sp.integrate.odeint(
             f, init, measurement_times,
             args=(pars["theta1"],pars["theta2"]))
    pause_ms = pow(10,-1*sleepvar)*np.random.lognormal(0,sleepvar)
    noise = np.random.normal(1,noisevar,len(sol[:,1]))
    noisysol = sol[:,1]*noise
    time.sleep(pause_ms)
    return {"X_2": noisysol} 

true_trajectory = model({"theta1": theta1_true,
                         "theta2": theta2_true})["X_2"]

def distance(simulation, data):
    return np.absolute(data["X_2"] - simulation["X_2"]).sum()
   
parameter_prior = pyabc.Distribution(theta1=pyabc.RV("uniform", 0, 1),
                                     theta2=pyabc.RV("uniform", 0, 1))
parameter_prior.get_parameter_names()


# OriSampling

redis_sampler = pyabc.sampler.RedisEvalParallelSampler(host=host,
                                                 port=port,
                                                 look_ahead = False)

for psize in pop_sizes:

        means1 = np.zeros(iters_ori)
        means2 = np.zeros(iters_ori)
        runtimes_original = np.zeros(iters_ori)
        totalstarttime =time.time()
        for i in range(0,iters_ori):

            abc = pyabc.ABCSMC(models=model,
                     parameter_priors=parameter_prior,
                     distance_function=distance,
                     population_size=psize,
                     sampler=redis_sampler,
                     eps=eps)
            
            db_path_ori = "sqlite:///" + os.path.join(basepath,
                                                      "results/Var"+str(sleepvar)+"/database",
                                                      "ORI"+str(nodes)+"_"+str(psize)+"_"+str(i)+".db")
            
            abc.new(db_path_ori, {"X_2": measurement_data});
            
            start = time.time()
            h = abc.run(minimum_epsilon=0.1, max_nr_populations=len(eps_list))
            end = time.time()
            
            kekse = h.get_distribution(m=0,t=h.max_t)
            npkekse1 = kekse[0][['theta1']].to_numpy()
            npkekse2 = kekse[0][['theta2']].to_numpy()
            means1[i] = npkekse1.mean()
            means2[i] = npkekse2.mean()
            runtimes_original[i] = end-start

        totalendtime = time.time()
        walltime_original = totalendtime-totalstarttime

        resultfile = open(resultfilepath, "a")
        resultfile.write(str(psize)+", ")
        resultfile.write("Ori, ")
        resultfile.write(str(iters_ori)+", ")
        resultfile.write(str(runtimes_original.mean())+", ")
        resultfile.write(str(runtimes_original.var())+", ")
        resultfile.write(str(walltime_original)+"\n")
        resultfile.close()


# LASampling

redis_sampler = pyabc.sampler.RedisEvalParallelSampler(host=host,
                                                 port=port,
                                                 look_ahead = True,
                                                 look_ahead_delay_evaluation=True)

for psize in pop_sizes:

        means1 = np.zeros(iters_PPP)
        means2 = np.zeros(iters_PPP)
        runtimes = np.zeros(iters_PPP)
        totalstarttime =time.time()

        for i in range(0,iters_PPP):

            abc = pyabc.ABCSMC(models=model,
                           parameter_priors=parameter_prior,
                           distance_function=distance,
                           population_size=psize,
                           sampler=redis_sampler,
                           eps=eps) 
            
            db_path = "sqlite:///" + os.path.join(basepath,
                                                  "results/Var"+str(sleepvar)+"/database",
                                                  "DYNLA"+str(nodes)+"_"+str(psize)+"_"+str(i)+".db")
            
            abc.new(db_path, {"X_2": measurement_data});
            
            start = time.time()
            h = abc.run(minimum_epsilon=0.1, max_nr_populations=len(eps_list))
            end = time.time()
            
            kekse = h.get_distribution(m=0,t=h.max_t)
            npkekse1 = kekse[0][['theta1']].to_numpy()
            npkekse2 = kekse[0][['theta2']].to_numpy()
            means1[i] = npkekse1.mean()
            means2[i] = npkekse2.mean()
            runtimes[i] = end-start

        totalendtime = time.time()
        walltime = totalendtime-totalstarttime

        resultfile = open(resultfilepath, "a")
        resultfile.write(str(psize)+", ")
        resultfile.write("PPP, ")
        resultfile.write(str(iters_PPP)+", ")
        resultfile.write(str(runtimes.mean())+", ")
        resultfile.write(str(runtimes.var())+", ")
        resultfile.write(str(walltime)+"\n")
        resultfile.close()


# Static Sampling

redis_sampler = pyabc.sampler.RedisStaticSampler(host=host, port=port)

for psize in pop_sizes:

        means1 = np.zeros(iters_stat)
        means2 = np.zeros(iters_stat)
        runtimes_static = np.zeros(iters_stat)
        totalstarttime =time.time()

        for i in range(0,iters_stat):

            abc = pyabc.ABCSMC(models=model,
                     parameter_priors=parameter_prior,
                     distance_function=distance,
                     population_size=psize,
                     sampler=redis_sampler,
                     eps=eps)
            
            db_path_stat = "sqlite:///" + os.path.join(basepath,
                                                       "results/Var"+str(sleepvar)+"/database",
                                                       "STATIC"+str(nodes)+"_"+str(psize)+"_"+str(i)+".db")

            abc.new(db_path_stat, {"X_2": measurement_data});
            
            start = time.time()
            h = abc.run(minimum_epsilon=0.1, max_nr_populations=len(eps_list))
            end = time.time()
           
            kekse = h.get_distribution(m=0,t=h.max_t)
            npkekse1 = kekse[0][['theta1']].to_numpy()
            npkekse2 = kekse[0][['theta2']].to_numpy()
            means1[i] = npkekse1.mean()
            means2[i] = npkekse2.mean()
            runtimes_static[i] = end-start

        totalendtime = time.time()
        walltime_static = totalendtime-totalstarttime

        resultfile = open(resultfilepath, "a")
        resultfile.write(str(psize)+", ")
        resultfile.write("Stat, ")
        resultfile.write(str(iters_stat)+", ")
        resultfile.write(str(runtimes_static.mean())+", ")
        resultfile.write(str(runtimes_static.var())+", ")
        resultfile.write(str(walltime_static)+"\n")
        resultfile.close()

