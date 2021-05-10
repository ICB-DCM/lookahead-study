import pyabc
import matplotlib.pyplot as plt
import os
import tempfile
import numpy as np

import seaborn as sns
import pandas as pd


font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 12}
from matplotlib import rc
rc('font', **font)


#Creates the plots for the MJP model (Fig.11&12) 
#based on the databases created in the Tests

#Works in the same way as the ODEBoxplotCreation.py file

#Fix some global variables
nodes = 4

basepath="/p/project/fitmulticell/felipe/scripts/Batch_pyABC/programs/MJP"

pop_sizes = [64, 256, 1024, 4096]

eps_list=[8, 4, 2, 1, 0.7,  0.5, 0.33, 0.25]
eps = pyabc.ListEpsilon(eps_list)

iters_LA=25
iters_ORI=25


# In[3]:

def load_histories(basepath, look_ahead, iters, popsize):
    """
    Reads in all databases as pyABC History objects
    """
    histories = []
    for i in range(iters):
        path=os.path.join(basepath,
                          "results/database", 
                          "MJP"+str(look_ahead)+"N"+str(nodes)+"DB"+str(i)+".db")
        print(path)
        histories.append(pyabc.History("sqlite:///" + path))
    return histories

histories_LA=load_histories(basepath, "LA", iters_LA, pop_sizes[-1])
histories_ORI=load_histories(basepath, "ORI", iters_ORI, pop_sizes[-1])


# In[4]:

col_names=["parameter_means", "parameter_stds", "effective_sample_size"]

def get_data(histories, iters):
    """
    Extracts the relevant data from the history objects
    """
    alldata=np.zeros((iters, len(col_names)))
    for i in range(iters):
        df,w= histories[i].get_distribution(m=0,t=histories[i].max_t)

        points = df['rate'].values
        alldata[i,0]=pyabc.weighted_statistics.weighted_mean(points,w)
        alldata[i,1]=pyabc.weighted_statistics.weighted_std(points,w)

        alldata[i,2]=pyabc.weighted_statistics.effective_sample_size(w)

    return np.array(alldata)

alldata_LA=get_data(histories_LA, iters_LA)
alldata_ORI=get_data(histories_ORI, iters_ORI)


# Creation of a histogram comparing the means, STDs and the effective sample size (Fig. 11)

fig = plt.figure(figsize=(12,4))
nx, ny=1, 3
binranges=[2,2.2,0.25,0.4,0,4096]
titles=["Mean","Standard Deviation","Effective Sample Size"]
for i in range(len(col_names)):
    ax = fig.add_subplot(nx,ny,i+1)
    ax.set_title(titles[i])

    bins=50
    binrange_min=min(alldata_LA[:,i].min(),alldata_ORI[:,i].min())
    binrange_max=max(alldata_LA[:,i].max(),alldata_ORI[:,i].max())
    sns.histplot(alldata_ORI[:,i], ax = ax, bins = bins,
                 binrange = (binranges[2*i],binranges[2*i+1]), common_bins = True,
                 kde=False, label="DYN", color = "blue", alpha = 0.2)
    sns.histplot(alldata_LA[:,i], ax = ax, bins = bins,
                 binrange = (binranges[2*i],binranges[2*i+1]), common_bins = True,
                 kde=False, label="LA", color = "orange", alpha = 0.3)
    
    ax.axvline(x = alldata_LA[:,i].mean(), label="means", color="blue", linestyle="dashed")
    ax.axvline(x = alldata_ORI[:,i].mean(), color="orange", linestyle="dashed")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend()
fig.tight_layout()
plt.savefig(os.path.join(basepath, "img", "MJPMeansIEcomparision.pdf"))


# Plots the posteriors and their development (Figure 12)

iter=3
max_t=histories_LA[iter].max_t
print(max_t)
fig,ax = plt.subplots(1,1, figsize=(8,5))
x_max = 7
for i in range(max_t):
    df,w = histories_LA[iter].get_distribution(m=0,t=i)
    pyabc.visualization.plot_kde_1d(df, w, x='rate',ax=ax, color='black', alpha=(i)/max_t, xmin = 0, xmax=x_max)

df,w = histories_LA[iter].get_distribution(m=0,t=max_t)
pyabc.visualization.plot_kde_1d(df, w, x='rate',ax=ax, color='blue', xmin = 0, xmax=x_max)

df_ori, w_ori = histories_ORI[iter].get_distribution(m=0,t=histories_ORI[1].max_t)
pyabc.visualization.plot_kde_1d(df_ori, w_ori, x='rate',ax=ax, color='orange', xmin = 0, xmax=x_max)

#ax.axvline(x = alldata_LA[:, 0].mean(), color="blue", linestyle="dashed")
#ax.axvline(x = alldata_ORI[:, 0].mean(), color="orange", linestyle="dashed")
ax.set_xlabel("Rate k")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ymin,ymax=ax.get_ylim()
ax.set_ylim(ymin=0,ymax=ymax)
ax.plot(0,0, color='black', alpha=0.5, label= "Generations 1-"+str(max_t)+" - LA")
ax.plot(0,0, color='blue', label='LA Posterior')
ax.plot(0,0, color='orange', label='DYN Posterior')
#ax.plot(0,0, linestyle='dashed', color='black', label="Means")
ax.set_yticks([])
ax.tick_params(axis='both', which='major', labelsize=10)

plt.legend()
fig.tight_layout()
plt.savefig(os.path.join(basepath, "img", "MJPPosteriorDevelopment.pdf"))
