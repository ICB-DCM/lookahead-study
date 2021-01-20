# In[1]:


import pyabc
import matplotlib.rc
import matplotlib.pyplot as plt
import os
import tempfile
import numpy as np

import seaborn as sns
import pandas as pd

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}

matplotlib.rc('font', **font)

# In[2]:


nodes = 8
noisevar = 0.03

sleepvars=[1,1.5,2]

basepath="/p/home/jusers/reck1/juwels/scripts/Batch_pyABC/programs/ODE"

pop_sizes = [64,256,1024,4096]
psize=pop_sizes[-1]


eps_list=[8, 4, 2, 1, 0.7,  0.5, 0.33, 0.25]
eps = pyabc.ListEpsilon(eps_list)

iters_LA=25
iters_ORI=25
iters_STAT=3


# In[3]:



histories_LA=[]
for i in range(iters_LA):
    histories_LA.append(pyabc.History("sqlite:///" +
                                      os.path.join(basepath,
                                                   "results/Var"+str(sleepvars[0])+"/database",
                                                   "DYNLA"+str(nodes)+"_"+str(psize)+"_"+str(i)+".db")))

histories_ORI=[]
for i in range(iters_ORI):
    histories_ORI.append(pyabc.History("sqlite:///" +
                                      os.path.join(basepath,
                                                   "results/Var"+str(sleepvars[0])+"/database",
                                                   "ORI"+str(nodes)+"_"+str(psize)+"_"+str(i)+".db")))

histories_STAT=[]
for i in range(iters_STAT):
    histories_STAT.append(pyabc.History("sqlite:///" +
                                      os.path.join(basepath,
                                                   "results/Var"+str(sleepvars[0])+"/database",
                                                   "STATIC"+str(nodes)+"_"+str(psize)+"_"+str(i)+".db")))


# In[4]:


allnames=[]
allnames.append("theta1_means")
allnames.append("theta1_stds")

allnames.append("effective_sample_size")

allnames.append("theta2_means")
allnames.append("theta2_stds")


# In[5]:


theta1_means_LA = np.zeros(iters_LA)
theta1_stds_LA = np.zeros(iters_LA)

theta2_means_LA = np.zeros(iters_LA)
theta2_stds_LA = np.zeros(iters_LA)

effective_sample_size_LA=np.zeros(iters_LA)

for i in range(iters_LA):
    
    df,w= histories_LA[i].get_distribution(m=0,t=histories_LA[i].max_t)
    
    points = df['theta1'].values
    theta1_means_LA[i]=pyabc.weighted_statistics.weighted_mean(points,w)
    theta1_stds_LA[i]=pyabc.weighted_statistics.weighted_std(points,w)
    
    points = df['theta2'].values
    theta2_means_LA[i]=pyabc.weighted_statistics.weighted_mean(points,w)
    theta2_stds_LA[i]=pyabc.weighted_statistics.weighted_std(points,w)
    
    effective_sample_size_LA[i]=pyabc.weighted_statistics.effective_sample_size(w)

alldata_LA = []
for i in range(len(allnames)):
    alldata_LA.append(eval(allnames[i]+"_LA"))


# In[6]:


theta1_means_ORI = np.zeros(iters_ORI)
theta1_stds_ORI = np.zeros(iters_ORI)

theta2_means_ORI = np.zeros(iters_ORI)
theta2_stds_ORI = np.zeros(iters_ORI)

effective_sample_size_ORI = np.zeros(iters_ORI)

for i in range(iters_ORI):
    
    df,w= histories_ORI[i].get_distribution(m=0,t=histories_ORI[i].max_t)
    
    points = df['theta1'].values
    theta1_means_ORI[i]=pyabc.weighted_statistics.weighted_mean(points,w)
    theta1_stds_ORI[i]=pyabc.weighted_statistics.weighted_std(points,w)
    
    points = df['theta2'].values
    theta2_means_ORI[i]=pyabc.weighted_statistics.weighted_mean(points,w)
    theta2_stds_ORI[i]=pyabc.weighted_statistics.weighted_std(points,w)
    
    effective_sample_size_ORI[i]=pyabc.weighted_statistics.effective_sample_size(w)

alldata_ORI = []
for i in range(len(allnames)):
    alldata_ORI.append(eval(allnames[i]+"_ORI"))


# In[7]:


fig = plt.figure(figsize=(15,8))
nx, ny=2, 3

for i in range(5):
    ax = fig.add_subplot(nx,ny,i+1)
    ax.set_title(allnames[i])

    bins=10
    binrange_min=min(alldata_LA[i].min(),alldata_ORI[i].min())
    binrange_max=max(alldata_LA[i].max(),alldata_ORI[i].max())
    sns.histplot(alldata_LA[i], ax = ax, bins = bins,
                 binrange = (binrange_min, binrange_max), common_bins = True,
                 kde=False, label="Look ahead", color = "orange", alpha = 0.3)
    sns.histplot(alldata_ORI[i], ax = ax, bins = bins,
                 binrange = (binrange_min, binrange_max), common_bins = True,
                 kde=False, label="Original", color = "blue", alpha = 0.2)
    ax.axvline(x = alldata_LA[i].mean(), label="means", color="blue", linestyle="dashed")
    ax.axvline(x = alldata_ORI[i].mean(),color="orange", linestyle="dashed")
        
    ax.legend()
fig.tight_layout()
plt.savefig(os.path.join(basepath, "img", "MeansIEcomparision.pdf"))


# In[8]:


df_LA=pd.DataFrame(np.array(alldata_LA).transpose(), columns=allnames)
df_LA['Look_ahead']=True
df_LA['LNVariance']=sleepvars[0]

df_ORI=pd.DataFrame(np.array(alldata_ORI).transpose(), columns=allnames)
df_ORI['Look_ahead']=False
df_ORI['LNVariance']=sleepvars[0]

if len(sleepvars)>1:
    for k in sleepvars:
        histories_LA_temp=[]
        for i in range(iters_LA):
            histories_LA_temp.append(pyabc.History("sqlite:///" +
                                                  os.path.join(basepath,
                                                               "results/Var"+str(k)+"/database",
                                                               "DYNLA"+str(nodes)+"_"+str(psize)+"_"+str(i)+".db")))
            
        histories_ORI_temp=[]
        for i in range(iters_ORI):
            histories_ORI_temp.append(pyabc.History("sqlite:///" +
                                                  os.path.join(basepath,
                                                               "results/Var"+str(k)+"/database",
                                                               "ORI"+str(nodes)+"_"+str(psize)+"_"+str(i)+".db")))
        
        theta1_means_LA = np.zeros(iters_LA)
        theta1_stds_LA = np.zeros(iters_LA)

        theta2_means_LA = np.zeros(iters_LA)
        theta2_stds_LA = np.zeros(iters_LA)

        effective_sample_size_LA=np.zeros(iters_LA)

        for i in range(iters_LA):

            df,w= histories_LA_temp[i].get_distribution(m=0,t=histories_LA[i].max_t)

            points = df['theta1'].values
            theta1_means_LA[i]=pyabc.weighted_statistics.weighted_mean(points,w)
            theta1_stds_LA[i]=pyabc.weighted_statistics.weighted_std(points,w)

            points = df['theta2'].values
            theta2_means_LA[i]=pyabc.weighted_statistics.weighted_mean(points,w)
            theta2_stds_LA[i]=pyabc.weighted_statistics.weighted_std(points,w)

            effective_sample_size_LA[i]=pyabc.weighted_statistics.effective_sample_size(w)

        alldata_LA_temp = []
        for i in range(len(allnames)):
            alldata_LA_temp.append(eval(allnames[i]+"_LA"))
        
        
        theta1_means_ORI = np.zeros(iters_ORI)
        theta1_stds_ORI = np.zeros(iters_ORI)

        theta2_means_ORI = np.zeros(iters_ORI)
        theta2_stds_ORI = np.zeros(iters_ORI)

        effective_sample_size_ORI = np.zeros(iters_ORI)

        for i in range(iters_ORI):

            df,w= histories_ORI_temp[i].get_distribution(m=0,t=histories_ORI[i].max_t)

            points = df['theta1'].values
            theta1_means_ORI[i]=pyabc.weighted_statistics.weighted_mean(points,w)
            theta1_stds_ORI[i]=pyabc.weighted_statistics.weighted_std(points,w)

            points = df['theta2'].values
            theta2_means_ORI[i]=pyabc.weighted_statistics.weighted_mean(points,w)
            theta2_stds_ORI[i]=pyabc.weighted_statistics.weighted_std(points,w)

            effective_sample_size_ORI[i]=pyabc.weighted_statistics.effective_sample_size(w)

        alldata_ORI_temp = []
        for i in range(len(allnames)):
            alldata_ORI_temp.append(eval(allnames[i]+"_ORI"))
        
        df_LA_temp=pd.DataFrame(np.array(alldata_LA_temp).transpose(), columns=allnames)
        df_LA_temp['Look_ahead']=True
        df_LA_temp['LNVariance']=k
        
        df_LA=pd.concat([df_LA, df_LA_temp])
        
        df_ORI_temp=pd.DataFrame(np.array(alldata_ORI_temp).transpose(), columns=allnames)
        df_ORI_temp['Look_ahead']=False
        df_ORI_temp['LNVariance']=k
        
        df_ORI=pd.concat([df_ORI, df_ORI_temp])

FullDFVar = pd.concat([df_ORI, df_LA])


# In[14]:


df_LA=pd.DataFrame(np.array(alldata_LA).transpose(), columns=allnames)
df_LA['Look_ahead']=True
df_LA['Pop_size']=psize

df_ORI=pd.DataFrame(np.array(alldata_ORI).transpose(), columns=allnames)
df_ORI['Look_ahead']=False
df_ORI['Pop_size']=psize

if len(pop_sizes)>1:
    for k in pop_sizes[:-1]:
        histories_LA_temp=[]
        for i in range(iters_LA):
            histories_LA_temp.append(pyabc.History("sqlite:///" +
                                                  os.path.join(basepath,
                                                               "results/Var"+str(sleepvars[0])+"/database",
                                                               "DYNLA"+str(nodes)+"_"+str(k)+"_"+str(i)+".db")))
            
        histories_ORI_temp=[]
        for i in range(iters_ORI):
            histories_ORI_temp.append(pyabc.History("sqlite:///" +
                                                  os.path.join(basepath,
                                                               "results/Var"+str(sleepvars[0])+"/database",
                                                               "ORI"+str(nodes)+"_"+str(k)+"_"+str(i)+".db")))
        
        theta1_means_LA = np.zeros(iters_LA)
        theta1_stds_LA = np.zeros(iters_LA)

        theta2_means_LA = np.zeros(iters_LA)
        theta2_stds_LA = np.zeros(iters_LA)

        effective_sample_size_LA=np.zeros(iters_LA)

        for i in range(iters_LA):

            df,w= histories_LA_temp[i].get_distribution(m=0,t=histories_LA[i].max_t)

            points = df['theta1'].values
            theta1_means_LA[i]=pyabc.weighted_statistics.weighted_mean(points,w)
            theta1_stds_LA[i]=pyabc.weighted_statistics.weighted_std(points,w)

            points = df['theta2'].values
            theta2_means_LA[i]=pyabc.weighted_statistics.weighted_mean(points,w)
            theta2_stds_LA[i]=pyabc.weighted_statistics.weighted_std(points,w)

            effective_sample_size_LA[i]=pyabc.weighted_statistics.effective_sample_size(w)

        alldata_LA_temp = []
        for i in range(len(allnames)):
            alldata_LA_temp.append(eval(allnames[i]+"_LA"))
        
        
        theta1_means_ORI = np.zeros(iters_ORI)
        theta1_stds_ORI = np.zeros(iters_ORI)

        theta2_means_ORI = np.zeros(iters_ORI)
        theta2_stds_ORI = np.zeros(iters_ORI)

        effective_sample_size_ORI = np.zeros(iters_ORI)

        for i in range(iters_ORI):

            df,w= histories_ORI_temp[i].get_distribution(m=0,t=histories_ORI[i].max_t)

            points = df['theta1'].values
            theta1_means_ORI[i]=pyabc.weighted_statistics.weighted_mean(points,w)
            theta1_stds_ORI[i]=pyabc.weighted_statistics.weighted_std(points,w)

            points = df['theta2'].values
            theta2_means_ORI[i]=pyabc.weighted_statistics.weighted_mean(points,w)
            theta2_stds_ORI[i]=pyabc.weighted_statistics.weighted_std(points,w)

            effective_sample_size_ORI[i]=pyabc.weighted_statistics.effective_sample_size(w)

        alldata_ORI_temp = []
        for i in range(len(allnames)):
            alldata_ORI_temp.append(eval(allnames[i]+"_ORI"))
        
        df_LA_temp=pd.DataFrame(np.array(alldata_LA_temp).transpose(), columns=allnames)
        df_LA_temp['Look_ahead']=True
        df_LA_temp['Pop_size']=k
        
        df_LA=pd.concat([df_LA, df_LA_temp])
        
        df_ORI_temp=pd.DataFrame(np.array(alldata_ORI_temp).transpose(), columns=allnames)
        df_ORI_temp['Look_ahead']=False
        df_ORI_temp['Pop_size']=k
        
        df_ORI=pd.concat([df_ORI, df_ORI_temp])

FullDFPop = pd.concat([df_ORI, df_LA])


# In[17]:


sns.set_theme(style="whitegrid")

plotnr=3

fig,axes=plt.subplots(1, plotnr, figsize=(plotnr*5,5))

sns.boxplot(x="Pop_size",y=allnames[0], hue="Look_ahead", data=FullDFPop, ax=axes[0])
    
axes[0].set_xlabel("Population size")


sns.boxplot(x="LNVariance",y=allnames[0], hue="Look_ahead", data=FullDFVar, ax=axes[1])
    
axes[1].set_xlabel("LogNormal Variance")


sns.boxplot(x="LNVariance",y=allnames[2], hue="Look_ahead", data=FullDFVar, ax=axes[2])
    
axes[2].set_xlabel("LogNormal Variance")

fig.tight_layout()

plt.savefig(os.path.join(basepath, "img", "BoxPlots.pdf"))


# In[18]:


max_t=histories_LA[-1].max_t

fig,ax = plt.subplots()

for i in range(max_t):
    df,w = histories_LA[0].get_distribution(m=0,t=i)
    pyabc.visualization.plot_kde_1d(df, w, x='theta1',ax=ax, color='black', alpha=(i+1)/max_t, label='Gen'+str(i+1))

df,w = histories_LA[0].get_distribution(m=0,t=max_t)
pyabc.visualization.plot_kde_1d(df, w, x='theta1',ax=ax, color='blue', alpha=0.8, label='LA Posterior')

df_ori, w_ori = histories_ORI[1].get_distribution(m=0,t=histories_ORI[1].max_t)
pyabc.visualization.plot_kde_1d(df_ori, w_ori, x='theta1',ax=ax, color='orange', alpha=0.8, label='Ori Posterior')

ax.axvline(x = alldata_LA[0].mean(), label="LA mean", color="blue", linestyle="dashed")
ax.axvline(x = alldata_ORI[0].mean(), label='Ori mean', color="orange", linestyle="dotted")

ax.set_yticks([])
ax.set_xlim(xmin=0, xmax=0.3)
plt.legend()

plt.savefig(os.path.join(basepath, "img", "PosteriorDevelopment.pdf"))


# In[21]:


fig,ax = plt.subplots()

for i in range(3):
    df,w = histories_LA[i].get_distribution(m=0,t=max_t)
    pyabc.visualization.plot_kde_1d(df, w, x='theta1',ax=ax, color='blue', alpha=0.5)

    df_ori, w_ori = histories_ORI[i].get_distribution(m=0,t=histories_ORI[1].max_t)
    pyabc.visualization.plot_kde_1d(df_ori, w_ori, x='theta1',ax=ax, color='orangered', alpha=0.5)

ax.plot(0,0, label = "LA", color = 'blue')
ax.plot(0,0, label = "Ori", color = 'orangered')


ax.axvline(x = alldata_LA[0].mean(), label="LA mean", color="blue", linestyle="dashed")
ax.axvline(x = alldata_ORI[0].mean(), label='Ori mean', color="orange", linestyle="dotted")

plt.legend()

plt.savefig(os.path.join(basepath, "img", "Posteriors.pdf"))

