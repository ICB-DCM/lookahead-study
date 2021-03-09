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


# In[2]:


nodes = 8
noisevar = 0.03

sleepvars=[0.5, 1.0, 1.5, 2.0]

basepath="/p/project/fitmulticell/felipe/scripts/Batch_pyABC/programs/ODE"

pop_sizes = [64, 256, 1024, 4096]

eps_list=[8, 4, 2, 1, 0.7,  0.5, 0.33, 0.25]
eps = pyabc.ListEpsilon(eps_list)

iters_LA=[10,10,10,10]
iters_ORI=[10,10,10,10]
iters_STAT=[3,3,1]


# In[3]:

def load_histories(basepath, look_ahead, iters, sleepvar, popsize):
    histories = []
    for i in range(iters):
        print(sleepvar, look_ahead, popsize, i)
        path=os.path.join(basepath,
                          "results/Var"+str(sleepvar)+"/database", 
                          str(look_ahead)+str(nodes)+"_"+str(popsize)+"_"+str(i)+".db")
        histories.append(pyabc.History("sqlite:///" + path))
    return histories

histories_LA=load_histories(basepath, "DYNLA", iters_LA[0], sleepvars[1], pop_sizes[-1])
histories_ORI=load_histories(basepath, "ORI", iters_ORI[0], sleepvars[1], pop_sizes[-1])
histories_STATIC=load_histories(basepath, "STATIC", iters_STAT[0], sleepvars[1], pop_sizes[-1])


# In[4]:


col_names=["theta1_means", "theta1_stds", "effective_sample_size", "theta2_means", "theta2_stds"]

def get_data(histories, iters):
    alldata=np.zeros((iters, len(col_names)))
    for i in range(iters):
        df,w= histories[i].get_distribution(m=0,t=histories[i].max_t)

        points = df['theta1'].values
        alldata[i,0]=pyabc.weighted_statistics.weighted_mean(points,w)
        alldata[i,1]=pyabc.weighted_statistics.weighted_std(points,w)

        points = df['theta2'].values
        alldata[i,3]=pyabc.weighted_statistics.weighted_mean(points,w)
        alldata[i,4]=pyabc.weighted_statistics.weighted_std(points,w)

        alldata[i,2]=pyabc.weighted_statistics.effective_sample_size(w)

    return np.array(alldata)

alldata_LA=get_data(histories_LA, iters_LA[0])
alldata_ORI=get_data(histories_ORI, iters_ORI[0])


# In[5]:




fig = plt.figure(figsize=(15,8))
nx, ny=2, 3

for i in range(len(col_names)):
    ax = fig.add_subplot(nx,ny,i+1)
    ax.set_title(col_names[i])

    bins=10
    binrange_min=min(alldata_LA[:,i].min(),alldata_ORI[:,i].min())
    binrange_max=max(alldata_LA[:,i].max(),alldata_ORI[:,i].max())
    sns.histplot(alldata_LA[:,i], ax = ax, bins = bins,
                 binrange = (binrange_min, binrange_max), common_bins = True,
                 kde=False, label="Look ahead", color = "orange", alpha = 0.3)
    sns.histplot(alldata_ORI[:,i], ax = ax, bins = bins,
                 binrange = (binrange_min, binrange_max), common_bins = True,
                 kde=False, label="Original", color = "blue", alpha = 0.2)
    ax.axvline(x = alldata_LA[:,i].mean(), label="means", color="blue", linestyle="dashed")
    ax.axvline(x = alldata_ORI[:,i].mean(),color="orange", linestyle="dashed")
        
    ax.legend()
fig.tight_layout()
plt.savefig(os.path.join(basepath, "img", "ODEMeansIEcomparision.pdf"))


# In[6]:



def data_to_frame(basepath, look_aheads, alliters, sleepvars, pop_sizes):
    allframes = []
    for i in range(len(alliters)):
        for k in pop_sizes:
            for l in range(len(sleepvars)):

                keks = load_histories(basepath, look_aheads[i], alliters[i][l] , sleepvars[l], k)
                keksarray = get_data(keks, alliters[i][l])
                keksframe = pd.DataFrame(keksarray, columns=col_names)
                if look_aheads[i]=='LA' or look_aheads[i]=='DYNLA' or look_aheads[i]=='PPP':
                    keksframe['Scheduling'] = "LA"
                elif look_aheads[i]=='Ori' or look_aheads[i]=='ORI':
                    keksframe['Scheduling'] = "DYN"
                else: keksframe['Scheduling'] = look_aheads[i]
                keksframe['LNVariance'] = sleepvars[l]
                keksframe['Pop_size'] = k
                allframes.append(keksframe)

    df = pd.concat(allframes)
    return df

FullDFVar=data_to_frame(basepath, ["DYNLA","ORI"], [iters_LA, iters_ORI], sleepvars, [pop_sizes[-1]])
FullDFVar_N256=data_to_frame(basepath, ["DYNLA","ORI"], [iters_LA, iters_ORI], sleepvars, [pop_sizes[0]])
FullDFPop=data_to_frame(basepath, ["DYNLA","ORI"], [iters_LA, iters_ORI], [1.0], pop_sizes)


# In[9]:


fig,axes=plt.subplots(2, 2, figsize=(10,10))

ax = axes[0][0]
sns.boxplot(x="Pop_size",y=col_names[0], hue="Scheduling",
            data=FullDFPop, ax=axes[0][0], hue_order=["DYN","LA"])
#l1=plt.bar([0],[0],0 ,label="A", color='white')
#axes[0][0].set_title("A: Means; z=1")
ymin,ymax=ax.get_ylim()
xmin,xmax=ax.get_xlim()
ax.text(xmin+0.02*(xmax-xmin),ymin+0.95*(ymax-ymin),'A')
ax.set_xlabel("Population size")
ax.set_ylabel("Mean")
ax.tick_params(axis='both', which='major', labelsize=10)
ax.legend([],[], frameon=False)#, loc='upper left', bbox_to_anchor = (-0.05,1))


ax = axes[0][1]
sns.boxplot(x="Pop_size",y=col_names[1], hue="Scheduling",
            data=FullDFPop, ax=ax, hue_order=["DYN","LA"])
ymin,ymax=ax.get_ylim()
xmin,xmax=ax.get_xlim()
ax.text(xmin+0.02*(xmax-xmin),ymin+0.95*(ymax-ymin),'B')
#axes[0][1].set_title("B: Standard deviation; z=1")
axes[0][1].set_xlabel("Population size")
axes[0][1].set_ylabel("Standard deviations")
axes[0][1].tick_params(axis='both', which='major', labelsize=10)
ax.legend([],[], frameon=False)


ax = axes[1][0]
sns.boxplot(x="LNVariance",y=col_names[0], hue="Scheduling",
            data=FullDFVar_N256, ax=axes[1][0], hue_order=["DYN","LA"])
ymin,ymax=ax.get_ylim()
xmin,xmax=ax.get_xlim()
ax.text(xmin+0.02*(xmax-xmin),ymin+0.95*(ymax-ymin),'C')
#axes[1][0].set_title("C: Means; N=4096")
axes[1][0].set_xlabel("LogNormal Variance")
axes[1][0].set_ylabel("Mean")
axes[1][0].tick_params(axis='both', which='major', labelsize=10)
ax.legend([],[], frameon=False)


ax = axes[1][1]
sns.boxplot(x="LNVariance",y=col_names[2], hue="Scheduling",
            data=FullDFVar_N256, ax=axes[1][1], hue_order=["DYN","LA"])
ymin,ymax=ax.get_ylim()
xmin,xmax=ax.get_xlim()
ax.text(xmin+0.02*(xmax-xmin),ymin+0.95*(ymax-ymin),'D')
#axes[1][1].set_title("D: Effective sample sizes; N=256")
axes[1][1].set_xlabel("STD of delay on log-scale")
axes[1][1].set_ylabel("Effective sample size")
axes[1][1].tick_params(axis='both', which='major', labelsize=10)
ax.legend([],[], frameon=False)


legend, axi = plt.subplots()
ax = axi
l1=ax.bar([xmax],[0.5*(ymax-ymin)],0, bottom = ymin)
l2=ax.bar([xmax],[0.5*(ymax-ymin)],0, bottom = ymin)
labels=["DYN","LA"]

axes[0][1].legend([l1,l2],labels, loc='lower right', bbox_to_anchor = (1.4,0.8), frameon=False)
fig.tight_layout()

fig.savefig(os.path.join(basepath, "img", "ODEBoxPlots.pdf"))


# In[10]:


max_t=histories_LA[-1].max_t

fig,ax = plt.subplots(1,1, figsize=(5,5))

for i in range(max_t):
    df,w = histories_LA[0].get_distribution(m=0,t=i)
    pyabc.visualization.plot_kde_1d(df, w, x='theta1',ax=ax, color='black', alpha=(i)/max_t)

df,w = histories_LA[0].get_distribution(m=0,t=max_t)
pyabc.visualization.plot_kde_1d(df, w, x='theta1',ax=ax, color='blue')

df_ori, w_ori = histories_ORI[1].get_distribution(m=0,t=histories_ORI[1].max_t)
pyabc.visualization.plot_kde_1d(df_ori, w_ori, x='theta1',ax=ax, color='orange')

ax.axvline(x = alldata_LA[:, 0].mean(), color="blue", linestyle="dashed")
ax.axvline(x = alldata_ORI[:, 0].mean(), color="orange", linestyle="dashed")


ax.plot(0,0, color='black', alpha=0.5, label= "Generations 1-"+str(max_t))
ax.plot(0,0, color='blue', label='LA Posterior')
ax.plot(0,0, color='orange', label='DYN Posterior')
ax.plot(0,0, linestyle='dashed', color='black', label="Means")
ax.set_yticks([])
ax.set_xlim(xmin=0, xmax=0.3)
ax.tick_params(axis='both', which='major', labelsize=10)

plt.legend()

plt.savefig(os.path.join(basepath, "img", "ODEPosteriorDevelopment.pdf"))


# In[11]:

def get_acceptance_rates(basepath, nodes, pop_sizes, sleepvar, iters):
    acceptance_rates_array=[]
    for psize in pop_sizes:
        for i in range(iters):
            path = os.path.join(basepath, "results/Var"+str(sleepvar)+"/logfiles", "Logs"+str(nodes)+"_"+str(psize)+"_"+str(i)+".csv")
            stat_df = pd.read_csv(path)[1:]
                                
            last_gen = stat_df.iloc[-1,:]
            last_gen["sleepvar"]=sleepvar
            last_gen["Pop_size"]=psize
            last_gen["LA fraction"]=min(last_gen['n_lookahead_accepted']/psize,1)
            last_gen_df = last_gen.to_frame().T
            acceptance_rates_array.append(last_gen_df)
    acceptance_rates_df = pd.concat(acceptance_rates_array)
    return acceptance_rates_df

all_acceptances_list=[]
for sleep_i in range(len(sleepvars)):
    acceptances_df = get_acceptance_rates(basepath, nodes, pop_sizes, sleepvars[sleep_i], iters_LA[sleep_i])
    all_acceptances_list.append(acceptances_df)
    
all_acceptances_df = pd.concat(all_acceptances_list)

#In[12]

fig,axes=plt.subplots(2,1,figsize=(5,10))

ax = axes[0]

sns.boxplot(x="Pop_size",y="LA fraction", data=all_acceptances_df[all_acceptances_df["sleepvar"]==1], ax=ax)
ax.text(xmin+0.02*(xmax-xmin),ymin+0.95*(ymax-ymin),'E')
ax.set_title(None)
ax.set_xlabel("Population size")
ax.set_ylabel("Samples from preliminary")

from matplotlib.ticker import PercentFormatter
ax.yaxis.set_major_formatter(PercentFormatter(1))
ax.tick_params(axis='both', which='major', labelsize=10)

ax = axes[1]

sns.boxplot(x="sleepvar",y="LA fraction", data=all_acceptances_df[all_acceptances_df["Pop_size"]==pop_sizes[-1]], ax=ax)
ax.text(xmin+0.02*(xmax-xmin),ymin+0.95*(ymax-ymin),'E')
ax.set_title(None)
ax.set_xlabel("STD of delay on log-scale")
ax.set_ylabel("Samples from preliminary")

from matplotlib.ticker import PercentFormatter
ax.yaxis.set_major_formatter(PercentFormatter(1))
ax.tick_params(axis='both', which='major', labelsize=10)

fig.tight_layout()
plt.savefig(os.path.join(basepath, "img", "ODEFinalPrelFraction.pdf"))
