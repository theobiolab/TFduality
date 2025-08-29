import numpy as np
import sys,os, re
import pandas as pd
from functools import partial
import itertools
import json, glob
repobase="/users/romartinez/romartinez/repos/MPRAfitting"
#repobase="/Users/rosamartinezcorral/OneDrive - CRG - Centre de Regulacio Genomica/papers/github_repos/MPRAfitting/"
sys.path.append(os.path.join(repobase,"bin"))
from functions_fractional_occupancy import fi_1sites, fi_2sites, fi_3sites, fi_4sites, fi_5sites, fi_6sites
from biophysmodel_utils import logfoldchange_from_avbinding_ratesdirectly_Hill, get_S_max
from utils_fit import run_genetic, get_parsetnames
from utils_fit import dataframe_from_bestpars
sys.path.append("../")
from error import newerror

fifunctions={1:fi_1sites, 2:fi_2sites, 3:fi_3sites,4:fi_4sites,5:fi_5sites,6:fi_6sites}



jid=int(sys.argv[1])-1 #launch 400 jobs
outfolder="2025_01_22_fit1"
outfname="./%s/best_%d.json"%(outfolder,jid)

affinities_tofit=["0.10","0.25","0.50","0.75","0.89","1.00"]
wlist=[5,1,1,5,1,5]
rates_affected_l=[[0,1,2,3],[1,2,3]]
relrange_l=[True,False]
seedlist=[np.arange(0,1000)[10*i:10*i+10] for i in range(100)]
combination=list(itertools.product(rates_affected_l,relrange_l,seedlist))[jid]
rates_affected,relrange,seeds=combination
print("seeds are", seeds)

TF="Sp1"
meandata=pd.read_csv("../data_for_fit_seq_to_loess.csv")
meandata["afcat"]=["%2.2f"%x for x in meandata["afcat"].values]
meandata_fit=meandata[meandata["afcat"].isin(affinities_tofit)]
afcat_idxs=[list(np.where(meandata_fit["afcat"]==i)[0]) for i in meandata_fit["afcat"].unique()]

#if not os.path.isdir(outfolder):
#    os.mkdir(outfolder)

# Genetic Algorithm parameters:
POPULATION_SIZE = 300 #300
MAX_GENERATIONS = 75 #75

MUTSIGMA=0.25 #sigma
INDPB=0.5

P_CROSSOVER = 0.25  # probability for crossover
P_MUTATION = 0.5  # probability for mutating an individ
ETACXBIN=0.25


S_max=get_S_max(TF,pwmpath=os.path.join(repobase,"PWMs"))    


parsetnames,idxs=get_parsetnames(rates_affected=rates_affected,nconditions=1)
parsetnames=parsetnames[:-1] #remove the last one, which is the concentration

idx1,idx2,idx3,idxs_conditions=idxs #we will not use idxs_conditions here


#get bounds
bounds=[]
for name in parsetnames:
   
    if name=="lambda_K":
        bounds.append([-4,4])
    elif name.startswith("n"):
        bounds.append([0,3])
    else:
        bounds.append([-7,7])


args_fitnessfunction={"logfcfunction":logfoldchange_from_avbinding_ratesdirectly_Hill,"fifunctions":fifunctions,
      "idx1":idx1, "idx2":idx2, "idx3":idx3,"S_max":S_max,"coop":False,"rates_affected":rates_affected,
      "data":meandata_fit,"afcat_idxs":afcat_idxs,
     "outerror":"inverse", "individualerror":False, "wlist":wlist, "relrange":relrange}
error_partial=partial(newerror,**args_fitnessfunction)
def fitnessfunc(pars):
    return error_partial(pars)[0],


pars_genetic={"POPULATION_SIZE":POPULATION_SIZE, "MAX_GENERATIONS": MAX_GENERATIONS, 
              "P_CROSSOVER": P_CROSSOVER,"P_MUTATION": P_MUTATION, "MUTSIGMA": MUTSIGMA, "INDPB": INDPB,"TOURNSIZE":10,
              "etacxbin":ETACXBIN,
              "HALL_OF_FAME_SIZE": 5,
             "fitnessfunc": fitnessfunc,
             "seeds":seeds,
              "plot_fitness_evo":False, 
              "plotintermediates":False, 
              "plotbest":False, "plotbestfunc":None, "bounds":bounds}
best=run_genetic(**pars_genetic)
fitness,parset,seed,refined=best
best_condition={"jid":jid,"seed":int(seed),"fitness":fitness,"parset":",".join(list(map(str,parset))),"refined":refined}
with open(outfname, "w") as outfile: 
    json.dump(best_condition, outfile)
