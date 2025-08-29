import glob,re,os,sys
import numpy as np
from functools import partial
import time,json
import itertools
sys.path.append('/home/rm335/repos/sharedposstpNov19/GeneRegulatoryFunctions/utilsGRF')
sys.path.append('/home/rm335/repos/myrepos/scriptsdu')
sys.path.append('../cfO2')
from utilsdu import return_fullpars_polcycle, get_constraints_npars_polcycle, get_score_up_down_v2, score

import BoundaryFinder as BF

import bPcycle_1bs_3states_rev1_TFA_prec50 as state3_rev1
import bPcycle_1bs_3states_rev2_TFA_prec50 as state3_rev2
import bPcycle_1bs_3states_rev1u2_TFA_prec50 as state3_rev12
import bPcycle_1bs_5states_rev1_TFA_prec50 as state5_rev1
import bPcycle_1bs_5states_rev1u2u3u4_TFA_prec50 as state5_rev1234




jid=int(sys.argv[1])-1
#first one is control, should not be non-monotonic
cases=[["3_1",["1f","2f"],"p_p"],
        ["3_1",["1f","2f"],"p_m"],
       ["3_1",["1f","2f"],"m_p"],
       ["3_1",["1r","2f"],"p_p"],
       ["3_1",["1r","2f"],"m_m"],
       ["3_1",["1r","3f"],"p_p"],
       ["3_1",["1r","3f"],"m_m"],
       ["3_2",["1f","2f"],"p_m"],
       ["3_2",["1f","2f"],"m_p"],
       ["3_2",["1f","2r"],"p_p"],
       ["3_2",["1f","2r"],"m_m"],
       ["3_2",["1f","3f"],"m_p"],
       ["3_2",["1f","3f"],"p_m"],
       ["3_12",["1f","2f"],"p_m"],
       ["3_12",["1f","2f"],"m_p"],
       ["3_12",["1f","3f"],"p_m"],
       ["3_12",["1f","3f"],"m_p"],
       ["5_1234",["1f","2f"],"p_m"],
       ["5_1234",["2f","4f"],"m_p"],
       ["5_1234",["2r","4f"],"p_p"],
       ["5_1234",["1r","3f"],"m_m"]]
#21 cases
#252 jobs

fcd=0.01
fcu=100
#parslimit=[[-2,2],[-3,3],[-4,4]]
extremesu=[[-2,2],[-1.5,1.5],[-1,1]]
prob_par=[0.2,0.5]
prob_replace=[0.2,0.6]

stepx=0.025
x_ar=np.arange(-3,15+stepx,stepx)
stepy=0.025
y_ar=np.arange(-15,15+stepy,stepy)

combination=list(itertools.product(cases,extremesu,prob_par,prob_replace))[jid]
case,extr_uniform,prob_par,prob_replace=combination
parslimit=[0,4] 
min_,max_=parslimit

model, transitions,direction=case

if model=="3_1":
    ccode=state3_rev1
elif model=="3_2":
    ccode=state3_rev2
elif model=="3_12":
    ccode=state3_rev12
elif model=="5_1":
    ccode=state5_rev1
elif model=="5_1234":
    ccode=state5_rev1234
else:
    print("model not understood")
    raise(ValueError)

constraints,npars=get_constraints_npars_polcycle(model,direction,fcd=fcd,fcu=fcu)
ssfunc=ccode.interfacess
return_fullpars=partial(return_fullpars_polcycle,model=model,transitions=transitions)

kwargs=dict(tol=1e-5,check_out0=False,out0_tol=0.001,fc_tol=0.070389327891398)
#tol is the tolerance to decide if at any given point, the function is increasing, decreasing or flat. 
#out0_tol is used to discard functions if they start with a value greater than out0_tol. This is only applied if check_out0 is set to True. It is also used to determine Amin and Amax in function score (named tol_A):
#if np.abs(np.log2(outA0/out0))<tol_A:Amin_found=True
#if np.abs(np.log2(outAm/outinf30))<tol_A:Amax_found=True
#fc_tol is set to 0.07039 which is log2(1.05/1). If the difference between maximum and minimum of the function (expressed as log2fc) is less than this, it is essentially flat, so discard.
  
myfunc=partial(score,return_fullpars=return_fullpars, ssfunc=ssfunc,scoref=get_score_up_down_v2,n=None,Amin=None,Amax=None,log2out=True,**kwargs)


name="model=%s_t=%s_d=%s"%(model,",".join(transitions),direction)
myfunc.__name__=name
settings={'pars_limit':[10**min_,10**max_],
          'constraints':constraints, 
          'compute_x_y_f':myfunc,
          'npars':npars,
          'row_ar':y_ar,
          'col_ar':x_ar,
          'seed':jid,
         'mat':None, #np.load('2019_08_06_matallrev.npy'),
         'mat_pars':None} #np.load('2019_08_06_matparsallrev.npy')}



niters_conv=25
niters=1000
L=10

name_save=name
#name_save="kinsyndifADsbnp"
dir_path = os.path.dirname(os.path.realpath(__file__))
dir_path2=dir_path.replace("/home","/n/scratch3/users/r")
outfolder=os.path.join(dir_path2,name_save+'_out_%d'%jid)

if not os.path.isdir(outfolder):
    os.mkdir(outfolder)


outfolder_final='final_results'

args={'niters':niters,
      'niters_conv':niters_conv,
      'niters_conv_points':100,
      'niters_save':1,
      'folder_save':outfolder,
       'name_save':name_save, 
      'prob_par':prob_par,
      'prob_replace':prob_replace,
      'extr_uniform':extr_uniform,
      'L_project':L,
      'plotting':False,
      'verbose':False,
     'dopulltangents':True} #

def function_tostring(x):
    if isinstance(x, np.ndarray):
        return ','.join(map(str,x))
    else:
        return x.__name__

outfnames=[os.path.join(outfolder,name_save+'_%d.sett'%jid),os.path.join(outfolder_final,name_save+'_%d.sett'%jid)]
for fname in outfnames:
    outf=open(fname,'w')
    #outf.write(time.ctime()+'\n')
    #with open(outf, 'w') as file:
    json.dump(dict({'time':time.ctime()},**settings),outf,default=function_tostring) # use `json.loads` to do the reverse
    outf.close()

    outf=open(fname.replace('.sett','.args'),'w')
    #outf.write(time.ctime()+'\n')
    #with open(outf, 'w') as file:
    json.dump(dict({'time':time.ctime()},**args),outf) # use `json.loads` to do the reverse
    outf.close()

BE=BF.BoundaryExplorer(**settings)
if settings['mat'] is None:
    BE.get_initial_points(10)
ti=time.time()
BE.extend_boundary(**args)
name='%s_%d_last'%(name_save,jid)
np.save(os.path.join(outfolder_final,'mat_'+name+'.npy'),BE.mat)
np.save(os.path.join(outfolder_final,'mat_pars_'+name+'.npy'),BE.mat_pars)
te=time.time()
print('time difference',te-ti)
print(BE.converged)
