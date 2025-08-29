import numpy as np
from scipy.special import logsumexp

def avb_1sites(Ks,x):
    K_1,=Ks
    #here directly we can work with the natural scale
    K_1=np.exp(K_1)
    x=np.exp(x)
    x1=x 
    mu_0=1
    mu_1=x1*K_1
    Z=mu_0+mu_1
    prob_ar=np.array([mu_0/Z,mu_1/Z])
    returnvalue=prob_ar[1]
    return returnvalue



def avb_2sites(Ks,x):
    K_1, K_2=Ks
    #here directly we can work with the natural scale
    K_1=np.exp(K_1)
    K_2=np.exp(K_2)
    x=np.exp(x)
    x1=x 
    x2=x**2
    mu_0=1
    mu_1=x1*K_1
    mu_2=x1*K_2
    mu_3=x2*K_1*K_2
    Z=mu_0+mu_1+mu_2+mu_3
    prob_ar=np.array([mu_0/Z,mu_1/Z,mu_2/Z,mu_3/Z])
    returnvalue=prob_ar[1]+prob_ar[2]+2*prob_ar[3] #average binding: i*p_i, where i is the number of bound molecules
    return returnvalue
def avb_3sites_nologsumexp(Ks,x):
    
    K_1,K_2,K_3=Ks #should be log
    K_1=np.exp(K_1)
    K_2=np.exp(K_2)
    K_3=np.exp(K_3)
    
    x=np.exp(x)

    x1=x 
    x2=x**2  
    x3=x**3
    mu_0=1
    mu_1=x1*K_1
    mu_2=x1*K_2
    mu_3=x1*K_3
    mu_4=x2*K_1*K_2
    mu_5=x2*K_1*K_3
    mu_6=x2*K_2*K_3
    mu_7=x3*K_1*K_2*K_3
        
        
    mu_ar=np.array([mu_0,mu_1,mu_2,mu_3,mu_4,mu_5,mu_6,mu_7])
    Z=np.sum(mu_ar)
    prob_ar=mu_ar/Z
    
    #output=="avbinding":
    idxs_1=[1,2,3]
    idxs_2=[4,5,6]
    idxs_3=[7]
    f1=np.sum(prob_ar[idxs_1])
    f2=np.sum(prob_ar[idxs_2])
    f3=np.sum(prob_ar[idxs_3])
    returnvalue=np.sum([f1,2*f2,3*f3])

    return returnvalue

def avb_3sites(Ks,x):
    
    K_1,K_2,K_3=Ks #should be log
    


    x1=x #log(x)
    x2=2*x  #log(x^2)
    x3=3*x
    mu_0=0
    mu_1=x1+K_1
    mu_2=x1+K_2
    mu_3=x1+K_3
    mu_4=x2+K_1+K_2
    mu_5=x2+K_1+K_3
    mu_6=x2+K_2+K_3
    mu_7=x3+K_1+K_2+K_3
        
        
    mu_ar=np.array([mu_0,mu_1,mu_2,mu_3,mu_4,mu_5,mu_6,mu_7])
    c=np.max(mu_ar)
    #print(mu_ar)
    #print(c)
    #print("mu_ar-c", mu_ar-c)
    y=c+np.log(np.sum(np.exp(mu_ar-c))) #log-sum-exp trick
    prob_ar=np.exp(mu_ar-y) #equilibrium probability of each state
    
    #output=="avbinding":
    idxs_1=[1,2,3]
    idxs_2=[4,5,6]
    idxs_3=[7]
    f1=np.sum(prob_ar[idxs_1])
    f2=np.sum(prob_ar[idxs_2])
    f3=np.sum(prob_ar[idxs_3])
    returnvalue=np.sum([f1,2*f2,3*f3])

    return returnvalue

def avb_4sites(Ks,x):
    
    K_1,K_2,K_3,K_4=Ks #should be log 
    
    
    x1=x #log(x)
    x2=2*x  #log(x^2)
    x3=3*x
    x4=4*x 
    mu_0=0
    mu_1=x1+K_1
    mu_2=x1+K_2
    mu_3=x1+K_3
    mu_4=x1+K_4
    mu_5=x2+K_1+K_2
    mu_6=x2+K_1+K_3
    mu_7=x2+K_1+K_4
    mu_8=x2+K_2+K_3
    mu_9=x2+K_2+K_4
    mu_10=x2+K_3+K_4
    mu_11=x3+K_1+K_2+K_3
    mu_12=x3+K_1+K_2+K_4
    mu_13=x3+K_1+K_3+K_4
    mu_14=x3+K_2+K_3+K_4
    mu_15=x4+K_1+K_2+K_3+K_4

    
    mu_ar=np.array([mu_0,mu_1,mu_2,mu_3,mu_4,mu_5,mu_6,mu_7,mu_8,mu_9,mu_10,mu_11,mu_12,mu_13,mu_14,mu_15])
    c=np.max(mu_ar)
    #print(mu_ar)
    #print(c)
    #print("mu_ar-c", mu_ar-c)
    y=c+np.log(np.sum(np.exp(mu_ar-c))) #log-sum-exp trick
    prob_ar=np.exp(mu_ar-y)

    #output=="avbinding":
    idxs_1=[1,2,3,4]
    idxs_2=[5,6,7,8,9,10]
    idxs_3=[11,12,13,14]
    idxs_4=[15]
    f1=np.sum(prob_ar[idxs_1])
    f2=np.sum(prob_ar[idxs_2])
    f3=np.sum(prob_ar[idxs_3])
    f4=np.sum(prob_ar[idxs_4])
    returnvalue=np.sum([f1,2*f2,3*f3,4*f4])

    return returnvalue

def avb_5sites(Ks,x):
    
    K_1,K_2,K_3,K_4,K_5=Ks
        
    x1=x
    x2=2*x
    x3=3*x
    x4=4*x
    x5=5*x
    mu_0=0 #log of mu
    mu_1=x1+K_1
    mu_2=x1+K_2
    mu_3=x1+K_3
    mu_4=x1+K_4
    mu_5=x1+K_5
    mu_6=x2+K_1+K_2
    mu_7=x2+K_1+K_3
    mu_8=x2+K_1+K_4
    mu_9=x2+K_1+K_5
    mu_10=x2+K_2+K_3
    mu_11=x2+K_2+K_4
    mu_12=x2+K_2+K_5
    mu_13=x2+K_3+K_4
    mu_14=x2+K_3+K_5
    mu_15=x2+K_4+K_5
    mu_16=x3+K_1+K_2+K_3
    mu_17=x3+K_1+K_2+K_4
    mu_18=x3+K_1+K_2+K_5
    mu_19=x3+K_1+K_3+K_4
    mu_20=x3+K_1+K_3+K_5
    mu_21=x3+K_1+K_4+K_5
    mu_22=x3+K_2+K_3+K_4
    mu_23=x3+K_2+K_3+K_5
    mu_24=x3+K_2+K_4+K_5
    mu_25=x3+K_3+K_4+K_5
    mu_26=x4+K_1+K_2+K_3+K_4
    mu_27=x4+K_1+K_2+K_3+K_5
    mu_28=x4+K_1+K_2+K_4+K_5
    mu_29=x4+K_1+K_3+K_4+K_5
    mu_30=x4+K_2+K_3+K_4+K_5
    mu_31=x5+K_1+K_2+K_3+K_4+K_5
    
    mu_ar=np.array([mu_0,mu_1,mu_2,mu_3,mu_4,mu_5,mu_6,mu_7,mu_8,mu_9,mu_10,mu_11,mu_12,mu_13,mu_14,mu_15,mu_16,mu_17,mu_18,mu_19,mu_20,mu_21,mu_22,mu_23,mu_24,mu_25,mu_26,mu_27,mu_28,mu_29,mu_30,mu_31])
    c=np.max(mu_ar)
    #print(mu_ar)
    #print(c)
    #print("mu_ar-c", mu_ar-c)
    y=c+np.log(np.sum(np.exp(mu_ar-c))) #log-sum-exp trick
    prob_ar=np.exp(mu_ar-y)
    
    idxs_1=[1,2,3,4,5]
    idxs_2=[6,7,8,9,10,11,12,13,14,15]
    idxs_3=[16,17,18,19,20,21,22,23,24,25]
    idxs_4=[26,27,28,29,30]
    idxs_5=[31]
    f1=np.sum(prob_ar[idxs_1])
    f2=np.sum(prob_ar[idxs_2])
    f3=np.sum(prob_ar[idxs_3])
    f4=np.sum(prob_ar[idxs_4])
    f5=np.sum(prob_ar[idxs_5])
    returnvalue=np.sum([f1,2*f2,3*f3,4*f4,5*f5])

    return returnvalue
        
    
def avb_6sites(Ks,x):
    
    K_1,K_2,K_3,K_4,K_5,K_6=Ks
    

    x1=x
    x2=2*x
    x3=3*x
    x4=4*x
    x5=5*x
    x6=6*x
    mu_0=0
    mu_1=x1+K_1
    mu_2=x1+K_2
    mu_3=x1+K_3
    mu_4=x1+K_4
    mu_5=x1+K_5
    mu_6=x1+K_6
    mu_7=x2+K_1+K_2
    mu_8=x2+K_1+K_3
    mu_9=x2+K_1+K_4
    mu_10=x2+K_1+K_5
    mu_11=x2+K_1+K_6
    mu_12=x2+K_2+K_3
    mu_13=x2+K_2+K_4
    mu_14=x2+K_2+K_5
    mu_15=x2+K_2+K_6
    mu_16=x2+K_3+K_4
    mu_17=x2+K_3+K_5
    mu_18=x2+K_3+K_6
    mu_19=x2+K_4+K_5
    mu_20=x2+K_4+K_6
    mu_21=x2+K_5+K_6
    mu_22=x3+K_1+K_2+K_3
    mu_23=x3+K_1+K_2+K_4
    mu_24=x3+K_1+K_2+K_5
    mu_25=x3+K_1+K_2+K_6
    mu_26=x3+K_1+K_3+K_4
    mu_27=x3+K_1+K_3+K_5
    mu_28=x3+K_1+K_3+K_6
    mu_29=x3+K_1+K_4+K_5
    mu_30=x3+K_1+K_4+K_6
    mu_31=x3+K_1+K_5+K_6
    mu_32=x3+K_2+K_3+K_4
    mu_33=x3+K_2+K_3+K_5
    mu_34=x3+K_2+K_3+K_6
    mu_35=x3+K_2+K_4+K_5
    mu_36=x3+K_2+K_4+K_6
    mu_37=x3+K_2+K_5+K_6
    mu_38=x3+K_3+K_4+K_5
    mu_39=x3+K_3+K_4+K_6
    mu_40=x3+K_3+K_5+K_6
    mu_41=x3+K_4+K_5+K_6
    mu_42=x4+K_1+K_2+K_3+K_4
    mu_43=x4+K_1+K_2+K_3+K_5
    mu_44=x4+K_1+K_2+K_3+K_6
    mu_45=x4+K_1+K_2+K_4+K_5
    mu_46=x4+K_1+K_2+K_4+K_6
    mu_47=x4+K_1+K_2+K_5+K_6
    mu_48=x4+K_1+K_3+K_4+K_5
    mu_49=x4+K_1+K_3+K_4+K_6
    mu_50=x4+K_1+K_3+K_5+K_6
    mu_51=x4+K_1+K_4+K_5+K_6
    mu_52=x4+K_2+K_3+K_4+K_5
    mu_53=x4+K_2+K_3+K_4+K_6
    mu_54=x4+K_2+K_3+K_5+K_6
    mu_55=x4+K_2+K_4+K_5+K_6
    mu_56=x4+K_3+K_4+K_5+K_6
    mu_57=x5+K_1+K_2+K_3+K_4+K_5
    mu_58=x5+K_1+K_2+K_3+K_4+K_6
    mu_59=x5+K_1+K_2+K_3+K_5+K_6
    mu_60=x5+K_1+K_2+K_4+K_5+K_6
    mu_61=x5+K_1+K_3+K_4+K_5+K_6
    mu_62=x5+K_2+K_3+K_4+K_5+K_6
    mu_63=x6+K_1+K_2+K_3+K_4+K_5+K_6
    mu_ar=np.array([mu_0,mu_1,mu_2,mu_3,mu_4,mu_5,mu_6,mu_7,mu_8,mu_9,mu_10,mu_11,mu_12,mu_13,mu_14,mu_15,mu_16,mu_17,mu_18,mu_19,mu_20,mu_21,mu_22,mu_23,mu_24,mu_25,mu_26,mu_27,mu_28,mu_29,mu_30,mu_31,mu_32,mu_33,mu_34,mu_35,mu_36,mu_37,mu_38,mu_39,mu_40,mu_41,mu_42,mu_43,mu_44,mu_45,mu_46,mu_47,mu_48,mu_49,mu_50,mu_51,mu_52,mu_53,mu_54,mu_55,mu_56,mu_57,mu_58,mu_59,mu_60,mu_61,mu_62,mu_63])
    c=np.max(mu_ar)
    #print(mu_ar)
    #print(c)
    #print("mu_ar-c", mu_ar-c)
    y=c+np.log(np.sum(np.exp(mu_ar-c))) #log-sum-exp trick
    prob_ar=np.exp(mu_ar-y)
       
    #elif output=="avbinding":
    idxs_1=[1, 2, 3, 4, 5, 6]
    idxs_2=[7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]
    idxs_3=[22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41]
    idxs_4= [42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56]
    idxs_5=[57,58,59,60,61,62]
    idxs_6=[63]
    f1=np.sum(prob_ar[idxs_1])
    f2=np.sum(prob_ar[idxs_2])
    f3=np.sum(prob_ar[idxs_3])
    f4=np.sum(prob_ar[idxs_4])
    f5=np.sum(prob_ar[idxs_5])
    f6=np.sum(prob_ar[idxs_6])
    returnvalue=np.sum([f1,2*f2,3*f3,4*f4,5*f5,6*f6])

    return returnvalue

def foldchange_from_avbinding_ratesdirectly_Hill(globalpars=[], TFpars=[], avbinding=1, rates_affected=[0,1,2,3],tolerance_low=1e-15, tolerance_high=1e15):
    
    #nonlinear function with saturation to a max or min value for the effect that the TF has on the activation energy of the transition.
    #directly consider the transition rates and the change in them 
    #length of globalpars: len(rates_affected)
    #TFpars are the parameters that determine the factor by which the bound TFs modulate each rate. 3 parameters per rate. Length=3*len(rates_affected)
    #avbinding is average number of TF molecules bound on DNA
    #rates_affected: should be indexes: 0: k1, 1: km1, 2: k2, 3:k3. Indexes of the rates that can be modified by the TFs
    
    
    
    
    nk=len(rates_affected)
    nratescycle=4 #we could explore other cycles but the 3-state cycle with the first transition reversible should be enough
    ksTF=[] #ln of rates of pol cycle when TF is present in the system
    ks0=globalpars #ln of basal rates of pol cycle
    j=0
    
    
    for i in range(nratescycle):#for each rate there is a c and theta parameter. This should be shared among all TFs
        log_k_i=globalpars[i] 
        k_i=np.exp(log_k_i)
        if i in rates_affected:
            #assume that the activation energy is a function of the TF bound as theta_TF=theta+ftheta;
            #ftheta=sign*alpha*avbinding/(K+avbinding). If alpha>0, increase in theta with more TF, if less than 0, decrease in theta. Small alpha: small increase/decrease. Large alpha: large increase/decrease
            #this function enables saturation
            K,n,foldchange=np.exp(TFpars[3*j:3*j+3])
            k_sat=k_i*foldchange
            avb_n=avbinding**n
            K_n=K**n
            
            #print(avbinding, K, n, avb_n, K_n)
            #if np.isnan(avb_n):
            #    print("nan in avbn", avbinding, n)
            #if np.isnan(K_n):
            #    print("nan in K_n", K, n)

            if avb_n<tolerance_low or K_n>tolerance_high: #no average binding, or too large a threshold that we will always have 0
                f0=0
            elif np.isinf(avb_n):
                f0=1
            else:
                f0=avb_n/(K_n+avb_n)

            if np.isnan(f0):
                print("nan", avbinding, K, n, avb_n, K_n)
                raise 

            #log_f=np.log(avbn)-np.log(K_n+avbn)
            
            if np.isinf(f0): #
                if f0>0:
                    print("inf f0, resetting to 1", f0, avbn, K_n)
                    f0=1.0
                else:
                    f0=0.0
                    print("inf f0, resetting to 0", f0, avbn, K_n)

            k_TF=k_i+(k_sat-k_i)*f0
            j+=1
            ksTF.append(np.log(k_TF))
        else:
            ksTF.append(log_k_i)
    #print("ks0", ks0)
    #print("ksTF", ksTF)
    k1,km1,k2,k3=ks0 #log of ks
    #(k3*k1*k2)/(km1*k3+k2*k3+k1*k3+k1*k2
    #output0=np.exp(k3+k1+k2)/(np.exp(km1+k3)+np.exp(k2+k3)+np.exp(k1+k3)+np.exp(k1+k2)) #division by 0 should not happen because we have bounds on c and theta
    log_den=logsumexp([km1+k3,k2+k3,k1+k3,k1+k2]) #log(np.exp(km1+k3)+np.exp(k2+k3)+np.exp(k1+k3)+np.exp(k1+k2))
    log_output0=k1+k2+k3-log_den #log_output0=log(k1*k2*k3)-log(np.exp(km1+k3)+np.exp(k2+k3)+np.exp(k1+k3)+np.exp(k1+k2))
    
    k1,km1,k2,k3=ksTF  #log of ks
    log_den=logsumexp([km1+k3,k2+k3,k1+k3,k1+k2])
    log_outputTF=k1+k2+k3-log_den
    #if output0<1e-15:
    #    print("small output 0", output0, km1)
    #if outputTF<1e-3:
    #    #print("output less 1e-3", outputTF)
    #    return -100 #avoid working with super small values which can lead to spurious results. 
    #else:
    return np.exp(log_outputTF-log_output0) #fold change        