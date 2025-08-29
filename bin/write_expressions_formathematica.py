import sympy
from sympy.printing import mathematica
import numpy as np

def write_expressions_formathematica(filename,combis_ij,variables,printi=7):
    math_output=open(filename,"w")
    all_terms=[]
    for tn,term in enumerate(combis_ij):
        termsc=[]
        print("term", tn)
        for monomial in sympy.Add.make_args(term):
            #print(monomial, sympy.Mul.make_args(monomial))
            noneps=[]
            eps=[]
            number=1
            eps_base=[] #bases, for ordering by letter
            noneps_base=[] #bases, for ordering by letter
            args= sympy.Mul.make_args(monomial) #separate each monomial into each variable (letters,number)
            for arg_ in args:
                if type(arg_)==sympy.core.power.Pow:
                    arg=sympy.core.power.Pow.as_base_exp(arg_)[0] #get the base of the power so that I can check if it is epsilon or not
                else:
                    arg=arg_
                if  arg in variables: #these parameters are generally called epsilons
                    eps.append(arg_)
                    eps_base.append(arg.name)
                else:
                    
                    #print(type(arg))
                    if type(arg)!=sympy.core.numbers.Integer and type(arg)!=sympy.core.numbers.NegativeOne:
                        noneps_base.append(arg.name)
                        noneps.append(arg_)
                    else:
                        number=arg
                        
            
            #print(number)
            noneps_o=np.array(noneps)[np.argsort(np.array(noneps_base))]
            eps_o=np.array(eps)[np.argsort(np.array(eps_base))]
            
            termsc.append([number,noneps_o,eps_o])
            if tn>printi:
                print(monomial)
                print(number,noneps_o,eps_o)
        print("--------")
        all_terms.append(termsc)

    for i in range(len(all_terms)):
        print("term ",i)
        coeffi=all_terms[i]
        #print(coeffi)
        base=[] #list of common factors that do not contain epsilons
        #find the unique common factors without epsilons
        for j in range(len(coeffi)):
            b=coeffi[j][1]
            present=False
            for k in range(len(base)): #check if it is already in the list base, and if it is not there, add it
                if np.all(base[k]==b):
                    present=True
            if not present:
                base.append(b)
        #find which term in the original coefficient expression correspond to each common factor
        common_bases_idx=[[] for i in range(len(base))] #which terms (indices) have those common factors
        for j in range(len(coeffi)):
            b=coeffi[j][1]
            for k in range(len(base)):
                if np.all(base[k]==b):
                    common_bases_idx[k].append(j)
        #now for each base, print the epsilon terms side by side to check which sign wins 
        print(len(base))
        if i>printi:
            print(base)
        for k in range(len(base)):
            if i>printi:
                print("epsilons for base", end=":")
                print(base[k])
            idxs=common_bases_idx[k]
            epsilons_terms=0
            for idx in idxs:
                terms=coeffi[idx]
                if i>printi:
                    print(terms)
                number=terms[0]
                epsilons=terms[2]
                prod=1
                for e in epsilons:
                    prod*=e

                #if np.abs(number)==1:
                epsilons_terms+=number*prod
                #else:
                #    if number<1:
                #        sign=-1
                #    else:
                #        sign=1
                #    for k in range(np.abs(number)):
                #        epsilons_terms+=sign*prod
            #display(epsilons_terms)
            #display(idx,sympy.factor(epsilons_terms),epsilons_terms)
            math_output.write(str(i)+";"+mathematica.mathematica_code(epsilons_terms).replace("_","")+"\n")
    math_output.close()
