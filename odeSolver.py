from scipy.integrate import solve_ivp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phaseportrait


class System:
    def __init__(self, species_dict, coefficients, exponents, kinetics):
        
        # save as attributes
        self.coeffs = coefficients
        self.kinetics = kinetics
        self.species_d = species_dict
        self.exps = exponents
        
    def __call__(self, t, C):
        variables = C.tolist()
        
        coefs = self.coeffs
        expos = self.exps
        kins = self.kinetics
        args = (coefs, expos, kins)
        return diffy(t, variables, *args)

        
        
def make_exps(stmat):
    ## takes in the st matrix to generate the exps
    ## STMAT SHOULD BE list of reaction lists
    exps = []
    
    # iterate across each reaction vector
    for reaction in stmat:
        r_exps = []
        for entry in reaction:
            if entry < 0:
                exp = -1 * entry
                r_exps.append(exp)
            else:
                r_exps.append(0)
        exps.append(r_exps)
    final_exps = tuple(exps)
    print(final_exps)
    return(final_exps)



def diffy(t, variables, coeffs, exps, kinetics):
    #print(variables)
    try:
        len(variables)
        n = len(variables)
    except TypeError:
        n = 1
    re = np.zeros(n)
    r_ct = 0
    # dCdt = big sum over reacs(c_a,c_b,..)
    for k, exp in zip(kinetics, exps):
        #print("16", kinetics, k, exp)
        weight_prod = float(k)
        weight = np.power(variables, exp)
        weight = weight.tolist()

        for w in weight:
            weight_prod = weight_prod*w
        for i in range(n):
            # n species

            coeff = coeffs[i][r_ct]
            re[i] += coeff * weight_prod

        r_ct += 1
    re = re.tolist()
    return re


def gen_coeffs(species_dict):
    """
    parameters:  # of reactions (int), species dict
    does: populates a list of m (reacs) lists of length n (sp) accordingly (the coeffs)
    also does: iterates through the reactants of each reaction to calculate the exps
    """
    while True:
        r = input("how many reactions?\n")
        try:
            int(r)
        except ValueError:
            pass
        else:
            r = int(r)
            break
    nums = len(species_dict)
    # get reactions
    # populate R matrix
    
    # generate 0 matrix
    array = np.zeros((nums,r), dtype=int)
    df = pd.DataFrame(array, columns =list(range(r)))
    
    # generate empty list for exps 
    # should be a list of m (reacts) lists of len n
    exps = []

    
    for i in range(r):
        reacts = input(f"reactants in r{i+1} ex. A,2*B \n")
        # set up default 0s for exps for reaction i
        k_array = np.zeros(nums)
        r_i_exps = list(k_array)
        
        #populate reactants
        for r in reacts.split(sep=","):
            if r == "":
                pass
            else:
                s, c = coeffs(r)
                row_indx = species_dict[s]
                curr = df.at[row_indx, i]
                df.at[row_indx, i] = curr - c
                r_i_exps[row_indx] = c
                #print(df)
        exps.append(tuple(r_i_exps))
        
        #pop prods
        prods = input(f"prods in r{i+1} ex. 2*B,3*C \n")
        for p in prods.split(sep=","):
            if p == "":
                pass
            else:
                sp, co = coeffs(p)
                row_idx = species_dict[sp]
                curr = df.at[row_idx, i]
                df.at[row_idx, i] = curr + co
                #print(df)
            
    coeffz = df.values.tolist()
    coeffz = tuple(coeffz)
    print("exps are", exps, "\n coeffs are", coeffz)

    return coeffz, exps

def coeffs(spec):
    ### takes a str of one species ex A, 2*B and returns the species and coeff
    if "*" in spec:
        id = spec.index("*")
        coeff = int(spec[:id])
        species = spec[id+1:]
    else:
        coeff = 1
        species = spec
    return species, coeff

def make_sys():
    
    ## gather species info
    species_dict = {}
    num_species = int(input("how many species?\n"))
    for i in range(num_species):
        num = i + 1
        specy = input(f"enter species #{num}\n")
        species_dict[specy] = i
        
    ## gather reaction info
    coeffz, exps = gen_coeffs(species_dict)
    while True:
        whatisk = input("values of k, separated by commas\n")
        whatisk = whatisk.split(sep= ",")
        if len(whatisk) == len(exps):
            break
    sys1 = System(species_dict, coeffz, exps, whatisk)
    
        
    return sys1

def solplot(sys_obj):
    
    while True:
        numys = input("how many y0s?")
        try:
            int(numys)
            break
        except ValueError:
            print("enter an integer!")
            
    ys = int(numys)
    for i in range(ys):
            
        while True:
            whatsy0 = input("values of y0, separated by commas\n")
            whatsy0 = whatsy0.split(sep= ",")
            coeffz = sys_obj.coeffs
            if len(whatsy0) == len(coeffz):
                break
        
        
        t_span = (0.0, 60.0)
        t = np.arange(0.0, 60, 0.1)
    
    
        sol = solve_ivp(sys_obj, t_span, whatsy0)
        
        ymax = 0
        for i in range(sol.y.shape[0]):
            plt.xlim(t_span)
    
            
            yvals = sol.y[i,:]
            floatvals = yvals.astype(np.float64)
            highesty = np.amax(floatvals)
            if highesty > ymax:
                ymax = highesty
            plt.plot(sol.t, floatvals, label=f'$C_{i}(t)$')
    
    plt.ylim(0,ymax+1)
    plt.xlabel('$t$') # the horizontal axis represents the time 
    plt.legend() # show how the colors correspond to the components of X
    plt.show()
    

    
def main():
    plt.clf()
    obj = make_sys()
    solplot(obj)
    
main()

