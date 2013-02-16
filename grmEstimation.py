#!/usr/bin/env python
''' Module for the estimation of a generalized Roy Model.

'''

# standard library
import sys
import numpy as np
import json
import os

from scipy.stats        import norm
from scipy.optimize     import fmin_bfgs

# edit PYTHONPATH
fileName = os.path.realpath(__file__).replace('/grmEstimation.py','')  
sys.path.insert(0, fileName)

# project library
import grmInitFileReader as auxInit

''' Maximization Functions.
'''
def _distributeEvaluationValues(x, numCovarsOut, isList = False):
    ''' Distribute the evaluation values.
    
    '''
    #Antibugging.
    assert (isList in [True, False])
    
    if(isList):
        
        pass
            
    else:
        
        assert (isinstance(x, np.ndarray))
        assert (np.all(np.isfinite(x)))
        assert (x.ndim == 1)
        
    assert (isinstance(numCovarsOut, int))
    assert (numCovarsOut > 0)
        
    # Distribution of parameters.
    rslt = {}

    rslt['Y1_beta'] = x[:numCovarsOut]
    rslt['Y0_beta'] = x[numCovarsOut:(2*numCovarsOut)]
    
    rslt['D_gamma'] = x[(2*numCovarsOut):(-4)]
    
    ''' Parameters with natural bounds. 
    '''    
    rslt['U1_var']  = np.exp(x[(-4)])
    rslt['U0_var']  = np.exp(x[(-3)])

    rslt['U1V_rho'] = -1.0 + 2.0/(1.0 + np.exp(-x[-2]))    
    rslt['U0V_rho'] = -1.0 + 2.0/(1.0 + np.exp(-x[-1])) 

    # Finishing.
    return rslt

def maxAlgorihtmInterface(x, Y, D, X, Z):
    ''' Interface to the SciPy maximization routines.
    
    '''
    
    # Collect maximization arguments.
    rslt = _distributeEvaluationValues(x, numCovarsOut)
    
    # Calculate likelihood.
    likl = negLogLiklContribution(rslt, Y, D, X, Z)
    
    # Finishing.
    return likl
     
def negLogLiklContribution(rslt, Y, D, X, Z):
    ''' Negative log-likelihood function of the Generalized Roy Model.
    
    '''
    # Distribute input.
    Y1_beta = rslt['Y1_beta']
    Y0_beta = rslt['Y0_beta']

    D_gamma = rslt['D_gamma']

    U1_var  = rslt['U1_var']
    U0_var  = rslt['U0_var']

    U1V_rho = rslt['U1V_rho']
    U0V_rho = rslt['U0V_rho']

    # Construct auxiliary objects.
    numAgents = D.shape[0]

    # Likelihood calculation.
    choiceCoeffs  = np.concatenate((Y1_beta  - Y0_beta, - D_gamma))
    choiceCovars  = np.concatenate((X, Z), axis = 1)
    choiceIndices = np.dot(choiceCoeffs, choiceCovars.T) 
    
    argOne = D * (Y - np.dot(Y1_beta, X.T))/np.sqrt(U1_var) + (1 - D)*(Y - np.dot(Y0_beta, X.T))/np.sqrt(U0_var)
    argTwo = D * (choiceIndices - U1V_rho*argOne)/np.sqrt(1.0 - U1V_rho**2) + \
            (1 - D) *  (choiceIndices - U0V_rho*argOne)/np.sqrt(1.0 - U0V_rho**2)
    
    cdfEvals = norm.cdf(argTwo)
    pdfEvals = norm.pdf(argOne)
    
    likl = D*(1.0/np.sqrt(U1_var))*pdfEvals*cdfEvals + \
                (1 - D) *(1.0/np.sqrt(U0_var))*pdfEvals*(1.0  - cdfEvals)
    
    # Transformations.
    likl =  np.clip(likl, 1e-20, 1.0)
    
    likl = -np.log(likl)
    
    likl = likl.sum()
    
    lik = (1.0/numAgents)*likl
        
    # Quality checks.
    assert (isinstance(lik, float))    
    assert (np.isfinite(lik))
    assert (lik > 0.0)
    
    #Finishing.        
    return lik

''' *.ini file processing.
'''
auxInit.readInitFile()

with open('init.json', 'r') as file_:
    
    initDict = json.load(file_)

''' Distribute useful covariates.
'''
numAgents = initDict['numAgents']
fileName  = initDict['fileName']

Y1_beta   = np.array(initDict['Y1_beta'])
Y0_beta   = np.array(initDict['Y0_beta'])

D_gamma   = np.array(initDict['D_gamma'])    

U1_var    = initDict['U1_var'] 
U0_var    = initDict['U0_var'] 
V_var     = initDict['V_var']

U1V_rho   = initDict['U1V_rho']  
U0V_rho   = initDict['U0V_rho']  

''' Construct auxiliary objects.
'''
numCovarsOut  = Y1_beta.shape[0]
numCovarsCost = D_gamma.shape[0]

''' Read and check dataset, distribute entries.
'''
data = np.genfromtxt(initDict['fileName'], dtype = 'float')

assert (data.shape == (numAgents, (1 + 1 + numCovarsOut +  numCovarsCost))) 
assert (np.all(np.isfinite(data)))
assert ((data[:,1].all() in [1.0, 0.0]))

Y = data[:,0]
D = data[:,1]

X = data[:,2:(numCovarsOut + 2)]
Z = data[:,-numCovarsCost:]

''' Maximization Script.
'''

startVals = np.concatenate((initDict['Y1_beta'], initDict['Y0_beta'], initDict['D_gamma'], \
                            np.array(initDict['U1_var'], ndmin = 1),\
                            np.array(initDict['U0_var'], ndmin = 1),\
                            np.array(initDict['U1V_rho'], ndmin = 1),\
                            np.array(initDict['U0V_rho'], ndmin = 1) ))

rslt = fmin_bfgs(maxAlgorihtmInterface, startVals, \
                 args = (Y, D, X, Z), maxiter = 10, \
                 full_output = True)

rslt = rslt[0].tolist()

rslt = _distributeEvaluationValues(rslt, numCovarsOut, isList = True)

''' Output results to *.json file.
'''
with open('rslt.json', 'w') as file_:
        
    json.dump(rslt, file_)