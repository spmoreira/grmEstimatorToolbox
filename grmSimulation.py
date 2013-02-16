#!/usr/bin/env python
''' Module for the simulation of a generalized Roy Model.



    Development Notes:
    
        The first element in X is always the intercept.
'''

# standard library
import numpy as np

import json
import os
import sys

from scipy.stats        import norm
from scipy.optimize     import fmin_bfgs

# edit PYTHONPATH
fileName = os.path.realpath(__file__).replace('/grmSimulation.py','')  
sys.path.insert(0, fileName)

# project library
import grmInitFileReader as auxInit

''' Set random seed
'''
np.random.seed(123)

''' Read specification from *.ini file.
'''
auxInit.readInitFile()

with open('init.json', 'r') as file_:
    
    initDict = json.load(file_)

''' Distribute parametrization and (limited) type conversions.
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

U1V_cov      = U1V_rho*np.sqrt(U1_var)*np.sqrt(V_var)
U0V_cov      = U0V_rho*np.sqrt(U0_var)*np.sqrt(V_var)

''' Simulate observable agent characteristics.
'''
means = np.tile(0.0, numCovarsOut)
covs  = np.identity(numCovarsOut)

X      = np.random.multivariate_normal(means, covs, numAgents)
X[:,0] = 1.0

means = np.tile(0.0, numCovarsCost)
covs  = np.identity(numCovarsCost)

Z     = np.random.multivariate_normal(means, covs, numAgents)

''' Construct level indicators for outcomes and choices. 
'''
Y1_level = np.dot(Y1_beta, X.T)
Y0_level = np.dot(Y0_beta, X.T)
D_level  = np.dot(D_gamma, Z.T)

''' Simulate unobservables from the model.
'''
means = np.tile(0.0, 3)
vars_ = [U1_var, U0_var, V_var]

covs  = np.diag(vars_)

covs[0,2] = U1V_cov 
covs[2,0] = covs[0,2]

covs[1,2] = U0V_cov
covs[2,1] = covs[1,2]

U = np.random.multivariate_normal(means, covs, numAgents)

''' Simulate individual outcomes and choices.
'''
Y1 = np.tile(np.nan, (numAgents))
Y0 = np.tile(np.nan, (numAgents))
Y  = np.tile(np.nan, (numAgents))

D  = np.tile(np.nan, (numAgents))

for i in range(numAgents):
    
    # Distribute unobservables.
    U1 = U[i,0]
    U0 = U[i,1]
    V  = U[i,2]

    # Decision Rule.
    expectedBenefits = Y1_level[i] - Y0_level[i]
    cost             = D_level[i]  + V 
    
    D[i] = np.float((expectedBenefits - cost > 0))
    
    # Potential outcomes.
    Y1[i] = Y1_level[i] + U1
    Y0[i] = Y0_level[i] + U0
    
    # Observed outcomes.
    Y[i]  = D[i]*Y1[i] + (1.0 - D[i])*Y0[i]
    
''' Check quality of simulated sample. 
'''
assert (np.all(np.isfinite(Y1)))
assert (np.all(np.isfinite(Y0)))

assert (np.all(np.isfinite(Y)))
assert (np.all(np.isfinite(D)))

assert (Y1.shape == (numAgents, ))
assert (Y0.shape == (numAgents, ))

assert (Y.shape  == (numAgents, ))
assert (D.shape  == (numAgents, ))

assert (Y1.dtype == 'float')
assert (Y0.dtype == 'float')

assert (Y.dtype == 'float')
assert (D.dtype == 'float')

assert ((D.all() in [1.0, 0.0]))
   
''' Export sample to *.txt file for further processing. 
'''
np.savetxt(fileName, np.column_stack((Y, D, X, Z)), fmt= '%8.3f')

print 'simulation request completed.'