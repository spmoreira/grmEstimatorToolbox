'''
    This module is an extension of the Generalized Roy Toolbox
    by  Philipp Eisenhauer and Stefano Mosso (Copyright 2013)
    
    It computes average treatment effects: ATE, TT, TUT
    
    No parallel version

'''

import os.path
import sys
import json

import  numpy as np

# project library
import grmReader

''' Public Interface.
'''
def treatmentEffects():
    ''' Public interface to request an estimation of TE using the Generalized
        Roy Model.
    '''
    
    # Checks.
    assert (os.path.exists('grmRslt.json')) 
    
    # Process initialization file.
    initDict = grmReader.read()
    
    # Checks.
    assert (os.path.exists(initDict['fileName'])) 
    
    # Process initFile.
    #_initializeLogging()
    
    ''' Distribute useful covariates.
    '''
    
    ''' Distribute parametrization and (limited) type conversions.
    '''
    numAgents  = initDict['numAgents']
    fileName   = initDict['fileName']
    
    Y1_beta    = np.array(initDict['Y1_beta'])
    Y0_beta    = np.array(initDict['Y0_beta'])
    
    D_gamma    = np.array(initDict['D_gamma'])
    
    U1_var     = initDict['U1_var'] 
    U0_var     = initDict['U0_var'] 
    V_var      = initDict['V_var']
    
    U1V_rho    = initDict['U1V_rho']  
    U0V_rho    = initDict['U0V_rho']  
    
    randomSeed = initDict['randomSeed']  
    
    ''' Set random seed
    '''
    np.random.seed(randomSeed)
  
    
    ''' Construct auxiliary objects.

    '''
    numCovarsOut  = Y1_beta.shape[0]
    numCovarsCost = D_gamma.shape[0]
    
    U1V_cov      = U1V_rho*np.sqrt(U1_var)*np.sqrt(V_var)
    U0V_cov      = U0V_rho*np.sqrt(U0_var)*np.sqrt(V_var)
    
    
    ''' Read and check dataset, distribute entries.
    '''
    data = np.genfromtxt(initDict['fileName'], dtype = 'float')

    assert (_checkData(data, numAgents, numCovarsOut, numCovarsCost) == True)
    
    #Y = data[:,0]
    #D = data[:,1]
    
    X = data[:,2:(numCovarsOut + 2)]
    Z = data[:,-numCovarsCost:]
    
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
    
    
    ''' Check quality of simulated sample. 
    '''
    assert (np.all(np.isfinite(Y1)))
    assert (np.all(np.isfinite(Y0)))
    
    assert (np.all(np.isfinite(D)))
    
    assert (Y1.shape == (numAgents, ))
    assert (Y0.shape == (numAgents, ))

    assert (D.shape  == (numAgents, ))
    
    assert (Y1.dtype == 'float')
    assert (Y0.dtype == 'float')
    
    assert ((D.all() in [1.0, 0.0]))    
    
    ''' Compute average treatment effects using the simulation approach
    '''
    
    ATE = np.mean(Y1-Y0)
    
    Dmean   = np.mean(D)
    
    ATT = np.average(a=Y1-Y0,weights= D)
    
    ATUT = np.average(a=Y1-Y0,weights= 1-D)
    
    print 'D', D
    
    print "The average treatment effect is", ATE
    print "The average treatment effect on the treated is", ATT
    print "The average treatment effect on the untreated is", ATUT
    
    '''Check
    '''
    check = ATE-(ATT*Dmean+ATUT*(1-Dmean))
    assert (check < 0.00001)
    

''' Private Functions
'''
    
def _checkData(data, numAgents, numCovarsOut, numCovarsCost):
    ''' Basic checks for the data.
    
    '''
    # Antibugging.
    assert (isinstance(data, np.ndarray))
    assert (isinstance(numAgents, int))
    assert (numAgents > 0)
    assert (isinstance(numCovarsOut, int))
    assert (numCovarsOut > 0)
    assert (isinstance(numCovarsCost, int))
    assert (numCovarsCost > 0)
        
    # Checks.
    assert (data.shape == (numAgents, (1 + 1 + numCovarsOut +  numCovarsCost))) 
    assert (np.all(np.isfinite(data)))
    assert ((data[:,1].all() in [1.0, 0.0]))
    
    # Finishing.
    return True
    
    
''' Executable.
'''
if __name__ == '__main__':
    
    treatmentEffects()
