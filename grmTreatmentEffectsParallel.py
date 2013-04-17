'''
    _____________________________________________________________
    
    This module is an extension of the Generalized Roy Toolbox
    by  Philipp Eisenhauer and Stefano Mosso (Copyright 2013)
    
    It was created to follow a coding structure similar to 
    the original toolbox.  
    
    The module generates simulated data (in parallel)
    and computes average treatment effects (ATE, TT, TUT)
    
    This is also the final exam in Computation Econometrics 
    by Sara Moreira
    ___________________________________________________________
    
    The public function treatmentEffects is Public interface to 
    request an estimation of TE using the Generalized Roy Model.

    To launch a MPI-parallel run with 10 processes:
         mpiexec -n 10 python grmTreatmentEffectsSM.py 

''' 

import os.path
import sys
import json
import random

import  numpy as np

from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
   
'''' Public Interface.
'''

def treatmentEffects(numAgentsSim=1000): 
    # "numAgentSim" is the number of observations simulated in each processor
    
    ''' Read estimated parameters
    '''
    assert (os.path.exists('grmRslt.json')) 
    resultsDict =  open('grmRslt.json').read()
    resultsDict =  json.loads(resultsDict)
       
    ''' Distribute parametrization and (limited) type conversions.
    '''
    numAgents  = resultsDict['numAgents']
    fileName   = resultsDict['fileName']
    
    Y1_beta    = np.array(resultsDict['Y1_beta'])
    Y0_beta    = np.array(resultsDict['Y0_beta'])
    
    D_gamma    = np.array(resultsDict['D_gamma'])
    
    U1_var     = resultsDict['U1_var'] 
    U0_var     = resultsDict['U0_var'] 
    V_var      = resultsDict['V_var']
    
    U1V_rho    = resultsDict['U1V_rho']  
    U0V_rho    = resultsDict['U0V_rho']  
    
    randomSeed = resultsDict['randomSeed']  
             
    ''' Set random seed.
    '''
    np.random.seed(randomSeed)
        
    '''Construct auxiliary objects.
    '''
    numCovarsOut  = Y1_beta.shape[0]
    numCovarsCost = D_gamma.shape[0]
        
    U1V_cov      = U1V_rho*np.sqrt(U1_var)*np.sqrt(V_var)
    U0V_cov      = U0V_rho*np.sqrt(U0_var)*np.sqrt(V_var)
        
    ''' Reading "observed" Data.
    '''
    data =  np.genfromtxt(resultsDict['fileName'], dtype = 'float')
    
    '''Debugging and Checks.
    '''
    assert (isinstance(data, np.ndarray))
    assert (isinstance(numAgents, int))
    assert (numAgents > 0)
    assert (isinstance(numCovarsOut, int))
    assert (numCovarsOut > 0)
    assert (isinstance(numCovarsCost, int))
    assert (numCovarsCost > 0)

    assert (data.shape == (numAgents, (1 + 1 + numCovarsOut +  numCovarsCost))) 
    assert (np.all(np.isfinite(data)))
    assert ((data[:,1].all() in [1.0, 0.0]))
    
         
    ''' Randomly draw individuals into different processors.
    '''
    processor  = random.choice  # choose a random element
        
    dataSim = [processor(data) for _ in xrange(numAgentsSim)]   # Create data
    dataSim = np.reshape(dataSim,(numAgentsSim,(1 + 1 + numCovarsOut +  numCovarsCost)))
    
    ''' Construct level indicators for outcomes and choices. 
    '''
    X = dataSim[:,2:(numCovarsOut + 2)]
    Z = dataSim[:,-numCovarsCost:]
    
    Y1_level = np.dot(Y1_beta, X.T)
    Y0_level = np.dot(Y0_beta, X.T)
    D_level  = np.dot(D_gamma, Z.T)
    
    ''' Simulate unobservables from the model using estimated parameters
    '''
    means = np.tile(0.0, 3)
    vars_ = [U1_var, U0_var, V_var]
    
    covs  = np.diag(vars_)
    
    covs[0,2] = U1V_cov 
    covs[2,0] = covs[0,2]
    
    covs[1,2] = U0V_cov
    covs[2,1] = covs[1,2]
    
    U = np.random.multivariate_normal(means, covs, numAgentsSim)
    
    ''' Simulate individual outcomes and choices.
    '''
    # Unobservables
    U1 = U[:,0]
    U0 = U[:,1]
    V  = U[:,2]
    
    # Potential outcomes.
    Y1 = Y1_level + U1
    Y0 = Y0_level + U0
    
    # Some calculations outside the loop
    expectedBenefits = Y1_level - Y0_level

    # Decision Rule.
    cost = D_level  + V     
    D = np.array((expectedBenefits - cost > 0), float)
        
    # Observed outcomes.
    Y  = D*Y1 + (1.0 - D)*Y0

    ''' Check quality of simulated sample.
    '''
    assert (np.all(np.isfinite(Y1)))
    assert (np.all(np.isfinite(Y0)))
    
    assert (np.all(np.isfinite(Y)))
    assert (np.all(np.isfinite(D)))
    
    assert (Y1.shape == (numAgentsSim, ))
    assert (Y0.shape == (numAgentsSim, ))
    
    assert (Y.shape  == (numAgentsSim, ))
    assert (D.shape  == (numAgentsSim, ))
    
    assert (Y1.dtype == 'float')
    assert (Y0.dtype == 'float')
    
    assert (Y.dtype == 'float')
    
    assert ((D.all() in [1.0, 0.0]))

    ''' Defining the average treatment effects
    '''
    Dmean = np.mean(D)
    D1sim    = np.int(np.sum(D))
    
    assert (0.0 <Dmean < 1.0) 
    
    ATEsim  = np.mean(a=Y1-Y0)
    ATTsim  = np.average(a=Y1-Y0, weights = D)
    ATUTsim = np.average(a=Y1-Y0, weights = 1-D)
    
    output = [numAgentsSim,D1sim,ATEsim,ATTsim,ATUTsim]
    
    ''' Consistency checks and quality checks
    '''
    check = ATEsim -(ATTsim*Dmean + ATUTsim*(1-Dmean))
    assert (check < 0.00001)
    
    assert(isinstance(D1sim,int))
    assert(isinstance(ATEsim,float))
    assert(isinstance(ATTsim,float))
    assert(isinstance(ATUTsim,float))
     
    ''' Exporting the individual estimates and weighting
    '''
    #results = comm.allgather(output)
    results = comm.gather(output, root=0)  
    
    if rank==0: # Processor 1 (rank 0) will compile the estimates
        
        results = np.asarray(results)
        np.savetxt('grmTESim.dat', np.column_stack((rank,numAgentsSim,D1sim,ATEsim,ATTsim,ATUTsim)), fmt= '%8.3f')
         
        ATE  = np.average(a=(results[:,2]), weights = (results[:,0]))
        ATT  = np.average(a=(results[:,3]), weights = (results[:,1]))
        ATUT = np.average(a=(results[:,4]), weights = ((results[:,0])-results[:,1]))
        
        Dict= {}
        Dict['ATE'] = ATE
        Dict['ATT']  = ATT
        Dict['ATUT'] = ATUT
        
        with open('grmTE.json', 'w') as file_:
           
            json.dump(Dict, file_)
            
           
''' Executable.
'''
if __name__ == '__main__':       

    treatmentEffects()
