''' This module contains the required functionality for the processing and checking of an *.ini 
    file that characterizes the simulation and estimation of the Generalized Roy Model.
    
    Public Function:
    
        readInitFile
        
    Private Function:
    
        _checkIntegrity
    
'''

# standard library
import json

''' Public Functions
'''
def readInitFile():
    ''' This function reads the *.ini file, collects the information in a dictionary and outputs
        it to a *.json file.
    '''
    # Initialization
    dict_ = {}
    
    with open('init.ini', 'r') as file_:
        
        for count in range(30):
            
            str_ = file_.readline()
            
            # Number of Agents.
            if(count == 6):
                
                str_ = str_.split('!')[0]
                
                dict_['numAgents'] = int(str_)

            # File Name.    
            if(count == 7):
                
                str_ = str_.split("'")[1]

                dict_['fileName'] = str_    
            
            # Coefficients for Treated State .   
            if(count == 11):
                
                str_ = str_.split('!')[0]
                str_ = str_.replace('\t', '')
                str_ = str_.replace('\n', '')
                str_ = str_.split(',')
                
                for i in range(len(str_)):
                
                    str_[i] = float(str_[i])

                
                dict_['Y1_beta'] = list(str_)       
            
            # Coefficients for Untreated State.   
            if(count == 12):
                
                str_ = str_.split('!')[0]
                str_ = str_.replace('\t', '')
                str_ = str_.replace('\n', '')
                str_ = str_.split(',')
                
                for i in range(len(str_)):
                
                    str_[i] = float(str_[i])
                
                dict_['Y0_beta'] = list(str_)  
            
            # Coefficients in Cost Equation. 
            if(count == 16):
                
                str_ = str_.split('!')[0]
                str_ = str_.replace('\t', '')
                str_ = str_.replace('\n', '')
                str_ = str_.split(',')
                
                for i in range(len(str_)):
                
                    str_[i] = float(str_[i])
                
                dict_['D_gamma'] = list(str_)  
            
            # Variance of Disturbance in Treated State. 
            if(count == 20):
                
                str_ = str_.split('!')[0]
                str_ = str_.replace('\t', '')
                str_ = str_.replace('\n', '')
                str_ = str_.split(',')
                
                dict_['U1_var'] = float(str_[0])  
            
            # Variance of Disturbance in Untreated State.
            if(count == 21):
                
                str_ = str_.split('!')[0]
                str_ = str_.replace('\t', '')
                str_ = str_.replace('\n', '')
                str_ = str_.split(',')
                
                dict_['U0_var'] = float(str_[0])  
            
            # Variance of Disturbance in Cost Equation.
            if(count == 22):
                
                str_ = str_.split('!')[0]
                str_ = str_.replace('\t', '')
                str_ = str_.replace('\n', '')
                str_ = str_.split(',')
                
                dict_['V_var'] = float(str_[0])  
                            
            # Correlation Coefficients (U1, V)            
            if(count == 23):
                
                str_ = str_.split('!')[0]
                str_ = str_.replace('\t', '')
                str_ = str_.replace('\n', '')
                str_ = str_.split(',')
                
                dict_['U1V_rho'] = float(str_[0])  
 
            # Correlation Coefficients (U1, V)            
            if(count == 24):
                
                str_ = str_.split('!')[0]
                str_ = str_.replace('\t', '')
                str_ = str_.replace('\n', '')
                str_ = str_.split(',')
                
                dict_['U0V_rho'] = float(str_[0]) 
                                                                                   
    # Check quality.
    _checkIntegrity(dict_)
    
    #  Write out the *.json file.
    with open('init.json', 'w') as file_:
        
        json.dump(dict_, file_)
    
    # Finishing.
    return dict_

''' Private Functions
'''
def _checkIntegrity(dict_):
    ''' Check integrity of initFile dict.
    '''
    # Antibugging.
    assert (isinstance(dict_, dict))
    
    # Check number of agents.
    assert (dict_['numAgents'] > 0)
    assert (isinstance(dict_['numAgents'], int))
    
    # Finishing.
    return True
    
''' initFile Reader Script.
'''     
initDict = readInitFile()

_checkIntegrity(initDict)




