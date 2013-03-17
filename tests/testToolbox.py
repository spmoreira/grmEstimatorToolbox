#!/usr/bin/env python
''' Testing module for the grmEstimatorToolbox.

'''
# standard library
import os
import sys

from    nose.core   import *
from    nose.tools  import *

# set working directory
dir_ = os.path.abspath(os.path.split(sys.argv[0])[0])
os.chdir(dir_)

# edit system path
sys.path.insert(0, dir_.replace('/tests', ''))

# project library
import grmToolbox

''' Auxiliary functions.
'''
class testingToolbox(object):
    ''' Testing the toolbox.
    '''
    def testOne(self):
        
        assert_true(True == True)

    def testTwo(self):
        
        assert_equal(0.00, 0.00)
                
if __name__ == '__main__':
    
    runmodule()   