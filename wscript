#!/usr/bin/env python
''' Example of a waf script, that executes the unit tests defined in the tests
    directory. 
    
'''

# standard library
import os
import shutil
import sys
import fnmatch

# edit pythonpath
project_root = os.getcwd()
sys.path.insert(0, os.path.join(project_root, 'dev/py'))

# directories
top = '.'
out = '.bld'

def configure(conf):
    ''' Configuration.
    '''
    conf.env.project_paths = {}
    
    conf.env.project_paths['GRM_TOOLBOX'] = os.getcwd()

    tools_dir = conf.env.project_paths['GRM_TOOLBOX'] + '/tools'

    conf.load('runPyScript', tooldir = tools_dir)

def build(bld):
    ''' Build.
    '''
    bld.env.PROJECT_PATHS = set_project_paths(bld)
    
    bld.recurse('tests')

def distclean(ctx):
    ''' Clean.
    '''
    remove_filetypes_distclean('.')
    
    remove_for_distclean('.waf-1.6.4-8c7ad4bb8e1ca65b04e5d8dd9d0dac54')

    remove_for_distclean('.bld')

''' Auxillary functions.

'''
def remove_for_distclean(path):
    ''' Remove path, where path can be either a directory or a file. The
        appropriate function is selected. Note, however, that if an 
        OSError occurs, the function will just pass.
    '''
    if os.path.isdir(path):

        shutil.rmtree(path)
    
    if os.path.isfile(path):

        os.remove(path)

def remove_filetypes_distclean(path):
    ''' Remove nuisance files from the directory tree.
    '''
    matches = []

    for root, _, filenames in os.walk('.'):

        for filetypes in ['*.aux','*.out','*.log','*.pyc', '*~', \
            '.waf*', '*lock*', '*.mod', '*.a', '*.txt', '*.json', '*.dat']:

                for filename in fnmatch.filter(filenames, filetypes):
                    
                    matches.append(os.path.join(root, filename))

    for files in matches:

        remove_for_distclean(files)
        
def set_project_paths(ctx):
    ''' Return a dictionary with project paths represented by Waf nodes. This is
        required such that the run_py_script works as the whole PROJECT_ROOT is
        added to the Python path during execution.
    ''' 
    pp = {}

    pp['PROJECT_ROOT'] = '.'
   
    for key, val in pp.items():

        pp[key] = ctx.path.make_node(val)
   
    return pp