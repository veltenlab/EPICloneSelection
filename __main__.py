#!/usr/bin/env python
# coding: utf-8

# In[4]:


import argparse as ap
from os import os.path
import sys


# In[2]:


pars = ap.ArgumentParser()
pars.add_argument('-s','--sample_annotation',type=str,help='Path to the sample annotation sheet')
pars.add_argument('-p','--path',type=str,help="Path to the data directory")

args = pars.parse_args()


# In[ ]:


if not os.path.isfile(args.sample_annotation):
    sys.exit('Sample annotation sheet at ' + args.sample_annotation + 'does not exist')    
    
if not os.path.isdir(args.path):
    sys.exit('Data directory at ' + args.sample_annotation + 'does not exist')        

