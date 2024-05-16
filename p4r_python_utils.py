#!/usr/bin/env python
# -*- coding: utf-8 -*-

## Import 
import logging
import os
import sys
logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
#handler.setFormatter(logging.Formatter('%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S'))
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)
pipe_replace = '-' # to replace | when writing files

def log_and_exit(code, rep): # temporary code to store the return status in file since it cannot be retrieve directly in shell Launch* when using the current version of ../bin/p4r (as of 05/14/2024)
	with open(os.path.join(rep, 'python_return_status'), 'w') as f:
		f.write(str(code))
	sys.exit(code)
    
