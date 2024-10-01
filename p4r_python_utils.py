#!/usr/bin/env python
# -*- coding: utf-8 -*-

## Import 
import logging
import os
import pandas as pd
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
	

def read_input_csv(cfg, file_name_key, **kwargs):
	file = os.path.join(cfg['inputpath'], cfg['csvfiles'][file_name_key])
	if not os.path.isfile(file):
		logger.error('File '+file+' does not exist. Use key inputpath in configuration file to specify input directory.')
		log_and_exit(2, cfg['path'])
	if 'csv_delim' in cfg.keys():
		kwargs.update({'sep' : cfg['csv_delim']})
	logger.info('Read file '+file)
	data = pd.read_csv(file, **kwargs)
	if data.index.has_duplicates:
		msg = 'Error: input csv '+file+' has some duplicated names '
		if 'index_col' in kwargs.keys():
			msg += 'in columns ['+', '.join(kwargs['index_col'])+']'
		msg += '. Duplicates are '+'\n\t'.join(str(i) for i in data.index[data.index.duplicated()])
		logger.error(msg)
		log_and_exit(3, cfg['path'])
	return data

def save_input_csv(cfg, file_name_key,data, index=False, **kwargs):
	file = os.path.join(cfg['inputpath'], cfg['csvfiles'][file_name_key])
	indexSave=0
	while os.path.isfile(os.path.join(cfg['inputpath'], cfg['csvfiles'][file_name_key]+'.save_'+str(indexSave)+'.csv')):
		indexSave=indexSave+1
	fileSave = os.path.join(cfg['inputpath'], cfg['csvfiles'][file_name_key]+'.save_'+str(indexSave)+'.csv')
	if 'csv_delim' in cfg.keys():
		kwargs.update({'sep' : cfg['csv_delim']})
	#data.to_csv(fileSave,**kwargs)
	data.to_csv(fileSave,index=index, **kwargs)
	
def write_input_csv(cfg, file_name_key,data, index=False, **kwargs):
	file = os.path.join(cfg['inputpath'], cfg['csvfiles'][file_name_key])
	indexSave=0
	if 'csv_delim' in cfg.keys():
		kwargs.update({'sep' : cfg['csv_delim']})
	#data.to_csv(file,**kwargs)
	data.to_csv(file,index=index, **kwargs)


def check_and_read_csv(cfg, file, **kwargs):
	if not os.path.isfile(file):
		logger.error('File '+file+' does not exist.')
		log_and_exit(2, cfg['path'])
	return pd.read_csv(file, **kwargs)