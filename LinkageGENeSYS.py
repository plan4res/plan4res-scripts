#!/usr/bin/env python
# -*- coding: utf-8 -*-

## Import packages
import pyam  
import pandas as pd ## necessary data analysis package
import numpy as np
import os
import os.path as osp
import sys
import yaml
#import nomenclature as nc
from math import ceil,fabs
from datetime import timedelta
from calendar import monthrange

from p4r_python_utils import *

path = get_path()
logger.info('path='+path)
p4rpath = os.environ.get("PLAN4RESROOT")
logger.info('p4rpath='+p4rpath)

def abspath_to_relpath(path, basepath):
	return osp.relpath(path, basepath) #if osp.abspath(path) else path

def replace_scenario_names(version, replacement_dict):
	for short_name, long_name in replacement_dict.items():
		if short_name in version:
			return long_name
	return version

def check_required_inputs(cfg):
	logger.info('\nCheck required inputs')
	required_input = ['Sets', 'Par_OutputActivityRatio', 'Par_TechnologyFromStorage', 'Par_TagTechnologyToSector']
	required_output = []
	missing = False
	for var in cfg['variables']:
		if 'source' in cfg['variables'][var]:
			if cfg['variables'][var]['source']!='internal':
				# get data
				if cfg['variables'][var]['source']!='input':
					required_output.append(cfg['variables'][var]['source'])
				else:
					required_input+=cfg['variables'][var]['sheets']
	for s in set(required_input):
		if s not in cfg['genesys_datafiles']['input']['Sheets']:
			missing = True
			logger.warning(f'/!\ Required input sheet {s} is missing from parameter genesys_datafiles>input>Sheets. Please add it in settings file settingsLinkageGENeSYS.yml.')
	for d in set(required_output):
		if d not in cfg['genesys_datafiles']['output'].keys():
			missing = True
			logger.warning(f'/!\ required_output output data {s} is missing from parameter genesys_datafiles>output. Please add it in settings file settingsLinkageGENeSYS.yml.')
	if missing:
		log_and_exit(1,os.getcwd()) 
	

def get_techno_power(data, for_storage=False):
	elements = data['input_Sets']['Storage' if for_storage else 'Technology'].dropna()
	technos_ouput_power = data['input_Par_OutputActivityRatio'].set_index('Technology')
	technos_ouput_power = technos_ouput_power[(technos_ouput_power['Fuel']=='Power')&technos_ouput_power['Value'].abs()>0]
	technos_is_power = technos_ouput_power.index.unique()
	if for_storage:
		storage_techno = data.loc['input_Par_TechnologyFromStorage'].set_index('Storage').loc[elements, 'Technology'].unique()
		technos_is_power = [data.loc['input_Par_TechnologyToStorage'].set_index('Technology').loc[t, 'Storage'].unique()[0] for t in storage_techno if t in technos_is_power] # we keep storage is associated  techno is outputting power
	return technos_is_power
		
def interactive_check_if_set_is_in_mapping(sets, set_col_name, mappings, mapping_name, data, cfg, filter_on_power=True):
	logger.info(f'\nCheck if elements in set {set_col_name} are in settings file mappings {mapping_name}.')
	elements = sets[set_col_name].dropna()
	technos_ouput_power = data['input_Par_OutputActivityRatio'].set_index('Technology')
	technos_ouput_power = technos_ouput_power[(technos_ouput_power['Fuel']=='Power')&technos_ouput_power['Value'].abs()>0]
	technos_is_power = get_techno_power(data, mapping_name=='storages')
	technos_to_sector = data['input_Par_TagTechnologyToSector'].set_index('Technology')
	missing_in_mapping = [val for val in elements.dropna() if (not filter_on_power or val in technos_is_power) and val not in mappings.loc[mapping_name].dropna().index]
	present_in_mapping = [val for val in elements.dropna() if (not filter_on_power or val in technos_is_power) and val in mappings.loc[mapping_name].dropna().index]
	def display_and_check_mapping(m, data):
		logger.info(f'{m}:')
		for k,d in data.items():
			logger.info(f'\t{k} : {data[k]}')
		if set_col_name in genesys_mappings.keys():
			techno_iamc_settings = data['TechnoIAMC']
			if m in genesys_mappings[set_col_name].index:
				techno_genmap = genesys_mappings[set_col_name][m]
				if techno_genmap != techno_iamc_settings:
					logger.warning(f'Mapping in setting file ({techno_iamc_settings}) is different from ({techno_genmap}) defined in {genesys_mappings_dir}.')
					query=input('Continue creation of IAMC data? [y]/n\n')
					if query=='n':
						log_and_exit(1, os.getcwd())
								
								
	def try_to_find_matching_IAMC_techno(m):
		tentative = [cfg[b+'_'.join(m[1:])]['TechnoIAMC'] for b in ['RES_', 'P_','D_', 'S_'] if b+'_'.join(m[1:]) in cfg.keys()] # we first try to match to another techno inside settings file ignoring the GENeSYS-MOD prefix
		if set_col_name in genesys_mappings.keys(): # if additional genesys-mod mappings were provided we also try to match against it
			tentative += [genesys_mappings[set_col_name][b+'_'.join(m[1:])] for b in ['RES_', 'P_','D_', 'S_'] if b+'_'.join(m[1:]) in genesys_mappings[set_col_name].index]
		return tentative[0] if len(tentative) > 0 else 'UNDEFINED'
		
	for m in present_in_mapping:
		display_and_check_mapping(m, mappings.loc[mapping_name].loc[m])
	if len(missing_in_mapping) > 0:
		logger.warning(f'The following items from set {set_col_name} are missing from mapping {mapping_name}: '+', '.join(missing_in_mapping))
		query=input('Display a proposal for additional mapping? [y]/n\n')
		if query!='n':
			for m in missing_in_mapping:
				logger.info(f'{m}:')
				if mapping_name == 'technos':
					logger.info(f'\tEmission: Yes')
					logger.info(f'\tSector: {technos_to_sector.loc[m,"Sector"]}')
				elif mapping_name == 'storages':
					logger.info(f'\tStorageRatio: 100')
				logger.info(f'\tTechnoIAMC: {try_to_find_matching_IAMC_techno(m)}')
			
		query=input('Continue creation of IAMC data? [y]/n\n')
		if query=='n':
			log_and_exit(1, os.getcwd())

def estimate_storage_energy_capacity(data, cfg, filter_on_power=True):
	logger.info('\nAttempt to check storage capacity')
	storages = data.loc['input_Sets']['Storage'].dropna()
	storages_in_settings_mappings = cfg['StorageMappings'].keys()
	storage_techno = data.loc['input_Par_TechnologyFromStorage'].set_index('Storage').loc[storages, 'Technology'].unique()
	storage_techno = pd.Index(data.loc['capacity']['Technology'].unique()).intersection(storage_techno)
	inst_capa_storages = data.loc['capacity'].set_index(['Type','Technology','Region','Year',])['Value'].loc[('TotalCapacity', storage_techno)].droplevel(0).unstack(level=['Region','Year'])
	inst_capa_storages.index = [data.loc['input_Par_TechnologyToStorage'].set_index('Technology').loc[i,'Storage'].unique()[0] for i in inst_capa_storages.index]
	inst_capa_storages = inst_capa_storages.reindex(storages).fillna(0.)
	capacity_to_energy = pd.Series(0, index=storages)
	technos_is_power = get_techno_power(data, True)
	if 'input_Par_StorageE2PRatio' in data.index:
		capacity_to_energy.update(data.loc['input_Par_StorageE2PRatio'].set_index('Technology').loc[storage_techno])
	elif'input_Par_ResidualStorageCapacity' in data.index:
		residualstoragecapacity = data['input_Par_ResidualStorageCapacity'].set_index(['Storage','Region','Year'])['Value'].unstack(['Region','Year'])
		if not residualstoragecapacity.empty:
			capacity_to_energy.update(residualstoragecapacity.divide(inst_capa_storages).groupby(level='Region',axis=1).max())
	for s in storages:
		if filter_on_power and s not in technos_is_power:
			continue
		logger.info(f'Check energy capacity for storage {s}')
		capacity_to_energy_settings = None
		capacity_to_energy_data = None
		if s in capacity_to_energy.index:
			capacity_to_energy_data = capacity_to_energy.loc[s]
		if s in cfg['StorageMappings'].keys():
			if 'StorageRatio' in cfg['StorageMappings'][s].keys():
				capacity_to_energy_settings = cfg['StorageMappings'][s]['StorageRatio']
		capacity_to_energy_data_invalid = (capacity_to_energy_data is None) or (fabs(capacity_to_energy_data)<=1e-6)
		capacity_to_energy_settings_invalid = (capacity_to_energy_settings is None) or (fabs(capacity_to_energy_settings)<=1e-6)
		if capacity_to_energy_data is None:
			logger.warning(f'/!\ Could not estimate capacity to energy ratio for storage {s}')
		if fabs(capacity_to_energy_data)<=1e-6:
			logger.warning(f'/!\ Capacity to energy ratio for storage computed from data is 0')
		if capacity_to_energy_data_invalid:
			if capacity_to_energy_settings is None:
				logger.warning(f'/!\ No value provided in setting file at StorageMappings>{s}>StorageRatio.')
				query=input('Continue creation of IAMC data? [y]/n\n')
				if query=='n':
					log_and_exit(1, os.getcwd())
			if fabs(capacity_to_energy_settings)<=1e-6:
				logger.warning(f'/!\ Value in setting file at StorageMappings>{s}>StorageRatio is 0.')
				query=input('Continue creation of IAMC data? [y]/n\n')
				if query=='n':
					log_and_exit(1, os.getcwd())
			else:
				logger.info(f'Value provided in setting file at StorageMappings>{s}>StorageRatio is {cfg["StorageMappings"][s]["StorageRatio"]}')
	if capacity_to_energy_data_invalid and capacity_to_energy_settings_invalid:
		logger.warning(f'/!\ Storage ratio infered from input GENeSYS-MOD data {gen_val} is different from value {settings_val} indicated in settings file (parameter StorageMappings>{s}>StorageRatio).')
		query=input('Continue creation of IAMC data? [y]/n\n')
		if query=='n':
			log_and_exit(1, os.getcwd())
		
nbargs=len(sys.argv)	
if nbargs>1: 
	settings=sys.argv[1]
else:
	settings="settingsLinkageGENeSYS.yml"
if osp.abspath(settings):
	settings = osp.relpath(settings, path)

cfg={}
with open(osp.join(path, settings),"r") as mysettings:
	cfg=yaml.load(mysettings,Loader=yaml.FullLoader)

# replace name of current dataset by name given as input
if nbargs>1:
	namedataset=sys.argv[2]
	if 'path' in cfg:
		cfg['path']=cfg['path'].replace(cfg['path'].split('/')[len(cfg['path'].split('/'))-1],namedataset)
	else:
		cfg['path']=osp.join(path, 'data', namedataset)
	if 'p4rpath' in cfg:
		cfg['p4rpath']=cfg['p4rpath'].replace(cfg['path'].split('/')[len(cfg['path'].split('/'))-1],namedataset)
	else:
		cfg['p4rpath']=p4rpath
		
if 'configDir' not in cfg: cfg['configDir']=osp.join(cfg['path'], 'settings/')
if 'genesys_inputpath' not in cfg: cfg['genesys_inputpath']=osp.join(cfg['path'], 'GENeSYS-MOD/inputs/')
if 'genesys_resultspath' not in cfg: cfg['genesys_resultspath']=osp.join(cfg['path'], 'GENeSYS-MOD/outputs/')
if 'genesys_git' not in cfg: cfg['genesys_gitpath']=osp.join(cfg['path'], 'GENeSYS-MOD/GENeSYS_MOD.data/Data/')
if 'timeseriespath' not in cfg: cfg['timeseriespath']=osp.join(cfg['path'], 'TimeSeries/')
if 'mappingspath' not in cfg: cfg['mappingspath']=osp.join(cfg['path'], 'settings/mappings_genesys/')
if 'outputpath' not in cfg: cfg['outputpath']=osp.join(cfg['path'], 'IAMC/')
if 'outputfile' not in cfg: cfg['outputfile']=namedataset+'.csv'
if 'pythonDir' not in cfg: cfg['pythonDir']=osp.join(cfg['p4rpath'],'scripts/python/plan4res-scripts/settings/')
if 'nomenclatureDir' not in cfg: cfg['nomenclatureDir']=osp.join(cfg['p4rpath'],'scripts/python/openentrance/definitions/')

if not osp.isdir(cfg['outputpath']): os.mkdir(cfg['outputpath'])
if not osp.isdir(cfg['timeseriespath']): os.mkdir(cfg['timeseriespath'])
if osp.exists(cfg['outputfile']):
	os.remove(cfg['outputfile'])

if 'treat' in cfg:
	if 'fixed_data' in cfg['treat']:
		treatFix=cfg['treat']['fixed_data']
	else:
		treatFix=True
	if 'hourly_data' in cfg['treat']:
		treatHourly=cfg['treat']['hourly_data']
	else:
		treatHourly=True
else:
	treatFix=True
	treatHourly=True

CreateInputFromGit=False
if 'input_from_git' in cfg:
	if cfg['input_from_git']: CreateInputFromGit=True

if CreateInputFromGit:
	cfg['genesys_datafiles']['input']='InputDataGenesysFromGit.xlsx'
	if not osp.isdir(cfg['genesys_inputpath']):os.mkdir(cfg['genesys_inputpath'])

# To help define missing mappings in settings if necessary
genesys_mappings = dict()
genesys_mappings_dir = None
if 'genesys_mappings' in cfg.keys():
	if not osp.isdir(cfg['genesys_mappings']):
		logger.warning('Error: key "genesys_mappings" in configuration fail does not point to a valid directory ('+cfg['genesys_mappings']+')')
	else:
		genesys_mappings_dir = cfg['genesys_mappings']
		logger.info('Read GENeSYS-MOD to IAMC mappings in '+cfg['genesys_mappings'])
		logger.info('NB: these mappings do not supersede mappings defined in the settings file. They are only used to propose new mappings in case of missing data')
		genesys_mappings['Technology'] = pd.read_csv(osp.join(cfg['genesys_mappings'], 'capacity_technologies.csv'),index_col=0,header=None)[1]
		genesys_mappings['Storage'] = pd.read_csv(osp.join(cfg['genesys_mappings'], 'storages.csv'),index_col=0,header=None)[1]
		genesys_mappings['Emissions'] = pd.read_csv(osp.join(cfg['genesys_mappings'], 'emissions_technologies.csv'),index_col=0,header=None)[1]


if treatFix:
	logger.info('create IAMC file for GENeSYS-MOD outputs in '+cfg['outputpath'])
	check_required_inputs(cfg)

	# loop on the different variables
	BigOut=pd.DataFrame()
	firstVar=True

	# create dictionnary for replacing scenario names
	scenarios_names_dict = {}
	for main_name, alt_names in cfg['Scenarios'].items():
		for alt in alt_names:
			scenarios_names_dict[alt] = main_name

	# open datafiles
	data=pd.Series(index=[f'input_{elem}' for elem in cfg['genesys_datafiles']['input']['Sheets']])
	file=cfg['genesys_datafiles']['input']['inputfile']
	logger.info('read '+file)
	xls=pd.ExcelFile(osp.join(cfg['genesys_inputpath'],file),engine='openpyxl')
	scen = None
	if 'Scenario' in cfg:
		scen=cfg['Scenario']
	for sheet in cfg['genesys_datafiles']['input']['Sheets']:
		logger.info('  treat sheet '+sheet)
		if not sheet in xls.sheet_names:
			logger.warning(f'/!\ sheet {sheet} listed in genesys_datafiles>input>Sheets is absent from input data {file}.')
			data.drop('input_'+sheet, inplace=True)
			query=input('Continue creation of IAMC data? [y]/n\n')
			if query=='n':
				log_and_exit(1, os.getcwd())
			continue
		df=pd.read_excel(xls,sheet_name=sheet)
		if pd.isna(data['input_'+sheet]):
			data['input_'+sheet]=df
		else:
			data['input_'+sheet]=pd.concat(data['input_'+sheet],df)
		if 'Scenario' not in data['input_'+sheet].columns and 'PathwayScenario' not in data['input_'+sheet].columns:
			data['input_'+sheet]['Scenario']= cfg['genesys_datafiles']['input']['scenario'] if 'scenario' in cfg['genesys_datafiles']['input'].keys() else scen
			
		if 'PathwayScenario' in data['input_'+sheet].columns:
			data['input_'+sheet].rename(columns={'PathwayScenario': 'Scenario'}, inplace=True)


	if 'output' in cfg['genesys_datafiles']:
		for file in cfg['genesys_datafiles']['output']:
			logger.info('read '+osp.join(cfg['genesys_resultspath'],cfg['genesys_datafiles']['output'][file]))
			data.loc[file]=pd.read_csv(osp.join(cfg['genesys_resultspath'],cfg['genesys_datafiles']['output'][file]))
			if 'Scenario' not in data[file].columns and 'PathwayScenario' not in data[file].columns:
				if 'Model Version' in data[file].columns:					
					data[file]['Model Version'] = data[file]['Model Version'].str.split('_').str[1]
					data[file]['Model Version'] = data[file]['Model Version'].apply(replace_scenario_names, args=(scenarios_names_dict,))
					data[file].rename(columns={'Model Version': 'Scenario'}, inplace=True)
				elif 'Scenario' in cfg['genesys_datafiles']:
					data[file]['Scenario']=cfg['Scenario']
					data[file]['Scenario']=data[file]['Scenario'].apply(replace_scenario_names, args=(scenarios_names_dict,))
				else:
					data[file]['Scenario']='No Scenario'
				
			if 'PathwayScenario' in data[file].columns:
				data[file].rename(columns={'PathwayScenario': 'Scenario'}, inplace=True)
				data[file]['Scenario']=data[file]['Scenario'].apply(replace_scenario_names, args=(scenarios_names_dict,))
				
			if 'Scenario' in data[file].columns:
				data[file]['Scenario']=data[file]['Scenario'].apply(replace_scenario_names, args=(scenarios_names_dict,))


	# read mappings
	mappings=pd.Series()
	logger.info('read mappings')
	mappings.loc['technos']=pd.DataFrame([(tech, details['TechnoIAMC']) for tech, details in cfg['TechnosMappings'].items()],columns=['Technology', 'TechnoIAMC']).set_index('Technology')
	interactive_check_if_set_is_in_mapping(data.loc['input_Sets'], 'Technology', mappings, 'technos', data, cfg['TechnosMappings'])
	mappings.loc['technosfuels']=pd.DataFrame([(tech, details['TechnoIAMC']) for tech, details in cfg['TechnosMappings'].items()],columns=['Technology', 'TechnoIAMC']).set_index('Technology')
	mappings.loc['finalenergy_sector']=pd.DataFrame([(tech, details['Sector']) for tech, details in cfg['TechnosMappings'].items()],columns=['Technology', 'Sector']).set_index('Technology')
	mappings.loc['emissions']=pd.DataFrame([(tech, details['Emission']) for tech, details in cfg['TechnosMappings'].items()],columns=['Technology', 'Emission']).set_index('Technology')
	mappings.loc['storages_ratios']=pd.DataFrame([(details['TechnoIAMC'], details['StorageRatio']) for tech, details in cfg['StorageMappings'].items() ] ,columns=['TechnoIAMC', 'StorageRatio']).set_index('TechnoIAMC')
	mappings.loc['storages']=pd.DataFrame([(tech, details['TechnoIAMC']) for tech, details in cfg['StorageMappings'].items() ] ,columns=['Technology', 'TechnoIAMC']).set_index('Technology')
	interactive_check_if_set_is_in_mapping(data.loc['input_Sets'], 'Storage', mappings, 'storages', data, cfg['StorageMappings'])
	# Check value of storage ratios against values found in input data
	estimate_storage_energy_capacity(data, cfg)

	out=pd.DataFrame()
	isFirst=True
	IAMCcols=['Model','Scenario','Region','Variable','Unit','Year','Value']
	colsAgg=['Region','PathwayScenario','Year','Unit']
	
	regions=[]
	regions_source=data.loc['input_Sets']['Region'].dropna()
	for reg in regions_source:
		if reg not in regions and reg!=0:
			regions.append(str(reg))
	regions_interco=[]
	for region1 in regions:
		if region1!=cfg['global_region']:
			for region2 in regions:
				if region2!=cfg['global_region']:
					reg=str(region1)+'>'+str(region2)
					if reg not in regions_interco and region2!=region1:
						regions_interco.append(reg)
	logger.info('\n')
	logger.info('regions in dataset '+str(regions))
	logger.info('interco in dataset '+str(regions_interco))

	Yearsdf=pd.Series(data.loc['input_Sets']['Year']).dropna()
	Yearsdf=Yearsdf.drop(Yearsdf.loc[Yearsdf ==0].index,axis=0 ).astype(int)
	Years=Yearsdf.to_list()
	logger.info('years in dataset '+', '.join([str(y) for y in Years]))
	
	for var in cfg['variables']:
		debug=False
		if 'debug' in cfg:
			if var in cfg['debug']:
				debug=True
		isInternal=False
		logger.info('treat '+var)
		if debug: 
			print('\n treat '+var)
		
		if 'source' in cfg['variables'][var]:
			if cfg['variables'][var]['source']=='internal':
				isInternal=True
			else:
				# get data
				if cfg['variables'][var]['source']!='input':
					vardata=pd.DataFrame(data=data.loc[cfg['variables'][var]['source']])
				else:
					firstSheet=True
					for sheet in cfg['variables'][var]['sheets']:
						vardatasheet=pd.DataFrame(data=data.loc[cfg['variables'][var]['source']+'_'+sheet])
						if firstSheet: 
							vardata=pd.DataFrame(data=vardatasheet)
							firstSheet=False
						else:
							vardata=pd.concat([vardata,vardatasheet],axis=0)
		elif 'sources' in cfg['variables'][var]:
			if cfg['variables'][var]['sources']=='input':
				logger.error('input cannot be in multiple source')
				log_and_exit(1, os.getcwd())
			else:
				firstFile=True
				for file in cfg['variables'][var]['sources']:
					logger.info(' read '+file)
					vardatafile=pd.DataFrame(data=data.loc[file])
					#if 'Unit' not in vardatafile.columns: vardatafile['Unit']=cfg['variables'][var]['unit']
					vardatafile['Unit']=cfg['variables'][var]['unit']
						
					if firstFile:
						vardata=pd.DataFrame(data=vardatafile)
						firstFile=False
					else:					
						vardata=pd.concat([vardata,vardatafile],axis=0)
					
		colsdata=[]
		
		# treat case with 2 columns Region instead of Region and Region2
		if 'Region.1' in vardata.columns and 'Region2' not in vardata.columns:
			# rename the second Region column in Region
			vardata.rename(columns={'Region.1': 'Region2'}, inplace=True)
			
		if 'Region' in vardata.columns:
			vardata=vardata[ vardata['Region'].isin(regions) ]
		if 'Region2' in vardata.columns:
			vardata=vardata[ vardata['Region'].isin(regions) ]
				
		#if 'Unit' not in vardata.columns:
		vardata['Unit']=cfg['variables'][var]['unit']		
				
		# replace scenario nameserie
		if 'PathwayScenario' in vardata.columns:
			vardata['PathwayScenario']=vardata['PathwayScenario'].replace(scenarios_names_dict)
		vardata.rename(columns={'PathwayScenario': 'Scenario'}, inplace=True)

		# treat column names with space
		vardata.columns=vardata.columns.str.rstrip() 

		for rulecat in cfg['variables'][var]['rules']:
			logger.info('   apply '+rulecat)
			if debug: 
				print('\n apply '+rulecat)
			if rulecat=='selectAndMap':
				if debug: 
					print('\n before selectandmap')
					print(vardata.columns)
					print(vardata)
				# select rows 
				colmap=cfg['variables'][var]['rules'][rulecat]['column']
				if colmap not in colsdata: colsdata.append(colmap)
				firstMap=True
				strmaps=""
				for map in cfg['variables'][var]['rules'][rulecat]['mappings']:
					strmaps=strmaps+'_'+str(map)
					mappingpart=mappings.loc[map]
					if firstMap:
						fullmapping=mappingpart
						firstMap=False
					else:
						fullmapping=pd.concat([fullmapping,mappingpart],axis=0)
				vardata=vardata[ vardata[colmap].isin(list(fullmapping.index)) ]
				non_mapped_values = vardata[~vardata[colmap].isin(list(fullmapping.index))][colmap].unique()
				if len(non_mapped_values) >0:
					print("############")
					print(" --- The following values of column ",colmap," are not present in mappings ",strmaps)
					print(non_mapped_values)
					print("############")
				if debug: 
					print('\n in select andmap after  mapping')
					print(vardata.columns)
					print(vardata)
				# create variable name
				dict={fullmapping.index[i]: fullmapping.iloc[i,0] for i in range(len(fullmapping.index))}			
				vardata['Variable']=vardata[colmap].map(lambda a: dict[a])
				if debug: 
					print('\n in select andmap after  name change')
					print(vardata.columns)
					print(vardata)
				# compute variable
				ruleagg=str(cfg['variables'][var]['rules'][rulecat]['rule'])
				
				colsToAggr=[]			
				for coldata in vardata.columns:
					if coldata != 'Value' and coldata not in colsToAggr:
						colsToAggr.append(coldata)
				if 'Year' in vardata.columns:
					vardata['Year']=vardata['Year'].astype(int)
				if debug: 
					print('\n in select andmap before groupby')
					print('colsToAggr:',colsToAggr)
					print(vardata.columns)
					print(vardata)
				vardata=pd.DataFrame(data=pd.DataFrame(data=vardata).groupby(colsToAggr).agg(ruleagg).reset_index())
				if debug: 
					print('\n after select andmap')
					print(vardata)

				colsKeep=[]
				for coldata in vardata.columns:
					if coldata in IAMCcols:
						colsKeep.append(coldata)
				if debug: 
					print('colsKeep:',colsKeep)
				vardata=vardata[ colsKeep ]
				if debug: 
					print('\n after select andmap')
					print(vardata)

			elif rulecat=='addyear':
				firstYear=True
				for year in Years:
					vardatayear=pd.DataFrame(data=vardata)
					vardatayear['Year']=year				
					if firstYear:
						vardataout=vardatayear
						firstYear=False
					else:
						vardataout=pd.concat([vardataout,vardatayear],axis=0)
				vardata=vardataout
			
			elif rulecat=='apply_abs':
				if debug: print(vardata)
				vardata['Value']=vardata['Value'].abs()
			
			elif rulecat=='selectFromMapping':
				# select rows 
				col=cfg['variables'][var]['rules'][rulecat]['column']
				firstMap=True
				for map in cfg['variables'][var]['rules'][rulecat]['mappings']:
					vardatamap=pd.DataFrame(data=vardata[ vardata[col].isin(list(mappings.loc[map].index)) ])
					if firstMap:
						vardataout=vardatamap
						firstMap=False
					else:
						vardataout=pd.concat([vardataout,vardatamap],axis=0)
				vardata=vardataout
			
			elif rulecat=='map':
				colmap=cfg['variables'][var]['rules'][rulecat]['column']
				map=cfg['variables'][var]['rules'][rulecat]['mapping']
				if colmap not in colsdata: colsdata.append(colmap)
				
				# map variable name
				dict={mappings.loc[map].index[i]: mappings.loc[map].iloc[i,0] for i in range(len(mappings.loc[map].index))}
				
				vardata['Variable']=vardata[colmap].map(lambda a: dict[a] if a in dict.keys() else 'None')
				vardata=vardata.drop( vardata[vardata['Variable']=='None'].index  )

				# compute variable
				ruleagg=str(cfg['variables'][var]['rules'][rulecat]['rule'])
				
				colsKeep=[]
				for col in vardata.columns:
					if col in IAMCcols+colsdata:
						colsKeep.append(col)
				vardata=vardata[ colsKeep ]

				colsToAggr=[]			
				for coldata in vardata.columns:
					if coldata != 'Value' and coldata not in colsToAggr:
						colsToAggr.append(coldata)
				vardata=vardata.groupby(colsToAggr).agg(ruleagg).reset_index()
			
			elif rulecat=='select':
				if debug: 
					print('\n before select')
					print(vardata)
				for colselect in cfg['variables'][var]['rules'][rulecat]:
					values=cfg['variables'][var]['rules'][rulecat][colselect]['values']
					vardata=vardata[ vardata[colselect].isin(values) ]
				if debug: 
					print('\n after select')
					print(vardata)
					
			elif rulecat=='group':
				if debug: 
					print('\n before group')
					print(vardata.columns)
					print(vardata)
				ruleagg=str(cfg['variables'][var]['rules'][rulecat]['rule'])
				colsKeep=[]
				for col in vardata.columns:
					if col in IAMCcols:
						colsKeep.append(col)
				vardata=vardata[ colsKeep ]
				colsToAggr=[]			
				for coldata in vardata.columns:
					if coldata != 'Value' and coldata not in colsToAggr:
						colsToAggr.append(coldata)
				vardata=vardata.groupby(colsToAggr).agg(ruleagg).reset_index()
				if debug: 
					print('\n after group')
					print(vardata.columns)
					print(vardata)
			
			elif rulecat=='addvariablecol':
				vardata['Variable']=var
			
			elif rulecat=='concatvariablename':
				vardata['startVar']=var
				vardata['Variable']=vardata['startVar'].str.cat(vardata['Variable'])
				vardata=vardata.drop(['startVar'],axis=1)
				
			elif rulecat=='complete_variable_name':
				completion=cfg['variables'][var]['rules'][rulecat]
				vardata['endVar']=completion
				vardata['Variable']=vardata['Variable'].str.cat(vardata['endVar'])
				vardata=vardata.drop(['endVar'],axis=1)
			
			elif rulecat=='combineWithOtherSources':
				for subrule in cfg['variables'][var]['rules'][rulecat]:
					logger.info('		apply '+subrule)
					if 'source' in cfg['variables'][var]['rules'][rulecat][subrule]:
						if cfg['variables'][var]['rules'][rulecat][subrule]['source']!='input':
							newdata=data.loc[cfg['variables'][var]['rules'][rulecat][subrule]['source']]
						else:
							newdata=data.loc[cfg['variables'][var]['rules'][rulecat][subrule]['source']+'_'+cfg['variables'][var]['rules'][rulecat][subrule]['sheet']]
					if 'select' in cfg['variables'][var]['rules'][rulecat][subrule]:
						for colselect in cfg['variables'][var]['rules'][rulecat][subrule]['select']:
							values=cfg['variables'][var]['rules'][rulecat][subrule]['select'][colselect]['values']
							newdata=newdata[ newdata[colselect].isin(values) ]
					if subrule=='SelectAndMap':
						# select rows 
						colmap=cfg['variables'][var]['rules'][rulecat][subrule]['column']
						firstMap=True
						strmaps=""
						for map in cfg['variables'][var]['rules'][rulecat][subrule]['mappings']:
							strmaps=strmaps+'_'+str(map)
							mappingpart=mappings.loc[map]
							if firstMap:
								fullmapping=mappingpart
								firstMap=False
							else:
								fullmapping=pd.concat([fullmapping,mappingpart],axis=0)
						newdata=newdata[ newdata[colmap].isin(list(fullmapping.index)) ]
						non_mapped_values = newdata[~newdata[colmap].isin(list(fullmapping.index))][colmap].unique()
						if len(non_mapped_values) >0:
							print(" --- The following values of column ",colmap," are not present in mappings ",strmaps)
							print(non_mapped_values)
						dict={fullmapping.index[i]: fullmapping.iloc[i,0] for i in range(len(fullmapping.index))}			
						newdata['Variable']=newdata[colmap].map(lambda a: dict[a])
						ruleagg=str(cfg['variables'][var]['rules'][rulecat][subrule]['rule'])
						colsToAggr=[]			
						for coldata in newdata.columns:
							if coldata != 'Value' and coldata not in colsToAggr:
								colsToAggr.append(coldata)
						if 'Year' in newdata.columns:
							newdata['Year']=newdata['Year'].astype(int)
						newdata=pd.DataFrame(data=pd.DataFrame(data=newdata).groupby(colsToAggr).agg(ruleagg).reset_index())
						if debug:
							print('compute from other sources, dubrule select and map')
							print(newdata)
					elif subrule=='multiply':
						if 'mapping' in cfg['variables'][var]['rules'][rulecat][subrule]:							
							map=cfg['variables'][var]['rules'][rulecat][subrule]['mapping']
							values_mult=mappings.loc[map]
							colval=cfg['variables'][var]['rules'][rulecat][subrule]['value']
							for index in newdata.index:
								newdata.loc[index,'Value']=newdata.loc[index,'Value']*values_mult.loc[newdata.loc[index,colmap],colval]
						if debug:
							print('compute from other sources, dubrule multiply')
							print(newdata)
					elif subrule=='mapAndAddCols':					
						colref=cfg['variables'][var]['rules'][rulecat][subrule]['column']
						if debug: print('colref:', colref)
						for newcol in cfg['variables'][var]['rules'][rulecat][subrule]['mappings']:
							colmap=cfg['variables'][var]['rules'][rulecat][subrule]['mappings'][newcol]
							combinedmap=newdata[[colref,colmap]].groupby([colref]).first().reset_index()
							combineddict={combinedmap.iloc[i,0]: combinedmap.iloc[i,1] for i in range(len(combinedmap.index))}
							vardata[newcol]=vardata[colref].map(lambda a: combineddict[a] if a in combineddict.keys() else 'None')
							if debug:
								print(' combineothersources/mapaddcols/mappings row:',newcol)
								print('combinedmap')
								print(combinedmap)
								print('combineddict')
								print(combineddict)
								print('vardata[',newcol,']')
								print(vardata[newcol])
						if 'product_cols' in cfg['variables'][var]['rules'][rulecat][subrule]:
							if debug: print(' apply product_cols')
							if debug: print(vardata)
							for col in cfg['variables'][var]['rules'][rulecat][subrule]['product_cols']: 
								col2=cfg['variables'][var]['rules'][rulecat][subrule]['product_cols'][col]
								if debug: 
									print( ' 	product by ',col2)
									print( ' vardata[',col,'] before product')
									print(vardata[col])
									print( 'multiplied by:')
									print( vardata[cfg['variables'][var]['rules'][rulecat][subrule]['product_cols'][col]])
									for i in vardata[col].index:
										print(i,vardata['Technology'].loc[i],vardata['Year'].loc[i],vardata[col].loc[i],vardata[cfg['variables'][var]['rules'][rulecat][subrule]['product_cols'][col]].loc[i])
									vardata[col]=vardata[col].astype(float)*vardata[cfg['variables'][var]['rules'][rulecat][subrule]['product_cols'][col]].astype(float)
									print(var)
									print(rulecat)
									print(subrule)
									print(col)
									print(cfg['variables'][var]['rules'][rulecat][subrule]['product_cols'][col])
									print(vardata[col])
									print(vardata[cfg['variables'][var]['rules'][rulecat][subrule]['product_cols'][col]])
								if debug: print(vardata[col])
					elif subrule=='changeValue':
						colref=cfg['variables'][var]['rules'][rulecat][subrule]['column']
						colval=cfg['variables'][var]['rules'][rulecat][subrule]['value']
						colmap=cfg['variables'][var]['rules'][rulecat][subrule]['map']
						newvalue=newdata[['Value',colmap]].groupby([colref]).first().reset_index()
						valuedict={newvalue.iloc[i,0]: newvalue.iloc[i,1] for i in range(len(newvalue.index))}
						rows_to_remove=[]
						if cfg['variables'][var]['rules'][rulecat][subrule]['rule']=='mult':
							for row in vardata.index:
								if vardata.loc[row,colmap] in valuedict.keys():
									vardata.loc[row,'Value']=vardata.loc[row,'Value']*valuedict[vardata.loc[row,colmap]]					
								else:
									# remove row
									rows_to_remove.append(row)
						vardata=vardata.drop(rows_to_remove,axis=0)
					elif subrule=='group':
						ruleagg=cfg['variables'][var]['rules'][rulecat][subrule]['rule']
						colsKeep=[]
						for col in vardata.columns:
							if col in IAMCcols:
								colsKeep.append(col)
						vardata=vardata[ colsKeep ]
						colsToAggr=[]			
						for coldata in vardata.columns:
							if coldata != 'Value' and coldata not in colsToAggr:
								colsToAggr.append(coldata)
						vardata=vardata.groupby(colsToAggr).agg(ruleagg).reset_index()
			
			elif rulecat=='convert_unit':
				vardata['Value']=vardata['Value']*cfg['variables'][var]['rules'][rulecat]['factor']
				vardata['Unit']=cfg['variables'][var]['rules'][rulecat]['to']
			elif rulecat=='compute':
				if isInternal:
					if 'mapping' in cfg['variables'][var]['rules'][rulecat]:
						map=mappings.loc[cfg['variables'][var]['rules'][rulecat]['mapping']]
						dict={map.index[i]: map.iloc[i,0] for i in range(len(map.index))}
						serieElems=pd.Series([[] for _ in range(len(map.index))], index=map.index)
					listComponents=[]
					isManyVar=False
					listElem=[]
					firstComponent=True
					
					for component in cfg['variables'][var]['rules'][rulecat]['from']:
						if component[-1]=='|':
							isManyVar=True
							# add mapping list to variable name
							for elem in map.index:
								if firstComponent: 
									if not elem in listElem: listElem.append(elem)
								if 'ruleaggr' in cfg['variables'][var]['rules']['compute']:
									if component+dict[elem] not in listComponents:
										listComponents.append(component+dict[elem])
									if component+dict[elem] not in serieElems[elem]:
										serieElems[elem].append(component+dict[elem])
								else:
									if component+elem not in listComponents:
										listComponents.append(component+elem)
						else:
							listComponents.append(component)
						firstComponent=False
					vardata=out[ out['Variable'].isin(listComponents) ]	
					colsKeep=[]
					for col in vardata.columns:
						if col in IAMCcols:
							colsKeep.append(col)
					vardata=vardata[ colsKeep ]
					if 'ruleaggr' in cfg['variables'][var]['rules']['compute']:
						ruleagg=cfg['variables'][var]['rules']['compute']['ruleaggr']
						if isManyVar:
							firstElem=True
							# for elem in listElem:
								# listpossible=[el+str(dict[elem]) for el in cfg['variables'][var]['rules'][rulecat]['from']]
								# vardataelem = pd.DataFrame(vardata[vardata['Variable'].isin(listpossible)]).reset_index().drop(columns='index')
								# colsToAggr=[]
								# vardataelem=vardataelem.drop(columns='Variable')
								# for col in vardataelem.columns:
									# if col != 'Value':
										# colsToAggr.append(col)
								# if len(vardataelem.index)>0:
									# vardataelem=vardataelem.groupby(colsToAggr).agg(ruleagg).reset_index()
								# vardataelem['Variable']=var+dict[elem]
								# if firstElem: 
									# if len(vardataelem.index)>0:
										# vardatanew=pd.DataFrame(vardataelem)
										# firstElem=False
								# else:
									# if len(vardataelem.index)>0:
										# vardatanew = pd.concat([vardatanew, vardataelem],ignore_index=True)
							# vardata=pd.DataFrame(vardatanew)							
							for elem in map.index:
								vardataelem = pd.DataFrame(vardata[vardata['Variable'].isin( serieElems[elem] )]).reset_index().drop(columns='index')
								colsToAggr=[]
								if 'PHS' in elem: 
									if debug:
										print(elem)
										print('phs 1')
										print(vardataelem)
								vardataelem=vardataelem.drop(columns='Variable')
								for col in vardataelem.columns:
									if col != 'Value':
										colsToAggr.append(col)
								if 'PHS' in elem: 
									if debug:
										print('phs 2')
										print(vardataelem)
								if len(vardataelem.index)>0:
									vardataelem=vardataelem.groupby(colsToAggr).agg(ruleagg).reset_index()
								if 'PHS' in elem: 
									if debug:
										print('phs 3')
										print(vardataelem)
								vardataelem['Variable']=var+dict[elem]
								if firstElem:
									if len(vardataelem.index)>0:
										vardatanew=pd.DataFrame(vardataelem)
										firstElem=False
								else:
									if len(vardataelem.index)>0:
										vardatanew = pd.concat([vardatanew, vardataelem],ignore_index=True)
							vardata=pd.DataFrame(vardatanew)
										
						else:
							colsToAggr=[] 
							if 'Variable' in vardata.columns: vardata=vardata.drop(columns='Variable')
							for col in vardata.columns:
								if col != 'Value':
									colsToAggr.append(col)
							vardata=vardata.groupby(colsToAggr).agg(ruleagg).reset_index()
							vardata['Variable']=var

					elif 'rulemap' in cfg['variables'][var]['rules']['compute']:
						for row in vardata.index:
							for componentfrom in cfg['variables'][var]['rules']['compute']['from']:
								if componentfrom in vardata.loc[row,'Variable']:
									if cfg['variables'][var]['rules']['compute']['rulemap']=='mult':
										vardata.loc[row,'Value']=vardata.loc[row,'Value']*dict[vardata.loc[row,'Variable'].replace(componentfrom,'')]
									vardata.loc[row,'Variable']=vardata.loc[row,'Variable'].replace(componentfrom,var)							
				if debug:
					print('end compute')
					print(vardata)
					print(vardata['Variable'].unique())
			elif rulecat=='create_interco':	
				vardata['>']='>'
				vardata['Region']=vardata['Region'].str.cat(vardata['>']).str.cat(vardata['Region2'])
			
			elif rulecat=='global':
				globaldata=pd.DataFrame(data=vardata)
				globalreg=cfg['global_region']
				isFirstRegion=True
				# case of interconnection variable
				regions_use=[globalreg]
				if 'Network' in var:				
					regions_use=regions_interco
				for region in regions_use:
					globaldata['Region']=region
					if isFirstRegion:
						vardataout=pd.DataFrame(data=globaldata)
						isFirstRegion=False
					else:
						vardataout=pd.concat([vardataout,globaldata],axis=0,ignore_index=True)
				vardata=vardataout
		
		if debug:
			print('after rules')
			print(vardata)

		if not vardata.empty:
			if 'Year' not in vardata.columns:
				firstYear=True
				for year in Years:
					if year in vardata.columns:
						if firstYear:
							vardata['Value']=vardata[year]
							firstYear=False
						else:
							vardatayear=pd.DataFrame(data=vardata)
							vardatayear['Value']=vardatayear[year]
							vardata=pd.concat([vardata,vardatayear],axis=0)


			if debug:
				print('Selected years = ', str(Years))
				print('after year')
				print(vardata)
				
			#fill scenario
			if 'PathwayScenario' in vardata.columns:
				vardata['Scenario']=vardata['PathwayScenario']
			
			if debug:
				print('after scenario')
				print(vardata)
				
			if 'Unit' not in vardata.columns:
				vardata['Unit']=cfg['variables'][var]['unit']
			if debug:
				print('after unit')
				print(vardata)
				
			#select columns to keep
			colsKeep=[]
			for col in vardata.columns:
				if col in IAMCcols:
					colsKeep.append(col)
			vardata=vardata[ colsKeep ]
			
			if debug:
				print('after colskeep')
				print(vardata)

			#fill missing columns
			for col in IAMCcols:
				if col not in colsKeep:
					if col == 'Model': 
						vardata[col]=cfg['Model']
					elif col == 'Scenario':
						vardata[col]=vardata[col].map(scenarios_names_dict)
					
			if debug:
				print('after misscols')
				print(vardata)

			vardata['Year']=vardata['Year'].astype(int)
			
			if isFirst:
				out=vardata
				isFirst=False
			else:
				out=pd.concat([out,vardata],axis=0,ignore_index=True)
		
		else: logger.info('empty data')
	
	# change country names from iso2 to real names
	with open(osp.join(cfg['nomenclatureDir'], "region/countries.yaml"),"r",encoding='UTF-8') as nutsreg:
		countries=yaml.safe_load(nutsreg)
	iso2_to_country = {}
	for country in countries[0]['Countries']:
		for country_name, details in country.items():
			iso2_to_country[details['iso2']] = country_name
			if 'iso2_alt' in details:
				iso2_to_country[details['iso2_alt']] = country_name
	for reg in cfg['Add_Map_Region']:
		iso2_to_country[reg]=cfg['Add_Map_Region'][reg]
	iso2_to_country2={}
	for reg1 in iso2_to_country:
		for reg2 in iso2_to_country:
			iso2_to_country2[reg1+'>'+reg2]=iso2_to_country[reg1]+'>'+iso2_to_country[reg2]
			iso2_to_country2[reg2+'>'+reg1]=iso2_to_country[reg2]+'>'+iso2_to_country[reg1]
	iso2_to_country={**iso2_to_country, **iso2_to_country2}
	out['Region'] = out['Region'].map(lambda _ : _ if _ not in iso2_to_country.keys() else _)
	out.to_csv(cfg['outputpath']+cfg['outputfile'],index=False)	

	# check for duplicated and output synthesis of data
	logger.info('scenarios in data '+', '.join([str(_) for _ in out['Scenario'].unique()]))
	logger.info('models in data '+', '.join([str(_) for _ in out['Model'].unique()]))
	logger.info('regions in data '+', '.join([str(_) for _ in out['Region'].unique()]))
	#logger.info('variables in data '+', '.join([str(_) for _ in out['Variable'].unique()]))
	logger.info('years in data '+', '.join([str(_) for _ in out['Year'].unique()]))


	duplicates=out.duplicated()
	duprows=[]
	for row in duplicates.index:
		if duplicates.loc[row]==True:
			duprows.append(row)
	if len(duprows)>0:
		logger.warning('there are duplicated rows')
		logger.warning(', '.join([str(_) for _ in duprows]))
		duplicated_rows=out.loc[duprows]
		logger.warning('  for variables'+', '.join([str(_) for _ in duplicated_rows['Variable'].unique()]))
	else:
		logger.info('no duplicated rows')
		
	df=pd.read_csv(cfg['outputpath']+cfg['outputfile'],index_col=0)
	df['Scenario']=df['Scenario'].apply(replace_scenario_names, args=(scenarios_names_dict,))

	 
	# conversion to IAMDataFrame
	BigIAM=pyam.IamDataFrame(df)
	 
	#filter on unwanted variables
	def filter_variable_list(variable_list):
		removed_variables = [] if 'removed_variables' not in cfg.keys() else cfg['removed_variables']
		filtered_list = []
		for variable in variable_list:
			exclude = False
			for removed_var in removed_variables:
				if removed_var.endswith('|'):
					if variable.startswith(removed_var[:-1]):
						exclude = True
						break
				elif variable == removed_var:
					exclude = True
					break
			if exclude:
				filtered_list.append(variable)
		return filtered_list
	
	logger.info('list of all variables')	
	variable_list=list(df['Variable'].unique())
	print(variable_list)
	#new_variable_list=[item for item in variable_list if item not in cfg['removed_variables']]
	removed_variable_list=filter_variable_list(variable_list)
	logger.info('filtering on variables')
	print(removed_variable_list)
	#logger.info('excluding: '.join([str(_) for _ in removed_variable_list]))
	BigIAM.to_excel(cfg['outputpath']+'allvar_'+cfg['outputfile'].replace('csv','xlsx'))

	BigIAM=BigIAM.filter(variable=removed_variable_list, keep=False)

	#filter on unwanted variables
	logger.info('validating')
	BigIAM.validate(exclude_on_fail=True)
	BigIAM.to_excel(cfg['outputpath']+cfg['outputfile'].replace('csv','xlsx'))

if treatHourly:
	logger.info('treat hourly data')
	# dates treatments
	dates=pd.Series()
	beginTS=pd.to_datetime(cfg['Calendar']['BeginTimeSeries'],dayfirst=cfg['Calendar']['dayfirst'])
	endTS=pd.to_datetime(cfg['Calendar']['EndTimeSeries'],dayfirst=cfg['Calendar']['dayfirst'])
	dates['BeginTS']=pd.Timestamp(year=beginTS.year,month=beginTS.month,day=beginTS.day,hour=beginTS.hour,minute=beginTS.minute)
	dates['EndTS']=pd.Timestamp(year=endTS.year,month=endTS.month,day=endTS.day,hour=endTS.hour,minute=endTS.minute)
	DurationTimeSeries=pd.Timedelta(dates['EndTS']-dates['BeginTS'])
	TimeStep=cfg['Calendar']['TimeStep']['Duration']
	UnitTimeStep=cfg['Calendar']['TimeStep']['Unit']
	if UnitTimeStep=='days': TimeStep=TimeStep*24
	if UnitTimeStep=='weeks': TimeStep=TimeStep*168
	NumberTimeSteps=int((DurationTimeSeries.days*24+DurationTimeSeries.seconds/3600)/TimeStep)
	durationTimeStep=pd.Timedelta(str(TimeStep)+' hours')
	logger.info('dates: timeseries start: '+str(dates['BeginTS'])+' end: '+str(dates['EndTS']))
	logger.info('Duration timeseries:'+str(DurationTimeSeries))
	logger.info('Number of time steps:'+str(NumberTimeSteps)+' of duration:'+str(durationTimeStep))
	datesTS=pd.DataFrame(index=list(range(NumberTimeSteps)),columns=['start','end'])
	start=dates['BeginTS']
	for i in range(NumberTimeSteps):
		datesTS.loc[i]=[start,start+durationTimeStep]
		start=start+durationTimeStep
	#TimeSeriesTemplate=pd.DataFrame(columns=['Timestamp [UTC]'])
	TimeSeriesTemplate=pd.read_csv(osp.join(cfg['timeseriespath'],'Example.csv'))
	# start=dates['BeginTS']
	# i=0
	# while start<=dates['EndTS']:
		# TimeSeriesTemplate.loc[i]=[start]
		# start=start+durationTimeStep
		# i=i+1
	
	# read genesys-mod timeseries and create plan4res timeseries
	NumberScenarios=1+len(cfg['AdditionnalScenarios'])
	AddScenarios=[elem for elem in cfg['AdditionnalScenarios']]
	Scenarios=['Base']+AddScenarios
	logger.info('Scenarios:')
	logger.info(Scenarios) 
	for sheet in cfg['genesys_datafiles']['timeseries']['sheets']:
		sheetname=cfg['genesys_datafiles']['timeseries']['sheets'][sheet]
		logger.info('  sheet '+sheetname)
		df=pd.read_excel(osp.join(cfg['genesys_inputpath'],cfg['genesys_datafiles']['timeseries']['xlsx']),sheet_name=sheetname,index_col=0).fillna(0)
		df=df.reset_index()
		if sheetname in cfg['TimeSeriesFactor']:
			multfactor=(1/cfg['TimeSeriesFactor'][sheetname])
		else:
			multfactor=1.0
		logger.info('  sheetname '+str(sheetname)+' mult '+str(multfactor)) 
		# create plan4res time series related to variable sheetname
		for region in df.columns:
			if not region=='HOUR':
				timeseries = pd.DataFrame({'Timestamp [UTC]': TimeSeriesTemplate['Timestamp [UTC]'],'Base':df[region]*multfactor})
				#timeseries=pd.DataFrame(TimeSeriesTemplate['Timestamp [UTC]'])
				#timeseries['Base']=df[region]
				for scenario in cfg['AdditionnalScenarios']:
					if sheet in cfg['AdditionnalScenarios'][scenario]:
						timeseries[scenario]=timeseries['Base']*cfg['AdditionnalScenarios'][scenario][sheet]
					else:
						timeseries[scenario]=timeseries['Base']
				nameserie=sheetname+'_'+region+'.csv'
				timeseries.to_csv(osp.join(cfg['timeseriespath'],nameserie),index=False)

logger.info('Completed')
log_and_exit(0, os.getcwd())