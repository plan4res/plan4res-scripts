#!/usr/bin/env python
# -*- coding: utf-8 -*-

## Import packages
import pyam
import pandas as pd ## necessary data analysis package
import numpy as np
import os
import yaml
import calendar
import math
from math import ceil
from datetime import timedelta
from calendar import monthrange
from itertools import product
import sys

from p4r_python_utils import *

def find_entry(variable_name, dictionary):
	for entry, variables in dictionary.items():
		for item in variables:
			if variable_name.startswith(item):
				return entry
	return None

def search_variable(variable_name, data):
	entry = find_entry(variable_name, data.get('coupling', {}))
	if entry:
		return entry
	for sub_dict in data.get('techno', {}).values():
		entry = find_entry(variable_name, sub_dict)
		if entry:
			return entry
	return None


path = get_path()
logger.info('path='+path)
p4rpath = os.environ.get("PLAN4RESROOT")
logger.info('p4rpath='+p4rpath)

def abspath_to_relpath(path, basepath):
	return os.path.relpath(path, basepath) #if os.path.abspath(path) else path
	
nbargs=len(sys.argv)
if nbargs>1: 
	settings_create=sys.argv[1]
else:
	settings_create="settingsCreateInputPlan4res.yml"
if os.path.abspath(settings_create):
	settings_create = os.path.relpath(settings_create, path)

cfg={}
# open the configuration file using the pathway defined below
with open(os.path.join(path, settings_create),"r") as mysettings:
	cfg=yaml.load(mysettings,Loader=yaml.FullLoader)
	
# replace name of current dataset by name given as input
if nbargs>2:
	namedataset=sys.argv[2]
	if 'path' in cfg:
		cfg['path']=cfg['path'].replace(cfg['path'].split('/')[len(cfg['path'].split('/'))-1],namedataset)
	else:
		cfg['path']=os.path.join(path, 'data', namedataset)
	if 'p4rpath' in cfg:
		cfg['p4rpath']=cfg['p4rpath'].replace(cfg['path'].split('/')[len(cfg['path'].split('/'))-1],namedataset)
	else:
		cfg['p4rpath']=p4rpath
if 'outputpath' not in cfg: 
	if cfg['ParametersCreate']['invest']:
		cfg['outputpath']=os.path.join(cfg['path'], 'csv_invest')
	else:
		cfg['outputpath']=os.path.join(cfg['path'], 'csv_simul')
else:
	if not os.path.isabs(cfg['outputpath']):
		cfg['outputpath'] = os.path.join(cfg['path'], cfg['outputpath'])
if 'dirTimeSeries' not in cfg: cfg['dirTimeSeries'] = os.path.join(cfg['path'], 'TimeSeries/')
if 'genesys_inputpath' not in cfg: cfg['genesys_inputpath'] = os.path.join(cfg['path'], 'genesys_inputs/')
if 'timeseriespath' not in cfg: cfg['timeseriespath'] = os.path.join(cfg['path'], 'TimeSeries/')
if 'configDir' not in cfg: cfg['configDir']=os.path.join(cfg['path'], 'settings/')
if 'pythonDir' not in cfg: cfg['pythonDir']=os.path.join(cfg['p4rpath'],'scripts/python/plan4res-scripts/settings/')
if 'nomenclatureDir' not in cfg: cfg['nomenclatureDir']=os.path.join(cfg['p4rpath'],'scripts/python/openentrance/definitions/')

for datagroup in cfg['datagroups']:
	if 'inputdatapath' not in cfg['datagroups'][datagroup]:
		cfg['datagroups'][datagroup]['inputdatapath']='IAMC'
	cfg['datagroups'][datagroup]['inputdatapath'] = os.path.join(cfg['path'], cfg['datagroups'][datagroup]['inputdatapath'])
	if 'inputdata' not in cfg['datagroups'][datagroup]:
		cfg['datagroups'][datagroup]['inputdata']=namedataset+'.xlsx'

cfg['dirTimeSeries']=cfg['timeseriespath']	

if not os.path.isdir(cfg['outputpath']):
	os.makedirs(cfg['outputpath'])
logger.info('results of this script will be available in: '+cfg['outputpath'])

isInertia= ( 'InertiaDemand' in cfg['CouplingConstraints'] )
isPrimary= ( 'PrimaryDemand' in cfg['CouplingConstraints'] )
isSecondary= ( 'SecondaryDemand' in cfg['CouplingConstraints'] )
isCO2=False
if 'PollutantBudget' in cfg['CouplingConstraints']:
	if 'CO2' in cfg['CouplingConstraints']['PollutantBudget']:
		isCO2=True

partitionDemand = cfg['CouplingConstraints']['ActivePowerDemand']['Partition']
if isInertia: partitionInertia=cfg['CouplingConstraints']['InertiaDemand']['Partition']
if isPrimary: partitionPrimary=cfg['CouplingConstraints']['PrimaryDemand']['Partition']
if isSecondary: partitionSecondary=cfg['CouplingConstraints']['SecondaryDemand']['Partition']

isInvest= cfg['ParametersCreate']['invest']


def log_debug(msg):
	if 'debug' in cfg['ParametersCreate']:
		if cfg['ParametersCreate']['debug']: logger.info(msg)
		
def filter_iamc_dataframe(df, filter_attribute, filter_set, display_detail=True, sep=', '):
	if not isinstance(filter_set, list):
		filter_set = [filter_set]
	excluded_data = sorted([d for d in getattr(df, filter_attribute) if d not in filter_set])
	if len(excluded_data) > 0:
		if display_detail:
			logger.warning(f'The following data for {filter_attribute} will be excluded:\n\t'+sep.join([str(_) for _ in excluded_data]))
		else:
			logger.warning(f'{len(excluded_data)} {filter_attribute}{("s" if len(excluded_data)>1 else "")} will be excluded from the dataset')
	return df.filter(**{filter_attribute:filter_set})

# connect to the openentrance scenario explorer (set credentials)
if cfg['ScenarioExplorer']['annual'] or cfg['ScenarioExplorer']['subannual']:
	pyam.iiasa.set_config(cfg['ScenarioExplorer']['user'],cfg['ScenarioExplorer']['password'])
	pyam.iiasa.Connection(cfg['ScenarioExplorer']['instance'])

# create the dictionnary of variables containing the correspondence between plan4res (SMS++) variable 
# names and openentrance nomenclature variable names
vardict={}
with open(os.path.join(path,cfg['configDir'], "VariablesDictionnary.yml"),"r") as myvardict:
	vardict=yaml.safe_load(myvardict)

# create the dictionnary of time series, containing the names of the timeseries to be included in 
# the dataset
timeseriesdict={}
timeseries_setting_file = os.path.join(path, cfg['configDir'], "DictTimeSeries.yml")
with open(timeseries_setting_file,"r") as mytimeseries:
	timeseriesdict=yaml.safe_load(mytimeseries)

# if only one scenario/year is defined in config file set the list of scenarios / years to 1 element
if 'scenarios' not in cfg: cfg['scenarios']= [ cfg['scenario'] ]
if 'years' not in cfg: cfg['years']= [ cfg['year'] ]

# create the list of options
if 'options' in cfg: 
	option_types=cfg['options']
	cfg['options']=[]
	for option_type in option_types:
		cfg['options'].append('WITH_'+option_type)
		cfg['options'].append('WITHOUT_'+option_type)
else: cfg['options']=['None']
cfg['treat']=cfg['csvfiles']

cfg['StochasticScenarios']=[str(x) for x in cfg['StochasticScenarios']]

# loop on scenarios, years, options => create one dataset per (scenario,year,option) triplet
for current_scenario, current_year, current_option in product(cfg['scenarios'],cfg['years'],cfg['options']):
	cfg['scenario']=current_scenario
	cfg['year']=current_year
	if current_option != "None":
		logger.info('\ncreate dataset for scenario '+current_scenario+', year '+str(current_year)+' and option '+current_option)
	else:
		logger.info('\ncreate dataset for scenario '+current_scenario+', year '+str(current_year))
	
	if len(cfg['scenarios'])==1 and len(cfg['years'])==1 and len(cfg['options'])<2:
		outputdir=cfg['outputpath']
		if not os.path.isdir(outputdir):
			os.makedirs(outputdir)
	elif not current_option=='None':
		outputdir=os.path.join(cfg['outputpath'], 'plan4res-'+cfg['scenario']+'-'+str(cfg['year'])+'-'+current_option)
		if not os.path.isdir(outputdir):
			os.makedirs(outputdir)
	else:
		outputdir=os.path.join(cfg['outputpath'], 'plan4res-'+cfg['scenario']+'-'+str(cfg['year']))
		if not os.path.isdir(outputdir):
			os.makedirs(outputdir)

	# upload of relevant Scenario data from platform
	# creation of a csv file and a pandas dataframe containing all necessary data
	####################################################################################################
	i=0

	ExistsAnnualData=False  # AnnualDataFrame is still empty
	ExistsSubAnnualData=False # SubAnnualDataFrame is still empty
	SubAnnualDataFrame=pyam.IamDataFrame
	AnnualDataFrame=pyam.IamDataFrame

	listLocalReg=[]
	listGlobalReg=[]

	# define list of aggregated regions
	if 'ExistsNuts' in cfg['ParametersCreate']:
		if cfg['ParametersCreate']['ExistsNuts']:
			with open(os.path.join(cfg['nomenclatureDir'], "region/nuts3.yaml"),"r",encoding='UTF-8') as nutsreg:
				nuts3=yaml.safe_load(nutsreg)
			with open(os.path.join(cfg['nomenclatureDir'], "region/nuts2.yaml"),"r",encoding='UTF-8') as nutsreg:
				nuts2=yaml.safe_load(nutsreg)
			with open(os.path.join(cfg['nomenclatureDir'], "region/nuts1.yaml"),"r",encoding='UTF-8') as nutsreg:
				nuts1=yaml.safe_load(nutsreg)
	with open(os.path.join(cfg['nomenclatureDir'], "region/countries.yaml"),"r",encoding='UTF-8') as nutsreg:
		countries=yaml.safe_load(nutsreg)
	with open(os.path.join(cfg['nomenclatureDir'], "region/ehighway.yaml"),"r",encoding='UTF-8') as nutsreg:
		subcountries=yaml.safe_load(nutsreg)
	with open(os.path.join(cfg['nomenclatureDir'], "region/european-regions.yaml"),"r",encoding='UTF-8') as nutsreg:
		aggregateregions=yaml.safe_load(nutsreg)

	# create table of correspondence for iso3 and iso2 names
	iso3=pd.Series(str)
	iso2=pd.Series(str)
	listcountries=[]
	for k in range(len(countries[0]['Countries'])):
		countryname=next(iter(countries[0]['Countries'][k]))
		iso=countries[0]['Countries'][k][countryname]['iso3']
		iso3[iso]=countryname
		iso=countries[0]['Countries'][k][countryname]['iso2']
		iso2[iso]=countryname
		listcountries.append(countryname)
		for k2 in range(len(countries[0]['Countries'])):
			countryname2=next(iter(countries[0]['Countries'][k2]))
			iso_2=countries[0]['Countries'][k2][countryname2]['iso3']
			iso3[iso+'>'+iso_2]=countryname+'>'+countryname2
			iso_2=countries[0]['Countries'][k2][countryname2]['iso2']
			iso2[iso+'>'+iso_2]=countryname+'>'+countryname2
	dict_iso3=iso3.to_dict()
	dict_iso2=iso2.to_dict()
	rev_dict_iso3={v:k for k,v in dict_iso3.items()}
	rev_dict_iso2={v:k for k,v in dict_iso2.items()}

	# create the list of regions to work on 
	for datagroup in cfg['listdatagroups']:
		if type(cfg['datagroups'][datagroup]['regions']['local'])==list: listLocalReg=listLocalReg+cfg['datagroups'][datagroup]['regions']['local']
		elif cfg['datagroups'][datagroup]['regions']['local']=='countries':listLocalReg=listLocalReg+countries.keys()
		elif cfg['datagroups'][datagroup]['regions']['local']=='countries_ISO3':listLocalReg=listLocalReg+iso3.index.tolist()
		elif cfg['datagroups'][datagroup]['regions']['local']=='countries_ISO2':listLocalReg=listLocalReg+iso2.index.tolist()
		
		if(cfg['datagroups'][datagroup]['regions']['global']): listGlobalReg.append(cfg['datagroups'][datagroup]['regions']['global'])
		
	# create list of lines
	lines=[]
	for region1 in cfg['listregionsGET']:
		for region2 in cfg['listregionsGET']:
			if (region1!=region2):
				lines.append(region1+'>'+region2)
				lines.append(region2+'>'+region1)
	
	# define list of regions including lines
	listLocalReg=listLocalReg+lines  
	for region in cfg['listregionsGET']:
		if region not in listGlobalReg: listLocalReg.append(region)
	listRegGet=listLocalReg+listGlobalReg
	
	# treat cases where disaggregation level is nuts1/2/3
	if 'ExistsNuts' in cfg['ParametersCreate']:
		if cfg['ParametersCreate']['ExistsNuts']:
			minaggr=''
			nuts3list= {}
			nuts2list= {}
			nuts1list= {}
			for k in nuts1[0]['NUTS1']: nuts1list.update(k)
			for k in nuts2[0]['NUTS2']: nuts2list.update(k)
			for k in nuts3[0]['NUTS3']: nuts3list.update(k)
			# create lists of nuts per countries
			countryNuts3={v:[] for v in listcountries}
			{countryNuts3[v['country']].append(k) for k,v in nuts3list.items()}
			countryNuts2={v:[] for v in listcountries}
			{countryNuts2[v['country']].append(k) for k,v in nuts2list.items()}
			countryNuts1={v:[] for v in listcountries}
			{countryNuts1[v['country']].append(k) for k,v in nuts1list.items()}
			
			# extend list of regions to get with regions from datagroups
			if 'nuts3' in listRegGet:
				while 'nuts3' in listRegGet: listRegGet.remove('nuts3')
				listRegGet=listRegGet+list(nuts3.keys())
				minaggr='nuts3'	
			if 'nuts2' in listRegGet:
				while 'nuts2' in listRegGet: listRegGet.remove('nuts2')
				listRegGet=listRegGet+list(nuts2.keys())
				minaggr='nuts2'	
			if 'nuts1' in listRegGet:
				while 'nuts1' in listRegGet: listRegGet.remove('nuts1')
				listRegGet=listRegGet+list(nuts1.keys())
				minaggr='nuts1'	
	if 'countries' in listRegGet:
		while 'countries' in listRegGet: listRegGet.remove('countries')
		listRegGet=listRegGet+list(countries.keys())

	# create list of variables
	for datagroup in cfg['listdatagroups']:
	# loop on all different data sources
		if 'scenario' in cfg['datagroups'][datagroup]: scenget=cfg['datagroups'][datagroup]['scenario']
		else: scenget=cfg['scenario']
		log_debug('reading '+datagroup)
		
		listvardatagroup=[]
		listlocalvar=[] # list of variables which are not 'global'
		listglobalvar=[]
		globalreg= cfg['datagroups'][datagroup]['regions']['global']
		
		# loop on category of variables : coupling or techno
		for varcat1 in cfg['datagroups'][datagroup]['listvariables']:
			# treat coupling variables (ie variables not depending on techno)
			if varcat1=='coupling':
				# loop on category of coupling variables: mean, add, flow or global
				for varcat2 in cfg['datagroups'][datagroup]['listvariables'][varcat1].keys():
					# loop on variables
					for var in cfg['datagroups'][datagroup]['listvariables'][varcat1][varcat2]:
						listvardatagroup.append(var)
					# specific treatment for global variables (meaning they do not depend on regions)
						if varcat2=='global':listglobalvar.append(var)
						else: listlocalvar.append(var)
							
			# treat variables depending on technos
			elif varcat1=='techno':
				# loop on category of techno variables: thermal, reservoir, .....
				for varcat2 in cfg['datagroups'][datagroup]['listvariables'][varcat1].keys():
					# loop on category of variables: add, mean, global
					for varcat3 in cfg['datagroups'][datagroup]['listvariables'][varcat1][varcat2].keys():
						# is there only a list of variables or subgroups of variables per fuels
						if(type(cfg['datagroups'][datagroup]['listvariables'][varcat1][varcat2][varcat3])==list):
							# loop on variables and fuels
							for var in cfg['datagroups'][datagroup]['listvariables'][varcat1][varcat2][varcat3]:
								if varcat2 not in cfg['technos'].keys():
									logger.warning(f'/!\ {varcat2} does not exist in technos listed in setting files')
									if 'InteractiveMode' in cfg and cfg['InteractiveMode']:
										query=input('Continue creation of plan4res input files? [y]/n\n')
										if query=='n':
											log_and_exit(1, os.getcwd())
								else:
									for fuel in cfg['technos'][varcat2]:
										newvar=var+fuel
										listvardatagroup.append(newvar)
										if varcat3=='global': listglobalvar.append(newvar) 
										else: listlocalvar.append(newvar)
						else: # there are subgroups of variables per fuels
							for subgroup in cfg['datagroups'][datagroup]['listvariables'][varcat1][varcat2][varcat3].keys():
								for var in cfg['datagroups'][datagroup]['listvariables'][varcat1][varcat2][varcat3][subgroup]['variables']:
									for fuel in cfg['datagroups'][datagroup]['listvariables'][varcat1][varcat2][varcat3][subgroup]['fuels']:
										newvar=var+fuel
										listvardatagroup.append(newvar)
										if varcat3=='global': listglobalvar.append(newvar) 
										else: listlocalvar.append(newvar)
		
		if ( cfg['datagroups'][datagroup]['subannual'] and cfg['ScenarioExplorer']['subannual'] ) or ( not cfg['datagroups'][datagroup]['subannual'] and cfg['ScenarioExplorer']['annual']):
			log_debug('download data from platform')
			
			groupdf=pyam.read_iiasa('openentrance',model=cfg['datagroups'][datagroup]['model'],
				variable=listvardatagroup,
				region=listRegGet,year=cfg['year'],
				scenario=scenget)
						
			# rename regions if necessary
			if cfg['datagroups'][datagroup]['regions']['rename']=='countries_ISO3':groupdf.rename(dict_iso3,inplace=True)
			if cfg['datagroups'][datagroup]['regions']['rename']=='countries_ISO2':groupdf.rename(dict_iso2,inplace=True)

			if cfg['datagroups'][datagroup]['subannual']:
				if not ExistsSubAnnualData:
					ExistsSubAnnualData=True
					SubAnnualDataFrame=groupdf
				else:
					SubAnnualDataFrame=SubAnnualDataFrame.append(groupdf)
			else:
				if not ExistsAnnualData:
					ExistsAnnualData=True
					AnnualDataFrame=groupdf
				else:
					AnnualDataFrame=AnnualDataFrame.append(groupdf)
				
		if not (cfg['ScenarioExplorer']['annual'] and cfg['ScenarioExplorer']['subannual']):
			# load data from files per data source (previously uploaded from platform)
			log_debug('open data files')
		
			if ( cfg['datagroups'][datagroup]['subannual'] and not cfg['ScenarioExplorer']['subannual']) or ( not cfg['datagroups'][datagroup]['subannual'] and not cfg['ScenarioExplorer']['annual']):
				logger.info('reading datagroup '+datagroup+' from '+os.path.join(cfg['genesys_inputpath'], cfg['datagroups'][datagroup]['inputdatapath'], cfg['datagroups'][datagroup]['inputdata']))
														
				if 'Start' in cfg['datagroups'][datagroup]['inputdata']:
					file=os.path.join(cfg['genesys_inputpath'], cfg['datagroups'][datagroup]['inputdatapath'], cfg['datagroups'][datagroup]['inputdata']['Start'], str(cfg['variant']), cfg['datagroups'][datagroup]['inputdata']['End'])
				else:
					file=os.path.join(cfg['genesys_inputpath'], cfg['datagroups'][datagroup]['inputdatapath'], cfg['datagroups'][datagroup]['inputdata'])
				
				# creation of empty df for storing annual and subannual data for the current group
				dfdatagroup=pyam.IamDataFrame(pd.DataFrame(columns=['model','scenario','region','variable','unit','subannual',str(cfg['year'])]))
				# filter on the listgetfegion regions
				SubAdfdatagroup=pyam.IamDataFrame(pd.DataFrame(columns=['model','scenario','region','variable','unit','subannual',str(cfg['year'])]))
				
				# read data as a IAMDataFrame
				log_debug('read file '+file)
				if not os.path.isfile(file):
					logger.error('\nError: '+file+' does not exist.') 
					logger.error('Check file '+settings_create+'. You can specify the input data repository using key inputdatapath (default=IAMC) and file name using key inputdata (default=name of the study+.xlsx).')
					log_and_exit(2, cfg['path'])
				if 'xlsx' in file:
					df=pd.read_excel(file,sheet_name='data')
				else:
					df=pd.read_csv(file)
				
				if 'Subannual' in df.columns:
					if len(df['Subannual'].unique()==1): df=df.drop(['Subannual'],axis=1)
				dfdatagroup=pyam.IamDataFrame(data=df)

				if 'rename' in cfg['datagroups'][datagroup]['regions'] and 'countries_ISO3' in cfg['datagroups'][datagroup]['regions']['rename']['dict']:				
					log_debug('renaming ISO3')
					dfdatagroup=dfdatagroup.rename(region=dict_iso3)
				if 'rename' in cfg['datagroups'][datagroup]['regions'] and 'countries_ISO2' in cfg['datagroups'][datagroup]['regions']['rename']['dict']:
					log_debug('renaming ISO2')
					dfdatagroup=dfdatagroup.rename(region=dict_iso2)
				if 'rename' in cfg['datagroups'][datagroup]['regions'] and 'added' in cfg['datagroups'][datagroup]['regions']['rename']:
					for reg in cfg['datagroups'][datagroup]['regions']['rename']['added']:
						dfdatagroup=dfdatagroup.rename(region={reg:cfg['datagroups'][datagroup]['regions']['rename']['added'][reg]})
				log_debug('change country names')
				excluded_data = [d for d in dfdatagroup.region if d not in listRegGet]
				if len(excluded_data)>0:
					logger.warning('/!\ data from the following regions will be excluded: '+', '.join(excluded_data))
				dfdatagroup=filter_iamc_dataframe(dfdatagroup, 'region', listRegGet)
				
				log_debug('filter countries')
				# if there are data at lower granularity than country or cluster (only country until now), aggregate
				firstcountry=1
				if 'ExistsNuts' in cfg['ParametersCreate']:
					if cfg['ParametersCreate']['ExistsNuts']:
						for country in listcountries:
							
							# create list of nuts of the country
							listNuts1=countryNuts1[country]
							listNuts2=countryNuts2[country]
							listNuts3=countryNuts3[country]
							listNuts=listNuts1+listNuts2+listNuts3
							NumberOfNutsLists=int(len(listNuts1)>0)+int(len(listNuts2)>0)+int(len(listNuts3)>0)
							# To be implemented: include weights to aggregation, weight=1/NumberOfNutsLists
							
							if (len(listNuts)>0 and ('NoNutsAggregation' not in cfg)): 
								log_debug('aggregating nuts')
								dfdatagroup.aggregate_region(dfdatagroup.variable,region=country, subregions=listNuts, append=True)

				log_debug(f'filter model '+cfg["datagroups"][datagroup]["model"])
				dfdatagroup=filter_iamc_dataframe(dfdatagroup, 'model', cfg['datagroups'][datagroup]['model'])
				log_debug(f'filter scenario {scenget}')
				dfdatagroup=filter_iamc_dataframe(dfdatagroup, 'scenario', scenget)
				log_debug(f'filter year {cfg["year"]}')
				dfdatagroup=filter_iamc_dataframe(dfdatagroup, 'year', cfg['year'])
				# treat case of variables which are only present with higher disaggregation	
				# eg. we need Capacity|Electricity|Gas and we have Capacity|Electricity|Gas|OCGT and Capacity|Electricity|Gas|CCGT
				# in that case we will compute Capacity|Electricity|Gas using the rules defined in the settings
				
				# get list of missing variables
				existing_variables = dfdatagroup.variable
				missing_variables = [var for var in listvardatagroup if var not in existing_variables]
				missing_variables = list(set(missing_variables))

				# Get for each of these missing variables the list of "sub-variables" in dfdatagroup
				sub_variables_dict = {}
				to_remove = []
				for var in missing_variables:
					sub_variables = [v for v in existing_variables if v.startswith(var + '|')]
					if len(sub_variables)>0:
						sub_variables_dict[var] = sub_variables
					else:
						to_remove.append(var)
				for var in to_remove:
					missing_variables.remove(var)
				
				# get aggregation method
				aggregation_methods = {}
				for var in missing_variables:
					method = search_variable(var, cfg['datagroups'][datagroup]['listvariables'])
					aggregation_methods[var] = method

				# Create rows for the missing variables
				for var, sub_vars in sub_variables_dict.items():
					if sub_vars:
						method = aggregation_methods[var]
						if method == 'global':
							method = 'mean'
						elif method in ['flow', 'add']:
							method = 'sum'
						dfdatagroup.aggregate(variable=var, components=sub_vars, method=method, append=True)

				# only keep variables in settings
				log_debug(f'filter on variables in settings: {listvardatagroup}')
				dfdatagroup=filter_iamc_dataframe(dfdatagroup, 'variable', listvardatagroup, display_detail=True, sep='\n\t')
				
				if cfg['datagroups'][datagroup]['subannual']: SubAdfdatagroup=dfdatagroup
				
				if cfg['datagroups'][datagroup]['subannual']:
					if not ExistsSubAnnualData:
						ExistsSubAnnualData=True
						SubAnnualDataFrame=SubAdfdatagroup
					else:
						SubAnnualDataFrame=SubAnnualDataFrame.append(SubAdfdatagroup)
				else:
					if not ExistsAnnualData:
						ExistsAnnualData=True
						AnnualDataFrame=dfdatagroup
					else:
						AnnualDataFrame=AnnualDataFrame.append(dfdatagroup)
	
	#conversion of units to plan4res usual units (MWh, MW, €/MWh, €/MW/yr, €/MW)
	if(ExistsAnnualData): #check if there exist annual data
		logger.info('converting units')
		for var_unit in cfg['ParametersCreate']['conversions']:
			if 'factor' in cfg['ParametersCreate']['conversions'][var_unit]:
				AnnualDataFrame=AnnualDataFrame.convert_unit(var_unit, to=cfg['ParametersCreate']['conversions'][var_unit]['to'], factor=float(cfg['ParametersCreate']['conversions'][var_unit]['factor'])) 
			else:
				AnnualDataFrame=AnnualDataFrame.convert_unit(var_unit, to=cfg['ParametersCreate']['conversions'][var_unit]['to']) 

		# validate the format of the data (prevents errors)
		AnnualDataFrame.validate(exclude_on_fail=True)
	
	if(ExistsSubAnnualData): # check if there exist subannual data
		for var_unit in cfg['ParametersCreate']['conversions']:
			if 'factor' in cfg['ParametersCreate']['conversions'][var_unit]:
				SubAnnualDataFrame=SubAnnualDataFrame.convert_unit(var_unit, to=cfg['ParametersCreate']['conversions'][var_unit]['to'], factor=cfg['ParametersCreate']['conversions'][var_unit]['factor']) 
			else:
				SubAnnualDataFrame=SubAnnualDataFrame.convert_unit(var_unit, to=cfg['ParametersCreate']['conversions'][var_unit]['to']) 

		SubAnnualDataFrame.validate(exclude_on_fail=True)
		
	#regional aggregations 
	logger.info('computing regional aggregations')
	
	# some variables are added (all energy, capacity), others are averageded, and finally for flow variable specific aggregations are done
	listvaradd=[]
	listvarmean=[]
	for datagroup in cfg['listdatagroups']:
		if 'coupling' in cfg['datagroups'][datagroup]['listvariables'].keys():
			for typeagr in cfg['datagroups'][datagroup]['listvariables']['coupling'].keys():
				for var in cfg['datagroups'][datagroup]['listvariables']['coupling'][typeagr]:
					variable=var
					if typeagr=='add': 
						if variable not in listvaradd: 
							listvaradd.append(variable)
					elif typeagr=='mean': 
						if variable not in listvarmean: 
							listvarmean.append(variable)
		if 'techno' in cfg['datagroups'][datagroup]['listvariables'].keys():
			for technogroup in cfg['datagroups'][datagroup]['listvariables']['techno'].keys():
				for typeagr in cfg['datagroups'][datagroup]['listvariables']['techno'][technogroup].keys():
					if type(cfg['datagroups'][datagroup]['listvariables']['techno'][technogroup][typeagr])==list:
						for var in cfg['datagroups'][datagroup]['listvariables']['techno'][technogroup][typeagr]:
							if technogroup not in cfg['technos'].keys():
								logger.warning(f'/!\ While creating data for {datagroup}>{typeagr}>{var} ! {technogroup} does not exist in technos listed in settings file.')
								if 'InteractiveMode' in cfg and cfg['InteractiveMode']:
									query=input('Continue creation of plan4res input files? [y]/n\n')
									if query=='n':
										log_and_exit(1, os.getcwd())
							else:
								for techno in cfg['technos'][technogroup]:
									newvar=var+techno
									if typeagr=='add': 
										if newvar not in listvaradd: 
											listvaradd.append(newvar)
									elif typeagr=='mean': 
										if newvar not in listvarmean: 
											listvarmean.append(newvar)
					else:
						for groupvar in cfg['datagroups'][datagroup]['listvariables']['techno'][technogroup][typeagr].keys():
							for fuel in cfg['datagroups'][datagroup]['listvariables']['techno'][technogroup][typeagr][groupvar]['fuels']:
								for var in cfg['datagroups'][datagroup]['listvariables']['techno'][technogroup][typeagr][groupvar]['variables']:
									newvar=var+fuel
									if typeagr=='add': 
										if newvar not in listvaradd: 
											listvaradd.append(newvar)
									elif typeagr=='mean': 
										if newvar not in listvarmean: 
											listvarmean.append(newvar)
	if(cfg['aggregateregions']!=None):
		for reg in cfg['aggregateregions'].keys():
			log_debug('aggregating ' +reg+' for subregions:')
			log_debug(cfg['aggregateregions'][reg])
			# creation of aggregated timeseries
			listTypes={'ZV': cfg['CouplingConstraints']['ActivePowerDemand']['SumOf'],'RES':cfg['technos']['res']+cfg['technos']['runofriver'],'SS':['Inflows']}
			for typeData in listTypes:
				for typeSerie in listTypes[typeData]:
					sumValTS=0.0
					firstSerie=True
					isSeries=False
					for region in cfg['aggregateregions'][reg]:
						if region in timeseriesdict[typeData][typeSerie].keys():
							isSeries=True
							file = os.path.join(cfg['dirTimeSeries'], timeseriesdict[typeData][typeSerie][region])
							if not os.path.isfile(file):
								logger.error('\nError: '+file+' does not exist.') 
								logger.error('Check file '+settings_create+'. You can specify the timeseries input repository using key dirTimeSeries (default=study path/TimeSeries).')
								logger.error('Also check the timeseries file name  for '+typeData+', '+typeSerie+', '+region+' in file '+timeseries_setting_file)
								log_and_exit(2, cfg['path'])
							timeserie=pd.read_csv(file,index_col=0)
							if len(timeserie.columns)>1:
								for col in timeserie.columns:
									if 'Unnamed' in col: timeserie.drop([col],axis=1)
								timeserie=timeserie[ cfg['StochasticScenarios'] ]
							else:
								timeserie.columns=['DET']
							if typeData=='ZV':
								varTS=vardict['Input']['Var'+typeData][typeSerie]
							elif typeData=='SS':
								varTS=vardict['Input']['Var'+typeData][typeSerie]+'Reservoir'
							else:
								varTS=vardict['Input']['Var'+typeData]['MaxPower']+typeSerie
							if varTS in AnnualDataFrame.variable:
								valTS_IAMdf=AnnualDataFrame.filter(variable=varTS,region=region,year=current_year).as_pandas()['value'].unique()
								if len(valTS_IAMdf)>0: valTS=valTS_IAMdf[0]
								else: valTS=0.0
							else:
								# the only variables which may not be in the data are the different parts of the ActiveDemand:
								if varTS in [timeseriesdict['ZV'][part] for part in cfg['CouplingConstraints']['ActivePowerDemand']['SumOf']]:
									# compute valTS as Total-sum of parts which are present
									valTS=AnnualDataFrame[ (AnnualDataFrame['Variable']==timeseriesdict['ZV']['Total']) & (df['Region']==region) ][str(current_year)].unique()[0]
									for part in cfg['CouplingConstraints']['ActivePowerDemand']:
										if vardict['Input']['Var'+typeData][part] in AnnualDataFrame['Variable'].unique():
											valTS=valTS-AnnualDataFrame[ (AnnualDataFrame['Variable']==timeseriesdict['ZV'][part]) & (df['Region']==region) ][str(current_year)].unique()[0]
							
							if 'valTS' not in locals():
								logger.error(f'No data for variable {varTS} region {region} year {current_year}')
								log_and_exit(2, cfg['path'])
							if valTS==0.0: valTS=cfg['ParametersCreate']['zerocapacity'] 
							
							if firstSerie: 
								newSerie=valTS*timeserie
								firstSerie=False
								sumValTS=valTS
							else:
								newSerie=newSerie+valTS*timeserie
								sumValTS=sumValTS+valTS
					if isSeries:
						if sumValTS>0: newSerie=(1/sumValTS)*newSerie
						else: newSerie=0.0*newSerie
						nameNewSerie='AggregatedTimeSerie_'+typeSerie+'_'+reg+'.csv'
						nameNewSerie=nameNewSerie.replace('|', pipe_replace)
						newSerie.to_csv(os.path.join(cfg['dirTimeSeries'], nameNewSerie))
						timeseriesdict[typeData][typeSerie][reg]=nameNewSerie
		
			# aggregation of variables
			for variable in listvaradd:
				if ExistsAnnualData: 
					if variable in AnnualDataFrame.variable: 
						AnnualDataFrame.aggregate_region(variable, region=reg, subregions=cfg['aggregateregions'][reg], append=True)
				if ExistsSubAnnualData: 
					if variable in SubAnnualDataFrame.variable: 
						SubAnnualDataFrame.aggregate_region(variable, region=reg, subregions=cfg['aggregateregions'][reg], append=True)
			for variable in listvarmean:
				if ExistsAnnualData: 
					if variable in AnnualDataFrame.variable: 
						AnnualDataFrame.aggregate_region(variable, region=reg, subregions=cfg['aggregateregions'][reg], method='mean',append=True)	
				if ExistsSubAnnualData: 
					if variable in SubAnnualDataFrame.variable: 
						SubAnnualDataFrame.aggregate_region(variable, region=reg, subregions=cfg['aggregateregions'][reg], method='mean',append=True)	
	#remove aggregated subregions
	listregion=listGlobalReg
	for partition in cfg['partition']:
		for region in cfg['partition'][partition]:
			if region not in listregion: listregion.append(region)

	log_debug('filter regions and links '+', '.join(listregion+lines))
	if ExistsAnnualData:
		AnnualDataFrame = filter_iamc_dataframe(AnnualDataFrame, 'region', listregion+lines)
		AnnualDataFrame.to_csv(os.path.join(outputdir, 'IAMC_annual_data.csv'))
		bigdata=AnnualDataFrame
	if ExistsSubAnnualData: 
		SubAnnualDataFrame = filter_iamc_dataframe(SubAnnualDataFrame, 'region', listregion+lines)
		SubAnnualDataFrame.to_csv(os.path.join(outputdir, 'IAMC_subannual_data.csv'))
		bigdata_SubAnnual=SubAnnualDataFrame
	

	
	# creation of plan4res dataset
	################################"

	# creating list of regions
	listregions=listGlobalReg
	for partition in list(cfg['partition'].keys()):
		listregions=listregions+cfg['partition'][partition]
	listregions = list(set(listregions))
	logger.info('regions in dataset:'+str(listregions))
	
	listregions_dataset=cfg['partition'][ partitionDemand  ]

	# create file ZP_ZonePartition
	#############################################################
	if cfg['csvfiles']['ZP_ZonePartition']:
		logger.info('\n\nCreating ZP_ZonePartition')
		nbreg=len(cfg['partition'][ partitionDemand  ])
		nbpartition=len(cfg['partition'])

		ZP = pd.DataFrame(columns=list(cfg['partition'].keys()),index=range(nbreg))
		ZP[ partitionDemand ]=pd.Series(cfg['partition'][ partitionDemand ],name=partition,index=range(nbreg))
		for partition in list(cfg['partition'].keys()):
			if not partition == cfg['CouplingConstraints']['ActivePowerDemand']['Partition']:
				if len(cfg['partition'][partition])==nbreg: ZP[partition]=pd.Series(cfg['partition'][partition],name=partition,index=range(nbreg))
				else: ZP[partition]=pd.Series(cfg['partition'][partition][0] ,name=partition,index=range(nbreg))
		ZP.to_csv(os.path.join(outputdir, cfg['csvfiles']['ZP_ZonePartition']), index=False)

	# create file IN_Interconnections
	###############################################################
	if cfg['csvfiles']['IN_Interconnections']:
		IN = pd.DataFrame()
		logger.info('\n\nCreating IN_Interconnections')
		for variable in vardict['Input']['VarIN']:
			varname=vardict['Input']['VarIN'][variable]
			vardf=bigdata.filter(variable=varname).as_pandas(meta_cols=False)
			vardf=vardf.set_index('region')
			vardf=vardf.rename(columns={"value":variable})
			dataIN=vardf[variable]
			IN=pd.concat([IN, dataIN], axis=1)	
		IN=IN.fillna(value=0.0)
		
		# delete lines which start/end in same aggregated
		log_debug('deleting lines which start and end in same aggregated region')
		IN['Name']=IN.index
		IN['StartLine']=IN['Name'].str.split('>',expand=True)[0]
		IN['EndLine']=IN['Name'].str.split('>',expand=True)[1]
		IN['AgrStart']=IN['StartLine']
		IN['AgrEnd']=IN['EndLine']
		for line in dataIN.index:
			# check if start / end is in an aggregated region
			regstart=line.split('>')[0]
			regend=line.split('>')[1]
			if(cfg['aggregateregions']!=None):
				for AggReg1 in cfg['aggregateregions'].keys():
					if (regstart in cfg['aggregateregions'][AggReg1]): IN.at[line,'AgrStart']=AggReg1
					if (regend in cfg['aggregateregions'][AggReg1]):  IN.at[line,'AgrEnd']=AggReg1
		# delete line with start and end in same aggregated region
		DeleteLines=IN[ IN.AgrStart == IN.AgrEnd ].index
		IN=IN.drop(DeleteLines)
		
		DeleteLines=[]
		# sum lines with start in same aggregated region AND end in same other aggregated region
		log_debug('aggregate lines which start or end in same aggregated region')
		log_debug(cfg['aggregateregions'])
		if(cfg['aggregateregions'] is not None and len(cfg['aggregateregions'])>0):
			for AggReg1 in cfg['partition'][ partitionDemand ]:
				if AggReg1 in cfg['aggregateregions'].keys(): # AggReg1 is an aggregated region
					for AggReg2 in cfg['partition'][partitionDemand]:
						if AggReg2 != AggReg1:
							if AggReg2 in cfg['aggregateregions'].keys(): # AggReg2 is another aggregated  region
								# sum all lines fitting this selection
								MaxPowerFlow=0.0
								InvCost=0.0
								LossFactor=0.0
								N=0
								for reg1 in cfg['aggregateregions'][AggReg1]:
									for reg2 in cfg['aggregateregions'][AggReg2]:
										if (reg1+'>'+reg2) in IN.index:
											N=N+1
											MaxPowerFlow=MaxPowerFlow+IN['MaxPowerFlow'][reg1+'>'+reg2]
											InvCost=InvCost+IN['InvestmentCost'][reg1+'>'+reg2]
											LossFactor=LossFactor+IN['LossFactor'][reg1+'>'+reg2]
											DeleteLines.append(reg1+'>'+reg2) # delete individual line
								if N==0: N=1
								InvCost=InvCost/N
								LossFactor=LossFactor/N
								IN = pd.concat([IN,pd.DataFrame(data=[[MaxPowerFlow,InvCost,LossFactor,AggReg1+'>'+AggReg2,AggReg1,AggReg2,AggReg1,AggReg2]],index=[AggReg1+'>'+AggReg2],columns=IN.columns)])
							else: 
								# AggReg2 is not an aggregated region
								MaxPowerFlow=0.0
								InvCost=0.0
								LossFactor=0.0
								N=0
								for reg1 in cfg['aggregateregions'][AggReg1]:
									if (reg1+'>'+AggReg2) in IN.index:
										N=N+1
										LossFactor=LossFactor+IN['LossFactor'][reg1+'>'+AggReg2]
										MaxPowerFlow=MaxPowerFlow+IN['MaxPowerFlow'][reg1+'>'+AggReg2]
										InvCost=InvCost+IN['InvestmentCost'][reg1+'>'+AggReg2]
										DeleteLines.append(reg1+'>'+AggReg2) # delete individual line
								if N==0: N=1
								InvCost=InvCost/N
								LossFactor=LossFactor/N
								IN = pd.concat([IN,pd.DataFrame(data=[[MaxPowerFlow,InvCost,LossFactor,AggReg1+'>'+AggReg2,AggReg1,AggReg2,AggReg1,AggReg2]],index=[AggReg1+'>'+AggReg2],columns=IN.columns)])
				else:
					for AggReg2 in cfg['partition'][partitionDemand]:
						if AggReg2 != AggReg1:
							if AggReg2 in cfg['aggregateregions'].keys(): # AggReg2 is an aggregated  region
								MaxPowerFlow=0.0
								InvCost=0.0
								LossFactor=0.0
								N=0
								for reg2 in cfg['aggregateregions'][AggReg2]:
									if (AggReg1+'>'+reg2) in IN.index:
										N=N+1
										InvCost=InvCost+IN['InvestmentCost'][AggReg1+'>'+reg2]
										LossFactor=LossFactor+IN['LossFactor'][AggReg1+'>'+reg2]
										MaxPowerFlow=MaxPowerFlow+IN['MaxPowerFlow'][AggReg1+'>'+reg2]
										DeleteLines.append(AggReg1+'>'+reg2) # delete individual line
								if N==0: N=1
								InvCost=InvCost/N
								LossFactor=LossFactor/N
								IN = pd.concat([IN,pd.DataFrame(data=[[MaxPowerFlow,InvCost,LossFactor,AggReg1+'>'+AggReg2,AggReg1,AggReg2,AggReg1,AggReg2]],index=[AggReg1+'>'+AggReg2],columns=IN.columns)])
		IN=IN.drop(DeleteLines )
		RowsToDelete = IN[ IN['MaxPowerFlow'] == 0 ].index
		IN=IN.drop(RowsToDelete )

		# merge lines reg1>reg2 reg2>reg1
		if 'debug' in cfg['ParametersCreate']: 
			if cfg['ParametersCreate']['debug']: logger.info('merging symetric lines')
		IN['MinPowerFlow']=0
		NewLines=[]
		LinesToDelete=[]
		for line in IN.index:
			regstart=line.split('>')[0]
			regend=line.split('>')[1]
			inverseline=regend+'>'+regstart
			if (line not in LinesToDelete):
				NewLines.append(line)
				if inverseline in IN.index: # NB : line are not symetrical in GENeSYS-MOD
					LinesToDelete.append(inverseline)
					IN.at[line,'MinPowerFlow']=-1.0*IN.loc[inverseline]['MaxPowerFlow']
		IN=IN.drop(LinesToDelete)

		IN['Name']=IN.index
		
		listcols=['Name','StartLine','EndLine','MaxPowerFlow','MinPowerFlow']
		if 'Impedance' in IN.columns: listcols.append('Impedance')
		if isInvest and 'interconnections' in cfg['ParametersCreate']['CapacityExpansion']:
			IN['MaxAddedCapacity']=0
			IN['MaxRetCapacity']=0
			if 'Share' in cfg['ParametersCreate']['CapacityExpansion']['interconnections']:
				# all lines can be invested
				IN['MaxAddedCapacity']=IN['MaxPowerFlow']*cfg['ParametersCreate']['CapacityExpansion']['interconnections']['Share']['MaxAdd']
				IN['MaxRetCapacity']=IN['MaxPowerFlow']*cfg['ParametersCreate']['CapacityExpansion']['interconnections']['Share']['MaxRet']
			else:
				for line in cfg['ParametersCreate']['CapacityExpansion']['interconnections']:
					if line != 'Share' and line !='InvestmentCost':
						IN.loc[line]['MaxAddedCapacity']=IN.loc[line]['MaxPowerFlow']*cfg['ParametersCreate']['CapacityExpansion']['interconnections'][line]['MaxAdd']
						IN.loc[line]['MaxRetCapacity']=IN.loc[line]['MaxPowerFlow']*cfg['ParametersCreate']['CapacityExpansion']['interconnections'][line]['MaxRet']
			if 'InvestmentCost' not in IN.columns:
				if 'InvestmentCost' in cfg['ParametersCreate']['CapacityExpansion']['interconnections']:
					IN['InvestmentCost']=cfg['ParametersCreate']['CapacityExpansion']['interconnections']['InvestmentCost']
				else:
					IN['InvestmentCost']=0
			listcols.append('MaxAddedCapacity')
			listcols.append('MaxRetCapacity')
			listcols.append('InvestmentCost')
		IN=IN[ listcols ]
		
		# delete lines which do not start and end in a zone in partition
		log_debug('delete lines which are not in partition')
		LinesToDelete=[]
		for line in IN.index:
			regstart=line.split('>')[0]
			regend=line.split('>')[1]
			if (regstart not in listregions_dataset):
				LinesToDelete.append(line)
			if (regend not in listregions_dataset):
				LinesToDelete.append(line)
		IN=IN.drop(LinesToDelete)
		IN.to_csv(os.path.join(outputdir, cfg['csvfiles']['IN_Interconnections']), index=False)
		
	# create file ZV_ZoneValues
	###############################################################
	numserie=0
	if cfg['csvfiles']['ZV_ZoneValues']:
		logger.info('\n\nCreating ZV_ZoneValues')
		
		ListTypesZV=[]
		for coupling_constraint in cfg['CouplingConstraints']:
			ListTypesZV=ListTypesZV+cfg['CouplingConstraints'][coupling_constraint]['SumOf']
		listvar=[]
		for var in vardict['Input']['VarZV'].keys(): listvar.append(vardict['Input']['VarZV'][var])
		datapartition=bigdata.filter(variable=listvar,region=listregions_dataset).as_pandas(meta_cols=False)
		datapartition=datapartition.rename(columns={"variable": "Type", "region": "Zone", "unit":"Unit"})

		# rename variables using Variables dictionnary
		# create reverse variables dictionnary
		dictZV={}
		dictZV=vardict['Input']['VarZV']
		reversedictZV={}
		for key, values in dictZV.items():
			for value in values:
				myvalue=dictZV[key]
				reversedictZV[myvalue]=key
		datapartition=datapartition.replace({"Type":reversedictZV})

		# include slack unit costs (slack means non served)
		datainertia=datapartition[ datapartition.Type == 'Inertia' ]
		Inertia=datainertia['value'].mean()
		MaxDemand=cfg['CouplingConstraints']['ActivePowerDemand']['MaxPower']
		CostDemand=cfg['CouplingConstraints']['ActivePowerDemand']['Cost']
		if 'PrimaryDemand' in cfg['CouplingConstraints']:
			isPrimary=True
			MaxPrimary=cfg['CouplingConstraints']['PrimaryDemand']['MaxPower']
			CostPrimary=cfg['CouplingConstraints']['PrimaryDemand']['Cost']
		else: isPrimary=False

		if 'SecondaryDemand' in cfg['CouplingConstraints']:
			isSecondary=True
			MaxSecondary=cfg['CouplingConstraints']['SecondaryDemand']['MaxPower']
			CostSecondary=cfg['CouplingConstraints']['SecondaryDemand']['Cost']
		else: isSecondary=False
			
		if 'InertiaDemand' in cfg['CouplingConstraints']:
			isInertia=True
			MaxInertia=cfg['CouplingConstraints']['InertiaDemand']['MaxPower']
			CostInertia=cfg['CouplingConstraints']['InertiaDemand']['Cost']
		else: isInertia=False
			
		# compute missing global values
		log_debug('compute missing Coupling constraints')
		isTotalEnergy=False		
		isHeating=False
		isTransport=False
		isCooling=False
		isShareCooling=False
		isElectrolyzer=False
		isCooking=False
		
		if 'Cooking' in datapartition.Type.unique(): isCooking=True
		if 'ElecVehicle' in datapartition.Type.unique(): isTransport=True
		if 'ElecHeating'  in datapartition.Type.unique(): isHeating=True
		if 'Total' in datapartition.Type.unique() : isTotalEnergy=True
		if 'AirCondition' in datapartition.Type.unique(): isCooling=True
		if 'Electrolyzer' in datapartition.Type.unique(): isElectrolyzer=True
		if 'Share|Final Energy|Electricity|Cooling' in bigdata.variable : isShareCooling=True
		for region in cfg['partition'][partitionDemand]:			
			if isTotalEnergy: 
				TotalEnergy=datapartition[ (datapartition.Zone==region) & (datapartition.Type == 'Total') ]['value'].mean()
				if math.isnan(TotalEnergy): TotalEnergy=0
			else:
				logger.info('missing variable for total electricity demand')
				exit(0)
			if isCooking: 
				CookingEnergy=datapartition[ (datapartition.Zone==region) & (datapartition.Type == 'Cooking') ]['value'].mean()
				if math.isnan(CookingEnergy): CookingEnergy=0
			else: CookingEnergy=0
			if isElectrolyzer: 
				ElectrolyzerEnergy=datapartition[ (datapartition.Zone==region) & (datapartition.Type == 'Electrolyzer') ]['value'].mean()
				if math.isnan(ElectrolyzerEnergy): ElectrolyzerEnergy=0
			else: ElectrolyzerEnergy=0
			if isHeating: 
				HeatingEnergy=datapartition[ (datapartition.Zone==region) & (datapartition.Type == 'ElecHeating') ]['value'].mean()
				if math.isnan(HeatingEnergy): HeatingEnergy=0
			else: HeatingEnergy=0
		
			
			if isTransport: 
				TransportEnergy=datapartition[ (datapartition.Zone==region) & (datapartition.Type == 'ElecVehicle') ]['value'].mean()
				if math.isnan(TransportEnergy): TransportEnergy=0
			else: TransportEnergy=0
					
			if isCooling: 
				CoolingEnergy=datapartition[ (datapartition.Zone==region) & (datapartition.Type == 'AirCondition') ]['value'].mean()
				if math.isnan(CoolingEnergy): CoolingEnergy=0 
			elif isShareCooling:
				ShareCooling=bigdata.filter(region=region,variable='Share|Final Energy|Electricity|Cooling').as_pandas(meta_cols=False)['value'].fillna(0).mean()
				if isTotalEnergy : 
					CoolingEnergy=TotalEnergy*(ShareCooling/100)
				else : CoolingEnergy=0
			else : CoolingEnergy=0
			
			OtherEnergy=TotalEnergy-CoolingEnergy-TransportEnergy-HeatingEnergy-ElectrolyzerEnergy-CookingEnergy
																									   					

			if (isShareCooling) and ('AirCondition' in cfg['CouplingConstraints']['ActivePowerDemand']['SumOf']):
				newdata=pd.DataFrame({'Type':['AirCondition'],'Unit':['MWh/yr'],\
					'Zone':[region],'model':['Recalculated'],'scenario':[cfg['listdatagroups'][0]],\
					'value':[CoolingEnergy],'year':[cfg['year']]})
				datapartition=pd.concat([datapartition,newdata],ignore_index=True)
			if ('Other' in cfg['CouplingConstraints']['ActivePowerDemand']['SumOf']):
				newdata=pd.DataFrame({'Type':['Other'],'Unit':['MWh/yr'],\
					'Zone':[region],'model':['Recalculated'],'scenario':[cfg['listdatagroups'][0]],\
					'value':[OtherEnergy],'year':[cfg['year']]})
				datapartition=pd.concat([datapartition,newdata],ignore_index=True)
			
			newdata=pd.DataFrame({'Type':['MaxActivePowerDemand','CostActivePowerDemand'],\
				'Unit':['MW','EUR/MWh'],\
				'Zone':[region, region],'model':['Recalculated','Recalculated'],\
				'scenario':[cfg['scenario'],cfg['scenario']],\
				'value':[MaxDemand,CostDemand],\
				'year':[cfg['year'],cfg['year']]})
			datapartition=pd.concat([datapartition,newdata],ignore_index=True)

			if isPrimary:
				Primary=datapartition[ (datapartition.Zone==region) & (datapartition.Type == 'PrimaryDemand') ]['value'].mean()
			if isSecondary:
				Secondary=datapartition[ (datapartition.Zone==region) & (datapartition.Type == 'SecondaryDemand') ]['value'].mean()			


			if isPrimary:
				newdata=pd.DataFrame({'Type':['MaxPrimaryDemand','CostPrimaryDemand'],\
					'Unit':['MW','EUR/MWh'],\
					'Zone':[region, region],'model':['Recalculated','Recalculated'],\
					'scenario':[cfg['scenario'],cfg['scenario']],\
					'value':[MaxDemand,CostDemand],\
					'year':[cfg['year'],cfg['year']]})
				datapartition=pd.concat([datapartition,newdata],ignore_index=True)

			if isSecondary:
				newdata=pd.DataFrame({'Type':['MaxSecondaryDemand','CostSecondaryDemand'],\
					'Unit':['MW','EUR/MWh'],\
					'Zone':[region, region],'model':['Recalculated','Recalculated'],\
					'scenario':[cfg['scenario'],cfg['scenario']],\
					'value':[MaxDemand,CostDemand],\
					'year':[cfg['year'],cfg['year']]})
				datapartition=pd.concat([datapartition,newdata],ignore_index=True)

			if isInertia:
				newdata=pd.DataFrame({'Type':['Inertia','MaxInertia','CostInertia'],'Unit':['MWs/MWA','MWs/MWA','EUR/MWs/MWA'],\
					'Zone':[cfg['partition'][partitionInertia], cfg['partition'][partitionInertia], cfg['partition'][partitionInertia]],\
					'model':['Parameter','Parameter','Parameter'],\
					'scenario':[cfg['listdatagroups'][0],cfg['listdatagroups'][0],cfg['listdatagroups'][0]],\
					'value':[Inertia,MaxInertia,CostInertia],'year':[cfg['year'],cfg['year'],cfg['year']]})

		# include timeseries names from TimeSeries dictionnary and compute scaling coefficient
		log_debug('include time series and compute scaling coefficients')
		for row in datapartition.index:
			mytype=datapartition.loc[row,'Type']
			if mytype in ListTypesZV:
				if mytype in timeseriesdict['ZV'].keys():
					filetimeserie=timeseriesdict['ZV'][datapartition.loc[row,'Type']][datapartition.loc[row,'Zone']]
					datapartition.loc[row, 'Profile_Timeserie']=filetimeserie

		columnsZV=['Type','Zone','value','Profile_Timeserie']
		if 'Profile_Timeserie' not in datapartition.columns:
			logger.warning('/!\ No profile time series provided')
			columnsZV.pop(-1)
		datapartition=datapartition[columnsZV]
		datapartition.to_csv(os.path.join(outputdir, cfg['csvfiles']['ZV_ZoneValues']), index=False)

	# treat global variables
	logger.info('\ntreat global variables ')
	globalvars=pd.Series(str)
	for datagroup in cfg['listdatagroups']:
		log_debug('\ntreat datagroup '+datagroup)
		if 'techno' in cfg['datagroups'][datagroup]['listvariables'].keys():
			for techno in cfg['datagroups'][datagroup]['listvariables']['techno'].keys():
				if 'global' in cfg['datagroups'][datagroup]['listvariables']['techno'][techno].keys():
					log_debug('there are global vars for:'+datagroup+',techno:'+techno)
					if type(cfg['datagroups'][datagroup]['listvariables']['techno'][techno]['global'])==list:
						for var in cfg['datagroups'][datagroup]['listvariables']['techno'][techno]['global']:
							for fuel in cfg['technos'][techno]:
								globalvars[var+fuel]=cfg['datagroups'][datagroup]['regions']['global'] 
					else:
						for key in cfg['datagroups'][datagroup]['listvariables']['techno'][techno]['global'].keys():
							for var in cfg['datagroups'][datagroup]['listvariables']['techno'][techno]['global'][key]['variables']:
								for fuel in cfg['datagroups'][datagroup]['listvariables']['techno'][techno]['global'][key]['fuels']:
									globalvars[var+fuel]=cfg['datagroups'][datagroup]['regions']['global'] 
		if 'coupling' in cfg['datagroups'][datagroup]['listvariables'].keys():
			if 'global' in cfg['datagroups'][datagroup]['listvariables']['coupling'].keys():
				for var in cfg['datagroups'][datagroup]['listvariables']['coupling']['global']:
					globalvars[var]=cfg['datagroups'][datagroup]['regions']['global']


	# create file TU_ThermalUnits
	###############################################################

	if cfg['csvfiles']['TU_ThermalUnits']:
		logger.info('\n\nCreating TU_ThermalUnits')
		listvar=[]
		isCO2=False
		isDynamic=cfg['ParametersCreate']['DynamicConstraints']
		isPrice = False
		if not 'thermal' in cfg['ParametersCreate'].keys() or not 'variablecost' in cfg['ParametersCreate']['thermal'].keys():
			logger.warning('/!\ No value was provided for parameter "ParametersCreate>thermal>variablecost" in settings file. Default value is False: data for variable cost is consider to be provided as a production cost and not a fuel price.')
		else:
			isPrice=(cfg['ParametersCreate']['thermal']['variablecost']=='Price')
		isMaintenance=False
		
		for variable in list(vardict['Input']['VarTU'].keys()):
			oevar=vardict['Input']['VarTU'][variable]
			if variable!='NumberUnits':
				listvar.append(oevar)
			
		v=0
		# loop on technos 
		for oetechno in cfg['technos']['thermal']:
			log_debug('\ntreat '+oetechno)
			TU=pd.DataFrame({'Name':oetechno,'region':listregions_dataset})
			TU=TU.set_index('region')
			for variable in vardict['Input']['VarTU']:
				log_debug('\ntreat variable '+variable)
				isFuel=False
				TreatVar=True
				if variable=='Price' and 'thermal' in cfg['ParametersCreate'] and cfg['ParametersCreate']['thermal']['variablecost']=='Price':
					# case where VariableCost is computed as Efficency*Price of fuel
					if oetechno in cfg['ParametersCreate']['thermal']['fuel']:
						fuel=cfg['ParametersCreate']['thermal']['fuel'][oetechno]
						varname=vardict['Input']['VarTU'][variable]+fuel
						isFuel=True
					else: 
						TreatVar=False
				else:	
					varname=vardict['Input']['VarTU'][variable]+oetechno
					
				
				if varname not in bigdata['variable'].unique():
					log_debug('variable '+varname+' not in dataset')
					if variable=='VariableCost':
						if 'thermal' in cfg['ParametersCreate'] and 'useVariableCost' in cfg['ParametersCreate']['thermal']:
							if oetechno in cfg['ParametersCreate']['thermal']['useVariableCost']:
								# the technology may does have a variable cost but can be replaced by the variable cost of another technology
								varname=vardict['Input']['VarTU'][variable]+cfg['ParametersCreate']['thermal']['useVariableCost'][oetechno]
								log_debug('using variable '+varname+' instead')
					else:
						TreatVar=False
					
				if TreatVar:
					log_debug('variable '+varname + ' is in dataset')
					
					isGlobal=False
					Global=0
					if vardict['Input']['VarTU'][variable]+oetechno in globalvars.index:
						log_debug('variable '+variable+' '+varname+' is global for region '+globalvars[vardict['Input']['VarTU'][variable]+oetechno]+' techno '+oetechno)
						isGlobal=True						
						
					
					if isGlobal:
						vardf=bigdata.filter(variable=varname,region=listGlobalReg).as_pandas(meta_cols=False)
						vardf=vardf.set_index('region')
						vardf=vardf.rename(columns={"value":variable})
						dataTU=vardf[variable]
					else:
						vardf=bigdata.filter(variable=varname,region=listregions_dataset).as_pandas(meta_cols=False)
						vardf=vardf.set_index('region')
						vardf=vardf.rename(columns={"value":variable})
						dataTU=vardf[variable]
					
					if isGlobal:
						Global=dataTU[globalvars[vardict['Input']['VarTU'][variable]+oetechno] ]
					
					if isFuel:
						if vardict['Input']['VarTU'][variable]+fuel in globalvars.index:
							log_debug('variable '+variable+' '+varname+' is global for region '+globalvars[vardict['Input']['VarTU'][variable]+fuel]+' fuel '+fuel)
							isGlobal=True
							Global=dataTU[globalvars[vardict['Input']['VarTU'][variable]+fuel] ]
					
					if isGlobal:
						new_dataTU=pd.Series(index=TU.index,name=dataTU.name)
						for reg in TU.index:
							if reg in dataTU.index: new_dataTU.loc[reg] = dataTU.loc[reg]
							else: new_dataTU.loc[reg] = Global
						TU=pd.concat([TU, new_dataTU], axis=1)
					else:
						TU=pd.concat([TU, dataTU], axis=1)

			if len(TU.index)>0:							
				if 'AvailabilityRate' in TU.columns:
					TU['AvailabilityRate'].fillna(1, inplace=True)
				else:
					TU['AvailabilityRate']=1
				TU=TU.fillna(value=0.0)
				
				if 'Capacity' not in TU.columns:
					TU['Capacity']=0
				if 'thermal' in cfg['ParametersCreate'] and 'CapacityFromEnergy' in cfg['ParametersCreate']['thermal'] and oetechno in cfg['ParametersCreate']['thermal']['CapacityFromEnergy']:
					if 'Energy' in TU.columns:
						TU['EnergyRate']=TU.apply(lambda row: row['Energy'] / (row['Capacity'] * 8760) if row['Capacity']>0 else 0, axis=1) 
					else:
						TU['EnergyRate']=0
				else:
					TU['EnergyRate']=1
				
				if 'VariableCost' not in TU.columns:
					if isInvest and ('thermal' in cfg['ParametersCreate']['CapacityExpansion']) and (oetechno in cfg['ParametersCreate']['CapacityExpansion']['thermal']) and ('VariableCost' in cfg['ParametersCreate']['CapacityExpansion']['thermal'][oetechno]):
						TU['VariableCost']=cfg['ParametersCreate']['CapacityExpansion']['thermal'][oetechno]['VariableCost']
					else:
						TU['VariableCost']=0
				
				# compute InvestmentCost as CapitalCost/Lifetime+FixedCost
				if 'CapitalCost' in TU.columns:
					if 'Lifetime' in TU.columns:
						TU['Lifetime']=TU['Lifetime'].apply(lambda x: x if x>1 else 1)
						TU['InvestmentCost']=TU['CapitalCost']/TU['Lifetime']
					else:
						TU['InvestmentCost']=TU['CapitalCost']
					if 'AnnualFixedCost' in TU.columns:
						TU['InvestmentCost']=TU['InvestmentCost']+TU['AnnualFixedCost']
							
				nbUnitsPerTechno = 1
				if not 'thermal' in cfg['ParametersCreate'].keys() or not 'NbUnitsPerTechno' in cfg['ParametersCreate']['thermal'].keys():
					logger.warning('/!\ No value was provided for parameter "ParametersCreate>thermal>NbUnitsPerTechno" in settings file. Default value is 1.')
				else:
					nbUnitsPerTechno = 	cfg['ParametersCreate']['thermal']['NbUnitsPerTechno']	
				if nbUnitsPerTechno==1:
					TU['MaxPower']=TU['Capacity']*TU['AvailabilityRate']*TU['EnergyRate']
				
				if not cfg['ParametersCreate']['DynamicConstraints']:
					TU['MinPower']=0.0			
			
				TU['NumberUnits']=1
				if 'thermal' in cfg['ParametersCreate']:
					if 'NbUnitsPerTechno' in cfg['ParametersCreate']['thermal']:
						if not cfg['ParametersCreate']['thermal']['NbUnitsPerTechno']==1:
							TU['NumberUnits']=np.ceil(TU['Capacity']*TU['AvailabilityRate']*TU['EnergyRate']/TU['MaxPower'])
				
				if 'CO2Rate' in TU.columns: 
					isCO2=True
				elif ('CO2Emission' in TU.columns) & ('Energy' in TU.columns): 
					TU['CO2Rate']=TU['CO2Emission']*TU['Energy'].apply(lambda x: 1/x if x>0 else 0)
					isCO2=True
				else: 
					TU['CO2Rate']=0.0
					
				# Case where Variable Cost is computed out of Efficiency and Price
				if isPrice and oetechno in cfg['ParametersCreate']['thermal']['fuel']: 
					if 'Efficiency' in TU.columns and 'Price' in TU.columns:
						TU['VariableCost']=TU['Efficiency']*TU['Price']
					elif 'Price' in TU.columns:
						TU['VariableCost']=TU['Price']

				TU.loc[ TU['MaxPower'] < cfg['ParametersCreate']['zerocapacity'], 'MaxPower' ]=0
				
				# case with investment optimisation
				if isInvest and ('thermal' in cfg['ParametersCreate']['CapacityExpansion']) and (oetechno in cfg['ParametersCreate']['CapacityExpansion']['thermal']):
					if 'InvestmentCost' not in TU.columns: # if no investment cost given use value from config file	
						if oetechno in cfg['ParametersCreate']['CapacityExpansion']['thermal']:
							if 'InvestmentCost' in cfg['ParametersCreate']['CapacityExpansion']['thermal'][oetechno]:
								TU['InvestmentCost']=cfg['ParametersCreate']['CapacityExpansion']['thermal'][oetechno]['InvestmentCost']
							else:
								TU['InvestmentCost']=0
						else:
							TU['InvestmentCost']=0
									
					# it is allowed to invest in oetechno
					isMaxAddConfig=False
					if 'MaxAdd' in cfg['ParametersCreate']['CapacityExpansion']['thermal'][oetechno]:
						# the user defined a maximum investment in the settings
						isMaxAddConfig=True
						maxaddconfig=cfg['ParametersCreate']['CapacityExpansion']['thermal'][oetechno]['MaxAdd']
					if 'MaxRet' in cfg['ParametersCreate']['CapacityExpansion']['thermal'][oetechno]:
						maxretconfig=cfg['ParametersCreate']['CapacityExpansion']['thermal'][oetechno]['MaxRet']
					else:
						maxretconfig=0
					
					if 'MaxCapacity' in TU.columns and isMaxAddConfig:
						# MaxCapacity, if >0 defines the existing investment potential
						TU['MaxAddedCapacity']= TU.apply(lambda row: row['MaxCapacity'] - row['Capacity'] if (row['Capacity'] < row['MaxCapacity'] and row['MaxCapacity'] < row['Capacity']+maxaddconfig) else maxaddconfig, axis=1)
					elif 'MaxCapacity' in TU.columns:
						TU['MaxAddedCapacity']= TU.apply(lambda row: row['MaxCapacity'] - row['Capacity'] if row['Capacity'] < row['MaxCapacity']  else 0, axis=1)
					elif isMaxAddConfig:
						TU['MaxAddedCapacity']= maxaddconfig
					else:
						TU['MaxAddedCapacity']=0
					TU['MaxRetCapacity']=maxretconfig
					
					newcapa=cfg['ParametersCreate']['zerocapacity']
					if 'NewCapacity' in cfg['ParametersCreate']['CapacityExpansion']['thermal'][oetechno]:
						newcapa=cfg['ParametersCreate']['CapacityExpansion']['thermal'][oetechno]['NewCapacity']
					
					TU.loc[ TU['MaxPower'] <= cfg['ParametersCreate']['zerocapacity'], 'MaxPower' ]=TU.apply(lambda row: newcapa if row['MaxAddedCapacity'] > 0  else 0, axis=1)
				
				elif isInvest:
					TU['MaxAddedCapacity']=0
					TU['MaxRetCapacity']=0
								
				RowsToDelete = TU[ TU['MaxPower'] == 0 ].index
				if len(RowsToDelete)>0: 
					log_debug('deleting rows '+', '.join(RowsToDelete))
				TU=TU.drop(RowsToDelete)
				if v==0:	
					BigTU=TU
					v=1
				else:
					BigTU=pd.concat([BigTU,TU])
		BigTU['Zone']=BigTU.index
		listTU=BigTU.columns.tolist()
		listcols= ['Zone','Name','NumberUnits','MaxPower','VariableCost']
		if 'Capacity' in listTU: listcols.append('Capacity')
		if 'debug' in cfg['ParametersCreate']: 
			if 'Energy' in listTU and cfg['ParametersCreate']['debug']: listcols.append('Energy')
		if 'CO2Rate' in listTU: listcols.append('CO2Rate')
		if 'Lifetime' in listTU: listcols.append('Lifetime')
		if 'CapitalCost' in listTU: listcols.append('CapitalCost')
		if 'AnnualFixedCost' in listTU: listcols.append('AnnualFixedCost')
		if 'FixedCost' in listTU: listcols.append('FixedCost')
		if 'MaxCapacity' in listTU: listcols.append('MaxCapacity')
		if 'AvailabilityRate' in listTU: listcols.append('AvailabilityRate')
		if 'EnergyRate' in listTU: listcols.append('EnergyRate')
		if isDynamic:
			if 'MinPower' in listTU: listcols.append('MinPower')
			if 'DeltaRampUp' in listTU: listcols.append('DeltaRampUp')
			if 'DeltaRampDown' in listTU: listcols.append('DeltaRampDown')
			if 'MinUpTime' in listTU: listcols.append('MinUpTime')
			if 'MinDownTime' in listTU: listcols.append('MinDownTime')
			if 'StartUpCost' in listTU: listcols.append('StartUpCost')
			if 'InitialPower' in listTU: listcols.append(['InitialPower'])
			if 'InitUpDownTime' in listTU: listcols.append(['InitUpDownTime'])

		if isPrice: 
			if 'Price' in listTU: listcols.append(['Price'])
			if 'Efficiency' in listTU : listcols.append(['Efficiency'])
		if isInertia and 'Inertia' in listTU: listcols.append(['Inertia'])
		if isPrimary and 'PrimaryRho' in listTU: listcols.append(['PrimaryRho'])
		if isSecondary and 'SecondaryRho' in listTU: listcols.append(['SecondaryRho'])
		if isMaintenance and 'MaxPowerProfile' in listTU: listcols.append(['MaxPowerProfile'])
		if 'Pauxiliary' in listTU: listcols.append(['Pauxiliary'])
		if 'QuadTerm' in listTU: listcols.append(['QuadTerm'])
		if isInvest and 'thermal' in cfg['ParametersCreate']['CapacityExpansion']:
			if 'InvestmentCost' in BigTU.columns: listcols.append('InvestmentCost')
			if 'MaxAddedCapacity' in BigTU.columns: listcols.append('MaxAddedCapacity')
			if 'MaxRetCapacity' in BigTU.columns: listcols.append('MaxRetCapacity')

		if len(BigTU)>0:
			BigTU=BigTU[ listcols ]
			
			BigTU=BigTU[ BigTU['Zone'].isin(cfg['partition'][partitionDemand]) ]
			BigTU=BigTU[BigTU.NumberUnits >0]
			BigTU=BigTU.fillna(0)
			BigTU.to_csv(os.path.join(outputdir, cfg['csvfiles']['TU_ThermalUnits']), index=False)

	# treat seasonal storage
	# filling sheet SS_SeasonalStorage and STS_ShortTermStorage
	############################################################
	AddedCapa=pd.DataFrame()
	if cfg['csvfiles']['SS_SeasonalStorage']:
		logger.info('\n\nCreating SS_SeasonalStorage')
		BigSS = pd.DataFrame()
		v=0
		for oetechno in	(cfg['technos']['reservoir'] if 'reservoir' in cfg['technos'].keys() else []):
			log_debug('\ntreat '+oetechno)
			SS=pd.DataFrame({'Name':oetechno,'region':listregions_dataset})
			SS=SS.set_index('region')
			
			for variable in vardict['Input']['VarSS']:
				varname=vardict['Input']['VarSS'][variable]+oetechno
				if varname in bigdata['variable'].unique():					
					isGlobal=False
					Global=0
					if varname in globalvars.index:
						log_debug('variable '+variable+' '+varname+' is global for region '+globalvars[varname]+' techno '+oetechno)
						isGlobal=True
					
					if isGlobal:
						vardf=bigdata.filter(variable=varname,region=listGlobalReg).as_pandas(meta_cols=False)
					else:
						vardf=bigdata.filter(variable=varname,region=listregions_dataset).as_pandas(meta_cols=False)
					vardf=vardf.set_index('region')					
					data=vardf[['value']]
					data=data.rename(columns={"value":variable})
					
					# treat global variables
					if isGlobal:
						Global=data[variable][globalvars[varname] ]
						
					SS=pd.concat([SS, data], axis=1)
					if isGlobal: SS[variable]=Global
							
			if len(SS)>0:
				# case when no maxvolume is provides
				isMaxVolume=False
				if 'MaxVolume' in SS.columns and not (SS==0).all()['MaxVolume']: 
					isMaxVolume=True
				else:
					multfactor = None
				if 'Volume2CapacityRatio' in cfg['ParametersCreate'].keys():
					multfactor=cfg['ParametersCreate']['Volume2CapacityRatio'][oetechno]
				if 'VolumeRate' not in SS.columns or SS['VolumeRate'].isna().any():
					if 'Volume2CapacityRatio' in cfg['ParametersCreate'].keys():
						logger.warning(f'Warning: no VolumeRate data for {oetechno}. Using value {multfactor} provided in setting file at ParametersCreate>Volume2CapacityRatio>{oetechno}.')
					else:
						logger.error(f'Error: no VolumeRate data for {oetechno}. Please provide a default value in setting file using ParametersCreate>Volume2CapacityRatio>{oetechno}.')
						log_and_exit(2, cfg['path'])
				if 'VolumeRate' in SS.columns:
					SS['VolumeRate'].fillna(multfactor/8760, inplace=True)
					SS['VolumeRate']=SS['VolumeRate']*8760
				else:
					SS['VolumeRate']=multfactor
				if not isMaxVolume and 'MaxPower' in SS.columns:				
					SS['MaxVolume']=SS['MaxPower']*SS['VolumeRate']
					logger.info('no Max Storage for '+oetechno+' replaced by MaxPower*VolumeRate')
				
				# replace low capacities with 0
				SS.loc[ SS['MaxPower'] < cfg['ParametersCreate']['zerocapacity'], 'MaxPower' ]=0
			
				# include inflows profiles: include timeseries names from TimeSeries dictionnary
				SS['InflowsProfile']=''
				SS['WaterValues']=''
				SS['Energy_Timeserie']=0
				if cfg['ParametersCreate']['reservoir']['coordinated']:
					SS['HydroSystem']=0
				else:
					hydrosystem=pd.Series([i for i in range(len(SS.index))])
					hydrosystem.index=SS.index
					SS['HydroSystem']=hydrosystem

				SS['Name']=oetechno
				SS['Zone']=SS.index
				SS['AddPumpedStorage']=0
				SS['AddPumpedStorageVolume']=0
				SS['InitialVolume']=0
				for row in SS.index:
					# treatment initial volume
					if cfg['aggregateregions']!=None:
						if row in cfg['aggregateregions']:
							if (row in cfg['aggregateregions'] and SS.loc[row,'MaxVolume']>0 and SS.loc[row,'MaxPower']>cfg['ParametersCreate']['reservoir']['minpowerMWh']):
								Rate=0
								N=0
								for reg in cfg['aggregateregions'][row]:
									if reg in cfg['ParametersCreate']['InitialFillingrate']:
										Rate=Rate+cfg['ParametersCreate']['InitialFillingrate'][reg]
										N=N+1
								if N>0:
									Rate=Rate/N
								SS.loc[row,'InitialVolume']=SS.loc[row,'MaxVolume']*Rate
						elif (SS.loc[row,'MaxVolume']>0 and SS.loc[row,'MaxPower']>cfg['ParametersCreate']['reservoir']['minpowerMWh']):
							SS.loc[row,'InitialVolume']=SS.loc[row,'MaxVolume']*cfg['ParametersCreate']['InitialFillingrate'][row]
					elif (SS.loc[row,'MaxVolume']>0 and SS.loc[row,'MaxPower']>cfg['ParametersCreate']['reservoir']['minpowerMWh']):
						SS.loc[row,'InitialVolume']=SS.loc[row,'MaxVolume']*cfg['ParametersCreate']['InitialFillingrate'][row]
							
					# treatment volume when no max storage is provided or MaxPower<minpowerMWh or no inflows
					if SS.loc[row,'MaxVolume']==0 or SS.loc[row,'MaxPower']<cfg['ParametersCreate']['reservoir']['minpowerMWh'] or SS.loc[row,'Inflows']==0 :
						# this storage is moved to Additionnal Pumped Storage
						SS.loc[row,'AddPumpedStorage']=SS.loc[row,'MaxPower']
						SS.loc[row,'AddPumpedStorageVolume']=SS.loc[row,'MaxVolume']
						if SS.loc[row,'MaxVolume']==0:
							logger.info('No volume for reservoir in region '+row+', adding capacity to Pumped Storage')
						elif SS.loc[row,'MaxPower']<cfg['ParametersCreate']['reservoir']['minpowerMWh']:
							logger.info('MaxPower in region '+row+' is very low, adding capacity to Pumped Storage')
						elif SS.loc[row,'Inflows']==0:
							logger.info('No inflows in region '+row+', adding capacity to Pumped Storage')
					# same treatment for regions excluded from seasonal storage
					if 'excluded_regions' in cfg['ParametersCreate']['reservoir']:
						if row in cfg['ParametersCreate']['reservoir']['excluded_regions']:
							SS.loc[row,'AddPumpedStorage']=SS.loc[row,'MaxPower']
							SS.loc[row,'AddPumpedStorageVolume']=SS.loc[row,'MaxVolume']
							logger.info('No volume for reservoir in region '+row+', adding capacity to Pumped Storage')
					# treatment inflows timeseries
					logger.info(f'Including inflow timeseries for {row}')
					if 'SS' not in timeseriesdict.keys() or timeseriesdict['SS'] is None or 'Inflows' not in timeseriesdict['SS'].keys() or row not in timeseriesdict['SS']['Inflows'].keys():
						logger.warning(f'/!\ No inflows provided for SS {row}')
					else:
						filetimeserie=timeseriesdict['SS']['Inflows'][row]
						SS.loc[row, 'InflowsProfile']=filetimeserie

				SS['NumberUnits']=SS.apply(lambda x: 1 if x['MaxVolume'] > 0 else 0, axis = 1)
				SS['Zone']=SS.index
				SS['MinVolume']=0
				SS['MinPower']=0
				if 'DischargingEfficiency' in SS.columns:
					SS['TurbineEfficiency']=SS['DischargingEfficiency']
				else:
					SS['TurbineEfficiency']=1.0
				SS['PumpingEfficiency']=0.0
				# create df for addedcapacity
				AddedCapa=SS[SS.AddPumpedStorage >0]
				
			# remove rows where MaxVolume=0 or where MaxPower < minPowerMWh	
			SS=SS.fillna(0)
			RowsToDelete = SS[ SS['MaxPower'] == 0 ].index
			SS=SS.drop(RowsToDelete)			
			SS = SS.drop(SS[SS.MaxVolume == 0].index)
			SS = SS.drop(SS[SS.AddPumpedStorage > 0].index)
			SS = SS.drop(SS[SS.MaxPower < cfg['ParametersCreate']['reservoir']['minpowerMWh']].index)
			
			# if SS is empty (=no seasonal storage) and requested in settings=> add a fake seasonal storage
			if len(SS)==0 and 'reservoir' in cfg['ParametersCreate'] and 'fakeunit' in cfg['ParametersCreate']['reservoir']:
				logger.info('Adding fictive reservoir')
				if 'region' in cfg['ParametersCreate']['reservoir']['fakeunit']:
					regfake=cfg['ParametersCreate']['reservoir']['fakeunit']['Region']
				else:
					regfake=listregions_dataset[0]
				if 'MaxPower' in cfg['ParametersCreate']['reservoir']['fakeunit']:
					maxpowerfake=cfg['ParametersCreate']['reservoir']['fakeunit']['MaxPower']
				else:
					maxpowerfake=100
				if 'MaxVolume' in cfg['ParametersCreate']['reservoir']['fakeunit']:
					maxvolfake=cfg['ParametersCreate']['reservoir']['fakeunit']['MaxVolume']
				else:
					maxvolfake=100
				FakeSS=pd.DataFrame(index=[fakereg],columns=['Maxpower,MaxVolume','NumberUnits','InflowsProfile'],data=[[maxpowerfake,maxvolumsfake,1,""]])
				SS=pd.concat([SS,FakeSS],axis=0)
			
			SS=SS.fillna(value=0.0)			
			BigSS=pd.concat([BigSS,SS])
		
		listSS=BigSS.columns.tolist()
		if len(listSS) == 0:
			logger.warning(f'No SS was found in data set. File {cfg["csvfiles"]["SS_SeasonalStorage"]} will not be created.')
		else:
			if 'debug' in cfg['ParametersCreate'] and cfg['ParametersCreate']['debug']:
				listcols= ['Name','Zone','HydroSystem','NumberUnits','MaxPower','MinPower',
					'MaxVolume','MinVolume','Inflows','InflowsProfile','InitialVolume', 
					'TurbineEfficiency','PumpingEfficiency','AddPumpedStorage']
			else:
				listcols= ['Name','Zone','HydroSystem','NumberUnits','MaxPower','MinPower',
					'MaxVolume','MinVolume','Inflows','InflowsProfile','InitialVolume', 
					'TurbineEfficiency','PumpingEfficiency']
	
			if isInertia and 'Inertia' in listSS: listcols.append('Inertia')
			if isPrimary and 'PrimaryRho' in listSS: listcols.append('PrimaryRho')
			if isSecondary and 'SecondaryRho' in listSS: listcols.append('SecondaryRho')
			if 'VolumeRate' in listSS: listcols.append('VolumeRate')

			BigSS=BigSS[ listcols ]
			BigSS=BigSS[ BigSS['Zone'].isin(cfg['partition'][partitionDemand]) ]
		if BigSS.empty:
			BigSS=pd.DataFrame(columns=['Name','Zone'])
		BigSS.to_csv(os.path.join(outputdir, cfg['csvfiles']['SS_SeasonalStorage']), index=False)

	# filling sheet STS_ShortTermStorage
	###############################################################	

	if cfg['csvfiles']['STS_ShortTermStorage']:	
		logger.info('\n\nCreating STS_ShortTermStorage')
		v=0
		STS=pd.DataFrame({'region':listregions_dataset})
		STS=STS.set_index('region')
		
		# treat short term hydro storage
		logger.info('	 Treating Hydro Pumped Storage')
		for oetechno in	cfg['technos']['hydrostorage']:
			log_debug('\ntreat '+oetechno)
			for variable in vardict['Input']['VarSTS|Hydro']:
				varname=vardict['Input']['VarSTS|Hydro'][variable]+oetechno
				
				isGlobal=False
				Global=0
				# treat global variables
				if varname in globalvars.index:
					log_debug('variable '+variable+' '+varname+' is global for region '+globalvars[varname]+' techno '+oetechno)
					isGlobal=True
								
				if isGlobal:
					vardf=bigdata.filter(variable=varname,region=listGlobalReg).as_pandas(meta_cols=False)
				else:
					vardf=bigdata.filter(variable=varname,region=listregions_dataset).as_pandas(meta_cols=False)
				vardf=vardf.set_index('region')
				data=vardf[['value']]
				data=data.rename(columns={"value":variable})
								
				if isGlobal:
					if len(data[variable])>0:
						Global=data[variable][globalvars[varname] ]

				STS=pd.concat([STS, data], axis=1)	

				if isGlobal: STS[variable]=Global
			
			if len(STS)>0:
				if 'AvailabilityRate' in STS.columns:
					STS['AvailabilityRate'].fillna(1, inplace=True)
				else:
					STS['AvailabilityRate']=1
							
				isPumpingEfficiency=False
				isChargingEfficiency=False
				isDischargingEfficiency=False
				if 'PumpingEfficiency' in STS.columns:
					if not STS['PumpingEfficiency'].isnull().all():
						isPumpingEfficiency=True
				if 'ChargingEfficiency' in STS.columns:
					if not STS['ChargingEfficiency'].isnull().all():
						isChargingEfficiency=True
				if 'DischargingEfficiency' in STS.columns:
					if not STS['DischargingEfficiency'].isnull().all():
						isDischargingEfficiency=True
				
				if (not isPumpingEfficiency) and isChargingEfficiency :
					STS['PumpingEfficiency']=STS['ChargingEfficiency']
					if isDischargingEfficiency :
						STS['TurbineEfficiency']=STS['DischargingEfficiency']
					else:
						STS['TurbineEfficiency']=1
				
				if (not isPumpingEfficiency) and (not isChargingEfficiency) :
					STS['PumpingEfficiency']=cfg['ParametersCreate']['PumpingEfficiency'][oetechno]
					STS['TurbineEfficiency']=1
								
				STS=STS.fillna(value=0.0)
				
				# replace low capacities with 0
				STS.loc[ STS['MaxPower'] < cfg['ParametersCreate']['zerocapacity'], 'MaxPower' ]=0
		
				# case with investment: replace 0 capacity with investment minimal capacity
				if cfg['ParametersCreate']['invest'] and 'hydrostorage' in cfg['ParametersCreate']['CapacityExpansion'] and oetechno in cfg['ParametersCreate']['CapacityExpansion']['hydrostorage']:
					newcapa=cfg['ParametersCreate']['zerocapacity']
					if 'NewCapacity' in cfg['ParametersCreate']['CapacityExpansion']['hydrostorage'][oetechno]:
						newcapa=cfg['ParametersCreate']['CapacityExpansion']['hydrostorage'][oetechno]['NewCapacity']
						STS.loc[ STS['MaxPower'] == 0, 'MaxPower' ]=newcapa
				
				# compute max storage or max power if not in data
				isMaxPower=False
				isMaxVolume=False
				if 'MaxPower' in STS.columns and not (STS==0).all()['MaxPower']: 
					isMaxPower=True
				if 'MaxVolume' in STS.columns and not (STS==0).all()['MaxVolume']: 
					isMaxVolume=True

				multfactor = None
				if 'Volume2CapacityRatio' in cfg['ParametersCreate'].keys() and oetechno in cfg['ParametersCreate']['Volume2CapacityRatio']:
					multfactor=cfg['ParametersCreate']['Volume2CapacityRatio'][oetechno]
				if ('Volume2CapacityRatio' not in STS.columns) or (STS['Volume2CapacityRatio'].isna().any()) or ((STS==0).all()['Volume2CapacityRatio']):
					if 'Volume2CapacityRatio' in cfg['ParametersCreate'].keys() and oetechno in cfg['ParametersCreate']['Volume2CapacityRatio']:
						logger.warning(f'Warning: no Volume2CapacityRatio data for {oetechno}. Using value {multfactor} provided in setting file at ParametersCreate>Volume2CapacityRatio>{oetechno}.')
					elif not ( ((STS==0).all()['MaxPower']) and ((STS==0).all()['MaxVolume'])):
						logger.error(f'Error: no Volume2CapacityRatio data for {oetechno}. Please provide a default value in setting file using ParametersCreate>Volume2CapacityRatio>{oetechno}.')
						log_and_exit(2, cfg['path'])
					else:
						multfactor=0

				if 'Volume2CapacityRatio' in STS.columns:
					STS['Volume2CapacityRatio'].fillna(multfactor, inplace=True)
				else:
					STS['Volume2CapacityRatio']=multfactor
				
				if isMaxVolume and (not isMaxPower): 
					STS['MaxPower']=STS['MaxVolume']/STS['Volume2CapacityRatio']	
					
				if isMaxPower and not isMaxVolume: 
					STS['MaxVolume']=STS['MaxPower']*STS['Volume2CapacityRatio']	
					
				STS['MaxPower']=STS['MaxPower']*STS['AvailabilityRate']
				# treat additionnal pumped storage from SS
				STS['AddPumpedStorage']=0
				STS['AddPumpedStorageVolume']=0
				STS['MinPower']=-1*STS['MaxPower']
				for row in STS.index:
					if row in AddedCapa.index:
						logger.info('Adding additionnal capacity: '+str(AddedCapa.loc[row,'AddPumpedStorage'])+' for Pumped Storage in region '+row)
						STS.loc[row,'MaxPower']=STS.loc[row,'MaxPower']+AddedCapa.loc[row,'AddPumpedStorage']
						STS.loc[row,'MaxVolume']=STS.loc[row,'MaxVolume']+AddedCapa.loc[row,'AddPumpedStorageVolume']
						STS.loc[row,'AddPumpedStorage']=AddedCapa.loc[row,'AddPumpedStorage']
						STS.loc[row,'AddPumpedStorageVolume']=AddedCapa.loc[row,'AddPumpedStorageVolume']
				RowsToDelete = STS[ STS['MaxPower'] == 0 ].index
				# Delete row with 0 capacity
				STS=STS.drop(RowsToDelete)
				STS['Name']=oetechno
				STS['Zone']=STS.index
				STS['NumberUnits']=1	
				STS['Inflows']=0	
				STS['MinVolume']=0	
				STS['InitialVolume']=0	
				STS['Zone']=STS.index
				if 'PrimaryRho' in STS.columns: STS['MaxPrimaryPower']=STS['MaxPower']*STS['PrimaryRho']
				if 'SecondaryRho' in STS.columns: STS['MaxSecondaryPower']=STS['MaxPower']*STS['SecondaryRho']
				STS['MinPowerCoef']=1.0
				STS['MaxPowerCoef']=1.0	

				if v==0:
					BigSTS=STS
					v=1
				else:
					BigSTS=pd.concat([BigSTS,STS])

		logger.info('	 Treating Batteries')

		# treat batteries
		for oetechno in	cfg['technos']['battery']:
			log_debug('\ntreat '+oetechno)
			BAT=pd.DataFrame({'region':listregions_dataset})
			BAT=BAT.set_index('region')
			for variable in vardict['Input']['VarSTS|Battery']:
				varname=vardict['Input']['VarSTS|Battery'][variable]+oetechno
				
				TreatVar=True
				if varname not in bigdata['variable'].unique():
					log_debug('variable '+varname+' not in dataset')
					TreatVar=False
				if TreatVar:				
					# treat global variables
					isGlobal=False
					Global=0
					if varname in globalvars.index:
						log_debug('variable '+variable+' '+varname+' is global for region '+globalvars[varname]+' techno '+oetechno)
						isGlobal=True
					
					if isGlobal:
						vardf=bigdata.filter(variable=varname,region=listGlobalReg).as_pandas(meta_cols=False)
					else:
						vardf=bigdata.filter(variable=varname,region=listregions_dataset).as_pandas(meta_cols=False)
					vardf=vardf.set_index('region')
					data=vardf[['value']]
					data=data.rename(columns={"value":variable})
					
					# treat global variables
					if isGlobal:
						if len(data[variable])>0:
							Global=data[variable][globalvars[varname] ]				
				
					BAT=pd.concat([BAT, data], axis=1)	
					if isGlobal: BAT[variable]=Global
			
			if len(BAT)>0:
				if 'AvailabilityRate' in BAT.columns:
					BAT['AvailabilityRate'].fillna(1, inplace=True)
				else:
					BAT['AvailabilityRate']=1
				isRoundTripEfficiency=False
				isChargingEfficiency=False
				isDischargingEfficiency=False
				if 'RoundTripEfficiency' in BAT.columns:
					if not BAT['RoundTripEfficiency'].isnull().all():
						isRoundTripEfficiency=True
				if 'PumpingEfficiency' in BAT.columns:
					if not BAT['PumpingEfficiency'].isnull().all():
						isChargingEfficiency=True
				if 'TurbineEfficiency' in BAT.columns:
					if not BAT['TurbineEfficiency'].isnull().all():
						isDischargingEfficiency=True
				if isRoundTripEfficiency:
					BAT['PumpingEfficiency']=BAT['RoundTripEfficiency']
					BAT['TurbineEfficiency']=BAT['RoundTripEfficiency']
				if (not isRoundTripEfficiency) and (not isDischargingEfficiency):
					if 'PumpingEfficiency' in cfg['ParametersCreate'] and oetechno in cfg['ParametersCreate']['PumpingEfficiency']:
						BAT['TurbineEfficiency']=cfg['ParametersCreate']['PumpingEfficiency'][oetechno]
					else:
						BAT['TurbineEfficiency']=1
				if (not isRoundTripEfficiency) and (not isChargingEfficiency):
					if 'PumpingEfficiency' in cfg['ParametersCreate'] and oetechno in cfg['ParametersCreate']['PumpingEfficiency']:
						BAT['PumpingEfficiency']=cfg['ParametersCreate']['PumpingEfficiency'][oetechno]
					else:
						BAT['PumpingEfficiency']=1
					
				BAT=BAT.fillna(value=0.0)
				isMaxPower=False
				isMaxVolume=False
				if 'MaxPower' in BAT.columns and not (BAT==0).all()['MaxPower']: 
					isMaxPower=True
				if 'MaxVolume' in BAT.columns and not (BAT==0).all()['MaxVolume']: 
					isMaxVolume=True
				if 'MaxPower' not in BAT.columns and 'MaxVolume' not in BAT.columns: 
					logger.info('no data for techno '+oetechno)
					BAT['Capacity']=0
					BAT['MaxPower']=0
					BAT['MaxVolume']=0
				

				multfactor = None
				if 'Volume2CapacityRatio' in cfg['ParametersCreate'].keys() and oetechno in cfg['ParametersCreate']['Volume2CapacityRatio']:
					multfactor=cfg['ParametersCreate']['Volume2CapacityRatio'][oetechno]
				if 'Volume2CapacityRatio' not in BAT.columns or ((BAT['Volume2CapacityRatio'].isna().any()) or ((BAT==0).all()['Volume2CapacityRatio'])):
					if 'Volume2CapacityRatio' in cfg['ParametersCreate'].keys() and oetechno in cfg['ParametersCreate']['Volume2CapacityRatio']:
						logger.warning(f'Warning: no Volume2CapacityRatio data for {oetechno}. Using value {multfactor} provided in setting file at ParametersCreate>Volume2CapacityRatio>{oetechno}.')
					elif not ( ((BAT==0).all()['MaxPower']) and ((BAT==0).all()['MaxVolume'])):
						logger.error(f'Error: no Volume2CapacityRatio data for {oetechno}. Please provide a default value in setting file using ParametersCreate>Volume2CapacityRatio>{oetechno}.')
						log_and_exit(2, cfg['path'])
					else:
						multfactor=0

				if 'Volume2CapacityRatio' in BAT.columns:
					BAT['Volume2CapacityRatio'].fillna(multfactor, inplace=True)
				else:
					BAT['Volume2CapacityRatio']=multfactor
				
				if isMaxVolume and not isMaxPower: 
					BAT['MaxPower']=BAT['MaxVolume']/BAT['Volume2CapacityRatio']	
					
				if isMaxPower and not isMaxVolume: 
					BAT['MaxVolume']=BAT['MaxPower']*BAT['Volume2CapacityRatio']	
				
				# replace low capacities with 0
				if isMaxPower: 
					BAT.loc[ BAT['MaxPower'] < cfg['ParametersCreate']['zerocapacity'], ['MaxPower','MaxVolume'] ]=0
					BAT['MaxPower']=BAT['MaxPower']*BAT['AvailabilityRate']
				elif isMaxVolume: 
					BAT.loc[ BAT['MaxVolume'] < cfg['ParametersCreate']['zerocapacity']*BAT['Volume2CapacityRatio'] , 'MaxVolume' ]=0
							
				# compute InvestmentCost as CapitalCost/Lifetime+AnnualFixedCost
				if 'CapitalCost' in BAT.columns:
					if 'Lifetime' in BAT.columns:
						BAT['Lifetime']=BAT['Lifetime'].apply(lambda x: x if x>1 else 1)
						BAT['InvestmentCost']=BAT['CapitalCost']/BAT['Lifetime']
					else:
						BAT['InvestmentCost']=BAT['CapitalCost']
					if 'AnnualFixedCost' in BAT.columns:
						BAT['InvestmentCost']=BAT['InvestmentCost']+BAT['AnnualFixedCost']
					
				# case with investment
				if isInvest and 'battery' in cfg['ParametersCreate']['CapacityExpansion']:
					if 'InvestmentCost' not in BAT.columns: # if no investment cost given use value from config file	
						if oetechno in cfg['ParametersCreate']['CapacityExpansion']['battery']:
							if 'InvestmentCost' in cfg['ParametersCreate']['CapacityExpansion']['battery'][oetechno]:
								BAT['InvestmentCost']=cfg['ParametersCreate']['CapacityExpansion']['battery'][oetechno]['InvestmentCost']
							else:
								BAT['InvestmentCost']=0
						else:
							BAT['InvestmentCost']=0
					if oetechno in cfg['ParametersCreate']['CapacityExpansion']['battery']:
						# replace 0 capacity with investment minimal capacity
						newcapa=cfg['ParametersCreate']['zerocapacity']
						if 'NewCapacity' in cfg['ParametersCreate']['CapacityExpansion']['battery'][oetechno]:
							newcapa=cfg['ParametersCreate']['CapacityExpansion']['battery'][oetechno]['NewCapacity']
						BAT.loc[ BAT['MaxPower'] == 0, 'MaxPower' ]=newcapa
						BAT.loc[ BAT['MaxVolume'] == 0, 'MaxVolume' ]=newcapa*cfg['ParametersCreate']['Volume2CapacityRatio'][oetechno]
						
						isMaxAddConfig=False
						if 'MaxAdd' in cfg['ParametersCreate']['CapacityExpansion']['battery'][oetechno]:
							# the user defined a maximum investment in the settings
							isMaxAddConfig=True
							maxaddconfig=cfg['ParametersCreate']['CapacityExpansion']['battery'][oetechno]['MaxAdd']
						if 'MaxRet' in cfg['ParametersCreate']['CapacityExpansion']['battery'][oetechno]:
							maxretconfig=cfg['ParametersCreate']['CapacityExpansion']['battery'][oetechno]['MaxRet']
						else:
							maxretconfig=0
						
						if 'MaxCapacity' in BAT.columns and isMaxAddConfig:						
							BAT['MaxAddedCapacity']= BAT.apply(lambda row: row['MaxCapacity'] - row['Capacity'] if row['Capacity'] < row['MaxCapacity'] and row['MaxCapacity'] < row['Capacity']+maxaddconfig else maxaddconfig, axis=1)
						elif 'MaxCapacity' in BAT.columns:
							BAT['MaxAddedCapacity']= BAT.apply(lambda row: row['MaxCapacity'] - row['Capacity'] if row['Capacity'] < row['MaxCapacity']  else 0, axis=1)
						elif isMaxAddConfig:
							BAT['MaxAddedCapacity']= maxaddconfig
						else:
							BAT['MaxAddedCapacity']=0
					else:
						BAT['MaxAddedCapacity']=0
						BAT['MaxRetCapacity']=0
									
				BAT['Name']=oetechno
				BAT['Zone']=BAT.index
				BAT['NumberUnits']=1	
				BAT['Inflows']=0	
									
				# case where Maximum Charge is not in the data
				if 'MinPower' not in BAT.columns or ('MinPower'  in BAT.columns and (BAT==0).all()['MinPower']): 
					BAT['MinPower']=-1*BAT['MaxPower']
				else:
					BAT['MinPower']=-1*BAT['MinPower']
					
				BAT['MinVolume']=0
				BAT['InitialVolume']=0			
				BAT['Zone']=BAT.index
				if 'PrimaryRho' in STS.columns: BAT['MaxPrimaryPower']=BAT['MaxPower']*BAT['PrimaryRho']
				if 'SecondaryRho' in STS.columns: BAT['MaxSecondaryPower']=BAT['MaxPower']*BAT['SecondaryRho']
				BAT['AddPumpedStorage']=0
				BAT['MinPowerCoef']=1.0
				BAT['MaxPowerCoef']=1.0
				
				RowsToDelete=[]
				RowsToDelete = BAT[ BAT['MaxPower'] == 0 ].index
				#if isMaxVolume: RowsToDelete = BAT[ BAT['MaxVolume'] == 0 ].index
				# Delete row with 0 capacity
				BAT=BAT.drop(RowsToDelete)
				
				if v==0:
					BigSTS=BAT
					v=1
				else:
					BigSTS=pd.concat([BigSTS,BAT])
				
		# create blank time series dataframe
		start2050=pd.to_datetime('2050-01-01T00:00+01:00')
		end2050=pd.to_datetime('2050-12-31T23:00+01:00')
		TimeIndex=pd.date_range(start=start2050,end=end2050, freq='1H')	
		
		# treat demand response 'load shifting'
		if 'demandresponseloadshifting' in cfg['technos'].keys():
			logger.info('	 Treating Load shifting')			
			DRTimeSeries=pd.DataFrame()
			
			ParticipationRate=pd.read_csv(cfg['ParametersCreate']['DemandResponseLoadShifting']['participationRateData'])
			ParticipationRate=ParticipationRate.groupby('countryname').mean()
				
			MaxDispatchName=vardict['Input']['VarSTS|DemandResponseLoadShifting']['MaxDispatch']
			MaxReductionName=vardict['Input']['VarSTS|DemandResponseLoadShifting']['MaxReduction']
			
			# loop on appliances
			for appliance in cfg['technos']['demandresponseloadshifting']:
				log_debug('\ntreat '+appliance)
				DRLS=pd.DataFrame({'region':cfg['partition'][partitionDemand]})
				DRLS=DRLS.set_index('region')
				EMminusE=pd.Series(index={'region':cfg['partition'][partitionDemand]})
				
				tshift=int(cfg['ParametersCreate']['DemandResponseLoadShifting']['tshift'][appliance])
				NumberBalancing=int(len(TimeIndex)/tshift)
				
				# loop on regions
				for reg in cfg['partition'][partitionDemand]:
					if cfg['ParametersCreate']['DemandResponseLoadShifting']['participationRate']:
						if cfg['aggregateregions']!=None:
							if reg in cfg['aggregateregions']:
								N=0.0
								PartRate=0.0
								for country in cfg['aggregateregions'][reg]:
									if country in ParticipationRate.index:
										N=N+1.0
										PartRate=PartRate+ParticipationRate.at[country,appliance]
								PartRate=PartRate/N
							else: 
								PartRate=ParticipationRate.at[reg,appliance]
						else: 
							PartRate=ParticipationRate.at[reg,appliance]
					else:
						PartRate=1.0
					MaxDispatchData=bigdata_SubAnnual.filter(variable=MaxDispatchName+appliance,region=reg).as_pandas(meta_cols=False).reset_index()
					MaxReductionData=bigdata_SubAnnual.filter(variable=MaxReductionName+appliance,region=reg).as_pandas(meta_cols=False).reset_index()
					
					# compute EM-E
					EMminusE[reg]=0
					Balancing=0
					NumberBalancingInData=int(12*24/tshift)
					for IndexBalancing in range(NumberBalancingInData):
						MaxDispatchData_IndexBalancing=MaxDispatchData[ (MaxDispatchData.index>=IndexBalancing*tshift) & (MaxDispatchData.index<=(IndexBalancing+1)*tshift-1) ]
						MaxReductionData_IndexBalancing=MaxReductionData[ (MaxReductionData.index>=IndexBalancing*tshift) & (MaxReductionData.index<=(IndexBalancing+1)*tshift-1) ]
						EMminusE_IndexBalancing=MaxDispatchData_IndexBalancing['value'].sum()  #-MaxReductionData_IndexBalancing['value'].sum()
						if EMminusE_IndexBalancing>	EMminusE[reg]: EMminusE[reg]=EMminusE_IndexBalancing
					EMminusE[reg]=EMminusE[reg]*PartRate
						
					# create MaxDispatch and MaxReduction time series
					# the data give a representative day for each month which has to be extended to all days
					MaxDispatchData['subannual']=pd.to_datetime(MaxDispatchData['subannual'],format='%m-%d %H:%S+01:00')+pd.DateOffset(years=150)
					MaxReductionData['subannual']=pd.to_datetime(MaxReductionData['subannual'],format='%m-%d %H:%S+01:00')+pd.DateOffset(years=150)
					MaxDispatch=pd.Series()
					MaxReduction=pd.Series()
					for month in range(12):
						numberdays=monthrange(int(cfg['year']),month+1)[1]

						# get data for the current month
						MaxDispatchDataFirstDay=MaxDispatchData[ MaxDispatchData.subannual.dt.month==month+1 ]
						MaxReductionDataFirstDay=MaxReductionData[ MaxReductionData.subannual.dt.month==month+1 ]
						
						# duplicate for all days of the month
						for day in range(numberdays):
							MaxDispatch=pd.concat([MaxDispatch,MaxDispatchDataFirstDay['value']])
							MaxReduction=pd.concat([MaxReduction,MaxReductionDataFirstDay['value']])
					
					MaxDispatch=np.maximum(MaxDispatch, 0)
					MaxDispatch=MaxDispatch*(1/cfg['ParametersCreate']['DemandResponseCoefficient'])
					NameSerie='DR__'+appliance+'__'+reg
					DRTimeSeries[NameSerie+'__MaxPower']=MaxReduction*PartRate
					DRTimeSeries[NameSerie+'__MinPower']=-MaxDispatch*PartRate  #-cfg['ParametersCreate']['DemandResponseCoefficient']
					DRTimeSeries[NameSerie+'__MaxDispatch']=MaxDispatch*PartRate
					DRTimeSeries[NameSerie+'__MinVolume']=0
					
					# include constraint MinVolume=EMminusE for all indexes =i*tshift
					DRTimeSeries[NameSerie+'__MinVolume'].mask( ((DRTimeSeries.index+1) % tshift) == 0, EMminusE[reg], inplace=True)

				DRLS['Name']=appliance
				DRLS['Zone']=DRLS.index
				DRLS['NumberUnits']=1	
				DRLS['Inflows']=0	
				DRLS['TurbineEfficiency']=1.0
				DRLS['PumpingEfficiency']=1.0
				DRLS['MinPower']='DR__'+appliance+'__'+DRLS.index+'__MinPower'	
				DRLS['MaxPower']='DR__'+appliance+'__'+DRLS.index+'__MaxPower'	
				DRLS['MinPowerCoef']=1.0
				DRLS['MaxPowerCoef']=1.0
				DRLS['MinVolume']='DR__'+appliance+'__'+DRLS.index+'__MinVolume'
				DRLS['MaxVolume']=EMminusE*cfg['ParametersCreate']['DemandResponseCoefficient']
				DRLS['VolumeLevelTarget']=EMminusE		
				DRLS['InitialVolume']=EMminusE		
				DRLS['Zone']=DRLS.index
				DRLS['MaxPrimaryPower']=0
				DRLS['MaxSecondaryPower']=0
				DRLS['Inertia']=0
				DRLS['AddPumpedStorage']=0
				
				if v==0:
					BigSTS=DRLS
					v=1
				else:
					BigSTS=pd.concat([BigSTS,DRLS])
					
			#####################################"
		
		if 'debug' in cfg['ParametersCreate']: 
			if cfg['ParametersCreate']['debug']:
				listcols= ['Name','Zone','NumberUnits','MaxPower','MaxVolume','TurbineEfficiency','PumpingEfficiency','MinPower','MinVolume','Inflows','InitialVolume','AddPumpedStorage']
			else:
				listcols= ['Name','Zone','NumberUnits','MaxPower','MaxVolume','TurbineEfficiency','PumpingEfficiency','MinPower','MinVolume']
		else:
			listcols= ['Name','Zone','NumberUnits','MaxPower','MaxVolume','TurbineEfficiency','PumpingEfficiency','MinPower','MinVolume']

		listSTS=BigSTS.columns.tolist()
		
		if isInertia and 'Inertia' in listSTS: listcols.append('Inertia')
		if isPrimary and 'MaxPrimaryPower' in listSTS: listcols.append('MaxPrimaryPower')
		if isSecondary and 'MaxSecondaryPower' in listSTS: listcols.append('MaxSecondaryPower')
		if 'VolumeLevelTarget' in listSTS: listcols.append('VolumeLevelTarget')
		if 'CO2Rate' in listSTS: listcols.append('CO2Rate')
		if 'Lifetime' in listSTS: listcols.append('Lifetime')
		if 'CapitalCost' in listSTS: listcols.append('CapitalCost')
		if 'Capacity' in listSTS: listcols.append('Capacity')
		if 'FixedCost' in listSTS: listcols.append('FixedCost')   
		if 'AnnualFixedCost' in listSTS: listcols.append('AnnualFixedCost')   
		if 'MaxCapacity' in listSTS: listcols.append('MaxCapacity')
		if 'AvailabilityRate' in listSTS: listcols.append('AvailabilityRate')
		if isInvest and 'battery' in cfg['ParametersCreate']['CapacityExpansion']:
			if 'MaxAddedCapacity' in BigSTS.columns: listcols.append('MaxAddedCapacity')
			if 'MaxRetCapacity' in BigSTS.columns: listcols.append('MaxRetCapacity')
			if 'InvestmentCost' in BigSTS.columns: listcols.append('InvestmentCost')

		if len(BigSTS)>0:
			BigSTS=BigSTS[ listcols ]
			BigSTS=BigSTS[ BigSTS['Zone'].isin(cfg['partition'][partitionDemand]) ]
			BigSTS=BigSTS.fillna(0)
			BigSTS.to_csv(os.path.join(outputdir, cfg['csvfiles']['STS_ShortTermStorage']), index=False)
		
	# treating res
	if cfg['csvfiles']['RES_RenewableUnits']:
		v=0
		logger.info('\n\nCreating RES_RenewableUnits')
		isVarRes=pd.Series()
		for oetechno in	cfg['technos']['res']+cfg['technos']['runofriver']:
			log_debug('\ntreat '+oetechno)
			RES=pd.DataFrame({'Name':oetechno,'region':listregions_dataset})
			RES=RES.set_index('region')
			for variable in vardict['Input']['VarRES']:
				varname=vardict['Input']['VarRES'][variable]+oetechno
				TreatVar=True
				if varname not in bigdata['variable'].unique():
					log_debug('variable '+varname+' not in dataset')
					TreatVar=False
				if TreatVar:					
					# treat global variables
					isGlobal=False
					Global=0
					if varname in globalvars.index:
						log_debug(str(varname)+' is global')
						isGlobal=True
					
					if isGlobal:
						vardf=bigdata.filter(variable=varname,region=listGlobalReg).as_pandas(meta_cols=False)
					else:
						vardf=bigdata.filter(variable=varname,region=listregions_dataset).as_pandas(meta_cols=False)
						
					if len(vardf.index)>0: 
						isVarRes[variable]=True
					else:
						isVarRes[variable]=False
					vardf=vardf.set_index('region')
					data=vardf[['value']]
					data=data.rename(columns={"value":variable})			
					
					if isGlobal:
						if len(data[variable])>0:
							Global=data[variable][globalvars[varname]]
					
					RES=pd.concat([RES, data], axis=1)	
					if isGlobal: RES[variable]=Global
				
			if len(RES)>0:
				if 'MaxPower' not in RES.columns:
					RES['MaxPower']=0
				
				RES=RES.fillna(value=0.0)
				
				# replace low capacities with 0
				RES.loc[ RES['MaxPower'] < cfg['ParametersCreate']['zerocapacity'], 'MaxPower' ]=0
				# case with investment
				# compute InvestmentCost as CapitalCost/Lifetime+AnnualFixedCost
				if 'CapitalCost' in RES.columns:
					if 'Lifetime' in RES.columns:
						RES['Lifetime']=RES['Lifetime'].apply(lambda x: x if x>1 else 1)
						RES['InvestmentCost']=RES['CapitalCost']/RES['Lifetime']
					else:
						RES['InvestmentCost']=RES['CapitalCost']
					if 'AnnualFixedCost' in RES.columns:
						RES['InvestmentCost']=RES['InvestmentCost']+RES['AnnualFixedCost']
						
				if isInvest and 'res' in cfg['ParametersCreate']['CapacityExpansion']:
					if 'InvestmentCost' not in RES.columns: # if no investment cost given use value from config file	
						if oetechno in cfg['ParametersCreate']['CapacityExpansion']['res']:
							if 'InvestmentCost' in cfg['ParametersCreate']['CapacityExpansion']['res'][oetechno]:
								RES['InvestmentCost']=cfg['ParametersCreate']['CapacityExpansion']['res'][oetechno]['InvestmentCost']
							else:
								RES['InvestmentCost']=0
						else:
							RES['InvestmentCost']=0
					if oetechno in cfg['ParametersCreate']['CapacityExpansion']['res']:
						# replace 0 capacity with investment minimal capacity
						newcapa=cfg['ParametersCreate']['zerocapacity']
						if 'NewCapacity' in cfg['ParametersCreate']['CapacityExpansion']['res'][oetechno]:
							newcapa=cfg['ParametersCreate']['CapacityExpansion']['res'][oetechno]['NewCapacity']
						RES.loc[ RES['MaxPower'] == 0, 'MaxPower' ]=newcapa
						maxaddconfig=cfg['ParametersCreate']['CapacityExpansion']['res'][oetechno]['MaxAdd']
						if 'Capacity' not in RES.columns:
							RES['Capacity']=0
						if 'MaxCapacity' in RES.columns:						
							RES['MaxAddedCapacity']= RES.apply(lambda row: row['MaxCapacity'] - row['Capacity'] if row['Capacity'] < row['MaxCapacity'] < row['Capacity']+maxaddconfig else maxaddconfig, axis=1)
						else:
							RES['MaxAddedCapacity']=maxaddconfig					
						RES['MaxRetCapacity']=cfg['ParametersCreate']['CapacityExpansion']['res'][oetechno]['MaxRet']
					else:
						RES['MaxAddedCapacity']=0
						RES['MaxRetCapacity']=0
				
				RowsToDelete = RES[ RES['MaxPower'] == 0 ].index
				# Delete row with 0 capacity
				RES=RES.drop(RowsToDelete)
				
				RES['Name']=oetechno
				RES['Zone']=RES.index
				RES['NumberUnits']=1	
				RES['MinPower']=0	
				RES['Gamma']=1	
				if 'PrimaryRho' in RES.columns and 'SecondaryRho' in RES.columns:
					RES['Gamma']=RES['PrimaryRho']+RES['SecondaryRho']
				
				# case of RoR: use energy instead of Capacity
				RES['Capacity']=RES['MaxPower']
				if 'MultFactor' in cfg['ParametersCreate'].keys() and oetechno in cfg['ParametersCreate']['MultFactor'].keys():
					RES['MaxPower']=RES['Energy']
					# The sum on 1 year of MaxPower*Profile should be equal to energy
					# Capacity*8760*r=energy				
					# MaxPower must be multiplied by Energy/(Capacity*8760)
				
				
				# include renewable potential profiles: include timeseries names from TimeSeries Disctionnary
				RES['MaxPowerProfile']=''
				RES['Energy_Timeserie']=0
				logger.info('\nIncluding RES time series')
				for row in RES.index:
					if 'RES' not in timeseriesdict.keys() or timeseriesdict['RES'] is None or oetechno not in timeseriesdict['RES']:
						logger.warning(f'/!\ No time series provided for RES {oetechno} for {row}')
						if 'InteractiveMode' in cfg and cfg['InteractiveMode']:
							query=input('Continue creation of plan4res input files? [y]/n\n')
							if query!='n':
								continue
					filetimeserie=timeseriesdict['RES'][oetechno][row]
					RES.loc[row, 'MaxPowerProfile']=filetimeserie
					logger.info(RES.loc[row, 'MaxPowerProfile'])
					
				if v==0:
					BigRES=RES
					v=1
				else:
					BigRES=pd.concat([BigRES,RES])
		
		listRES=BigRES.columns.tolist()
		listcols= ['Name','Zone','NumberUnits','MaxPower','MinPower','MaxPowerProfile','Energy','Capacity']
		if (isPrimary or isSecondary) and 'Gamma' in listRES and isInvest: listcols.append('Gamma')
		if isInertia and 'Inertia' in listRES: listcols.append('Inertia')
		if 'Lifetime' in listRES: listcols.append('Lifetime')
		if 'CapitalCost' in listRES: listcols.append('CapitalCost')
		if 'FixedCost' in listRES: listcols.append('FixedCost') 
		if 'AnnualFixedCost' in listRES: listcols.append('AnnualFixedCost') 
		if 'MaxCapacity' in listRES: listcols.append('MaxCapacity')
		if isInvest and 'res' in cfg['ParametersCreate']['CapacityExpansion']:
			if 'MaxAddedCapacity' in BigRES.columns: listcols.append('MaxAddedCapacity')
			if 'MaxRetCapacity' in BigRES.columns: listcols.append('MaxRetCapacity')
			if 'InvestmentCost' in BigRES.columns: listcols.append('InvestmentCost')

		if len(BigRES)>0:
			BigRES=BigRES[ listcols ]
			BigRES=BigRES[ BigRES['Zone'].isin(cfg['partition'][partitionDemand]) ]
			BigRES=BigRES.fillna(0)
			BigRES.to_csv(os.path.join(outputdir, cfg['csvfiles']['RES_RenewableUnits']), index=False)

log_and_exit(0, cfg['path'])
