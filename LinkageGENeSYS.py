#!/usr/bin/env python
# -*- coding: utf-8 -*-

## Import packages
import pyam
import pandas as pd ## necessary data analysis package
import numpy as np
import os
import sys
import yaml
#import nomenclature as nc
from math import ceil
from datetime import timedelta
from calendar import monthrange

from p4r_python_utils import *

path = os.environ.get("PLAN4RESROOT")
print('path ',path)
nbargs=len(sys.argv)
if nbargs>1: 
	settings=sys.argv[1]
else:
	settings="settingsLinkageGENeSYS.yml"
if os.path.abspath(settings):
	settings = os.path.relpath(settings, path)

cfg={}
with open(os.path.join(path, settings),"r") as mysettings:
	cfg=yaml.load(mysettings,Loader=yaml.FullLoader)

# replace name of current dataset by name given as input
if nbargs>1:
	namedataset=sys.argv[2]
	if cfg['USEPLAN4RESROOT']: 
		cfg['path']=os.path.join(path, 'data/local', namedataset)
	else: cfg['path']=cfg['path'].replace(cfg['path'].split('/')[len(cfg['path'].split('/'))-2],namedataset)
	
if 'configDir' not in cfg: cfg['configDir']=os.path.join(cfg['path'], 'settings/')
if 'genesys_inputpath' not in cfg: cfg['genesys_inputpath']=os.path.join(cfg['path'], 'genesys_inputs/')
if 'genesys_resultspath' not in cfg: cfg['genesys_resultspath']=os.path.join(cfg['path'], 'genesys_inputs/')
if 'timeseriespath' not in cfg: cfg['timeseriespath']=os.path.join(cfg['path'], 'TimeSeries/')
if 'mappingspath' not in cfg: cfg['mappingspath']=os.path.join(cfg['path'], 'settings/mappings_genesys/')
if 'outputpath' not in cfg: cfg['outputpath']=os.path.join(cfg['path'], 'IAMC/')
if 'outputfile' not in cfg: cfg['outputfile']=namedataset+'.csv'

cfg['outputpath']=os.path.join(path,cfg['outputpath'])
cfg['timeseriespath']=os.path.join(path,cfg['timeseriespath'])
cfg['configDir']=os.path.join(path,cfg['configDir'])
cfg['genesys_inputpath']=os.path.join(path,cfg['genesys_inputpath'])
cfg['genesys_resultspath']=os.path.join(path,cfg['genesys_resultspath'])
cfg['mappingspath']=os.path.join(path,cfg['mappingspath'])
cfg['outputfile']=os.path.join(cfg['outputpath'],cfg['outputfile'])
if not os.path.isdir(cfg['outputpath']):os.mkdir(cfg['outputpath'])
if not os.path.isdir(cfg['timeseriespath']):os.mkdir(cfg['timeseriespath'])
if os.path.exists(cfg['outputfile']):
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

if treatFix:
	logger.info('create IAMC file for GENeSYS-MOD outputs in '+cfg['genesys_inputpath'])

	# loop on the different variables
	BigOut=pd.DataFrame()
	firstVar=True

	# open datafiles
	data=pd.Series()
	logger.info('read '+os.path.join(cfg['genesys_inputpath'],cfg['genesys_datafiles']['input']))
	for sheet in cfg['genesys_datafiles']['input_sheets']:
		logger.info('  sheet '+sheet)
		data.loc['input_'+sheet]=pd.read_excel(os.path.join(cfg['genesys_inputpath'],cfg['genesys_datafiles']['input']),sheet_name=sheet).fillna(0)

	for file in cfg['genesys_datafiles']:
		if (file!='input') and (file!='input_sheets'):
			logger.info('read '+os.path.join(cfg['genesys_inputpath'],cfg['genesys_datafiles'][file]))
			data.loc[file]=pd.read_csv(os.path.join(cfg['genesys_inputpath'],cfg['genesys_datafiles'][file]))

	# read mappings
	mappings=pd.Series()
	logger.info('read mappings')
	for mapping in cfg['mappings']:
		logger.info('   mapping '+mapping+' in '+os.path.join(cfg['mappingspath'],cfg['mappings'][mapping]))
		mappings.loc[mapping]=pd.read_csv(os.path.join(cfg['mappingspath'],cfg['mappings'][mapping]),index_col=0,header=None)
		rows_to_remove=[elem for elem in mappings.loc[mapping].index if str(elem)[0]=='#']
		mappings.loc[mapping].drop(rows_to_remove,inplace=True)
	out=pd.DataFrame()
	isFirst=True
	IAMCcols=['Model','Scenario','Region','Variable','Unit','Year','Value']
	colsAgg=['Region','PathwayScenario','Year','Unit']

	regions=[]
	regions_source=data.loc['input_Sets']['Region']
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
	logger.info('regions in dataset '+str(regions))
	logger.info('interco in dataset '+str(regions_interco))

	Yearsdf=pd.Series(data.loc['input_Sets']['Year'])
	Yearsdf=Yearsdf.drop(Yearsdf.loc[Yearsdf ==0].index,axis=0 ).astype(int)
	Years=Yearsdf.to_list()
	logger.info('years in dataset '+', '.join([str(y) for y in Years]))
	for var in cfg['variables']:
		isInternal=False
		logger.info('treat '+var)
		
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
				log_and_exit(1, cfg['path'])
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
		
		if 'Region' in vardata.columns:
			vardata=vardata[ vardata['Region'].isin(regions) ]
		if 'Region2' in vardata.columns:
			vardata=vardata[ vardata['Region'].isin(regions) ]
		
		
		#if 'Unit' not in vardata.columns:
		vardata['Unit']=cfg['variables'][var]['unit']

		for rulecat in cfg['variables'][var]['rules']:
			logger.info('   apply '+rulecat)
			if rulecat=='selectAndMap':
				# select rows 
				colmap=cfg['variables'][var]['rules'][rulecat]['column']
				if colmap not in colsdata: colsdata.append(colmap)
				firstMap=True
				for map in cfg['variables'][var]['rules'][rulecat]['mappings']:
					mappingpart=mappings.loc[map]
					if firstMap:
						fullmapping=mappingpart
						firstMap=False
					else:
						fullmapping=pd.concat([fullmapping,mappingpart],axis=0)
				vardata=vardata[ vardata[colmap].isin(list(fullmapping.index)) ]
				
				# create variable name
				dict={fullmapping.index[i]: fullmapping.iloc[i,0] for i in range(len(fullmapping.index))}			
				vardata['Variable']=vardata[colmap].map(lambda a: dict[a])
				
				# compute variable
				ruleagg=str(cfg['variables'][var]['rules'][rulecat]['rule'])
				colsKeep=[]
				for coldata in vardata.columns:
					if coldata in IAMCcols:
						colsKeep.append(coldata)
				vardata=vardata[ colsKeep ]
				
				colsToAggr=[]			
				for coldata in vardata.columns:
					if coldata != 'Value' and coldata not in colsToAggr:
						colsToAggr.append(coldata)
				if 'Year' in vardata.columns:
					vardata['Year']=vardata['Year'].astype(int)
				vardata=pd.DataFrame(data=pd.DataFrame(data=vardata).groupby(colsToAggr).agg(ruleagg).reset_index())

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
				for colselect in cfg['variables'][var]['rules'][rulecat]:
					values=cfg['variables'][var]['rules'][rulecat][colselect]['values']
					vardata=vardata[ vardata[colselect].isin(values) ]
		
			elif rulecat=='group':
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
					logger.info('        apply '+subrule)
					if 'source' in cfg['variables'][var]['rules'][rulecat][subrule]:
						if cfg['variables'][var]['rules'][rulecat][subrule]['source']!='input':
							newdata=data.loc[cfg['variables'][var]['rules'][rulecat][subrule]['source']]
						else:
							newdata=data.loc[cfg['variables'][var]['rules'][rulecat][subrule]['source']+'_'+cfg['variables'][var]['rules'][rulecat][subrule]['sheet']]
					if 'select' in cfg['variables'][var]['rules'][rulecat][subrule]:
						for colselect in cfg['variables'][var]['rules'][rulecat][subrule]['select']:
							values=cfg['variables'][var]['rules'][rulecat][subrule]['select'][colselect]['values']
							newdata=newdata[ newdata[colselect].isin(values) ]
					if subrule=='mapAndAddCols':					
						colref=cfg['variables'][var]['rules'][rulecat][subrule]['column']
					
						for newcol in cfg['variables'][var]['rules'][rulecat][subrule]['mappings']:
							colmap=cfg['variables'][var]['rules'][rulecat][subrule]['mappings'][newcol]
							combinedmap=newdata[[colref,colmap]].groupby([colref]).first().reset_index()
							combineddict={combinedmap.iloc[i,0]: combinedmap.iloc[i,1] for i in range(len(combinedmap.index))}
							vardata[newcol]=vardata[colref].map(lambda a: combineddict[a] if a in combineddict.keys() else 'None')
						if 'product_cols' in cfg['variables'][var]['rules'][rulecat][subrule]:
							for col in cfg['variables'][var]['rules'][rulecat][subrule]['product_cols']: 
								col2=cfg['variables'][var]['rules'][rulecat][subrule]['product_cols'][col]
								vardata[col]=vardata[col]*vardata[cfg['variables'][var]['rules'][rulecat][subrule]['product_cols'][col]]
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
			elif rulecat=='compute':
				if isInternal:
					if 'mapping' in cfg['variables'][var]['rules'][rulecat]:
						map=mappings.loc[cfg['variables'][var]['rules'][rulecat]['mapping']]
						dict={map.index[i]: map.iloc[i,0] for i in range(len(map.index))}
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
								if component[-1]+elem not in listComponents:
									#listComponents.append(component+elem)
									if 'ruleaggr' in cfg['variables'][var]['rules']['compute']:
										listComponents.append(component+dict[elem])
									else:
										listComponents.append(component+elem)
						else:
							listComponents.append(component)
						firstComponent=False
					vardata=out[ out['Variable'].isin(listComponents) ]	
					namezz='ZZZ_'+var.replace('|','_').replace('/','_')+'vardata.csv'
					if 'Capacity' in var: vardata.to_csv(namezz)
					#if 'Capacity' in var: vardataelem.to_csv(namezz)
					colsKeep=[]
					for col in vardata.columns:
						if col in IAMCcols:
							colsKeep.append(col)
					vardata=vardata[ colsKeep ]
					if 'ruleaggr' in cfg['variables'][var]['rules']['compute']:
						ruleagg=cfg['variables'][var]['rules']['compute']['ruleaggr']
						if isManyVar:
							firstElem=True
							for elem in listElem:
								listpossible=[el+str(dict[elem]) for el in cfg['variables'][var]['rules'][rulecat]['from']]
								vardataelem = pd.DataFrame(vardata[vardata['Variable'].isin(listpossible)]).reset_index().drop(columns='index')
								colsToAggr=[]
								vardataelem=vardataelem.drop(columns='Variable')
								for col in vardataelem.columns:
									if col != 'Value':
										colsToAggr.append(col)
								if len(vardataelem.index)>0:
									vardataelem=vardataelem.groupby(colsToAggr).agg(ruleagg).reset_index()
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


			#fill scenario
			if 'PathwayScenario' in vardata.columns:
				vardata['Scenario']=vardata['PathwayScenario']
				
			vardata['Unit']=cfg['variables'][var]['unit']
			
			#select columns to keep
			colsKeep=[]
			for col in vardata.columns:
				if col in IAMCcols:
					colsKeep.append(col)
			vardata=vardata[ colsKeep ]
			
			#fill missing columns
			for col in IAMCcols:
				if col not in colsKeep:
					vardata[col]=cfg['defaultvalues'][col]

			vardata['Year']=vardata['Year'].astype(int)
			
			if isFirst:
				out=vardata
				isFirst=False
			else:
				out=pd.concat([out,vardata],axis=0,ignore_index=True)
		
		else: logger.info('empty data')
		
	out.to_csv(cfg['outputfile'])	

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
		
	df=pd.read_csv(cfg['outputfile'],index_col=0)
	 
	# conversion to IAMDataFrame
	BigIAM=pyam.IamDataFrame(df)
	logger.info('converting units')
	for var_unit in cfg['conversions']:
		logger.info('converting '+str(var_unit)+' to '+str(cfg['conversions'][var_unit]['to']))
		if 'factor' in cfg['conversions'][var_unit]:
			BigIAM=BigIAM.convert_unit(var_unit, to=cfg['conversions'][var_unit]['to'], factor=float(cfg['conversions'][var_unit]['factor']),inplace=False) 
		else:
			BigIAM=BigIAM.convert_unit(var_unit, to=cfg['conversions'][var_unit]['to'],inplace=False) 

	#filter on unwanted variables
	logger.info('filtering on variables')
	logger.warning('excluding: '.join([str(_) for _ in cfg['removed_variables']]))
	variable_list=list(df['Variable'].unique())
	new_variable_list=[item for item in variable_list if item not in cfg['removed_variables']]
	BigIAM=BigIAM.filter(variable=cfg['removed_variables'], keep=False)

	#filter on unwanted variables
	logger.info('validating')
	BigIAM.validate(exclude_on_fail=True)
	BigIAM.to_excel(cfg['outputfile'].replace('csv','xlsx'))

if treatHourly:
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
	TimeSeriesTemplate=pd.read_csv(os.path.join(cfg['timeseriespath'],'Example.csv'))
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
	for sheet in cfg['genesys_timeseriesfiles']['timeseries_sheets']:
		sheetname=cfg['genesys_timeseriesfiles']['timeseries_sheets'][sheet]
		logger.info('  sheet '+sheet)
		df=pd.read_excel(os.path.join(cfg['genesys_inputpath'],cfg['genesys_timeseriesfiles']['xlsx']),sheet_name=sheetname,index_col=0).fillna(0)
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
				timeseries.to_csv(os.path.join(cfg['timeseriespath'],nameserie),index=False)
