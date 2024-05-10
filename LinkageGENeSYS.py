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


path = os.environ.get("PLAN4RESROOT")
nbargs=len(sys.argv)
if nbargs>1: 
	settings=sys.argv[1]
else:
	settings="settingsLinkageGENeSYS.yml"

cfg={}
with open(path+settings,"r") as mysettings:
	cfg=yaml.load(mysettings,Loader=yaml.FullLoader)

# replace name of current dataset by name given as input
if nbargs>1:
	namedataset=sys.argv[2]
	if cfg['USEPLAN4RESROOT']: cfg['path']='/data/local/'+namedataset+'/'
	else: cfg['path']=cfg['path'].replace(cfg['path'].split('/')[len(cfg['path'].split('/'))-2],namedataset)
	
if 'configDir' not in cfg: cfg['configDir']=cfg['path']+'settings/'
if 'genesys_inputpath' not in cfg: cfg['genesys_inputpath']=cfg['path']+'genesys_inputs/'
if 'genesys_resultspath' not in cfg: cfg['genesys_resultspath']=cfg['path']+'genesys_inputs/'
if 'timeseriespath' not in cfg: cfg['timeseriespath']=cfg['path']+'TimeSeries/'
if 'mappingspath' not in cfg: cfg['mappingspath']=cfg['path']+'settings/mappings_genesys/'
if 'outputpath' not in cfg: cfg['outputpath']=cfg['path']+'IAMC/'
if 'outputfile' not in cfg: cfg['outputfile']=namedataset+'.csv'

cfg['outputpath']=path+cfg['outputpath']
cfg['timeseriespath']=path+cfg['timeseriespath']
cfg['configDir']=path+cfg['configDir']
cfg['genesys_inputpath']=path+cfg['genesys_inputpath']
cfg['genesys_resultspath']=path+cfg['genesys_resultspath']
cfg['mappingspath']=path+cfg['mappingspath']
outputfile=cfg['outputpath']+cfg['outputfile']
if not os.path.isdir(cfg['outputpath']):os.mkdir(cfg['outputpath'])
if not os.path.isdir(cfg['timeseriespath']):os.mkdir(cfg['timeseriespath'])
if os.path.exists(outputfile):
	os.remove(outputfile)

print('create IAMC file for GENeSYS-MOD outputs in ',cfg['genesys_inputpath'])

# loop on the different variables
BigOut=pd.DataFrame()
firstVar=True

# open datafiles
data=pd.Series()
print('read ',cfg['genesys_inputpath']+cfg['genesys_datafiles']['input'])
for sheet in cfg['genesys_datafiles']['input_sheets']:
	print('  sheet ',sheet)
	data.loc['input_'+sheet]=pd.read_excel(cfg['genesys_inputpath']+cfg['genesys_datafiles']['input'],sheet_name=sheet).fillna(0)

for file in cfg['genesys_datafiles']:
	if (file!='input') and (file!='input_sheets'):
		print('read ',cfg['genesys_inputpath']+cfg['genesys_datafiles'][file])
		data.loc[file]=pd.read_csv(cfg['genesys_inputpath']+cfg['genesys_datafiles'][file])

# read mappings
mappings=pd.Series()
print('read mappings')
for mapping in cfg['mappings']:
	print('   mapping ',mapping,' in ',cfg['mappingspath']+cfg['mappings'][mapping])
	mappings.loc[mapping]=pd.read_csv(cfg['mappingspath']+cfg['mappings'][mapping],index_col=0,header=None)
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
print('regions in dataset ',regions)
print('interco in dataset ',regions_interco)

Yearsdf=pd.Series(data.loc['input_Sets']['Year'])
Yearsdf=Yearsdf.drop(Yearsdf.loc[Yearsdf ==0].index,axis=0 ).astype(int)
Years=Yearsdf.to_list()
print('years in dataset ',Years)
for var in cfg['variables']:
	isInternal=False
	print('treat ',var)
	
	if 'Network' in var: print(cfg['variables'][var]['unit'])
	if 'source' in cfg['variables'][var]:
		if cfg['variables'][var]['source']=='internal':
			isInternal=True
		else:
			# get data
			if cfg['variables'][var]['source']!='input':
				vardata=pd.DataFrame(data.loc[cfg['variables'][var]['source']])
			else:
				firstSheet=True
				for sheet in cfg['variables'][var]['sheets']:
					vardatasheet=pd.DataFrame(pd.DataFrame(data.loc[cfg['variables'][var]['source']+'_'+sheet]))
					if firstSheet: 
						vardata=pd.DataFrame(vardatasheet)
						firstSheet=False
					else:
						vardata=pd.concat([vardata,vardatasheet],axis=0)
	elif 'sources' in cfg['variables'][var]:
		if cfg['variables'][var]['sources']=='input':
			print('input cannot be in multiple source')
			exit()
		else:
			firstFile=True
			for file in cfg['variables'][var]['sources']:
				print(' read ',file)
				if firstFile:
					vardata=pd.DataFrame(data.loc[file])
	
					if 'Unit' not in vardata.columns: vardata['Unit']=cfg['variables'][var]['unit']
					firstFile=False
				else:
					vardata=pd.concat([vardata,pd.DataFrame(data.loc[file])],axis=0)
	
	colsdata=[]
	if 'Unit' not in vardata.columns:
		vardata['Unit']=cfg['variables'][var]['unit']

	for rulecat in cfg['variables'][var]['rules']:
		print('   apply ',rulecat)
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
			#dict={k : v for k,v in fullmapping.values}
			
			vardata['Variable']=vardata[colmap].map(lambda a: dict[a])
	
			# compute variable
			ruleagg=str(cfg['variables'][var]['rules'][rulecat]['rule'])
			colsKeep=[]
			for coldata in vardata.columns:
				#if coldata in IAMCcols+colsdata:
				if coldata in IAMCcols:
					colsKeep.append(coldata)
			vardata=vardata[ colsKeep ]

			colsToAggr=[]			
			#vardata=vardata.drop(columns=colmap)
			for coldata in vardata.columns:
				if coldata != 'Value' and coldata not in colsToAggr:
					colsToAggr.append(coldata)
			vardata['Unit']=vardata['Unit'].fillna(cfg['variables'][var]['unit'])
			vardata=vardata.groupby(colsToAggr).agg(ruleagg).reset_index()
	
		elif rulecat=='addyear':
			firstYear=True
			for year in Years:
				vardatayear=pd.DataFrame(vardata)
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
				vardatamap=pd.DataFrame(vardata[ vardata[col].isin(list(mappings.loc[map].index)) ])
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
			#dict=mappings.loc[map].transpose().to_dict()
			#dict={k : v for k,v in mappings.loc[map].values}
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
			#if colmap in vardata.columns: vardata=vardata.drop(columns=colmap)
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
				print('        apply ',subrule)
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
					#valuedict={k : v for k,v in newvalue.values}
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
				for component in cfg['variables'][var]['rules'][rulecat]['from']:
					if component[-1]=='|':
						# add mapping list to variable name
						for elem in map.index:
							if component[-1]+elem not in listComponents:
								listComponents.append(component+elem)
					else:
						listComponents.append(component)						
				vardata=out[ out['Variable'].isin(listComponents) ]	
				colsKeep=[]
				for col in vardata.columns:
					if col in IAMCcols:
						colsKeep.append(col)
				vardata=vardata[ colsKeep ]
				if 'ruleaggr' in cfg['variables'][var]['rules']['compute']:
					ruleagg=cfg['variables'][var]['rules']['compute']['ruleaggr']
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
			globaldata=pd.DataFrame(vardata)
			globalreg=cfg['global_region']
			isFirstRegion=True
			# case of interconnection variable
			regions_use=[globalreg]
			if 'Network' in var:				
				regions_use=regions_interco
			for region in regions_use:
				globaldata['Region']=region
				if isFirstRegion:
					vardataout=pd.DataFrame(globaldata)
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
						vardatayear=pd.DataFrame(vardata)
						vardatayear['Value']=vardatayear[year]
						vardata=pd.concat([vardata,vardatayear],axis=0)

		#fill scenario
		if 'PathwayScenario' in vardata.columns:
			vardata['Scenario']=vardata['PathwayScenario']
		
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
	else: print('empty data')
 
	#vardata.to_csv(cfg['outputpath']+var.replace('|','_').replace(',','').replace('.','').replace(' ','')+'_'+cfg['outputfile'])
	
	#IAM=pyam.IamDataFrame(vardata)
	#IAM.to_csv(cfg['outputpath']+'IAMC_'+var.replace('|','_').replace(',','').replace('.','').replace(' ','')+'_'+cfg['outputfile'])

out.to_csv(cfg['outputpath']+'csv_'+cfg['outputfile'])	

# check for duplicated and output synthesis of data
print('scenarios in data ',out['Scenario'].unique())
print('models in data ',out['Model'].unique())
print('regions in data ',out['Region'].unique())
print('variables in data ',out['Variable'].unique())
print('years in data ',out['Year'].unique())


duplicates=out.duplicated()
duprows=[]
for row in duplicates.index:
	if duplicates.loc[row]==True:
		duprows.append(row)
if len(duprows)>0:
	print('there are duplicated rows')
	print(duprows)
	duplicated_rows=out.loc[duprows]
	print('  for variables',duplicated_rows['Variable'].unique())
	#out_duplicated=out.loc[
else:
	print('no duplicated rows')
	
df=pd.read_csv(cfg['outputpath']+'csv_'+cfg['outputfile'],index_col=0)
 
# conversion to IAMDataFrame
BigIAM=pyam.IamDataFrame(df)

BigIAM.validate(exclude_on_fail=True)
BigIAM.to_excel(cfg['outputpath']+cfg['outputfile'].replace('csv','xlsx'))





