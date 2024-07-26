#!/usr/bin/env python
# -*- coding: utf-8 -*-

## Import packages
import pandas as pd ## necessary data analysis package
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
import numpy as np
import yaml
import os
import math
import os.path
import pyam
from itertools import product
import sys
from shapely import wkt
from p4r_python_utils import *

OE_SUBANNUAL_FORMAT = lambda x: x.strftime("%m-%d %H:%M%z").replace("+0100", "+01:00") 



# definition of seasons - this is used to aggregate variables per seasons
spring = range(80, 172)
summer = range(172, 264)
fall = range(264, 355)
def season(x):
	if x in spring:
		return 'Spring'
	if x in summer:
		return 'Summer'
	if x in fall:
		return 'Autumn'
	else :
		return 'Winter'
path = os.environ.get("PLAN4RESROOT")

def abspath_to_relpath(path, basepath):
		return os.path.relpath(path, basepath) if os.path.abspath(path) else path

nbargs=len(sys.argv)
if nbargs>1: 
	settings_posttreat=abspath_to_relpath(sys.argv[1], path)
	if nbargs>2:
		settings_format=abspath_to_relpath(sys.argv[2], path)
		if nbargs>3:
			settings_create=abspath_to_relpath(sys.argv[3], path)
		else:
			settings_create="settingsCreateInputPlan4res.yml"
	else:
		settings_format="settings_format.yml"
else:
	settings_posttreat="settingsPostTreatPlan4res.yml"

# read config file
cfg={}
# open the configuration files 
with open(os.path.join(path, settings_posttreat),"r") as myyaml:
	cfg1=yaml.load(myyaml,Loader=yaml.FullLoader)
with open(os.path.join(path, settings_create),"r") as mysettings:
	cfg2=yaml.load(mysettings,Loader=yaml.FullLoader)
with open(os.path.join(path, settings_format),"r") as mysettings:
	cfg3=yaml.load(mysettings,Loader=yaml.FullLoader)

cfg = {**cfg1, **cfg2, **cfg3}
if 'partition' in cfg.keys():
	if 'Continent' in cfg['partition'].keys():
		cfg['partition']['continent'] = cfg['partition']['Continent']


# replace name of current dataset by name given as input
if nbargs>4:
	namedataset=sys.argv[4]
	if 'path' in cfg:
		cfg['path']=cfg['path'].replace(cfg['path'].split('/')[len(cfg['path'].split('/'))-2],namedataset)
	else:
		cfg['path']='/data/local/'+namedataset+'/'
	logger.info('posttreat '+namedataset)

TimeStepHours=cfg['Calendar']['TimeStep']['Duration']

cfg['dir']=cfg['path']+cfg['Resultsdir']
cfg['inputpath']=cfg['path']+cfg['inputDir']
if 'configDir' not in cfg: cfg['configDir']=cfg['path']+'settings/'
if cfg['USEPLAN4RESROOT']:
	path = os.environ.get("PLAN4RESROOT")
	cfg['dir']=path+cfg['dir']
	cfg['inputpath']=path+cfg['inputpath']
	cfg['path']=path+cfg['path']
cfg['dayfirst']=cfg['Calendar']['dayfirst']
cfg['BeginDataset']=cfg['Calendar']['BeginDataset']
if 'pythonDir' not in cfg: 
	if cfg['USEPLAN4RESROOT']: 
		cfg['pythonDir']='/scripts/python/plan4res-scripts/settings/'
	else:
		logger.error('pythonDir missing in settingsCreateInputPlan4res')
		log_and_exit(1, cfg['path'])

logger.info('results in:'+cfg['dir'])
logger.info('dataset in:'+cfg['inputpath'])
# define latex functions
############################
isLatex=cfg['PostTreat']['Volume']['latex']+cfg['PostTreat']['Flows']['latex']\
	+cfg['PostTreat']['Power']['latex']+cfg['PostTreat']['Demand']['latex']\
	+cfg['PostTreat']['MarginalCost']['latex']+cfg['PostTreat']['MarginalCostFlows']['latex']\
	+cfg['PostTreat']['InstalledCapacity']['latex']+cfg['PostTreat']['SpecificPeriods']['latex']

if isLatex:
	# include list of packages
	def packages(*packages):
		p= ""
		for i in packages:
			p= p+"\\usepackage{"+i+"}\n"
		return p

	# include list of ticklibraries
	def tikzlibraries(*tikzlibrary):
		p=""
		for i in tikzlibrary:
			p= p+"\\usetikzlibrary{"+i+"}\n"	
		return p

	# define list of colors
	def definecolors(*color):
		p=""
		for i in color:
			p= p+"\\definecolor"+i+"\n"	
		return p

	# include an image
	def figure(fig,caption,figlabel):
		p="\\begin{figure}[H]\n\\centering\n\\includegraphics[width=\\textwidth,keepaspectratio]{"+fig+"}\n\\caption{"+caption+"}\n\\label{fig:"+figlabel+"}\n\\end{figure}\n"
		return p

	# include a table
	def tablelatex(df,label,column_format=None,header=True,index=True):
		p="\\begin{center}"
		p=p+"\\small"
		p=p+df.to_latex(header=header,index=index,label=label,caption=label,column_format=column_format)
		p=p+"\\end{center}"
		# the following is necessary if booktabs is not in the latex distribution used
		p=p.replace("\\toprule","\\hline")
		p=p.replace("\\midrule","\\hline")
		p=p.replace("\\bottomrule","\\hline")
		p=p.replace(">","$\\Rightarrow$ ")
		p=p+"\n"
		return p
	############################################################################

# create the dictionnary of variables containing the correspondence between plan4res (SMS++) variable 
# names and openentrance nomenclature variable names
vardict={}
with open(path+cfg['pythonDir']+"VariablesDictionnary.yml","r") as myvardict:
	vardict=yaml.safe_load(myvardict)

# it may be difficult to install geopandas, this allows to skip it (not install, not use)
if cfg['geopandas']:
	import geopandas as gpd
	from shapely.geometry import Polygon, LineString, Point

# loop on scenarios, years, options => post-treat and create outputs for each (scenario,year,option) triplet
if 'variants' not in cfg: cfg['variants']=['']
if 'option' not in cfg: cfg['option']=['']
for variant,option,year in product(cfg['variants'],cfg['option'],cfg['years']):
	if len(variant)>0 or len(option)>0: 
		logger.info('treat '+variant+' '+option)
	
	# create the paths to the different directories
	# dirSto is the main directory for one (scenario,year,option)
	# dirL is the same but with name of path written such that latex understands
	if 'dir' not in cfg: cfg['dir']=cfg['outputpath']
	cfg['dirL']=cfg['dir'].replace("/","\\")

	if len(option)>0: 
		cfg['dirSto']=cfg['dir']+str(variant)+'-'+str(option)+'/'
	elif len(variant)>0: 
		cfg['dirSto']=cfg['dir']+str(variant)+'/'
	else:
		cfg['dirSto']=cfg['dir']
	# dirOUT is used to store the output files (plan4res format)
	# dir IAMC is used to store the results in Open ENTRANCE format (IAMC)
	# dirOUTScen is used to store outputs when post-treating individual results for one climatic scenario
	# dirIMG stores the images created by this script
	cfg['dirOUT']=cfg['dirSto']+'OUT/'
	cfg['dirIAMC']=cfg['dirSto']+'IAMC/'
	if 'SpecificPeriods' in cfg['PostTreat']: 
		cfg['dirOUTScenario']=[cfg['dirSto']+'Scenario'+str(i)+'/' for i in cfg['PostTreat']['SpecificPeriods']['scenarios']]
	cfg['dirIMG']=cfg['dirSto']+'IMG/'
	cfg['dirIMGLatex']=cfg['dirIMG'].replace("/","\\")
	if not os.path.isdir(cfg['dirOUT']): os.mkdir(cfg['dirOUT'])
	if not os.path.isdir(cfg['dirIAMC']): os.mkdir(cfg['dirIAMC'])
	if not os.path.isdir(cfg['dirIMG']): os.mkdir(cfg['dirIMG'])
	if cfg['PostTreat']['SpecificPeriods']: 
		for dir in cfg['dirOUTScenario']:
			if not os.path.isdir(dir): os.mkdir(dir)

	# these variables are used to check wether it is the first IAMC format (annual or subannual) output
	# if not, the IAMC output will be merged with the already existing one
	firstAnnualIAMC=True
	firstSubAnnualIAMC=True
	writeIAMC=cfg['PostTreat']['Volume']['iamc']+cfg['PostTreat']['Power']['iamc']+cfg['PostTreat']['MarginalCost']['iamc']+cfg['PostTreat']['MarginalCostFlows']['iamc']+cfg['PostTreat']['Demand']['iamc']+cfg['PostTreat']['InstalledCapacity']['iamc']
	if 'SpecificPeriods' in cfg['PostTreat']: writeIAMC=writeIAMC+cfg['PostTreat']['SpecificPeriods']['iamc']

	# list of storage assets with/without pumping
	cfg['turbine']=cfg['nopumping']+cfg['pumping']
	cfg['pump']=[]
	for elem in cfg['pumping']: cfg['pump'].append(elem+'_PUMP')

	# preamble of latex file
	if isLatex:
		startlatex = "\\documentclass[10pt]{report}\n"
		startlatex = startlatex+packages('float','amsfonts','amssymb','amsmath','makeidx','amsgen','epsf','fancyhdr','acro','etoolbox')
		startlatex=startlatex+"\\usepackage[hidelinks]{hyperref}\n"
		startlatex = startlatex+packages('subfigure','lastpage','ltablex','graphicx','color','url')
		startlatex=startlatex+"\\usepackage[table]{xcolor}\n"
		startlatex = startlatex+packages('multirow','enumitem','cite','tikz')
		startlatex=startlatex+"\\usepackage[latin1]{inputenc}\n"
		startlatex=startlatex+"\\usepackage[official]{eurosym}\n"
		startlatex=startlatex+"\\setlength{\\hoffset}{0pt}\n"
		startlatex=startlatex+"\\setlength{\\oddsidemargin}{0pt}\n\\setlength{\\evensidemargin}{9pt}\n"
		startlatex=startlatex+"\\setlength{\\marginparwidth}{54pt}\n"
		startlatex=startlatex+"\\setlength{\\textwidth}{481pt}\n"
		startlatex=startlatex+"\\setlength{\\voffset}{-18pt}\n"
		startlatex=startlatex+"\\setlength{\\marginparsep}{7pt}\n"
		startlatex=startlatex+"\\setlength{\\topmargin}{0pt}\n"
		startlatex=startlatex+"\\setlength{\\headheight}{13pt}\n"
		startlatex=startlatex+"\\setlength{\\headsep}{10pt}\n"
		startlatex=startlatex+"\\setlength{\\footskip}{27pt}\n"
		startlatex=startlatex+"\\setlength{\\textheight}{708pt}\n"
		startlatex = startlatex+"\\begin{document}\n"
		startlatex=startlatex+"\\title{"+cfg['titlereport']+"}\n\\maketitle\n"
		startlatex=startlatex+"\\newpage\n\\listoffigures\n\\newpage\n\\listoftables\n\\newpage\n\\tableofcontents\n\\newpage\n"
		endlatex = "\\end{document}"
		bodylatex_IC=""
		bodylatex_Sto=""
		bodylatex_Det=""
		writelatex=False

	# create list of regions
	##########################################################
	inputdata_save = dict()
	logger.info('read regions')
	if os.path.isfile(os.path.join(cfg['inputpath'], cfg['csvfiles']['ZP_ZonePartition'])):
		InputData=read_input_csv(cfg, 'ZP_ZonePartition')
		colregion=cfg['CouplingConstraints']['ActivePowerDemand']['Partition']
		list_regions=InputData[colregion].unique().tolist()
	else:
		logger.error('no file ZP_ZonePartition in dataset')
		log_and_exit(1, cfg['path'])
		
	# create list of technos, installed capacities and costs
	##########################################################
	logger.info('read installed capacity and costs')
	InputVarCost=pd.DataFrame(index=list_regions)
	InputStartUpCost=pd.DataFrame(index=list_regions)
	InputFixedCost=pd.DataFrame(index=list_regions)
	InstalledCapacity=pd.DataFrame(index=list_regions)
	IsInvestedTechno=pd.DataFrame(index=list_regions)
	isInvest=cfg['ParametersCreate']['invest']
	listTechnosInDataset=[]
	listLinesInDataset=[]
	listInvestedAssets=[]				  
	
	# seasonal storage mix
	if os.path.isfile(os.path.join(cfg['inputpath'], cfg['csvfiles']['SS_SeasonalStorage'])):
		InputData=read_input_csv(cfg,'SS_SeasonalStorage')
		inputdata_save['SS_SeasonalStorage'] = InputData
		listTechnosSS=InputData.Name.unique().tolist()
		for techno in listTechnosSS:
			if techno not in listTechnosInDataset:
				listTechnosInDataset.append(techno)
			df=InputData[ InputData.Name==techno ]
			df.index=df.Zone
			InstalledCapacity[techno]=df['MaxPower']*df['NumberUnits']
			InputVarCost[techno]=0
			InputStartUpCost[techno]=0
			InputStartUpCost[techno]=0
		if isInvest:
			IsInvestedTechno[techno]=False
	else: listTechnosSS=[]	
	
	# thermal mix
	if os.path.isfile(os.path.join(cfg['inputpath'], cfg['csvfiles']['TU_ThermalUnits'])):
		InputTUData=read_input_csv(cfg, 'TU_ThermalUnits')
		inputdata_save['TU_ThermalUnits'] = InputTUData
		listTechnosTU=InputTUData.Name.unique().tolist()
		if isInvest:
			for row in InputTUData.index:
				if (InputTUData.loc[row,'MaxAddedCapacity']>0)+(InputTUData.loc[row,'MaxRetCapacity']>0):
					listInvestedAssets.append( (InputTUData.loc[row,'Zone'],InputTUData.loc[row,'Name']) )	  
		for techno in listTechnosTU:
			if techno not in listTechnosInDataset:
				listTechnosInDataset.append(techno)
			df=InputTUData[ InputTUData.Name==techno ]
			df.index=df.Zone
			InputVarCost[techno]=df['VariableCost']
			InstalledCapacity[techno]=df['MaxPower']*df['NumberUnits']
			if 'VariableCost' in df.columns:
				InputStartUpCost[techno]=df['VariableCost']
			else:
				InputStartUpCost[techno]=0
			if 'FixedCost' in df.columns:
				InputStartUpCost[techno]=df['VariableCost']/8760*TimeStepHours
			else:
				InputStartUpCost[techno]=0
			if isInvest:
				IsInvestedTechno[techno]=(df['MaxAddedCapacity']>0)+(df['MaxRetCapacity']>0)
	else: listTechnosTU=[]
		
	# intermittent generation mix
	if os.path.isfile(os.path.join(cfg['inputpath'], cfg['csvfiles']['RES_RenewableUnits'])):
		InputData=read_input_csv(cfg, 'RES_RenewableUnits')
		inputdata_save['RES_RenewableUnits'] = InputData
		listTechnosRES=InputData.Name.unique().tolist()
		if isInvest:
			for row in InputData.index:
				if (InputData.loc[row,'MaxAddedCapacity']>0)+(InputData.loc[row,'MaxRetCapacity']>0):
					listInvestedAssets.append( (InputData.loc[row,'Zone'],InputData.loc[row,'Name']) )	  
		for techno in listTechnosRES:
			if techno not in listTechnosInDataset:
				listTechnosInDataset.append(techno)
			df=InputData[ InputData.Name==techno ]
			df.index=df.Zone
			InstalledCapacity[techno]=df['Capacity']
			InputVarCost[techno]=0
			InputStartUpCost[techno]=0
			InputStartUpCost[techno]=0
			if isInvest:
				IsInvestedTechno[techno]=(df['MaxAddedCapacity']>0)+(df['MaxRetCapacity']>0)
	else: listTechnosRES=[]	

	# short term storage mix
	if os.path.isfile(os.path.join(cfg['inputpath'], cfg['csvfiles']['STS_ShortTermStorage'])):
		InputData=read_input_csv(cfg, 'STS_ShortTermStorage')
		inputdata_save['STS_ShortTermStorage'] = InputData
		listTechnosSTS=InputData.Name.unique().tolist()
		if isInvest:
			for row in InputData.index:
				if (InputData.loc[row,'MaxAddedCapacity']>0)+(InputData.loc[row,'MaxRetCapacity']>0):
					listInvestedAssets.append( (InputData.loc[row,'Zone'],InputData.loc[row,'Name']) )	  
		for techno in listTechnosSTS:
			if techno not in listTechnosInDataset:
				listTechnosInDataset.append(techno)
				listTechnosInDataset.append(techno+'_PUMP')
			df=InputData[ InputData.Name==techno ]
			df.index=df.Zone
			InstalledCapacity[techno]=df['MaxPower']*df['NumberUnits']
			InputVarCost[techno]=0
			InputStartUpCost[techno]=0
			InputStartUpCost[techno]=0
			if isInvest:
				IsInvestedTechno[techno]=(df['MaxAddedCapacity']>0)+(df['MaxRetCapacity']>0)
	else: listTechnosSTS=[]

	listTechnosInDataset.append('SlackUnit')
	
	# create list of aggregated technos
	listaggrTechnosInDataset=[]
	for aggrTech in cfg['technosAggr']:
		for tech in listTechnosInDataset:
			if tech in cfg['technosAggr'][aggrTech]['technos'] and aggrTech not in listaggrTechnosInDataset:
				listaggrTechnosInDataset.append(aggrTech)

	# interconnections
	if os.path.isfile(os.path.join(cfg['inputpath'], cfg['csvfiles']['IN_Interconnections'])):
		InputData=read_input_csv(cfg, 'IN_Interconnections')
		inputdata_save['IN_Interconnections'] = InputData
		if 'Name' in InputData.columns:
			InputData=InputData.set_index('Name')
		listLinesInDataset=InputData.index.tolist()
		cfg['lines']=listLinesInDataset
		LineCapacity=InputData[ ['MaxPowerFlow', 'MinPowerFlow'] ]
		if isInvest:
			IsInvestedLine=(InputData['MaxAddedCapacity']>0)+(InputData['MaxRetCapacity']>0)
	else:
		listLinesInDataset=[]
		cfg['lines']=[]

	InstalledCapacity=InstalledCapacity.fillna(0)
	
	LineCapacity=LineCapacity.fillna(0)
	invest_factor_by_asset = dict()
	if isInvest:
		InstalledCapacity.to_csv(cfg['dirOUT']+'InitialInstalledCapacity.csv')
		LineCapacity.to_csv(cfg['dirOUT']+'InitialLineCapacity.csv')
		IsInvestedTechno=IsInvestedTechno.fillna(False)
		IsInvestedLine=IsInvestedLine.fillna(False)
		
		InvestedCapacity=pd.DataFrame(index=InstalledCapacity.index,columns=InstalledCapacity.columns,data=0.0)
		InvestedLine=pd.DataFrame(index=LineCapacity.index,columns=LineCapacity.columns,data=0.0)
		
		# get solution of capacity expansion
		sol=check_and_read_csv(cfg, cfg['dir']+'Solution_OUT.csv',header=None)
		
		indexSol=0
		for asset in listInvestedAssets:
			techno=asset[1]
			region=asset[0]
			invest_factor_by_asset[(region, techno)] = sol[0].loc[indexSol]
			print('region ',region,' techno ',techno,' indexsol ',indexSol,' sol ',sol[0].loc[indexSol],' added ',InstalledCapacity[techno].loc[region]*sol[0].loc[indexSol]-InstalledCapacity[techno].loc[region])
			InvestedCapacity[techno].loc[region]=np.round(InstalledCapacity[techno].loc[region]*sol[0].loc[indexSol]-InstalledCapacity[techno].loc[region],decimals=cfg['arrondi'])
			indexSol=indexSol+1
		InvestedCapacity=InvestedCapacity.fillna(0.0)
		for line in listLinesInDataset:
			if IsInvestedLine.loc[line]:
				invest_factor_by_asset[line] = sol[0].loc[indexSol]
				InvestedLine.loc[line]=np.round(LineCapacity.loc[line]*sol[0].loc[indexSol],decimals=cfg['arrondi'])																										
				indexSol=indexSol+1
		InvestedLine=InvestedLine.fillna(0.0)
		InvestedCapacity.to_csv(cfg['dirOUT']+'InvestedCapacity.csv')
		InvestedLine.to_csv(cfg['dirOUT']+'InvestedLines.csv')
		LineCapacity=LineCapacity+InvestedLine
		InstalledCapacity=InstalledCapacity+InvestedCapacity

	InputVarCost.to_csv(cfg['dirOUT']+'InputVariableCost.csv')
	InputStartUpCost.to_csv(cfg['dirOUT']+'InputStartUpCost.csv')
	InputFixedCost.to_csv(cfg['dirOUT']+'InputFixedCost.csv')
	InstalledCapacity.to_csv(cfg['dirOUT']+'InstalledCapacity.csv')
	LineCapacity.to_csv(cfg['dirOUT']+'LineCapacity.csv')

	# write csv_after_invest
	if isInvest:
		dir_csv_after_invest = os.path.join(cfg['dirOUT'], 'csv_after_invest')
		if not os.path.isdir(dir_csv_after_invest):
			os.makedirs(dir_csv_after_invest)
		installed_capacity = InstalledCapacity.stack().swaplevel().sort_index()
		for k in inputdata_save.keys():
			path = os.path.join(dir_csv_after_invest, cfg['csvfiles'][k])
			logger.info('Writing '+path)
			for c in ['MaxPower', 'MinPower', 'Capacity']:
				if c in inputdata_save[k].columns:
					for i in inputdata_save[k].index:
						asset = tuple(inputdata_save[k].loc[i, ['Zone', 'Name']]) if 'Zone' in inputdata_save[k].columns else inputdata_save[k].loc[i, 'Name']
						if asset in invest_factor_by_asset.keys():
							inputdata_save[k][c].loc[i] = np.round(inputdata_save[k][c].loc[i], invest_factor_by_asset[asset],decimals=cfg['arrondi'])
			inputdata_save[k].to_csv(path, index=None)

	# compute aggregated Installed capacity
	AggrInstalledCapacity=pd.DataFrame(index=list_regions)
	for aggrtechno in cfg['technosAggr']:
		AggrCols=[]
		for techno in cfg['technosAggr'][aggrtechno]['technos']:
			if techno in InstalledCapacity.columns: AggrCols=AggrCols+[ techno ]
		if len(AggrCols)>0: 
			AggrInstalledCapacity[aggrtechno]=InstalledCapacity[ AggrCols ].sum(axis=1)
	AggrInstalledCapacity.to_csv(cfg['dirOUT']+'AggrInstalledCapacity.csv')

	# treat hex color codes
	for techno in cfg['Technos']:
		cfg['Technos'][techno]['color']=mcolors.hex2color(cfg['Technos'][techno]['color'])

	# build tables of colors
	TechnoColors=pd.DataFrame(index=listTechnosInDataset,columns=['color'])
	TechnoAggrColors=pd.DataFrame(index=listaggrTechnosInDataset,columns=['color'])
	ChloroColors=pd.DataFrame(index=listaggrTechnosInDataset,columns=['color'])

	for techno in listTechnosInDataset:
		TechnoColors.loc[techno,'color']=cfg['Technos'][techno]['color']
	for techno in listaggrTechnosInDataset:
		TechnoAggrColors.loc[techno,'color']=cfg['technosAggr'][techno]['color']
		ChloroColors.loc[techno,'colormap']=cfg['technosAggr'][techno]['colors']
	

	# treat lists of regions
	#################################################
	
	# get list of scenarios and timeframe in dataset from file Demand0.csv
	if os.path.isdir(cfg['dirSto']+cfg['PostTreat']['Demand']['Dir']):
		listFiles=os.listdir(cfg['dirSto']+cfg['PostTreat']['Demand']['Dir'])
		listscen=[]
		for elem in listFiles:
			if '-' not in elem: listscen.append(int(elem.split('.')[0].replace("Demand","")))
		listscen=list(range(len(listscen)))
	else:
		logger.error('missing Demand in results')
		log_and_exit(2, cfg['path'])
		
	logger.info('scenarios in dataset: '+str(cfg['ParametersFormat']['Scenarios']))	
	logger.info('  indexed: '+', '.join([str(s) for s in listscen]))
	
	#create shortened names of regions, used in the latex report (for graphics to be readable)
	cfg['regions_short']=[]
	for region in list_regions:
		cfg['regions_short'].append(region[0:4])
		
	# list of regions to account for in visualisation
	if 'regionsANA' not in cfg:
		cfg['regionsANA']=list_regions
	
	logger.info('Regions in dataset: '+', '.join(list_regions))
	logger.info('   regions shortened names: '+str(cfg['regions_short']))
	logger.info('   regions analysed: '+str(cfg['regionsANA']))
	
	# update list of regions with reservoirs
	toremove=[]
	cfg['ReservoirRegions']=[]
	if os.path.isfile(os.path.join(cfg['dirSto'], cfg['PostTreat']['Volume']['Dir'], 'Volume'+str(listscen[0])+'.csv')):
		df=pd.read_csv(cfg['dirSto']+cfg['PostTreat']['Volume']['Dir']+'Volume'+str(listscen[0])+'.csv',index_col=0)
		for elem in df.columns:
			if 'Reservoir' in elem:
				cfg['ReservoirRegions'].append(elem.split('_')[1])
	
	for region in cfg['ReservoirRegions']:
		if region not in cfg['regionsANA']:
			toremove.append(region)
	for region in toremove:
		cfg['ReservoirRegions'].remove(region)
	logger.info('   regions with reservoirs: '+str(cfg['ReservoirRegions']))
	logger.info('Lines in dataset: '+str(cfg['lines']))
	
	# treat dates
	number_timesteps=len(df.index)
	BeginDataset=pd.Timestamp(pd.to_datetime(cfg['BeginDataset'],dayfirst=cfg['dayfirst']))
	if cfg['Calendar']['TimeStep']['Unit']=='days': TimeStepHours=TimeStepHours*24
	elif cfg['Calendar']['TimeStep']['Unit']=='weeks': TimeStepHours=TimeStepHours*168
	elif cfg['Calendar']['TimeStep']['Unit']!='hours': 
		logger.error('only hours, days, weeks possible as timestep unit')
		log_and_exit(2, cfg['path'])
	EndDataset=BeginDataset+number_timesteps*pd.Timedelta(str(TimeStepHours)+' hours')-pd.Timedelta('1 hours')
	BeginTreat=pd.Timestamp(pd.to_datetime(cfg['BeginTreatData'],dayfirst=cfg['dayfirst']))
	EndTreat=pd.Timestamp(pd.to_datetime(cfg['EndTreatData'],dayfirst=cfg['dayfirst']))
	
	if BeginTreat>EndDataset:
		logger.error('Treatment start date is after end of available data')
		log_and_exit(2, cfg['path'])
	if EndTreat<BeginDataset:
		logger.error('Treatment end date is before start of available data')
		log_and_exit(2, cfg['path'])
	if BeginTreat<BeginDataset: BeginTreat=BeginDataset
	if EndTreat>EndDataset: EndTreat=EndDataset
	if EndTreat<BeginTreat:
		logger.error('Treatment end date is before treatment start date')
		log_and_exit(2, cfg['path'])
	cfg['p4r_start']=BeginDataset
	cfg['p4r_end']=EndDataset
	cfg['plot_start']=BeginTreat
	cfg['plot_end']=EndTreat
	
	logger.info('data available between '+str(BeginDataset)+' and '+str(EndDataset))
	logger.info('PostTreating between '+str(BeginTreat)+' and '+str(EndTreat))

	# compute date ranges
	cfg['Tp4r'] = pd.date_range(start=pd.to_datetime(cfg['p4r_start'],dayfirst=cfg['dayfirst']),end=pd.to_datetime(cfg['p4r_end'],dayfirst=cfg['dayfirst']), freq=str(TimeStepHours)+'h')	
	cfg['PlotDates'] = pd.date_range(start=pd.to_datetime(cfg['plot_start'],dayfirst=cfg['dayfirst']),end=pd.to_datetime(cfg['plot_end'],dayfirst=cfg['dayfirst']), freq=str(TimeStepHours)+'h')
	NbTimeStepsToRemoveBefore=len(pd.date_range(start=pd.to_datetime(cfg['p4r_start'],dayfirst=cfg['dayfirst']),end=pd.to_datetime(cfg['plot_start'],dayfirst=cfg['dayfirst']), freq=str(TimeStepHours)+'h'))-1
	NbTimeStepsToRemoveAfter=len(pd.date_range(start=pd.to_datetime(cfg['plot_end'],dayfirst=cfg['dayfirst']),end=pd.to_datetime(cfg['p4r_end'],dayfirst=cfg['dayfirst']), freq=str(TimeStepHours)+'h'))-1
	
	# date ranges for graph outputs
	datestart=cfg['plot_start']
	dateend=cfg['plot_end']
	start=pd.to_datetime(datestart,dayfirst=cfg['dayfirst'])
	end=pd.to_datetime(dateend,dayfirst=cfg['dayfirst'])
	TimeIndex=pd.date_range(start=start,end=end, freq=str(TimeStepHours)+'h')	
	MonthIndex=pd.to_datetime(TimeIndex,format ='%Y-%m-%d %H%M%s+01:00').strftime('%Y-%m')
	NbTimeSteps=len(TimeIndex)	
	
	ScenarioIndex=[]
	for i in range(listscen[0],listscen[0]+len(listscen)): 
		ScenarioIndex=ScenarioIndex+[i]
		
	FailureCost=cfg['CouplingConstraints']['ActivePowerDemand']['Cost']
	
	world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
	world.to_csv('world.csv')
	if cfg['map']:
		####################################################################################################
		# Create map of selected regions
		####################################################################################################
		
		logger.info('create map')
		listcountrymap=list_regions
		if 'nameregions' in cfg.keys(): listcountrymap=cfg['nameregions']
		
		if 'private_map' in cfg:
			logger.info('use private map')
			df_continent = pd.read_csv(cfg['path']+cfg['private_map'],index_col=0)
			df_continent['geometry'] = df_continent['geometry'].apply(wkt.loads)
			continent=gpd.GeoDataFrame(df_continent, crs='epsg:4326')
		else:
			world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
			continent = world [ world['continent'].isin(cfg['partition']['continent'])]
		continent['Aggr']=continent['name']
		continent['country']=continent['name']
		continent['index']=continent['name']
		continent=continent.set_index('index')
		if cfg['aggregateregions']:
			for Aggr in cfg['aggregateregions']:
				for country in cfg['aggregateregions'][Aggr]:
					continent.loc[country,'Aggr']=Aggr
		else: 
			logger.info('no aggregated regions')

		# create list of countries and sample of values (to delete)
		mycountries=[]
		myvalue=[]
		i=0
		for country in listcountrymap:
			if cfg['aggregateregions']:
				if country in cfg['aggregateregions']: 
					#if not (cfg['countryaggregates'] and country in cfg['countryaggregates']):
					for smallcountry in cfg['aggregateregions'][country]:
						mycountries.append(smallcountry)
					else:
						mycountries.append(country)
				else: mycountries.append(country)
			else:
				mycountries.append(country)
		# keep only country names and geometry
		continent['name']=continent['Aggr']
		continent = continent[ ['continent','name','Aggr','country','geometry','pop_est'] ]
		continent['isin']=0
		for i in continent.index:
			if continent.loc[i,'country'] in mycountries: 
				continent.at[i,'isin']=1
		mycontinent = continent [ continent['isin'] ==1 ]
		restofcontinent = continent [ continent['isin'] == 0 ]
		del mycontinent['isin']
		del restofcontinent['isin']

		# aggregate countries in the map
		if cfg['aggregateregions']: mycontinent = mycontinent.dissolve(by='Aggr')

		# plot continent avec frontieres et representative point
		figmapcontinent, axmapcontinent = plt.subplots(1,1,figsize=(10,10))
		
		mycontinent.boundary.plot(ax=axmapcontinent,figsize=(10,10))
		axmapcontinent.set_xticklabels([])
		axmapcontinent.set_yticklabels([])
		
		# compute representative point of countries and get their coordinates
		centers=mycontinent.representative_point()
		numcenters=pd.Series(mycontinent.index,index=mycontinent['name'])
		mycontinent['x']=centers.x
		mycontinent['y']=centers.y
		mycontinent['centers']=centers
		mycontinent=mycontinent.set_index('name')
		plt.savefig(cfg['dirIMG']+'mycontinent.jpeg')
		plt.close()

	####################################################################################################
	# plot funtions
	####################################################################################################

	# functions for drawing continent maps with pies or bars
	# df: dataframe with cfg['partition']['Level1'] as index
	# NameFile: name of file => 'file.jpeg'
	# create one pie per region and draw it at the right place over the continent map
	# returns a jpg file cfg['dirIMG']+NameFile
	####################################################################################################
	def FlowsMap(df,NameFile): 
	# draws a map of continent with arrows representing the mean annual flows
	########################################################### 
		namefigpng=cfg['dirIMG']+NameFile
	
		if cfg['geopandas']:							   		 
			figflows, axflows = plt.subplots(1,1,figsize=(10,10))
			mycontinent.boundary.plot(ax=axflows,figsize=(10,10))
			centers.plot(color='r',ax=axflows,markersize=1)

			# fill data with flows and start/end
			lines=df.columns
			data={'start':[ x.split('>')[0] for x in lines ],'end':[ x.split('>')[1] for x in lines ],'flow':df.sum(axis=0).transpose()} 
			flows = pd.DataFrame(data, index = df.columns )
				
			# compute width of arrows
			if flows['flow'].abs().max()>0: flows['width']=(10.0*flows['flow']/(flows['flow'].abs().max())).abs()
			else: flows['width']=flows['flow'].abs()
			
			flows['newstart']=flows['start']
			flows['newend']=flows['end']
			flows['newstart']=np.where(flows['flow']<0,flows['end'],flows['start'])
			flows['newend']=np.where(flows['flow']<0,flows['start'],flows['end'])
			flows['newflow']=flows['flow'].abs()
			
			flows['startpoint']=flows['newstart'].apply(lambda x: centers[numcenters[x]])
			flows['endpoint']=flows['newend'].apply(lambda x: centers[numcenters[x]])
			
			
			flows['line']=flows.apply(lambda x: LineString([x['startpoint'], x['endpoint']]),axis=1)
			geoflows=gpd.GeoDataFrame(flows,geometry=flows['line'])

			# reduce length of arrows
			scaledgeometry=geoflows.scale(xfact=0.7,yfact=0.7,zfact=1.0,origin='center')
			geoflows.geometry=scaledgeometry

			# plot lines
			geoflows.plot(ax=axflows,column='newflow',linewidth=geoflows.width,cmap='Reds',vmin=-100,vmax=500)

			# plot arrows
			for line in geoflows.index:
				if geoflows['flow'][line]!=0:
					plt.arrow(list(geoflows['geometry'][line].coords)[0][0],list(geoflows['geometry'][line].coords)[0][1],
						list(geoflows['geometry'][line].coords)[1][0]-list(geoflows['geometry'][line].coords)[0][0],
						list(geoflows['geometry'][line].coords)[1][1]-list(geoflows['geometry'][line].coords)[0][1],
						head_width=1,head_length=0.5,color='black',linewidth=0,zorder=2)

			axflows.set_title("Import/Exports (MWh)",fontsize=10)
			axflows.tick_params(labelbottom=False,labelleft=False)
			plt.savefig(namefigpng)
			figflows.clf()
			plt.close('all')
		return namefigpng
		
	def EuropePieMap(df,NameFile,Colors):  
	# draw a map of continent with pies in each region
		namefigpng=cfg['dirIMG']+NameFile
		TechnosInGraph=[]
		for tech in df.columns:
			if ~df[tech].isin([0,cfg['ParametersCreate']['zerocapacity']]).all() :
				if tech not in TechnosInGraph:
					TechnosInGraph.append(tech)
		ColorsInGraph=pd.DataFrame(data=Colors.loc[ TechnosInGraph ])		
		if cfg['geopandas']:			 
			figmapcontinent, axmapcontinent = plt.subplots(1,1,figsize=(10,10))
			mycontinent.boundary.plot(ax=axmapcontinent,figsize=(10,10))
			axmapcontinent.set_xticklabels([])
			axmapcontinent.set_yticklabels([])
			xbounds=axmapcontinent.get_xbound()
			ybounds=axmapcontinent.get_ybound()
			xmin,xmax=axmapcontinent.get_xbound()
			ymin,ymax=axmapcontinent.get_ybound()
			
			for country in cfg['regionsANA']:
				todraw=df.loc[country]
				ColorsKept=Colors['color'].loc[ todraw.index ]
				x=mycontinent['x'][country]
				y=mycontinent['y'][country]
				axmapcontinent.pie(todraw,colors=ColorsKept,center=(x,y),radius=1.5)
			step=1/len(df.columns)
			X=np.arange(xmin,xmin+1,step)			
			patches=axmapcontinent.bar(X,todraw,bottom=ymax,width=0,color=ColorsInGraph['color'])
			axmapcontinent.legend(patches,ColorsInGraph.index,loc="lower left",bbox_to_anchor=(0, -0.09),facecolor='white',ncol=2)
			axmapcontinent.set_xbound(xbounds)
			axmapcontinent.set_ybound(ybounds)

			namefigpng=cfg['dirIMG']+NameFile
			plt.savefig(namefigpng)
			figmapcontinent.clf()
			plt.close('all')
		return namefigpng

	def EuropeBarMap(df,NameFile,Colors):
	# draw a map of continent with a bargraph in each region
		namefigpng=cfg['dirIMG']+NameFile
		TechnosInGraph=[]
		for tech in df.columns:
			if ~df[tech].isin([0,cfg['ParametersCreate']['zerocapacity']]).all() :
				if tech not in TechnosInGraph:
					TechnosInGraph.append(tech)
		ColorsInGraph=pd.DataFrame(data=Colors.loc[ TechnosInGraph ])
		if cfg['geopandas']:			 
			figmapcontinent, axmapcontinent = plt.subplots(1,1,figsize=(10,10))
			mycontinent.boundary.plot(ax=axmapcontinent,figsize=(10,10))
			axmapcontinent.set_xticklabels([])
			axmapcontinent.set_yticklabels([])
			xbounds=axmapcontinent.get_xbound()
			ybounds=axmapcontinent.get_ybound()	
			barsize=3
			for country in cfg['regionsANA']:
				todraw=df.loc[country]
				x=mycontinent['x'][country]
				y=mycontinent['y'][country]
				step=barsize/len(df.columns)
				X=np.arange(x,x+barsize,step)
				maxtodraw=todraw.max(axis=0)
				todraw=todraw*barsize/maxtodraw
				patches=axmapcontinent.bar(X,todraw,bottom=y,width=0.9*step,color=ColorsInGraph['color'])

			axmapcontinent.legend(patches,ColorsInGraph.index,loc="upper left",facecolor='white',ncol=2)
				
			namefigpng=cfg['dirIMG']+NameFile
			plt.savefig(namefigpng)
			figmapcontinent.clf()
			plt.close('all')
		return namefigpng

	def Chloromap(NbCols,NbRows,List,mytable,SizeCol,SizeRow,TitleSize,LabelSize,name,dpi=80):
	# draw a serie of maps of continent with each region lighter if the value of the table is lower
		namefigpng=cfg['dirIMG']+'ChloroMap-'+name+'.jpeg'
	
		if cfg['geopandas']:																		   
			fig, axes = plt.subplots(figsize=(SizeCol*NbCols,SizeRow*NbRows),nrows=NbRows, ncols=NbCols)
			x=0
			y=0
			for item in List:
				if item in mytable.columns:			   
					chloroeurope=mycontinent
					chloroeurope[item]=mytable[item]
					mintech=mytable[item].min()
					maxtech=mytable[item].max()
					chloroeurope=mycontinent.fillna(0)
					ax=chloroeurope.plot(column=item,ax=axes[y][x],cmap=ChloroColors.loc[item,'colormap'],vmin=mintech,vmax=maxtech,legend=True,\
						legend_kwds={'label':"energy (MWh)",'orientation':"vertical"})
					fig=ax.figure
					cb_ax=fig.axes[1]
					cb_ax.tick_params(labelsize=40)
					axes[y][x].set_title(item,fontsize=TitleSize)

					if x<NbCols-1: x=x+1
					else:
						x=0
						y=y+1
			namefigpng=cfg['dirIMG']+'ChloroMap-'+name+'.jpeg'
			plt.savefig(namefigpng,dpi=dpi)
			fig.clf()
			plt.close('all')	
		return namefigpng
		
	def OneChloromap(mycol,TitleSize,LabelSize,name,dpi=80,mytitle=''):
	# draw a unique map of continent with each region lighter if the value of the table is lower

		namefigpng=cfg['dirIMG']+'ChloroMap-'+name+'.jpeg'
	
		if cfg['geopandas']:	

			chloroeurope=mycontinent
			chloroeurope[mytitle]=mycol
			mintech=0
			maxtech=1
			chloroeurope=chloroeurope.fillna(0)
			figchloroeurope, axchloroeurope = plt.subplots(1,1)
			chloroeurope.plot(column=mytitle,ax=axchloroeurope,cmap='Greens',legend=True,vmin=0.2,vmax=1,legend_kwds={'label':"Decarbonized Energy Share (%)",'orientation':"horizontal"})
			namefigpng='ChloroMap-'+mytitle+'.jpeg'
			plt.savefig(cfg['dirIMG']+namefigpng)
			plt.close()

		return namefigpng
		
	def StackedBar(dfin,NameFile,Colors,drawScale=True):
	# draw a graphic with technos stacked (in the order of the df)
		df=pd.DataFrame(data=dfin)
		df = df.drop(columns=df.columns[df.le(cfg['ParametersCreate']['zerocapacity']).all()])
		
		x=list_regions
		df=df.fillna(0.0)
		Bottom=pd.Series(index=x,data=0.0)
		fig, axes = plt.subplots(figsize=(10,5),nrows=1, ncols=1)
		for techno in df.columns:
			plt.bar(x,df[techno],color=Colors.loc[techno,'color'],bottom=Bottom)
			Bottom=Bottom+df[techno]
		plt.legend(df.columns,bbox_to_anchor=(1.1, 1.05),loc='upper left')
		plt.xticks(rotation=45,fontsize=14)
		plt.tight_layout()
		
		if not drawScale:
			plt.yticks([])  
			plt.ylabel(None)  # remove the y-axis label
			plt.tick_params(left=False) 
		namefigpng=cfg['dirIMG']+NameFile
		plt.savefig(namefigpng)
		fig.clf()
		plt.close('all')
		return namefigpng		

	# function for drawing stochastic graph
	def StochasticGraphs(NbCols,NbRows,List,What,Dir,SizeCol,SizeRow,TitleSize,LabelSize,DrawMean=False,max=0,drawScale=True,isTimeIndex=True):
	# draw one graphic per item in the List, organised as a table with NbCols and NbRows, ...
	# each graphic includes all scenarios plus the mean of scenarios
		fig, axes = plt.subplots(figsize=(SizeCol*NbCols,SizeRow*NbRows),nrows=NbRows, ncols=NbCols)
		indexres=0
		x=0
		y=0
		i=0
		for item in List:
			Data=pd.read_csv(cfg['dirSto']+Dir+'/'+What+'-'+item+'.csv',nrows=NbTimeSteps,index_col=0).fillna(0.0)
			if 'time' in Data.columns:Data=Data.drop('time',axis=1)
			maxData=Data.max(axis=1).max()
			minY=Data.min(axis=1).min()
			if maxData<max: 
				maxY=maxData 
			else: maxY=max
			if isTimeIndex: Data.index=pd.to_datetime(Data.index)
			axes[y][x].plot(Data.index.to_list(),Data)
			axes[y][x].set_title(item,fontsize=TitleSize)
			
			if isTimeIndex:
				axes[y][x].tick_params(axis='x',  labelsize =LabelSize)
				period=pd.to_datetime(Data.index[-1])-pd.to_datetime(Data.index[0])
				nb_seconds=period.total_seconds()
				nb_hours=nb_seconds/3600
				nb_days=nb_hours/24
				if period> pd.Timedelta('180 days'):
					locator=mdates.MonthLocator(bymonth=range(1,13))
					axes[y][x].xaxis.set_major_locator(locator)
					axes[y][x].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
				elif period> pd.Timedelta('90 days'):
					locator=mdates.WeekdayLocator(byweekday=1,interval=2)
					axes[y][x].xaxis.set_major_locator(locator)
					axes[y][x].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
				elif period> pd.Timedelta('30 days'):
					locator=mdates.WeekdayLocator(byweekday=1)
					axes[y][x].xaxis.set_major_locator(locator)
					axes[y][x].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
				elif period> pd.Timedelta('14 days'):
					locator=mdates.DayLocator(bymonthday=range(1,32))
					axes[y][x].xaxis.set_major_locator(locator)
					axes[y][x].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
				elif period> pd.Timedelta('2 days'):
					locator=mdates.HourLocator(byhour=range(24),interval=12)
					axes[y][x].xaxis.set_major_locator(locator)
					axes[y][x].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d:%H'))
				else:
					locator=mdates.HourLocator(byhour=range(24))
					axes[y][x].xaxis.set_major_locator(locator)
					axes[y][x].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d:%H'))
			for label in axes[y][x].get_xticklabels(which='major'):
				label.set(rotation=30, horizontalalignment='right')
			
			if drawScale:
				axes[y][x].tick_params(axis='y',  labelsize =LabelSize)
			else:
				axes[y][x].set_yticks([])
			if max>0: axes[y][x].set_ylim(ymax=maxY,ymin=minY)
			mean=Data.mean(axis=1)
			del Data
			if i==0:
				meanScen=pd.DataFrame(mean,columns=[item])
				i=1
			else:
				meanScen[item]=mean
			if DrawMean: axes[y][x].plot(mean,linewidth=5,color="k")
			#axes[y][x].set_xticklabels(MonthIndex)
			if x<NbCols-1: x=x+1
			else:
				x=0
				y=y+1
			indexres=indexres+1
		namefig=cfg['dirIMG']+What+'.jpeg'
		meanScen.to_csv(cfg['dirOUT']+What+'.jpeg')
		fig.savefig(namefig,dpi=80)
		fig.clf()
		del fig
		plt.close('all')
		return namefig

	def DeterministicGraph(What,Index,nbRows,SizeCol,SizeRow,LabelSize,isTimeIndex=True):
	# draw one graphic per item in the List, organised as a table with NbCols and NbRows, ...
		fig, axes = plt.subplots(figsize=(SizeCol,SizeRow),nrows=1, ncols=1)
		Data=pd.read_csv(cfg['dirOUT']+What+'.csv',nrows=nbRows,index_col=0).fillna(0.0)
		#axes.plot(Data)
		####################
		if isTimeIndex:
			Data.index=pd.to_datetime(Data.index)
		else:
			Data.index=Index
		
		axes.plot(Data.index.to_list(),Data)
			
		if isTimeIndex:
			period=pd.to_datetime(Data.index[-1])-pd.to_datetime(Data.index[0])
			nb_seconds=period.total_seconds()
			nb_hours=nb_seconds/3600
			nb_days=nb_hours/24
			if period> pd.Timedelta('180 days'):
				locator=mdates.MonthLocator(bymonth=range(1,13))
				axes.xaxis.set_major_locator(locator)
				axes.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
			elif period> pd.Timedelta('90 days'):
				locator=mdates.WeekdayLocator(byweekday=1,interval=2)
				axes.xaxis.set_major_locator(locator)
				axes.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
			elif period> pd.Timedelta('30 days'):
				locator=mdates.WeekdayLocator(byweekday=1)
				axes.xaxis.set_major_locator(locator)
				axes.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
			elif period> pd.Timedelta('14 days'):
				locator=mdates.DayLocator(bymonthday=range(1,32))
				axes.xaxis.set_major_locator(locator)
				axes.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
			elif period> pd.Timedelta('2 days'):
				locator=mdates.HourLocator(byhour=range(24),interval=12)
				axes.xaxis.set_major_locator(locator)
				axes.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d:%H'))
			else:
				locator=mdates.HourLocator(byhour=range(24))
				axes.xaxis.set_major_locator(locator)
				axes.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d:%H'))
		for label in axes.get_xticklabels(which='major'):
			label.set(rotation=30, horizontalalignment='right')
				
		##############
		axes.tick_params(axis='x',  labelsize =LabelSize/2)
		axes.tick_params(axis='y',  labelsize =LabelSize/2)
		axes.legend(Data.columns,bbox_to_anchor=(0.95, 1.0),loc='upper left')
		namefig=cfg['dirIMG']+What+'.jpeg'
		plt.savefig(namefig)
		fig.clf()
		plt.close('all')
		return namefig

	def StackedGraph(dfin,serie,title,nameserie,Colors,Labels,namefigpng):
		df=pd.DataFrame(data=dfin)
		df = df.drop(columns=df.columns[df.le(cfg['ParametersCreate']['zerocapacity']).all()])
		
		ColorsInGraph=[]
		for col in df.columns:
			ColorsInGraph.append(TechnoColors['color'].loc[col])
		
		fig, axes = plt.subplots(figsize=(10,5),nrows=1, ncols=1)
		axes.set_title(title,fontsize=12)
		patches =df.plot.area(color=ColorsInGraph,ax=axes,legend=False,linewidth=0,fontsize=10)
		l11=serie.plot(color="black", linewidth=2.0, linestyle="-",ax=axes,legend=False)
		labels=Labels.append(nameserie)
		axes.tick_params(axis="x",labelsize=6)
		axes.tick_params(axis="y",labelsize=6)
		axes.legend(labels=labels,bbox_to_anchor=(1,1),loc="upper left",fontsize=4,ncol=2)
		plt.tight_layout()
		plt.savefig(cfg['dirIMG']+namefigpng)
		fig.clf()
		plt.close('all')
		
	def Curve(Serie,Y_min,Y_max,title,ylabel,namefigpng):
		fig, axes = plt.subplots(figsize=(10,5),nrows=1, ncols=1)				
		l1=Serie.plot(kind='line',ax=axes,legend=False)
		axes.legend()
		axes.set_ylabel(ylabel)
		axes.set_ylim([Y_min,MaxCmar+Y_max])
		axes.set_title(title,fontsize=12)			
		plt.tight_layout()
		plt.savefig(cfg['dirIMG']+namefigpng)
		fig.clf()
		plt.close('all')
	####################################################################################################
		
	####################################################################################################
	# treat Installed capacity
	####################################################################################################
		
	if(cfg['PostTreat']['InstalledCapacity']['draw']):
	# create graphs for installed capacity
		logger.info('Draw Installed capacity')
		InstalledCapacity=pd.read_csv(cfg['dirOUT']+'InstalledCapacity.csv',index_col=0).fillna(0.0)
		AggrInstalledCapacity=pd.read_csv(cfg['dirOUT']+'AggrInstalledCapacity.csv',index_col=0).fillna(0.0)
		namefigpng=StackedBar(InstalledCapacity,'InstalledCapacityBar.jpeg',TechnoColors)
		namefigpng=StackedBar(AggrInstalledCapacity,'AggrInstalledCapacityBar.jpeg',TechnoAggrColors)
		if isInvest:
			InitialInstalledCapacity=pd.read_csv(cfg['dirOUT']+'InitialInstalledCapacity.csv',index_col=0).fillna(0.0)
			InvestedInstalledCapacity=pd.read_csv(cfg['dirOUT']+'InvestedCapacity.csv',index_col=0).fillna(0.0)
			namefigpng=StackedBar(InitialInstalledCapacity,'InitialInstalledCapacityBar.jpeg',TechnoColors)
			namefigpng=StackedBar(InvestedInstalledCapacity,'InvestedInstalledCapacityBar.jpeg',TechnoColors)
			
		if cfg['map']: 
			namefigpng=EuropePieMap(InstalledCapacity,'InstalledCapacityMapPieEurope.jpeg',TechnoColors)
			if isInvest:
				namefigpng=EuropePieMap(InitialInstalledCapacity,'InitialInstalledCapacityMapPieEurope.jpeg',TechnoColors)

		del namefigpng

	if(cfg['PostTreat']['InstalledCapacity']['latex']):
	# create latex chapter about installed capacity
		logger.info('Include Installed capacity in latex report')
		writelatex=True
		namefigpng='InstalledCapacityMapPieEurope.jpeg'
		namefigpng2='InstalledCapacityMapBarEurope.jpeg'
		namefigpng3='InstalledCapacityBar.jpeg'
		namefigpng4='AggrInstalledCapacityBar.jpeg'
		if isInvest:
			namefigpng5='InitialInstalledCapacityMapPieEurope.jpeg'
			namefigpng6='InvestedInstalledCapacityMapPieEurope.jpeg'
		bodylatex_IC=bodylatex_IC+"\\chapter{Installed Capacity}\n" \
			+figure(cfg['dirIMGLatex']+namefigpng,'Map of Installed Capacity (Pies)',namefigpng) \
			+figure(cfg['dirIMGLatex']+namefigpng2,'Map of Installed Capacity (MW)',namefigpng2) \
			+figure(cfg['dirIMGLatex']+namefigpng3,'Installed Capacity (MW)',namefigpng3) \
			+figure(cfg['dirIMGLatex']+namefigpng4,'Installed Aggregated Capacity (MW)',namefigpng4)
		if isInvest:
			bodylatex_IC=bodylatex_IC \
			+figure(cfg['dirIMGLatex']+namefigpng5,'Initial Installed Aggregated Capacity (MW)',namefigpng5)\
			+figure(cfg['dirIMGLatex']+namefigpng6,'Invested Installed Aggregated Capacity (MW)',namefigpng6)

	if(cfg['PostTreat']['InstalledCapacity']['iamc']):
	# create Open ENTRANCE format files with Installed Capacity and Variable Costs
		logger.info('Create Installed capacity Iamc Files')
		firstAnnualIAMC==True
		InstalledCapacity=pd.read_csv(cfg['dirOUT']+'InstalledCapacity.csv',index_col=0)
		VariableCost=pd.read_csv(cfg['dirOUT']+'InputVariableCost.csv',index_col=0)
		for InstTech in InstalledCapacity.columns:
			#loop on technologies
			treat=True
			if InstTech in cfg['technos']['reservoir']: 
				var=vardict['Output']['Capacity']+'|'+vardict['Output']['TechnoLongNames']['reservoir']+InstTech
			elif InstTech in cfg['technos']['hydrostorage']:
				var=vardict['Output']['Capacity']+'|'+vardict['Output']['TechnoLongNames']['hydrostorage']+InstTech
			elif InstTech in cfg['technos']['battery'] or InstTech=="Import" or InstTech=="Export" or "_PUMP" in InstTech:
				treat=False
			else: var=vardict['Output']['Capacity']+'|'+InstTech
			if treat:
				data={'region':InstalledCapacity.index,str(BeginDataset.year):InstalledCapacity[InstTech]}
				df=pd.DataFrame(data)
				Annualdf=pyam.IamDataFrame(data=pd.DataFrame(data),model=cfg['model'],scenario=cfg['scenario'],unit='MW',variable=var)
				Annualdf=Annualdf.convert_unit('MW', to='GW')
				# if 'regions' in cfg['InstalledCapacity']:
					# Annualdf=Annualdf.filter(region=cfg['InstalledCapacity']['regions'])		
				if firstAnnualIAMC==True: 
					firstAnnualIAMC=False
					IAMCAnnualDf=Annualdf
				else: 
					IAMCAnnualDf=IAMCAnnualDf.append(Annualdf)
		for InstTech in VariableCost.columns:
			treat=True
			#loop on technologies
			if InstTech in cfg['technos']['reservoir']: 
				var=vardict['Output']['VariableCost']+'|'+vardict['Output']['TechnoLongNames']['reservoir']+InstTech
			elif InstTech in cfg['technos']['hydrostorage']:
				var=vardict['Output']['VariableCost']+'|'+vardict['Output']['TechnoLongNames']['hydrostorage']+InstTech
			elif InstTech in cfg['technos']['battery'] or InstTech=="SlackUnit" or InstTech=="Import" or InstTech=="Export" or "_PUMP" in InstTech:
				treat=False
			else: var=vardict['Output']['VariableCost']+'|'+InstTech
			
			if treat:
				data={'region':VariableCost.index,str(BeginDataset.year):VariableCost[InstTech]}
				df=pd.DataFrame(data)
				Annualdf=pyam.IamDataFrame(data=pd.DataFrame(data),model=cfg['model'],scenario=cfg['scenario'],unit='EUR/MWh',variable=var)
				# if 'regions' in cfg['InstalledCapacity']:
					# Annualdf=Annualdf.filter(region=cfg['InstalledCapacity']['regions'])
				if firstAnnualIAMC==True: 
					firstAnnualIAMC=False
					IAMCAnnualDf=Annualdf
				else: 
					IAMCAnnualDf=IAMCAnnualDf.append(Annualdf)
		 
		IAMCAnnualDf.to_excel(cfg['dirIAMC']+'InstalledCapacity&VariableCosts.xlsx', sheet_name='data', iamc_index=False, include_meta=True)


	####################################################################################################
	#Post treatment of stochastic results
	####################################################################################################
	nbscen=len(listscen)
	scen0=listscen[0]
	nbscengroup=nbscen # number of scenarised results
	if 'max' in cfg['PostTreat']['MarginalCost'].keys(): 
		cfg['marginalcostlimits']={'max':cfg['PostTreat']['MarginalCost']['max'],'min':(-1)*cfg['PostTreat']['MarginalCost']['max']}
	else:
		cfg['marginalcostlimits']={'max':10000,'min':-10000}
	logger.info('Treat Stochastic results')
	if isLatex: bodylatex_Sto=bodylatex_Sto+"\\chapter{Stochastic Results}\n" 
	
	# read VolumeOUT.csv per scenario and create 1 file per reservoir with all scenario data
	if (cfg['PostTreat']['Volume']['read']):
		logger.info('   Read Stochastic results for Volumes')
		reservoirs=cfg['ReservoirRegions']
		nbreservoirs=len(reservoirs)
		Volumes=[]

		if cfg['usevu']: 
			SMSVB=pd.read_csv(cfg['dirSto'] +'BellmanValuesOUT.csv')
			listres=cfg['ReservoirRegions']
			cols=['Timestep']+listres+['b']
			SMSVB.columns=cols
		BV=pd.DataFrame(index=listscen)		
		i=0
		
		usecols=[i for i in range(0,len(cfg['ReservoirRegions'])+1)]
		
		# treat scenarised results for Volume
		for scen in listscen:
			logger.info('     opening file Volumes for scenario '+str(scen))
			df=pd.read_csv(cfg['dirSto']+cfg['PostTreat']['Volume']['Dir']+'Volume'+str(scen)+'.csv',skiprows = lambda x : x > 0 and x <= NbTimeStepsToRemoveBefore,skipfooter=NbTimeStepsToRemoveAfter,engine='python',usecols=usecols,index_col=0)
			df.index=TimeIndex
			
			if scen==listscen[0]:
				reservoirs=[]
				for col in df.columns:
					if 'Reservoir' in col:
						reg=col.split('_')[1]
						reservoirs.append(reg)
			newreservoirs=[]
			for reservoir in reservoirs:
				newreservoirs.append(reservoir+'-'+str(scen))			
				
			df.columns=newreservoirs
			if i==0:
				for indexres in range(nbreservoirs):
					Volumes.append(df[newreservoirs[indexres]])
				i=1
			else:
				for indexres in range(nbreservoirs):
					Volumes[indexres]=pd.concat([Volumes[indexres],df[newreservoirs[indexres]]], axis=1)
			
			if cfg['usevu']: 
				# compute Bellman value at begin and end of time period
				VolSMS=pd.DataFrame(columns=['Name','VolumeIni','VolumeFin']) 
				for col in df.columns:
					name=col.split('-')[0]
					VolSMS.loc[name]=[ name, df.at[df.index[0],col], df.at[df.index[-1],col] ]
				value=[0]
				SMSVBIni=SMSVB[ SMSVB.Timestep.isin(value) ]
				value=[cfg['timestepvusms']]
				SMSVBFin=SMSVB[ SMSVB.Timestep.isin(value) ]
				aFin=SMSVBFin[ cfg['ReservoirRegions'] ]
				bFin=SMSVBFin['b']
				aIni=SMSVBIni[ cfg['ReservoirRegions'] ]
				bIni=SMSVBIni['b']

				valuesmsIni=-10000000000000
				valuesmsFin=-10000000000000
				for row in aIni.index:
					valIni=bIni.loc[row]+(aIni.loc[row]*VolSMS['VolumeIni']).sum()
					if valIni > valuesmsIni: 
						valuesmsIni=valIni
				for row in aFin.index:
					valFin=bFin.loc[row]+(aFin.loc[row]*VolSMS['VolumeFin']).sum()
					if valFin > valuesmsFin: 
						valuesmsFin=valFin	
				
				BV.at[scen,'InitialBellmanValue']=valuesmsIni
				BV.at[scen,'FinalBellmanValue']=valuesmsFin
				for reg in reservoirs:
					BV.at[scen,'InitialVolume-'+reg]=VolSMS.at[reg,'VolumeIni']
					BV.at[scen,'FinalVolume-'+reg]=VolSMS.at[reg,'VolumeFin']
					
		for indexres in range(nbreservoirs):
			reservoir=reservoirs[indexres]
			Volumes[indexres].to_csv(cfg['dirSto'] +cfg['PostTreat']['Volume']['Dir']+'Volume-Reservoir-'+reservoir+'.csv')
		if cfg['usevu']: BV.to_csv(cfg['dirOUT'] +'BellmanValuesPerScenario.csv')
		del Volumes, BV, df, reservoirs

	# write volumes in IAMC data format
	if (cfg['PostTreat']['Volume']['iamc']):
		logger.info('   creating pyam for Volume')
		reservoirs=cfg['ReservoirRegions']
		nbreservoirs=len(reservoirs)
		var=vardict['Output']['Volume']
		for indexres in range(nbreservoirs):
			firstIAMC=True
			reservoir=reservoirs[indexres]
			Volume=pd.read_csv(cfg['dirSto']+cfg['PostTreat']['Volume']['Dir']+'Volume-Reservoir-'+reservoir+'.csv')
			for scen in listscen:
				df=pd.DataFrame()
				df['time']=TimeIndex
				df['value']=Volume[reservoir+'-'+str(scen)]
				mydf=pyam.IamDataFrame(data=df,model=cfg['model'],scenario=cfg['scenario']+'|'+str(scen),region=reservoir,unit='MWh',variable=var) 
				Subannualdf=mydf.swap_time_for_year(subannual=OE_SUBANNUAL_FORMAT)
				Subannualdf=Subannualdf.convert_unit('MWh', to='GWh') 
				# if 'regions' in cfg['PostTreat']['Volume']:
					# Subannualdf=Subannualdf.filter(region=cfg['InstalledCapacity']['regions'])
				if firstIAMC==True: 
					firstIAMC=False
					IAMCSubAnnualDf=Subannualdf
				else: 
					IAMCSubAnnualDf=IAMCSubAnnualDf.append(Subannualdf)
			IAMCSubAnnualDf.to_excel(cfg['dirIAMC']+reservoir+'_Volume.xlsx', sheet_name='data', iamc_index=False, include_meta=True)
						
	# postreat electricity demand per scenarised timeserie
	# read DemandOUT.csv per scenario and create 1 file per region with all scenario data
	if (cfg['PostTreat']['Demand']['read']):
		logger.info('   Read Stochastic results for Demand')
		first=1
		Demand=[]
		i=0
		for scen in listscen:
			logger.info('     opening file Demand for scenario '+str(scen))
			regions=[]
			for reg in list_regions:
				regions.append(reg+'-'+str(scen))
			df=pd.read_csv(cfg['dirSto']+cfg['PostTreat']['Demand']['Dir']+'Demand'+str(scen)+'.csv',skiprows = lambda x : x > 0 and x <= NbTimeStepsToRemoveBefore,skipfooter=NbTimeStepsToRemoveAfter,engine='python',index_col=0)
			df.index=TimeIndex
			
			df.columns=regions
			if i==0:
				for index in range(len(cfg['regionsANA'])):
					Demand.append(df[regions[index]])
				i=1
			else:
				for index in range(len(cfg['regionsANA'])):
					Demand[index]=pd.concat([Demand[index],df[regions[index]]], axis=1)
		print(list_regions)
		for index in range(len(list_regions)):
			reg=list_regions[index]
			logger.info('     writing file for '+reg)
			Demand[index].to_csv(cfg['dirSto'] +cfg['PostTreat']['Demand']['Dir']+'Demand-'+reg+'.csv')
			
		del df, Demand
			
	# write demands in Open ENTRANCE data format
	if (cfg['PostTreat']['Demand']['iamc']):
		logger.info('   creating pyam for demand')
		var=vardict['Output']['Demand']
		for reg in cfg['regionsANA']:
			firstIAMC=True
			logger.info('     creating pyam for demand in region '+reg)
			treat=True
			if 'regions' in cfg['PostTreat']['Demand']:
				if reg not in cfg['PostTreat']['Demand']['regions']: 
					treat=False
			if treat:
				Demand=pd.read_csv(cfg['dirSto'] +cfg['PostTreat']['Demand']['Dir']+'Demand-'+reg+'.csv')
				for scen in listscen:
					df=pd.DataFrame()
					df['time']=TimeIndex
					df['value']=Demand[reg+'-'+str(scen)]
					mydf=pyam.IamDataFrame(data=df,model=cfg['model'],scenario=cfg['scenario']+'|'+str(scen),region=reg,unit='MWh',variable=var) 
					Subannualdf=mydf.swap_time_for_year(subannual=OE_SUBANNUAL_FORMAT)
					
					if 'regions' in cfg['PostTreat']['Demand']:
						Subannualdf=Subannualdf.filter(region=cfg['InstalledCapacity']['regions'])
					if firstIAMC==True: 
						firstIAMC=False
						IAMCSubAnnualDf=Subannualdf
					else: 
						IAMCSubAnnualDf=IAMCSubAnnualDf.append(Subannualdf)
			IAMCSubAnnualDf.to_excel(cfg['dirIAMC']+reg+'_Demand.xlsx', sheet_name='data', iamc_index=False, include_meta=True)
			del mydf, Subannualdf, IAMCSubAnnualDf,df,Demand

	# postreat marginal costs
	# read MarginalCostActivePowerDemandOUT.csv per scenario and create 1 file per region with all scenario data
	if (cfg['PostTreat']['MarginalCost']['read']):
		logger.info('   Read Stochastic results for Marginal Costs')
		first=1
		MargCosts=[]
		i=0
		for scen in listscen:
			logger.info('     opening file MarginalCost for scenario '+str(scen))
			regions=[]
			for reg in cfg['regionsANA']:
				regions.append(reg+'-'+str(scen))
			df=pd.read_csv(cfg['dirSto']+cfg['PostTreat']['MarginalCost']['Dir']+'MarginalCostActivePowerDemand'+str(scen)+'.csv',index_col=0,skiprows = lambda x : x > 0 and x <= NbTimeStepsToRemoveBefore,skipfooter=NbTimeStepsToRemoveAfter,engine='python')
			df.index=TimeIndex
			df.columns=regions
			if i==0:
				for index in range(len(cfg['regionsANA'])):
					MargCosts.append(df[regions[index]])
				i=1
			else:
				for index in range(len(cfg['regionsANA'])):
					MargCosts[index]=pd.concat([MargCosts[index],df[regions[index]]], axis=1)
			del df
		
		meanTime=pd.DataFrame(columns=cfg['regionsANA'],index=listscen)
		#meanScen=pd.DataFrame(columns=cfg['regionsANA'],index=range(NbTimeSteps))
		meanScen=pd.DataFrame(columns=cfg['regionsANA'],index=TimeIndex)
		meanScenMonotone=pd.DataFrame(columns=cfg['regionsANA'],index=range(NbTimeSteps))
		SlackCmar=pd.DataFrame(columns=listscen)
		nbHoursSlack=pd.DataFrame(columns=listscen)									 
		i=0
		for index in range(len(cfg['regionsANA'])):
			reg=cfg['regionsANA'][index]
			logger.info('     writing file for '+reg)
			MargCosts[index].to_csv(cfg['dirSto'] +cfg['PostTreat']['MarginalCost']['Dir']+'MarginalCostActivePowerDemand-'+reg+'.csv')
			#meanScenReg=MargCosts[index].mean(axis=1).reset_index()
			if len(listscen)>1:		  
				meanScenReg=MargCosts[index].mean(axis=1)			
				meanScen[reg]=MargCosts[index].mean(axis=1)							   
				meanTimeReg=MargCosts[index].mean(axis=0).transpose().reset_index().drop('index',axis=1)
				meanTime[reg]=meanTimeReg
				SortedSlack=pd.DataFrame(columns=MargCosts[index].columns, index=MargCosts[index].index)
				List=[]
				for col in MargCosts[index].columns:
					data=MargCosts[index][col].tolist()
					data.sort(reverse=True)
					List.append(data)
			else:
				meanScenReg=MargCosts[index]
				meanScen[reg]=MargCosts[index]
				meanTimeReg=MargCosts[index].mean(axis=0)
				meanTime[reg]=meanTimeReg
				SortedSlack=pd.Series(index=MargCosts[index].index)
				List=[]
				data=MargCosts[index].tolist()
				data.sort(reverse=True)
				List.append(data)
			SortedSlack=pd.DataFrame(data=List).transpose()
			if len(listscen)>1: SortedSlack.columns=MargCosts[index].columns
			SortedSlack.to_csv(cfg['dirSto'] +cfg['PostTreat']['MarginalCost']['Dir']+'HistCmar-'+reg+'.csv')
			del data, SortedSlack, List
						
			SlackCmarReg=(MargCosts[index].value_counts().sort_index()).tail(1).fillna(0.0)
			   
			SlackCmarReg.columns=listscen
			if not(FailureCost in SlackCmarReg.index): 
				logger.info('       no non-served energy for zone '+reg)
				nbHoursSlack.loc[reg]=0
			else:
				nbHoursSlack.loc[reg]=SlackCmarReg.iloc[0]
			del SlackCmarReg
			
		nbHoursSlack.to_csv(cfg['dirOUT'] +'nbHoursSlack.csv',index=True)
		meanTime.index=cfg['ParametersFormat']['Scenarios']
		meanTime.to_csv(cfg['dirOUT'] +'meanTimeCmar.csv',index=True)
		meanScen.to_csv(cfg['dirOUT'] +'meanScenCmar.csv',index=True)
		SlackCmardf=pd.DataFrame(index=list_regions,data=SlackCmar)
		SlackCmardf=SlackCmardf.fillna(0)
		SlackCmardf.to_csv(cfg['dirOUT'] +'SlackCmar.csv',index=True)
		meanScenMonotone=meanScen
		List=[]
		for col in meanScen.columns: 
			data=meanScen[col].tolist()
			data.sort(reverse=True)
			List.append(data)
		meanScenMonotone=pd.DataFrame(data=List).transpose()
		meanScenMonotone.columns=meanScen.columns
		meanScenMonotone.to_csv(cfg['dirOUT'] +'MonotoneCmar.csv',index=True)	
			
		del MargCosts, meanTime, meanScen, meanScenMonotone, SlackCmar, SlackCmardf, data, List

	# create Open ENTRANCE data format files for marginal costs
	if (cfg['PostTreat']['MarginalCost']['iamc']):
		logger.info('   creating pyam for marginal costs')
		firstIAMC=True
		for reg in cfg['regionsANA']:
			
			logger.info('     creating pyam for marginal costs in region '+reg)
			treat=True
			if 'regions' in cfg['PostTreat']['MarginalCost']:
				if reg not in cfg['PostTreat']['MarginalCost']['regions']: 
					treat=False
			if treat:
				logger.info('     treat '+reg)
				MargCosts=pd.read_csv(cfg['dirSto'] +cfg['PostTreat']['MarginalCost']['Dir']+'MarginalCostActivePowerDemand-'+reg+'.csv',index_col=0)
				# aggregate per season
				
				MargCosts.index=pd.to_datetime(MargCosts.index)
				MargCosts['year']=MargCosts.index.year
				MargCosts['year']=MargCosts['year'].apply(np.int64)
				MargCosts['DOY'] = MargCosts.index.dayofyear
				MargCosts['day'] = MargCosts.index.weekday
				MargCosts['WeekHour']= 24*MargCosts['day']+MargCosts.index.hour
				MargCosts['season']= MargCosts['DOY'].map(lambda x: season(x))
				MargCosts['subannual']='Average Week|'+MargCosts['season']+'|Hour '+MargCosts['WeekHour'].astype(str)
				MargCosts[reg]=MargCosts[ [reg+'-'+str(s) for s in listscen] ].mean(axis=1)
				MargCosts.to_csv(cfg['dirSto'] +cfg['PostTreat']['MarginalCost']['Dir']+'XMarginalCostActivePowerDemand-'+reg+'.csv')
				MargCosts=MargCosts.drop(['DOY','day','WeekHour'],axis=1)
				# create IAMC with mean margcosts per scenarios
				MargCostsScenario=MargCosts.drop(['season','subannual'],axis=1)
				MargCostsScenario['time']=MargCostsScenario.index
				MargCostsScenario=MargCostsScenario.drop([reg+'-'+str(s) for s in listscen],axis=1)
				yearsdata=MargCosts['year'].unique().tolist()
				firstY=True
			
				MargCostsScenario=MargCostsScenario.rename(columns={reg:'value'})
				varscenario=cfg['scenario']+variant
				if len(option)>0: 
					varscenario=cfg['scenario']+variant+'|'+option
				IAMC=pyam.IamDataFrame(MargCostsScenario.drop(['year'],axis=1),region=reg,model=cfg['model'],scenario=varscenario,unit='EUR_2020/MWh',variable='Marginal Cost|Final Energy|Electricity')
				IAMC.swap_time_for_year(subannual=OE_SUBANNUAL_FORMAT).to_excel(cfg['dirIAMC']+reg+'-MarginalCost.xlsx', sheet_name='data', iamc_index=False, include_meta=True)			

				# create IAMC with mean marg costs per season
				MargCostsSeason=MargCosts.drop(['subannual'],axis=1).rename(columns={'season': 'subannual'}).groupby('subannual').mean()
				
				for scen in listscen:
					data={'region':reg,str(BeginDataset.year):MargCostsSeason[reg+'-'+str(scen)]}
					varscenario=cfg['scenario']+variant+'|'+str(scen)
					if len(option)>0: varscenario=cfg['scenario']+variant+'|'+option+'|'+str(scen)
					Subannualdf=pyam.IamDataFrame(data=pd.DataFrame(data),model=cfg['model'],scenario=varscenario,unit='EUR_2020/MWh',variable='Marginal Cost|Final Energy|Electricity') 

					if firstIAMC==True: 
						firstIAMC=False
						IAMCSubAnnualDf=Subannualdf
					else: 
						IAMCSubAnnualDf=IAMCSubAnnualDf.append(Subannualdf)

				# create IAMC with mean marg costs per average week
				MargCostsWeek=MargCosts.drop(['season'],axis=1).groupby('subannual').mean()
				MargCostsWeek['year']=MargCostsWeek['year'].apply(np.int64)
				
				firstIAMCReg=True
				for scen in listscen:
					df=pd.DataFrame()
					
					df['time']=TimeIndex
					df['value']=MargCostsWeek[reg+'-'+str(scen)]
					data={'value':MargCostsWeek[reg+'-'+str(scen)],'year':MargCostsWeek['year']}
					mydf=pyam.IamDataFrame(data=pd.DataFrame(data),model=cfg['model'],scenario=cfg['scenario']+'|'+str(scen),region=reg,unit='EUR_2020/MWh',variable='Marginal Cost|Final Energy|Electricity') 
					if firstIAMCReg==True: 
						firstIAMCReg=False
						IAMCSubAnnualRegDf=mydf
					else: 
						IAMCSubAnnualRegDf=IAMCSubAnnualRegDf.append(mydf)
				IAMCSubAnnualRegDf.to_excel(cfg['dirIAMC']+reg+'_MarginalCost.xlsx', sheet_name='data', iamc_index=False, include_meta=True)
		IAMCSubAnnualDf.to_excel(cfg['dirIAMC']+'MarginalCost.xlsx', sheet_name='data', iamc_index=False, include_meta=True)			
		del Subannualdf, IAMCSubAnnualDf,MargCosts,MargCostsWeek,MargCostsSeason

	# post-treat marginal costs for flows in interconnections
	# read MarginalCostFlowsOUT.csv per scenario and create 1 file per line with all scenario data
	if (cfg['PostTreat']['MarginalCostFlows']['read']):
		logger.info('   Reading marginal costs of interconnections')
		MargCosts=pd.Series()
		i=0
		for scen in listscen:
			logger.info('     opening file MarginalCost for scenario '+str(scen))
			
			df=pd.read_csv(cfg['dirSto']+cfg['PostTreat']['MarginalCost']['Dir']+'MarginalCostFlows'+str(scen)+'.csv',index_col=0,skiprows = lambda x : x > 0 and x <= NbTimeStepsToRemoveBefore,skipfooter=NbTimeStepsToRemoveAfter,engine='python')
			df.index=TimeIndex
			lines=df.columns.tolist()
			for line in lines:
				df=df.rename(columns={line:line+'-'+str(scen)})
			if i==0:
				for line in lines:
					MargCosts.loc[line]=pd.DataFrame(df[line+'-'+str(scen)])
				i=1
			else:
				for line in lines:
					MargCosts.loc[line][line+'-'+str(scen)]=df[line+'-'+str(scen)]
		for line in lines:
			logger.info('     writing file '+line)
			line2=line.replace('>','2')
			MargCosts[line].to_csv(cfg['dirSto'] +cfg['PostTreat']['MarginalCost']['Dir']+'MarginalCostFlows-'+line2+'.csv')
		del df, MargCosts

	# create open ENTRANCE data format output files for Marginal costs of flows
	if (cfg['PostTreat']['MarginalCostFlows']['iamc']):
		logger.info('   creating pyam for marginal costs of interconnections')
		
		for reg in cfg['lines']:
			firstIAMC=True
			logger.info('     creating pyam for marginal costs of flows for line '+reg)
		 
			treat_reg=True
			reg1=reg.split('>')[0]
			reg2=reg.split('>')[1]
			if 'regions' in cfg['PostTreat']['MarginalCostFlows'] and reg1 not in cfg['PostTreat']['MarginalCostFlows']['regions'] and reg2 not in cfg['PostTreat']['MarginalCostFlows']['regions']:
				treat_reg=False
			
			if treat_reg:
				reg2=reg.replace('>','2')
				MargCosts=pd.read_csv(cfg['dirSto']+cfg['PostTreat']['MarginalCost']['Dir'] +'MarginalCostFlows-'+reg2+'.csv')
				for scen in listscen:
					df=pd.DataFrame()
					df['time']=TimeIndex
					df['value']=MargCosts[reg+'-'+str(scen)]
					mydf=pyam.IamDataFrame(data=df,model=cfg['model'],scenario=cfg['scenario']+'|'+str(scen),region=reg,unit='EUR_2020/MWh',variable='Marginal Cost|Maximum Flow|Electricity|Transmission') 
					Subannualdf=mydf.swap_time_for_year(subannual=OE_SUBANNUAL_FORMAT)
					if firstIAMC: 
						firstIAMC=False
						IAMCSubAnnualDf=Subannualdf
					else: 
						IAMCSubAnnualDf=IAMCSubAnnualDf.append(Subannualdf)
			IAMCSubAnnualDf.to_excel(cfg['dirIAMC']+reg2+'_MarginalCostFlows.xlsx', sheet_name='data', iamc_index=False, include_meta=True)
			del mydf, Subannualdf, IAMCSubAnnualDf,df,MargCosts

	# post-treat flows in interconnections
	# read flowsOUT.csv per scenario and create 1 file per region with imports and exports with all scenario data
	if (cfg['PostTreat']['Flows']['read']):
		logger.info(' reading flows')
		CreateMean=True
		MeanFlows=[]
		MeanImportExport=pd.DataFrame(columns=cfg['lines'],index=TimeIndex,data=0.0)
		for indexzone in range(len(cfg['regionsANA'])):
			MeanFlows.append(pd.DataFrame(columns=['Export','Import'],index=TimeIndex,data=0.0))
		for scen in listscen:
			logger.info('     opening file Flows for scenario '+str(scen))
			newzones=[]
			for zone in cfg['regionsANA']:
				newzones.append(zone+'-'+str(scen))
			Fulldf=pd.read_csv(cfg['dirSto']+cfg['PostTreat']['Flows']['Dir']+'Flows'+str(scen)+'.csv',skiprows = lambda x : x > 0 and x <= NbTimeStepsToRemoveBefore,skipfooter=NbTimeStepsToRemoveAfter,engine='python')
			Fulldf.index=TimeIndex
			i=0
			
			# replace colums names (from Line_i to reg1>reg2)
			#newcols=['Timestep']+cfg['lines']
			#Fulldf.columns=newcols
			Fulldf=Fulldf.drop(['Timestep'],axis=1)
			Fulldf.to_csv(cfg['dirSto']+cfg['PostTreat']['Flows']['Dir']+'ImportExport'+str(scen)+'.csv')
			MeanImportExport=MeanImportExport+Fulldf/nbscen
			leftcols=[]
			rightcols=[]
			indexregion=0
			for region in cfg['regionsANA']:
				leftcols=[elem for elem in Fulldf.columns if region in elem.split('>')[0]]
				rightcols=[elem for elem in Fulldf.columns if region in elem.split('>')[1]]
				leftdf=Fulldf[ leftcols ]
				rightdf=Fulldf[ rightcols ]
				left=leftdf.sum(axis=1)
				right=rightdf.sum(axis=1)
				total=left-right
			
				# separate directionnal flows
				exp=0.5* ( total+ total.abs())
				imp=-0.5* ( total- total.abs())

				TotalFlows=pd.DataFrame(index=TimeIndex)
				TotalFlows['Export']=exp
				TotalFlows['Import']=imp
				TotalFlows.to_csv(cfg['dirSto']+cfg['PostTreat']['Flows']['Dir']+'ImportExport-'+region+'-'+str(scen)+'.csv')
				MeanFlows[indexregion]['Export']=MeanFlows[indexregion]['Export']+TotalFlows['Export']/nbscen
				MeanFlows[indexregion]['Import']=MeanFlows[indexregion]['Import']+TotalFlows['Import']/nbscen
				indexregion=indexregion+1

		for indexregion in range(len(cfg['regionsANA'])):
			region=cfg['regionsANA'][indexregion]
			MeanFlows[indexregion].to_csv(cfg['dirOUT'] + 'MeanImportExport-'+region+'.csv')
			indexregion=indexregion+1
		MeanImportExport.to_csv(cfg['dirOUT'] + 'MeanImportExport.csv')
		del Fulldf, MeanFlows, MeanImportExport, leftdf, rightdf, TotalFlows

	# post-treat generation schedules
	# read ActivePowerOUT.csv per scenario and create 1 file per region with all scenario data, for selected technos
	if (cfg['PostTreat']['Power']['read']):	
		logger.info('   Read Stochastic results for Power')
		#listTechnosAggr=cfg['technosAggr'].keys()
		#listTechnos=cfg['Technos'].keys()
		listTechnosAggr=listaggrTechnosInDataset
		listTechnos=listTechnosInDataset
		
		if "treat" in cfg['PostTreat']['Power']: 
			listTechnosAggr=[]
			for aggTech in cfg['PostTreat']['Power']['treat']:
				if aggTech in listaggrTechnosInDataset and aggTech not in listTechnosAggr:
					listTechnosAggr.append(aggTech)
			listTechnos=[]
			for technoAggr in listTechnosAggr:
				for tech in cfg['technosAggr'][technoAggr]['technos']:
					if tech in listTechnosInDataset and tech not in listTechnos:
						listTechnos.append(tech)
			listTechnos.append('SlackUnit')
		MeanTechnoValuesAggr=[]
		TechnoValues=[]
		NonServed=[]
		MeanTechnoValues=[]
		listScenRegions=[]
		for region in cfg['regionsANA']:
			for scen in listscen:
				listScenRegions.append(region + '-' +str(scen))
		VarCost=pd.DataFrame(columns=listScenRegions,index=listTechnos,data=0.0)
		
		for AggrTechno in listTechnosAggr: 
			MeanTechnoValuesAggr.append([])

		for Techno in listTechnos: 
			TechnoValues.append([])
			MeanTechnoValues.append([])
		CreateMeanTechnoValues=True
		
		for scen in listscen:
			logger.info('     opening file Activepower for scenario '+str(scen))
			Fulldf=pd.read_csv(cfg['dirSto']+cfg['PostTreat']['Power']['Dir']+'ActivePower'+str(scen)+'.csv',skiprows = lambda x : x > 0 and x <= NbTimeStepsToRemoveBefore,skipfooter=NbTimeStepsToRemoveAfter,engine='python')
			
			# separate pumping from turbining for _PUMP technos
			for col in Fulldf.columns:
				if '_' in col and col.split('_')[0] in cfg['pumping'] and '_PUMP' not in col and '_TURB' not in col:
					Fulldf.loc[:,col.split('_')[0]+'_PUMP_'+col.split('_')[1]+'_'+col.split('_')[2]]=0.5*(Fulldf[col]-Fulldf[col].abs())
					Fulldf.loc[:,col.split('_')[0]+'_TURB_'+col.split('_')[1]+'_'+col.split('_')[2]]=0.5*(Fulldf[col]+Fulldf[col].abs())

			Fulldf['time']=TimeIndex
			Fulldf=Fulldf.set_index('time')
			IndexAggrTechno=0
			for techno in listTechnosAggr:
				i=0
				slackcols=[]
				for tech in cfg['technosAggr'][techno]['technos']:
					if tech in listTechnos:
						if '_PUMP' not in tech:
							if tech in cfg['pumping']:
								slackcols=slackcols+[elem for elem in Fulldf.columns if (tech in elem and '_TURB' in elem) ]
							else:
								slackcols=slackcols+[elem for elem in Fulldf.columns if (tech in elem ) ]
						else:
							slackcols=slackcols+[elem for elem in Fulldf.columns if (tech in elem and '_PUMP' in elem) ]
				df=Fulldf [ slackcols ]
				
				# aggreger les colonnes par pays
				for reg in cfg['regionsANA']:
					regcols=[elem for elem in df.columns if reg in elem]
					df.loc[:,reg]=df[ regcols ].sum(axis=1)
				df=df[ cfg['regionsANA'] ]
				newzones=[]

				for col in df.columns:
					newzones.append(col+'-'+str(scen))
				df.columns=newzones
					
				if CreateMeanTechnoValues:
					for indexzone in range(len(cfg['regionsANA'])):
						MeanTechnoValuesAggr[IndexAggrTechno].append(pd.DataFrame(columns=['Mean'],index=TimeIndex,data=0.0))
					i=1
				for indexzone in range(len(cfg['regionsANA'])):
					MeanTechnoValuesAggr[IndexAggrTechno][indexzone]['Mean']=MeanTechnoValuesAggr[IndexAggrTechno][indexzone]['Mean']+df[newzones[indexzone]]
				IndexAggrTechno=IndexAggrTechno+1
			IndexTechno=0
			del df
			for region in cfg['regionsANA']:
				regcols=[elem for elem in Fulldf.columns if region in elem]
				dfregion=Fulldf[ regcols ]
				colstechno=[]
				for techno in listTechnos:
					if '_PUMP' not in techno:
						if techno in cfg['pumping']:
							technocols=[elem for elem in dfregion.columns if (techno in elem and '_TURB' in elem)]
						else:
							technocols=[elem for elem in dfregion.columns if (techno in elem)]
					else:
						technocols=[elem for elem in dfregion.columns if (techno in elem and '_PUMP' in elem)]
					dfregion.loc[:,techno]=dfregion[ technocols ].sum(axis=1)
					colstechno=colstechno+ [techno ]
					#colstechno=colstechno+ technocols
				dfregion=dfregion[ colstechno ]
				dfregion.to_csv(cfg['dirSto']+cfg['PostTreat']['Power']['Dir']+'Generation-'+region+'-'+str(scen)+'.csv')
			del dfregion
			
			for techno in listTechnos:
				i=0
				slackcols=[]
				if '_PUMP' not in techno:
					if techno in cfg['pumping']:
						slackcols=[elem for elem in Fulldf.columns if (techno in elem and '_TURB' in elem)]
					else:
						slackcols=[elem for elem in Fulldf.columns if (techno in elem)]
				else:
					slackcols=[elem for elem in Fulldf.columns if (techno in elem and '_PUMP' in elem)]

				df=Fulldf [ slackcols ]
				# aggreger les colonnes par pays
				for reg in list_regions:
					regcols=[elem for elem in df.columns if reg in elem]
					df.loc[:,reg]=df[ regcols ].sum(axis=1)
				df=df[ cfg['regionsANA'] ]
				newzones=[]
				for col in df.columns:
					newzones.append(col+'-'+str(scen))
				df.columns=newzones
				#if i==0:
				if CreateMeanTechnoValues:
					for indexzone in range(len(cfg['regionsANA'])):
						MeanTechnoValues[IndexTechno].append(pd.DataFrame(columns=['Mean'],index=TimeIndex,data=0.0))
					i=1
				
				for indexzone in range(len(cfg['regionsANA'])):
					MeanTechnoValues[IndexTechno][indexzone]['Mean']=MeanTechnoValues[IndexTechno][indexzone]['Mean']+df[newzones[indexzone]]
					if techno=='SlackUnit':
						if scen==listscen[0]:
							colsnonserved=[cfg['regionsANA'][indexzone]+'-'+str(s) for s in listscen]
							NonServed.append(pd.DataFrame(columns=colsnonserved,index=TimeIndex,data=0.0))
						NonServed[indexzone][cfg['regionsANA'][indexzone]+'-'+str(scen)]=df[newzones[indexzone]]
				IndexTechno=IndexTechno+1
				
				
				# Compute Costs
				###############################
				# test if the techno is present in the data
				IsTechno=False
				if techno in InputVarCost.columns:
					IsTechno=True
					vcost=InputVarCost[techno]
				
				if IsTechno:
					# compute variable cost
					cols=df.columns
					newcols=[elem.split('-')[0] for elem in cols]
					df.columns=newcols
					varcost=(df*vcost).sum(axis=0)
					renamedict=dict(zip(newcols,cols))
					inv_renamedict=dict(zip(cols,newcols))
					varcost=varcost.rename(renamedict)
					for reg in cfg['regionsANA']: 
						scenreg=reg + '-' + str(scen)
						VarCost.loc[techno,scenreg]=varcost[scenreg]
			del df
			CreateMeanTechnoValues=False
		del Fulldf
		
		# compute mean varcost 
		VarCost.to_csv(cfg['dirOUT'] + 'VariableCostPerScenario.csv')
		for reg in cfg['regionsANA']:
			cols = [elem for elem in VarCost.columns if reg in elem]
			VarCost[reg]=VarCost[ cols ].mean(axis=1)
		VarCost=VarCost[ cfg['regionsANA'] ]
		VarCost.to_csv(cfg['dirOUT'] + 'MeanVariableCost.csv')
		
		MeanGenerationAggr=pd.DataFrame(columns=listTechnosAggr,index=TimeIndex)
		MeanGeneration=pd.DataFrame(columns=listTechnos,index=TimeIndex)
		MeanEnergyAggr=pd.DataFrame(columns=listTechnosAggr,index=list_regions)
		MeanEnergy=pd.DataFrame(columns=listTechnos,index=list_regions)
		for indexzone in range(len(cfg['regionsANA'])):
			region=cfg['regionsANA'][indexzone]
			IndexAggrTechno=0
			for techno in listTechnosAggr:
				MeanGenerationAggr[techno]=MeanTechnoValuesAggr[IndexAggrTechno][indexzone]/nbscen
				IndexAggrTechno=IndexAggrTechno+1
			MeanGenerationAggr.to_csv(cfg['dirOUT'] + 'AggrGeneration-'+region+'.csv')
			
			MeanEnergyAggr.loc[region]=MeanGenerationAggr.sum(axis=0).transpose()
			IndexAggrTechno=IndexAggrTechno+1
			
			IndexTechno=0
			for techno in listTechnos:
				if 'Import' not in techno or 'Export' not in techno or '_PUMP' not in techno:
					MeanGeneration[techno]=MeanTechnoValues[IndexTechno][indexzone]/nbscen
					IndexTechno=IndexTechno+1
			MeanGeneration.to_csv(cfg['dirOUT'] + 'Generation-'+region+'.csv')
			MeanEnergy.loc[region]=MeanGeneration.sum(axis=0).transpose()
			IndexTechno=IndexTechno+1
			NonServed[indexzone].to_csv(cfg['dirOUT'] + 'Slack-'+region+'.csv')
		MeanEnergyAggr.to_csv(cfg['dirOUT'] + 'AggrGeneration.csv',index=True)
		MeanEnergy.to_csv(cfg['dirOUT'] + 'Generation.csv',index=True)
		
		del MeanTechnoValues, MeanTechnoValuesAggr, NonServed, VarCost, InputVarCost
		del MeanEnergy, MeanEnergyAggr
		del MeanGeneration, MeanGenerationAggr
	
	# compute variable costs of technologies
	if (cfg['PostTreat']['Power']['computecost']):
		logger.info('   Compute costs')
		listTechnosAggr=listaggrTechnosInDataset
		listTechnos=listTechnosInDataset
		
		if "treat" in cfg['PostTreat']['Power']: 
			listTechnosAggr=[]
			for aggTech in cfg['PostTreat']['Power']['treat']:
				if aggTech in listaggrTechnosInDataset and aggTech not in listTechnosAggr:
					listTechnosAggr.append(aggTech)
			listTechnos=[]
			for technoAggr in listTechnosAggr:
				for tech in cfg['technosAggr'][technoAggr]['technos']:
					if tech in listTechnosInDataset and tech not in listTechnos:
						listTechnos.append(tech)
		
		TechnoValues=[]
		NonServed=[]
		MeanTechnoValues=[]
		listScenRegions=[]
		
		for region in cfg['regionsANA']:
			for scen in listscen:
				listScenRegions.append(region + '-' +str(scen))
		VarCost=pd.DataFrame(columns=listScenRegions,index=listTechnos,data=0.0)
		
		# open costs file
		InputVarCost=pd.read_csv(cfg['dirOUT']+'InputVariableCost.csv',index_col=0)
		
		for scen in listscen:
			logger.info('     compute cost for scenario '+str(scen))
			for region in cfg['regionsANA']:
				scenreg=region + '-' + str(scen)
				dfregion=pd.read_csv(cfg['dirSto']+cfg['PostTreat']['Power']['Dir']+'Generation-'+region+'-'+str(scen)+'.csv')
			
				for techno in listTechnos:
					slackcols=[]
					if '_PUMP' in techno:
						slackcols=[elem for elem in dfregion.columns if (techno in elem and '_PUMP' in elem)]
					else:
						slackcols=[elem for elem in dfregion.columns if (techno in elem and '_PUMP' not in elem)]

					df=dfregion [ slackcols ]
					# aggreger les colonnes par pays
					newzones=[]
					for col in df.columns:
						newzones.append(col+'-'+str(scen))
					df.columns=newzones
					
					# Compute Costs
					###############################
					vcost=0.0
					if techno in InputVarCost.columns:
						vcost=InputVarCost[techno][region]
						varcost=(df*vcost).sum(axis=0).sum()
						VarCost.loc[techno,scenreg]=varcost
			
		# compute mean varcost 
		VarCost.to_csv(cfg['dirOUT'] + 'VariableCostPerScenario.csv')
		for reg in cfg['regionsANA']:
			cols = [elem for elem in VarCost.columns if reg in elem]
			VarCost[reg]=VarCost[ cols ].mean(axis=1)
		VarCost=VarCost[ cfg['regionsANA'] ]
		VarCost.to_csv(cfg['dirOUT'] + 'MeanVariableCost.csv')
		
	# create demand graphs
	if(cfg['PostTreat']['Demand']['draw']):
		logger.info('   Create graphs for demand')
		drawMean=False
		if 'DrawMean' in cfg['PostTreat']['Demand']: 
			drawMean=cfg['PostTreat']['Demand']['DrawMean']
		elif 'DrawMean' in cfg['PostTreat']:
			drawMean=cfg['PostTreat']['DrawMean']
		namefigpng1=StochasticGraphs(cfg['Graphs']['Demand']['nbcols'],cfg['Graphs']['Demand']['nblines'],\
				cfg['regionsANA'],'Demand',cfg['PostTreat']['Demand']['Dir'],cfg['Graphs']['Demand']['SizeCol'],\
				cfg['Graphs']['Demand']['SizeRow'],cfg['Graphs']['Demand']['TitleSize'],cfg['Graphs']['Demand']['LabelSize'],drawMean,0,drawScale=True)
		del namefigpng1
	
	# create demand chapter in latex report
	if(cfg['PostTreat']['Demand']['latex']):
		writelatex=True
		bodylatex_Sto=bodylatex_Sto+"\\section{Demand}\n"
		namefigpng='Demand.jpeg'
		bodylatex_Sto=bodylatex_Sto+figure(cfg['dirIMGLatex']+namefigpng,'Demands for all scenarios (MWh)',namefigpng) 
	
	# create marginal costs graphs
	if(cfg['PostTreat']['MarginalCost']['draw']):
		logger.info('   Create graphs for marginal costs')
		if 'max' in cfg['PostTreat']['MarginalCost'].keys(): 
			maxCmarSto=cfg['PostTreat']['MarginalCost']['max']
		else:
			maxCmarSto=cfg['marginalcostlimits']['max']
		maxCmar=cfg['marginalcostlimits']['max']
		drawMean=False
		if 'DrawMean' in cfg['PostTreat']['MarginalCost']: 
			drawMean=cfg['PostTreat']['Demand']['MarginalCost']
		elif 'DrawMean' in cfg['PostTreat']:
			drawMean=cfg['PostTreat']['DrawMean']
		namefigpng=StochasticGraphs(cfg['Graphs']['MarginalCost']['nbcols'],cfg['Graphs']['MarginalCost']['nblines'],\
				cfg['regionsANA'],'MarginalCostActivePowerDemand',cfg['PostTreat']['MarginalCost']['Dir'], \
				cfg['Graphs']['MarginalCost']['SizeCol'],cfg['Graphs']['MarginalCost']['SizeRow'], \
				cfg['Graphs']['MarginalCost']['TitleSize'], cfg['Graphs']['MarginalCost']['LabelSize'],max=maxCmarSto,DrawMean=drawMean,drawScale=True)
		namefigpng=StochasticGraphs(cfg['Graphs']['MarginalCost']['nbcols'],cfg['Graphs']['MarginalCost']['nblines'],\
				cfg['regionsANA'],'HistCmar',cfg['PostTreat']['MarginalCost']['Dir'],\
				cfg['Graphs']['MarginalCost']['SizeCol'],cfg['Graphs']['MarginalCost']['SizeRow'], \
				cfg['Graphs']['MarginalCost']['TitleSize'], cfg['Graphs']['MarginalCost']['LabelSize'],DrawMean=False,max=maxCmarSto,isTimeIndex=False,drawScale=True)
		namefigpng=DeterministicGraph('meanScenCmar',TimeIndex,NbTimeSteps,10,5,16)
		namefigpng=DeterministicGraph('MonotoneCmar',list(range(NbTimeSteps)),NbTimeSteps,10,5,16,isTimeIndex=False)
		namefigpng=DeterministicGraph('meanTimeCmar',cfg['ParametersFormat']['Scenarios'],nbscen,10,5,16,isTimeIndex=False)
		del namefigpng
	
	# create marginal costs chapter un latex report
	if(cfg['PostTreat']['MarginalCost']['latex']):
		writelatex=True
		bodylatex_Sto=bodylatex_Sto+"\\section{Marginal Costs}\n"
		namefigpng1='MarginalCostActivePowerDemand.jpeg'
		namefigpng2='HistCmar.jpeg'
		namefigpng3='meanScenCmar.jpeg'
		namefigpng4='MonotoneCmar.jpeg'
		namefigpng5='meanTimeCmar.jpeg'

		bodylatex_Sto=bodylatex_Sto+figure(cfg['dirIMGLatex']+namefigpng1,'Marginal Costs for all scenarios (Euro/MWh)',namefigpng1) \
			+ figure(cfg['dirIMGLatex']+namefigpng2,'Histograms of Marginal Costs for all scenarios (Euro/MWh)',namefigpng2) \
			+ figure(cfg['dirIMGLatex']+namefigpng3,'Mean Marginal Costs (Euro/MWh)',namefigpng3) \
			+ figure(cfg['dirIMGLatex']+namefigpng4,'Histogram of Mean Marginal Costs (Euro/MWh)',namefigpng4) \
			+ figure(cfg['dirIMGLatex']+namefigpng5,'Mean on all Timesteps of Marginal Costs per Scenario (Euro/MWh)',namefigpng5) 	
	if(cfg['PostTreat']['Volume']['draw']):
		logger.info('   create stochastic volumes graphs')
		drawMean=False
		if 'DrawMean' in cfg['PostTreat']['Volume']: 
			drawMean=cfg['PostTreat']['Volume']['DrawMean']
		elif 'DrawMean' in cfg['PostTreat']:
			drawMean=cfg['PostTreat']['DrawMean']
		namefigpng=StochasticGraphs(cfg['Graphs']['Volume']['nbcols'],cfg['Graphs']['Volume']['nblines'],\
			cfg['ReservoirRegions'],'Volume-Reservoir',cfg['PostTreat']['Volume']['Dir'],\
				cfg['Graphs']['Volume']['SizeCol'],cfg['Graphs']['Volume']['SizeRow'], \
				cfg['Graphs']['Volume']['TitleSize'], cfg['Graphs']['Volume']['LabelSize'],DrawMean=drawMean,drawScale=True)
		del namefigpng

	# create volume chapter in latex report
	if(cfg['PostTreat']['Volume']['latex']):
		writelatex=True
		namefigpng='Volume-Reservoir.jpeg'
		bodylatex_Sto=bodylatex_Sto+"\\section{Volumes in seasonal storages}\n" \
			+figure(cfg['dirIMGLatex']+namefigpng,'Volumes of seasonal storages (Resservoirs) (MWh)',namefigpng)
	
	# create graphs for generation
	if(cfg['PostTreat']['Power']['draw']):
		NbRows=len(TimeIndex)
		nbres=len(list_regions)
		logger.info('   create graphs for generation')
		MeanEnergyAggr=pd.read_csv(cfg['dirOUT'] + 'AggrGeneration.csv',index_col=0).fillna(0.0)
		MeanEnergy=pd.read_csv(cfg['dirOUT'] + 'Generation.csv',index_col=0).fillna(0.0)

		cols=[elem for elem in MeanEnergy.columns]
		remove=[]
		for col in cols:
			if '_PUMP' in col: remove.append(col)
			if col in cfg['technos']['battery']: remove.append(col)
			if col in ['SlackUnit','Export','Import']: remove.append(col)
		for elem in remove: cols.remove(elem)
		
		colsaggr=[elem for elem in MeanEnergyAggr.columns]
		removeaggr=[]
		for col in colsaggr:
			if '_PUMP' in col: removeaggr.append(col)
			if col in ['Non Served','Battery']: removeaggr.append(col)
		for elem in removeaggr: colsaggr.remove(elem)
		
		MeanEnergyProduced=MeanEnergy[ [elem for elem in MeanEnergy.columns if elem in cols] ]
		MeanEnergyAggrProduced=MeanEnergyAggr[ [elem for elem in MeanEnergyAggr.columns if elem in colsaggr] ]
		TechnoColorsProduced=TechnoColors.loc[ [elem for elem in MeanEnergy.columns if elem in cols] ] 
		TechnoColorsAggrProduced=TechnoAggrColors.loc[ [elem for elem in MeanEnergyAggr.columns if elem in colsaggr] ] 
		
		namefigpng=StackedBar(MeanEnergyProduced,'MeanEnergyProducedBar.jpeg',TechnoColorsProduced,False)
		namefigpng=StackedBar(MeanEnergyAggrProduced,'MeanEnergyAggrProducedBar.jpeg',TechnoColorsAggrProduced,False)
		
		listTechnosAggr=cfg['technosAggr']
		listTechnos=cfg['Technos']
		if "treat" in cfg['Graphs']['Power']: 
			listTechnosAggr=cfg['Graphs']['Power']['treat']
			listTechnos=[]
			for technoAggr in listTechnosAggr: 
				listTechnos=listTechnos+cfg['technosAggr'][technoAggr]['technos']
				
		# plot chloromap of continent 
		if cfg['map']: namefigpng=Chloromap(cfg['Graphs']['Power']['ChloroGraph']['nbcols'],cfg['Graphs']['Power']['ChloroGraph']['nblines'],\
			listTechnosAggr,MeanEnergyAggr,cfg['Graphs']['Power']['SizeCol'],\
			cfg['Graphs']['Power']['SizeRow'],cfg['Graphs']['Power']['TitleSize'],cfg['Graphs']['Power']['LabelSize'],\
			'MeanEnergyAggr')
		
			
		# for techno in listTechnosAggr:
		logger.info('     create stochastic graph for non served energy')
		drawMean=False
		if 'DrawMean' in cfg['PostTreat']['Demand']: 
			drawMean=cfg['PostTreat']['Demand']['MarginalCost']
		elif 'DrawMean' in cfg['PostTreat']:
			drawMean=cfg['PostTreat']['DrawMean']
		namefigpng=StochasticGraphs(cfg['Graphs']['Demand']['nbcols'],cfg['Graphs']['Demand']['nblines'],\
			cfg['regionsANA'],'Slack','OUT',\
			cfg['Graphs']['Demand']['SizeCol'],cfg['Graphs']['Demand']['SizeRow'], \
			cfg['Graphs']['Demand']['TitleSize'], cfg['Graphs']['Demand']['LabelSize'],DrawMean=drawMean,max=0,drawScale=True)
		del namefigpng
		del MeanEnergy, MeanEnergyAggr
		
	# create output files in OpenENtrance data format for generation
	if(cfg['PostTreat']['Power']['iamc']):
		firstIAMC=True
		firstSubAnnualIAMC=True
		
		logger.info('   Create Energy Iamc Files')
		Generation=pd.read_csv(cfg['dirOUT']+'Generation.csv',index_col=0)
		for Tech in Generation.columns:
			if ('Import' not in Tech and 'Export' not in Tech and 'SlackUnit' not in Tech and Tech not in cfg['technos']['battery']):
				TechName=Tech
				BeginVar=vardict['Output']['Energy']
				if '_PUMP' in Tech:
					TechName=Tech[0:len(Tech)-5]
					BeginVar=vardict['Output']['Pump']
				if Tech in cfg['technos']['reservoir']: 
					var=BeginVar+'|'+vardict['Output']['TechnoLongNames']['reservoir']+TechName
				elif Tech in cfg['technos']['hydrostorage']:
					var=BeginVar+'|'+vardict['Output']['TechnoLongNames']['hydrostorage']+TechName
				elif Tech in cfg['technos']['battery']:
					var=vardict['Output']['EnergyFromPumping']+'|'+vardict['Output']['TechnoLongNames']['battery']+TechName
				else: var=BeginVar+'|'+TechName
				data={'region':Generation.index,str(BeginDataset.year):Generation[Tech]}
				df=pd.DataFrame(data)
				if '_PUMP' in Tech:
					data={'region':Generation.index,str(BeginDataset.year):(-1)*Generation[Tech]}
				varscenario=cfg['scenario']+variant
				if len(option)>0: varscenario=cfg['scenario']+variant+'|'+option
				Annualdf=pyam.IamDataFrame(data=pd.DataFrame(data),model=cfg['model'],scenario=varscenario,unit='MWh',variable=var)
				Annualdf=Annualdf.convert_unit('MWh', to='GWh') 
				if 'regions' in cfg['PostTreat']['Power']:
					Annualdf=Annualdf.filter(region=cfg['InstalledCapacity']['regions'])
					
				if firstIAMC==True: 
					firstIAMC=False
					IAMCAnnualDf=Annualdf
				else: 
					IAMCAnnualDf=IAMCAnnualDf.append(Annualdf)
					
		for reg in cfg['regionsANA']:
			firstscen=True
			
			Generation=pd.read_csv(cfg['dirOUT']+'Generation-'+reg+'.csv',index_col=0)
			Generation.index=pd.to_datetime(Generation.index)
			Generation['DOY'] = Generation.index.dayofyear
			Generation['day'] = Generation.index.weekday
			Generation['hour'] = Generation.index.hour
			
			######################################################
			logger.info('     Create Energy Iamc Files for region '+reg)
			# aggregate per season
			Generation['subannual']=Generation['DOY'].map(lambda x: season(x))
			Generation=Generation.groupby('subannual').sum()
			for Tech in Generation.columns:
				if ('Import' not in Tech and 'Export' not in Tech and 'SlackUnit' not in Tech and Tech not in cfg['technos']['battery']):
					TechName=Tech
					BeginVar=vardict['Output']['Energy']
					if '_PUMP' in Tech:
						TechName=Tech[0:len(Tech)-5]
						BeginVar=vardict['Output']['Pump']
					if Tech in cfg['technos']['reservoir']: 
						var=BeginVar+'|'+vardict['Output']['TechnoLongNames']['reservoir']+TechName
					elif Tech in cfg['technos']['hydrostorage']:
						var=BeginVar+'|'+vardict['Output']['TechnoLongNames']['hydrostorage']+TechName
					elif Tech in cfg['technos']['battery']:
						var=vardict['Output']['EnergyFromPumping']+'|'+vardict['Output']['TechnoLongNames']['battery']+TechName
					else: var=BeginVar+'|'+TechName
					data={'region':reg,str(BeginDataset.year):Generation[Tech]}
					df=pd.DataFrame(data)
					if '_PUMP' in Tech:
						data={'region':reg,str(BeginDataset.year):(-1)*Generation[Tech]}
					varscenario=cfg['scenario']+variant
					if len(option)>0: varscenario=cfg['scenario']+variant+'|'+option
					SubAnnualdf=pyam.IamDataFrame(data=pd.DataFrame(data),model=cfg['model'],scenario=varscenario,unit='MWh',variable=var)
					SubAnnualdf=SubAnnualdf.convert_unit('MWh', to='GWh') 
					if 'regions' in cfg['PostTreat']['Power']:
						SubAnnualdf=SubAnnualdf.filter(region=cfg['InstalledCapacity']['regions'])
						
					if firstSubAnnualIAMC==True: 
						firstSubAnnualIAMC=False
						IAMCSubAnnualDf=SubAnnualdf
					else: 
						IAMCSubAnnualDf=IAMCSubAnnualDf.append(SubAnnualdf)
		
		IAMCAnnualDf.to_excel(cfg['dirIAMC']+'Energy.xlsx', sheet_name='data', iamc_index=False, include_meta=True)
		IAMCSubAnnualDf.to_excel(cfg['dirIAMC']+'Energy_SubAnnual.xlsx', sheet_name='data', iamc_index=False, include_meta=True)

		VarCost=pd.read_csv(cfg['dirOUT']+'MeanVariableCost.csv',index_col=0)
		VarCost=VarCost.transpose()
		firstIAMC=True
		for Tech in VarCost.columns:
			if ('_PUMP' not in Tech) and ('SlackUnit' not in Tech) and ('Import' not in Tech) and ('Export' not in Tech) and (Tech not in cfg['technos']['reservoir']+cfg['technos']['hydrostorage']+cfg['technos']['battery']+cfg['technos']['res']+cfg['technos']['runofriver']):
				var=vardict['Output']['OperationCost']+'|'+Tech
				data={'region':VarCost.index,str(BeginDataset.year):VarCost[Tech]}
				df=pd.DataFrame(data)
				varscenario=cfg['scenario']+variant
				if len(option)>0: varscenario=cfg['scenario']+variant+'|'+option
				Annualdf=pyam.IamDataFrame(data=pd.DataFrame(data),model=cfg['model'],scenario=varscenario,unit='EUR_2020',variable=var)
		
				if firstIAMC==True: 
					firstIAMC=False
					IAMCAnnualDf=Annualdf
				else: 
					IAMCAnnualDf=IAMCAnnualDf.append(Annualdf)
		IAMCAnnualDf.to_excel(cfg['dirIAMC']+'Cost.xlsx', sheet_name='data', iamc_index=False, include_meta=True)

	# create generation chapter in latex report 
	if(cfg['PostTreat']['Power']['latex']):
		writelatex=True
		bodylatex_Sto=bodylatex_Sto+"\\section{Mean Energy Generated}\n"
		namefigpng1='MeanEnergyMapPieEurope.jpeg'
		namefigpng2='MeanEnergyBar.jpeg'
		namefigpng3='ChloroMap-MeanEnergyAggr.jpeg'
		bodylatex_Sto=bodylatex_Sto+figure(cfg['dirIMGLatex']+namefigpng1,'Map of Mean Power Generation (MWh)',namefigpng1)
		bodylatex_Sto=bodylatex_Sto+figure(cfg['dirIMGLatex']+namefigpng2,'Mean Power Generation (MWh)',namefigpng2)
		bodylatex_Sto=bodylatex_Sto+figure(cfg['dirIMGLatex']+namefigpng3,'Maps of Mean Generation per technology (MWh)',namefigpng3)
		
		listTechnos=cfg['technosAggr']
		if "treat" in cfg['Graphs']['Power']: listTechnos=cfg['Graphs']['Power']['treat']
		
		bodylatex_Sto=bodylatex_Sto+"\\section{Non Served Energy}\n"
		techno="SlackUnit"
		namefigpng='Slack.jpeg'
		bodylatex_Sto=bodylatex_Sto+figure(cfg['dirIMGLatex']+namefigpng,'non served energy (MWh)',namefigpng)
		column_format="l"
		for reg in list_regions: column_format=column_format+"m{1cm}"
	
		bodylatex_Sto=bodylatex_Sto+"\\section{Mean Costs}\n"
		VarCost=pd.read_csv(cfg['dirOUT'] + 'MeanVariableCost.csv',index_col=0)
		#bodylatex_Sto=bodylatex_Sto+tablelatex(VarCost,'Mean Variable Costs (Euro)')+"\n"
		if cfg['usevu']:
			BV=pd.read_csv(cfg['dirOUT'] + 'BellmanValuesPerScenario.csv',index_col=0)
			BV=BV.mean(axis=0)
			#bodylatex_Sto=bodylatex_Sto+tablelatex(BV,'Mean BellmanValue of Seasonal Storages (Euro)')+"\n"

	# create graphs for flows
	if(cfg['PostTreat']['Flows']['draw']):
		logger.info('   Create graphs for Flows')
		MeanFlows=pd.read_csv(cfg['dirOUT'] + 'MeanImportExport.csv',index_col=0).fillna(0.0)
		if cfg['map']: namefigpng=FlowsMap(MeanFlows,'MeanFlows.jpeg')
		
	# create flows chapter in latex report
	if(cfg['PostTreat']['Flows']['latex']):
		writelatex=True
		bodylatex_Sto=bodylatex_Sto+"\\section{Mean Flows}\n"
		namefigpng1='MeanFlows.jpeg'
		bodylatex_Sto=bodylatex_Sto+figure(cfg['dirIMGLatex']+namefigpng1,'Map of Mean Flows (MWh)',namefigpng1)
		
		MeanFlows=pd.read_csv(cfg['dirOUT'] + 'MeanImportExport.csv',index_col=0)
		MeanFlowsEnergy=MeanFlows.sum(axis=0)
		#bodylatex_Sto=bodylatex_Sto+tablelatex(MeanFlowsEnergy,'Mean Energy Flows (MWh)')+"\n"
	
	# create output files in openentrance data format for flows
	if(cfg['PostTreat']['Flows']['iamc']):
		logger.info('   Create Flows Iamc Files')
		Flows=pd.read_csv(cfg['dirOUT']+'MeanImportExport.csv',index_col=0)
		var=vardict['Output']['Flow']
		Flows=Flows.sum()
		Flows=Flows.transpose()
		data={'region':Flows.index,str(BeginDataset.year):Flows}
		varscenario=cfg['scenario']+variant
		if len(option)>0: varscenario=cfg['scenario']+variant+'|'+option
		Annualdf=pyam.IamDataFrame(data=pd.DataFrame(data),model=cfg['model'],scenario=varscenario,unit='MWh',variable=var)
		Annualdf=Annualdf.convert_unit('MWh', to='GWh') 
				
		IAMCAnnualDf=Annualdf
		IAMCAnnualDf.to_excel(cfg['dirIAMC']+'Flows.xlsx', sheet_name='data', iamc_index=False, include_meta=True)

		
	####################################################################################################	
	# treat results of chosen scenarios and periods for visualisation
	####################################################################################################

	# preamble
	colsTablebilans=[]
	for tec in cfg['technosAggr']:
		colsTablebilans.append(tec)
	colsTablebilans.append('CO2')

	# declare table of energy balance 
	tablebilansNRJ=pd.DataFrame(index= list_regions, columns=colsTablebilans)
	
	if(cfg['PostTreat']['SpecificPeriods']['read']):	
		####################################################################################################	
		# read plan4res résults 
		indexScen=0
		for NumScen in cfg['PostTreat']['SpecificPeriods']['scenarios']:
			logger.info('Treat scenario '+str(NumScen))
			useNumScen=True
			cfg['dirOUTScen']=cfg['dirOUTScenario'][indexScen]
			if(cfg['PostTreat']['SpecificPeriods']['latex']): bodylatex_Det=bodylatex_Det+"\\chapter{Detailed Results on Scenario "+str(NumScen)+"}\n"
		
			# Read Generation schedule
			##########################
			logger.info('   reading ActivePowerOUT.csv')
			SMSPower=pd.read_csv(cfg['dirSto']+cfg['PostTreat']['Power']['Dir'] +'ActivePower'+str(NumScen)+'.csv')
			SMSPower.drop(['Timestep'], axis='columns', inplace=True)

			SMSPower['Time']=cfg['Tp4r']
			SMSPower=SMSPower.set_index('Time')
			SMSPower=np.round(SMSPower,decimals=cfg['arrondi'])
			SMSPower=SMSPower[ (SMSPower.index>=datestart) & (SMSPower.index<=dateend) ]

			# reading Demand
			logger.info('   reading DemandOUT.csv')
			SMSDemand=pd.read_csv(cfg['dirSto']+cfg['PostTreat']['Demand']['Dir'] +'Demand' +str(NumScen)+'.csv')
			SMSDemand.drop(['Timestep'], axis='columns', inplace=True)
			SMSDemand['Time']=cfg['Tp4r']
			SMSDemand=SMSDemand.set_index('Time')
			SMSDemand=np.round(SMSDemand,decimals=cfg['arrondi'])
			SMSDemand=SMSDemand[ (SMSDemand.index>=datestart) & (SMSDemand.index<=dateend) ]

			# reading marginal costs
			logger.info('   reading MarginalCostActivePowerDemandOUT.csv')
			SMSCMarAPD=pd.read_csv(cfg['dirSto'] +cfg['PostTreat']['MarginalCost']['Dir']+'MarginalCostActivePowerDemand'+str(NumScen)+'.csv')
			SMSCMarAPD.drop(['Timestep'], axis='columns', inplace=True)
			SMSCMarAPD['Time']=cfg['Tp4r']
			SMSCMarAPD=SMSCMarAPD.set_index('Time')
			SMSCMarAPD=np.round(SMSCMarAPD,decimals=cfg['arrondi'])
			SMSCMarAPD=SMSCMarAPD[ (SMSCMarAPD.index>=datestart) & (SMSCMarAPD.index<=dateend) ]

			logger.info('   reading MarginalCostFlowsOUT.csv')
			SMSCMarFlo=pd.read_csv(cfg['dirSto'] +cfg['PostTreat']['MarginalCost']['Dir']+'MarginalCostFlows'+str(NumScen)+'.csv')
			SMSCMarFlo.drop(['Timestep'], axis='columns', inplace=True)
			SMSCMarFlo['Time']=cfg['Tp4r']
			SMSCMarFlo=SMSCMarFlo.set_index('Time')
			SMSCMarFlo=np.round(SMSCMarFlo,decimals=cfg['arrondi'])
			SMSCMarFlo=SMSCMarFlo[ (SMSCMarFlo.index>=datestart) & (SMSCMarFlo.index<=dateend) ]

			# rename columns for marg costs and demand
			SMSCMarAPD.columns=list_regions
			SMSDemand.columns=list_regions
			SMSCMarFlo.columns=cfg['lines']

			# reading Volumes
			logger.info('   reading VolumeOUT.csv')
			SMSVolume=pd.read_csv(cfg['dirSto'] +cfg['PostTreat']['Volume']['Dir']+'Volume'+str(NumScen)+'.csv')
			SMSVolume.drop(['Timestep'], axis='columns', inplace=True)
			SMSVolume['Time']=cfg['Tp4r']
			SMSVolume=SMSVolume.set_index('Time')
			SMSVolume=np.round(SMSVolume,decimals=cfg['arrondi'])
			SMSVolume=SMSVolume[ (SMSVolume.index>=datestart) & (SMSVolume.index<=dateend) ]

			# reading Flows
			logger.info('   reading FlowsOUT.csv')
			SMSFlows=pd.read_csv(cfg['dirSto']+cfg['PostTreat']['Flows']['Dir'] +'Flows'+str(NumScen)+'.csv')
			SMSFlows.drop(['Timestep'], axis='columns', inplace=True)
			SMSFlows['Time']=cfg['Tp4r']
			SMSFlows=SMSFlows.set_index('Time')
			SMSFlows=np.round(SMSFlows,decimals=cfg['arrondi'])
			SMSFlows=SMSFlows[ (SMSFlows.index>=datestart) & (SMSFlows.index<=dateend) ]
			SMSFlowsPerCountry=pd.DataFrame(index=SMSFlows.index)
			CountryFlowsCols={}

			# separate directionnal flows
			i=0
			for line in cfg['lines']:
				FlowsCol=SMSFlows[line]
				inverseline=line.split('>')[1]+'>'+line.split('>')[0]
				SMSFlows.loc[:,line]=0.5* ( FlowsCol+ FlowsCol.abs())
				SMSFlows[inverseline]=-0.5* ( FlowsCol- FlowsCol.abs())
				i=i+1
					
			# initialilisation table for volume at last time step and bellman value
			VolSMS=pd.DataFrame(columns=['Name','VolumeIni','VolumeFin']) 


			####################################################################################################	
			# Loop on countries
			####################################################################################################	
			firstcountry=1
			CoutTotalAllCountriesSMShorsVU=0.0
			if cfg['PostTreat']['SpecificPeriods']['read']:
				for country in cfg['regionsANA']:
					logger.info('   Treat '+country)
					if cfg['PostTreat']['SpecificPeriods']['latex']: bodylatex_Det=bodylatex_Det+"\\newpage\\section{"+country+"}\n"
				
					################################################################################################
					# treat active power and compute active power stackings
					################################################################################################
					logger.info('     treating Active power')
					SMSPowercols=SMSPower.columns
					# sélect country
					SMSPowerCountryCols=[elem for elem in SMSPowercols if country in elem]
					SMSPowerCountry=SMSPower[SMSPowerCountryCols]
					
					# remove spillage columns
					spillagecols=[]
					for col in SMSPowerCountryCols:
						for technopump in cfg['nopumping']:
							# columns _1 are spillage
							if (technopump in col) & (col[-1:]=="1"):
								del SMSPowerCountry[col]
						for technopump in cfg['pumping']:
							# columns _2 are spillage
							if (technopump in col) & (col[-1:]=="2"):
								del SMSPowerCountry[col]
								
					# case of Pumped Storage
					for col in SMSPowerCountryCols:
						for technopump in cfg['pumping']:
							if technopump in col:
								# separate turbining from pumping
								colturb=col+'_0'
								colpump=col+'_1'
								SMSPowerCountry[colturb]= SMSPowerCountry[col].where(SMSPowerCountry[col] > 0, 0)
								SMSPowerCountry[colpump]= SMSPowerCountry[col].where(SMSPowerCountry[col] < 0, 0)
								
								del SMSPowerCountry[col]
						
					# rename columns
					SMSPowerCountryColsWOCountryASupprimer=[]
					SMSPowerCountrycols=[]
					for col in SMSPowerCountry.columns:
						newcol=col
						newcolwithoutcountrynorindex=col
					
						# is it a réservoir?
						techno=col.split('_')[0]
						if techno in cfg['pumping']:
							# it is a reservoir with pumping
							if col[-1:] == "0":
								# turbine
								newcol=col.split('_')[0]
								newcolwithoutcountrynorindex=col.split('_')[0]
							else:
								# pump
								newcol=col.split('_')[0]+'_PUMP'
								newcolwithoutcountrynorindex=col.split('_')[0]+'_PUMP'
						elif techno in cfg['nopumping']:
							# it is a resevroir without pumping: one column only
							newcol=col.split('_')[0]
							newcolwithoutcountrynorindex=col.split('_')[0]
						elif col.count('_')==2:
							# other kind of unit
							newcol=col.split('_')[0]+'_'+col.split('_')[2]
							newcolwithoutcountrynorindex=col.split('_')[0]
						else:
							# slack has no name
							newcol=col.split('_')[0]+'_0'
							newcolwithoutcountrynorindex=col.split('_')[0]
					
						SMSPowerCountryColsWOCountryASupprimer.append(newcol)
						SMSPowerCountrycols.append(newcolwithoutcountrynorindex)
					SMSPowerCountry.columns=SMSPowerCountryColsWOCountryASupprimer

					################################################################################################
					# treat flows
					################################################################################################
					logger.info('     treating Flows')
					SMSFlowscols=SMSFlows.columns
					
					Importcols=[elem for elem in SMSFlowscols if country in elem.split('>')[1]]
					Exportcols=[elem for elem in SMSFlowscols if country in elem.split('>')[0]]

					SMSFlowsCountry=SMSFlows[ Importcols+Exportcols ]
					PrintFlowsCountrycols=SMSFlowsCountry.columns

					SMSFlowsPerCountry=pd.concat([SMSFlowsPerCountry,SMSFlowsCountry[ PrintFlowsCountrycols ]],axis=1)
					CountryFlowsCols[country]=PrintFlowsCountrycols
							
					# residual demand
					SMSFlowsCountry['TotalExportsCountriesIn']=SMSFlowsCountry[ Exportcols ].sum(axis=1)
					SMSFlowsCountry['TotalImportsCountriesIn']=SMSFlowsCountry[ Importcols ].sum(axis=1)
					SMSPowerCountry['Export']=-1*SMSFlowsCountry['TotalExportsCountriesIn']
					SMSPowerCountry['Import']=SMSFlowsCountry['TotalImportsCountriesIn']
					SMSPowerCountrycols.append('Export')
					SMSPowerCountrycols.append('Import')
								
					#################################################################################################
					# sort columns per techno costs (from settings)
					# and only keep technos which are in the considered country 
					#################################################################################################
					CountryTechnos=[]
					CountryColors=[]
					for techno in cfg['Technos']:
						if str(techno) in SMSPowerCountrycols:
							# techno is in the country
							CountryTechnos.append(str(techno))
							CountryColors.append(cfg['Technos'][techno]['color'])
					SMSAggrPowerCountry = pd.DataFrame({}, index = SMSPowerCountry.index)
					
					# aggregation per techno/country
					for techno in CountryTechnos:
						if techno in listTechnosInDataset:
							if techno.count('PUMP'):
								techno=techno.split('_')[0]+'_'+techno.split('_')[1]
							else:
								techno=techno.split('_')[0]
							listcoltechno=[ ]
							for unit in SMSPowerCountry.columns:
								if 'PUMP' in techno:
									if techno in unit:
										listcoltechno.append(unit)
								else:
									if techno in unit and 'PUMP' not in unit:
										listcoltechno.append(unit)
							newdf=SMSPowerCountry[ listcoltechno ]
							SMSAggrPowerCountry[techno]= newdf.sum(axis=1)
							
							# case of pumping
							if techno in cfg['turbine']:
								SMSAggrPowerCountry[techno]=SMSAggrPowerCountry[techno].abs()
							if techno in cfg['pump']:
								SMSAggrPowerCountry[techno]=-1*SMSAggrPowerCountry[techno].abs()
						

					if 'periods' in cfg['PostTreat']['SpecificPeriods']:
					################################################################################################
					# Analysis n specific weeks: stacked and marg costs graphs
					################################################################################################
						for week in cfg['PostTreat']['SpecificPeriods']['periods']:
							start_week=pd.to_datetime(cfg['PostTreat']['SpecificPeriods']['periods'][week]['begin'],dayfirst=cfg['Calendar']['dayfirst'])
							end_week = pd.to_datetime(cfg['PostTreat']['SpecificPeriods']['periods'][week]['end'],dayfirst=cfg['Calendar']['dayfirst'])
							start_short=str(start_week)[0:10]
							end_short= str(end_week)[0:10]
							logger.info('     treat specific period '+str(week)+' from '+str(start_short)+' to '+str(end_short))
							WeekTimeIndex=pd.date_range(start=start_week,end=end_week, freq='1H')	
							NbTimeStepsWeek=len(WeekTimeIndex)
								
							if cfg['PostTreat']['SpecificPeriods']['latex']: bodylatex_Det=bodylatex_Det+"\\subsection{Results on week "+start_short+" to "+end_short+"}\n"
							
							# compute and draw stack for activePower
							namefigpng='Scenario_'+str(NumScen)+'_StackedActivePower_'+country+start_week.strftime('%Y-%m-%d')+'.jpeg'

							SMSAggrPowerCountryWeek=SMSAggrPowerCountry[ (SMSAggrPowerCountry.index>=start_week) & (SMSAggrPowerCountry.index<=end_week) ]
							SMSDemandWeek=SMSDemand[ (SMSDemand.index>=start_week) & (SMSDemand.index<=end_week) ]
							SMSCmarAPDCountry=SMSCMarAPD[country]
							SMSCMarAPDCountryWeek=SMSCmarAPDCountry[ (SMSCmarAPDCountry.index>=start_week) & (SMSCmarAPDCountry.index<=end_week) ]
															
							MaxCmar=SMSCMarAPDCountryWeek[ (SMSCMarAPDCountryWeek<cfg['marginalcostlimits']['max']) & (SMSCMarAPDCountryWeek>cfg['marginalcostlimits']['min']) ].max()
							MinCmar=SMSCMarAPDCountryWeek[ (SMSCMarAPDCountryWeek<cfg['marginalcostlimits']['max']) & (SMSCMarAPDCountryWeek>cfg['marginalcostlimits']['min']) ].min()
							if not np.isfinite(MaxCmar) : MaxCmar=cfg['marginalcostlimits']['max']
							if not np.isfinite(MinCmar) : MinCmar=cfg['marginalcostlimits']['min']
							
							# stacked graphs
							P4Rtitle=country + ' ' + start_short + ' to ' + end_short
							StackedGraph(SMSAggrPowerCountryWeek,SMSDemandWeek[country],P4Rtitle,'Demand',CountryColors,CountryTechnos, namefigpng)
							
							if cfg['PostTreat']['SpecificPeriods']['latex']:
								bodylatex_Det=bodylatex_Det+"\\subsubsection{Electricity Generation}\n"
								bodylatex_Det=bodylatex_Det+figure(cfg['dirIMGLatex']+namefigpng,'Electricity Generation in '+country+' (MWh) from '+str(start_short)+' to '+str(end_short),namefigpng)

							# compute and draw marginal costs of activePower
							namefigpng='Scenario_'+str(NumScen)+'_CmarAPD_'+country+start_week.strftime('%Y-%m-%d')+'.jpeg'
						
							# Marginal costs
							title='Marginal Costs ' + country + ' from ' + start_short + ' to ' + end_short	
							Curve(SMSCMarAPDCountryWeek,MinCmar-0.1*MinCmar,MaxCmar+0.1*MaxCmar,title,'Marginal Cost Euro/MWh',namefigpng)

							if cfg['PostTreat']['SpecificPeriods']['latex']:
								bodylatex_Det=bodylatex_Det+"\\subsubsection{Marginal Cost}\n"+figure(cfg['dirIMGLatex']+namefigpng,'Marginal Costs in'+country+' (Euro/MWh) from '+str(start_short)+' to '+str(end_short),namefigpng)
								
							# compute and draw stacks of primary reserve
							namefigpng='Scenario_'+str(NumScen)+'_StackedPrimary_'+country+start_week.strftime('%Y-%m-%d')+'.jpeg'


					################################################################################################
					# Marginal costs for demand 
					################################################################################################
					
					# reinit start week and end_week to start and end of full period 
					start_week=pd.to_datetime(cfg['p4r_start'])
					end_week = pd.to_datetime(cfg['p4r_end'])
					start_short=str(start_week)[0:10]
					end_short= str(end_week)[0:10]
					logger.info('     treat marginal costs whole period')
						
					# treat and draw Marginal Costs Active Power whole period
					namefigpng='Scenario_'+str(NumScen)+'_MarginalCostDemand_'+country+start_week.strftime('%Y-%m-%d')+'.jpeg'
					logger.info('        for Active Power')
					# graphs for demand marginal costs
					MaxCmar=SMSCMarAPD[ (SMSCMarAPD<cfg['marginalcostlimits']['max']) & (SMSCMarAPD>cfg['marginalcostlimits']['min']) ].max().iloc[0]	
					MinCmar=SMSCMarAPD[ (SMSCMarAPD<cfg['marginalcostlimits']['max']) & (SMSCMarAPD>cfg['marginalcostlimits']['min']) ].min().iloc[0]

					title='Demand Marginal Cost ' + country  
					Curve(SMSCMarAPD[country],MinCmar,MaxCmar,title,'Marginal Cost Euro/MWh',namefigpng)
					
					if cfg['PostTreat']['SpecificPeriods']['latex']: 
						bodylatex_Det=bodylatex_Det+"\\subsection{Global Results (from"+ start_short +" to "+ end_short+")}\n"
						bodylatex_Det=bodylatex_Det+"\\subsubsection{Demand Marginal Cost}\n"+figure(cfg['dirIMGLatex']+namefigpng,'Marginal Cost Demand (Euro/MWh)',namefigpng)
					
					################################################################################################
					# treatment of volumes
					################################################################################################
					logger.info('     treating Volumes')
					SMSVolcols=SMSVolume.columns

					# select country
					SMSVolCountryCols=[elem for elem in SMSVolcols if country in elem]
					if len(SMSVolCountryCols)>0:
						SMSVolCountry=SMSVolume[SMSVolCountryCols]
						
						# rename columns
						SMSVolCountryNewCols=[]
						deletecols=[]
						for col in SMSVolCountry.columns:
							# is it a fictive reservoir?
							if col[-4:].count('_')>1 and int(col[-4:].split('_')[2])==1: 
								logger.info(        col+' is a fictive resevoir')
								deletecols.append(col)
							else:
								# replace techno names
								newcol=col.split('_')[0]+'_Volume'
								SMSVolCountryNewCols.append(newcol)
						for col in deletecols: del SMSVolCountry[col]
						SMSVolCountry.columns=SMSVolCountryNewCols
						for col in SMSVolCountry.columns:
							name=country+'_'+col
							VolSMS.loc[name]=[ name, SMSVolCountry.at[SMSVolCountry.index[0],col], SMSVolCountry.at[SMSVolCountry.index[-1],col] ]
						# affect colors
						CountryColors=[]
						for techno in cfg['Technos']:
							if str(techno)+'_Volume' in SMSVolCountryNewCols:
								# techno is in current country
								CountryColors.append(cfg['Technos'][techno]['color'])
						
						# volumes graphs
						# seasonal storage, pumped storage , battery and demand response
						if len(SMSVolCountry.columns)>0:
							for hydrotech in cfg['graphVolumes']:
								logger.info('        creating volume graph for '+hydrotech)
							
								P4Rtitle='Volumes '+ hydrotech + ' '+country  
								nbtechs=0
								for tech in cfg['graphVolumes'][hydrotech]['Technos']:
									for elem in SMSVolCountry.columns: 
										if tech in elem: 
											nbtechs=nbtechs+1
								
								if nbtechs>0 and nbtechs<4 : 
									nbcols=1
									nbrows=nbtechs
									fig, axes = plt.subplots(figsize=(10,5),nrows=nbrows, ncols=nbcols)
								elif nbtechs>0:
									nbcols=2
									nbrows=math.ceil(nbtechs/nbcols)
									fig, axes = plt.subplots(figsize=(10,3*nbrows),nrows=nbrows, ncols=nbcols)
								x=0
								y=0
								for tech in cfg['graphVolumes'][hydrotech]['Technos']:
									Cols=[elem for elem in SMSVolCountry if tech in elem]
									SMSVolCountryHydroTech=SMSVolCountry[ Cols ]
									if len(Cols)>0:
										MaxVol=SMSVolCountryHydroTech.max().iloc[0]
										MinVol=SMSVolCountryHydroTech.min().iloc[0]
										if nbcols==1 and nbrows==1:
											axes.plot(SMSVolCountryHydroTech)
											axes.set_title('Storage '+tech,fontsize=20)
											axes.set_xticklabels([])
											axes.set_yticklabels([])
											axes.set_ylim([MinVol-0.1*MinVol,MaxVol+0.1*MaxVol])
										elif nbcols==1 or nbrows==1:
											axes[x].plot(SMSVolCountryHydroTech)
											axes[x].set_title('Storage '+tech,fontsize=20)
											axes[x].set_xticklabels([])
											axes[x].set_yticklabels([])
											axes[x].set_ylim([MinVol-0.1*MinVol,MaxVol+0.1*MaxVol])
											x=x+1
										else:
											axes[y][x].plot(SMSVolCountryHydroTech)
											axes[y][x].set_title('Storage '+tech,fontsize=20)
											axes[y][x].set_xticklabels([])
											axes[y][x].set_yticklabels([])
											axes[y][x].set_ylim([MinVol-0.1*MinVol,MaxVol+0.1*MaxVol])
											if x<nbcols-1: x=x+1
											else:
												x=0
												y=y+1
									
								fig.tight_layout()
								namefigpng='Scenario_'+str(NumScen)+'_Volume_'+hydrotech+'_'+country+start_week.strftime('%Y-%m-%d')+'.jpeg'
								plt.savefig(cfg['dirIMG']+namefigpng)
								plt.close()
								
					if cfg['PostTreat']['SpecificPeriods']['latex']:				
						for hydrotech in cfg['graphVolumes']:
							namefigpng='Scenario_'+str(NumScen)+'_Volume_'+hydrotech+'_'+country+start_week.strftime('%Y-%m-%d')+'.jpeg'
							if (hydrotech == 'Reservoir' and country in cfg['ReservoirRegions']) or (not hydrotech == 'Reservoir'):
								bodylatex_Det=bodylatex_Det+"\\subsubsection{"+cfg['graphVolumes'][hydrotech]['Name']+ \
									" Storages}\n"+figure(cfg['dirIMGLatex']+namefigpng,hydrotech+' storage level in '+country+' (MWh)',namefigpng)
						
					################################################################################################
					# creation of indicators
					################################################################################################
					logger.info('     computing indicators')
					# compute annual averages
					BilansSMS=SMSAggrPowerCountry.sum(axis=0)
					# aggregation per group of technologies
					for grouptechno in cfg['technosAggr']:
						tablebilansNRJ.loc[country][grouptechno]=BilansSMS.filter(items=cfg['technosAggr'][grouptechno]['technos'],axis=0).sum()
					
					# compute variable costs 
					CoutsVariables=pd.Series(index=BilansSMS.index,dtype='float64')
					CoutPropSMS=pd.Series(index=BilansSMS.index,dtype='float64')
					for techno in SMSAggrPowerCountry.columns:
						if 'varcost' in cfg['Technos'][techno].keys(): 
							CoutsVariables[techno]=cfg['Technos'][techno]['varcost']
						else: 
							CoutsVariables[techno]=0
						CoutPropSMS[techno]=(SMSAggrPowerCountry[techno]*CoutsVariables[techno]).sum()
					
					# compute nb of hours with non served energy
					NbHeuresDefSMS=SMSAggrPowerCountry['SlackUnit'].astype('bool').sum(axis=0)
					CoutTotalSMShorsVU=CoutPropSMS.sum()	
					CoutThSMS=CoutPropSMS.sum()
					if firstcountry==1:
						CoutTotalAllCountriesSMShorsVU=CoutTotalSMShorsVU
					else:
						CoutTotalAllCountriesSMShorsVU=CoutTotalAllCountriesSMShorsVU+CoutTotalSMShorsVU


					BilanNRJ=pd.DataFrame({'Energy (MWh)': BilansSMS})
					BilanCouthTH=pd.DataFrame({'Variable Cost (Euro)': CoutPropSMS })
					NbHeuresDef=pd.DataFrame([['Number Hours Loss of Load ',NbHeuresDefSMS]])
					CoutTotal=pd.DataFrame([['Total Cost w/o Hydro(Euro)', CoutTotalSMShorsVU]])
					#BilanNRJ['CO2 (t)']=BilanNRJ['Energy (MWh)']*CO2
					
					#tablebilansNRJ.loc[country]['CO2']=BilanNRJ['CO2 (t)'].sum()
					BilanNRJ=BilanNRJ.dropna()
					BilanCouthTH=BilanCouthTH.dropna()
					blankline=pd.DataFrame([[]])
					namecountry=pd.DataFrame([[country]])
					
					mytitle='Synthesis from '+ start_short +' to '+ end_short
					title=pd.DataFrame([[mytitle]])
					
					SMSAggrPowerCountry.to_csv(cfg['dirOUTScen']+'Prod_'+country+'.csv',sep=',',mode='w')
					SMSFlowsCountry.to_csv(cfg['dirOUTScen']+'Flows_'+country+'.csv',sep=',',mode='w')
					BilanNRJ.to_csv(cfg['dirOUTScen']+'Balance-'+country+'.csv',sep=',',mode='w')
					
					NbHeuresDef.to_csv(cfg['dirOUTScen']+'NbHoursNonServed-'+country+'.csv',sep=',',mode='w',header=False,index=False)
					BilanCouthTH.to_csv(cfg['dirOUTScen']+'ThermalCost-'+country+'.csv',sep=',',mode='w')
					CoutTotal.to_csv(cfg['dirOUTScen']+'TotalCost-'+country+'.csv',sep=',',mode='w',header=False,index=False)
					
					# graphs with average energy and ssy  
					P4RNRJtitle='Energy generated from '  + start_week.strftime('%Y-%m-%d') + ' to ' + end_week.strftime('%Y-%m-%d') 
					
					LinesPie=[i for i in BilanNRJ.index if 'PUMP' in i]
					BilanNRJ.drop(labels=LinesPie,axis=0,inplace=True)

					fig, axes = plt.subplots(nrows=1, ncols=1,figsize=(5,5))
					limitsms=BilanNRJ['Energy (MWh)'].sum(axis=0)*cfg['pielimit']
					BilanNRJ=BilanNRJ[ (BilanNRJ['Energy (MWh)'] > limitsms) ]
					if len(BilanNRJ.index)>0:
						PieCountryColors=[]
						labels=[]
						for techno in cfg['Technos']:
							if str(techno) in BilanNRJ.index:
								# techno is in the country
								PieCountryColors.append(cfg['Technos'][techno]['color'])
								labels.append('')
						BilanNRJ['Energy (MWh)'].plot.pie(ax=axes,colors=PieCountryColors,startangle=0,labels=labels,legend=False,autopct = lambda x: str(round(x, 2)) + '%',textprops=dict(color="w"))	
							
						axes.set_title(P4RNRJtitle,fontsize=10)
					
					fig.tight_layout()
					namefigpng='Scenario_'+str(NumScen)+'_Pie_'+country+start_week.strftime('%Y-%m-%d')+'.jpeg'
					plt.savefig(cfg['dirIMG']+namefigpng)
					plt.close()

					firstcountry=0
					
					if cfg['PostTreat']['SpecificPeriods']['latex']:
						BilanNRJ=pd.read_csv(cfg['dirOUTScen']+'Balance-'+country+'.csv',index_col=0)
						#bodylatex_Det=bodylatex_Det+"\\subsection{Energy Generated}\n"+tablelatex(BilanNRJ,'Energy Balance(MWh) in '+country)+"\n"
						namefigpng='Scenario_'+str(NumScen)+'_Pie_'+country+start_week.strftime('%Y-%m-%d')+'.jpeg'
						bodylatex_Det=bodylatex_Det+figure(cfg['dirIMGLatex']+namefigpng,'Energy Pie in '+country,namefigpng)
						NbHeuresDef=pd.read_csv(cfg['dirOUTScen']+'NbHoursNonServed-'+country+'.csv',header=None)
						bodylatex_Det=bodylatex_Det+str(NbHeuresDef[1][0])+" Hours with Loss of Load\n\n"
						
						BilanCouthTH=pd.read_csv(cfg['dirOUTScen']+'ThermalCost-'+country+'.csv',index_col=0)
						CoutTotal=pd.read_csv(cfg['dirOUTScen']+'TotalCost-'+country+'.csv',header=None)
						bodylatex_Det=bodylatex_Det+"\\subsection{Costs}\n"+"\nTotal Cost without Bellman Value="+str(CoutTotal[1][0])+"\n"
						#bodylatex_Det=bodylatex_Det+tablelatex(BilanCouthTH,'Thermal costs in '+country,column_format='|p{3cm}|p{3cm}|p{3cm}|p{3cm}|')+"\n"
					
				################################################################################################
				# end loop on countries
				################################################################################################

			if cfg['usevu']:
				logger.info('compute SMS VB')
				SMSVB=pd.read_csv(cfg['dirScen'] +'BellmanValuesOUT.csv')
				cols=['Timestep']+cfg['ReservoirRegions']+['b']
				SMSVB.columns=cols
				VolSMS=pd.DataFrame(columns=['Name','VolumeIni','VolumeFin']) 
				
				# keep only index corresponding to last period
				value=[0]
				SMSVBIni=SMSVB[ SMSVB.Timestep.isin(value) ]
				value=[cfg['timestepvusms']]
				SMSVBFin=SMSVB[ SMSVB.Timestep.isin(value) ]
				aFin=SMSVBFin[ cfg['ReservoirRegions'] ]
				bFin=SMSVBFin['b']
				aIni=SMSVBIni[ cfg['ReservoirRegions'] ]
				bIni=SMSVBIni['b']

				valuesmsIni=-10000000000000
				valuesmsFin=-10000000000000
				for row in aIni.index:
					valIni=bIni.iloc[row]+(aIni.iloc[row]*VolSMS['VolumeIni']).sum()
					if valIni > valuesmsIni: valuesmsIni=valIni
					valFin=bFin.iloc[row]+(aFin.iloc[row]*VolSMS['VolumeFin']).sum()
					if valFin > valuesmsFin: valuesmsFin=valFin

			SMSFlowsPerCountry=SMSFlowsPerCountry.drop_duplicates()
			SMSFlowsPerCountry = SMSFlowsPerCountry.loc[:,~SMSFlowsPerCountry.columns.duplicated()]
			BilansSMSFlows=SMSFlowsPerCountry.sum(axis=0).transpose()
			#CostFlowsSMS=pd.Series(index=SMSFlowsPerCountry.columns,dtype='float64')
						
			#for col in SMSFlowsPerCountry.columns:
			#	CostFlowsSMS[col]=cfg['flow_cost']*SMSFlowsPerCountry[col].abs().sum()

			if cfg['usevu']: CoutTotalAllCountriesSMS=CoutTotalAllCountriesSMShorsVU+valuesmsIni-valuesmsFin
			else: CoutTotalAllCountriesSMS=CoutTotalAllCountriesSMShorsVU
			if cfg['usevu'] : VBSMS=pd.DataFrame([['BV begin (Euro)', valuesmsIni, 'BV end (Euro)', valuesmsFin]])

			CoutTotalAllCountries=pd.DataFrame([['Total Cost (Euro)', CoutTotalAllCountriesSMS]])
			BilanFlows=pd.DataFrame({'Import/Export (MWh)': BilansSMSFlows})
			#CostFlows=pd.DataFrame({'Import/Export Cost (MWh)': CostFlowsSMS})

			BilanFlows=BilanFlows.dropna()
			BilanFlows=BilanFlows.drop_duplicates()
			BilanFlows.to_csv(cfg['dirOUTScen']+'BalanceFlows.csv',sep=',',mode='w')
			CoutTotalAllCountries.to_csv(cfg['dirOUTScen']+'TotalCostWithoutBV.csv',sep=',',mode='w',header=False,index=False)
			if cfg['usevu'] : VBSMS.to_csv(cfg['dirOUTScen']+'BV.csv',sep=',',mode='w',header=False,index=False)
			tablebilansNRJ.to_csv(cfg['dirOUTScen']+'Balance.csv',sep=',',mode='w')
				
			logger.info('Treat Import Exports')
			#CostFlows.to_csv(cfg['dirOUTScen']+'CostFlows'+'.csv',sep=',',mode='w',header=True,index=True)
			
			# build flows dataframe (import/export in same col)
			SMS_IE=pd.DataFrame(columns=cfg['lines'],index=SMSFlowsPerCountry.index)
			for i in range(len(cfg['lines'])):
				line=cfg['lines'][i]
				inverseline=line.split('>')[1]+'>'+line.split('>')[0]
				SMS_IE[line]=SMSFlowsPerCountry[line]-SMSFlowsPerCountry[inverseline]
			# graphs exports/imports
			if len(list_regions)>1:
				P4Rtitle='Imports/Exports ' 
				fig, axes = plt.subplots(figsize=(10,2*len(list_regions)),nrows=len(list_regions), ncols=1)
				i=0
				for country in list_regions:
					IEcols=[elem for elem in SMS_IE.columns if country in elem]
					if (len(IEcols)>0): SMS_IE[ IEcols ].plot(kind='line',ax=axes[i],legend=False)	
					axes[i].set_ylabel('Imports/Exports MWh',fontsize=8)
					axes[i].set_xlabel('')
					axes[i].tick_params(labelsize=6)
					axes[i].legend(labels=SMSFlowsPerCountry[ CountryFlowsCols[country] ].columns,bbox_to_anchor=(1,1),loc="upper left",fontsize=8)
					axes[i].set_title(country,fontsize=8)
					i=i+1
				fig.suptitle(P4Rtitle,fontsize=10)
				namefigpng='Scenario_'+str(NumScen)+'_ImportsExports'+start_week.strftime('%Y-%m-%d')+'.jpeg'
				i=0
				fig.tight_layout()
				plt.savefig(cfg['dirIMG']+namefigpng)
				plt.close()

			if cfg['PostTreat']['SpecificPeriods']['latex']:
				writelatex=True
				BilanFlows=pd.read_csv(cfg['dirOUTScen']+'BalanceFlows.csv',index_col=0)
				CoutTotalAllCountries=pd.read_csv(cfg['dirOUTScen']+'TotalCostWithoutBV.csv',header=None)
				if cfg['usevu'] : BV=pd.read_csv(cfg['dirOUTScen']+'BV.csv',header=None)
				
				bodylatex_Det=bodylatex_Det+"\\newpage\\section{Synthesis All Countries}\n"
				TotalCost=CoutTotalAllCountries[1][0]
				bodylatex_Det=bodylatex_Det+"\subsection{Costs}\n"+"Total Cost for all countries, without Bellman Values="+str(TotalCost)+" Euros\n\n"
				if cfg['usevu'] : TotalCost=TotalCost+BV[1][0]-BV[3][0]
				bodylatex_Det=bodylatex_Det+"Total Cost for all countries, with Bellman Values="+str(TotalCost)+" Euros\n\n"

				if cfg['usevu'] : 
					bodylatex_Det=bodylatex_Det+"Initial Bellman Value (Euro) ="+str(BV[1][0])+"\n"
					bodylatex_Det=bodylatex_Det+"Final Bellman Value (Euro) ="+str(BV[3][0])+"\n"
				
				#CostFlows=pd.read_csv(cfg['dirOUTScen']+'CostFlows.csv')
				if len(list_regions)>1:
					namefigpng='ImportsExports'+start_week.strftime('%Y-%m-%d')+'.jpeg'
					#bodylatex_Det=bodylatex_Det+"\\subsection{Imports/Exports}\n"+tablelatex(BilanFlows,'Imports/Exports (MWh)')+figure(cfg['dirIMGLatex']+namefigpng,'Imports/Exports (MWh)',namefigpng)


			###########################################################		
			# graphs continent maps	
			###########################################################		
				
			# draw continent with flows
			##################################################################
			if cfg['PostTreat']['SpecificPeriods']['draw'] and cfg['map']:
				# draw map
				figflows, axflows = plt.subplots(1,1,figsize=(10,10))
				mycontinent.boundary.plot(ax=axflows,figsize=(10,10))
				centers.plot(color='r',ax=axflows,markersize=1)

				# fill data with flows and start/end
				lines=SMSFlows.columns
				data={'start':[ x.split('>')[0] for x in lines ],'end':[ x.split('>')[1] for x in lines ],'flow':SMSFlows.sum(axis=0).transpose()} 
				flows = pd.DataFrame(data, index = SMSFlows.columns )

				# add reverseflows and delete
				for i in range(len(cfg['lines'])):
					line=cfg['lines'][i]
					inverseline=  line.split('>')[1]+'>'+line.split('>')[0]
					flow=flows.loc[line,'flow']
					inverseflow=flows.loc[inverseline,'flow']
					if(flow>inverseflow):
						flows.loc[line,'flow']=flow-inverseflow
						flows=flows.drop(inverseline,axis=0)
					else:
						flows.loc[inverseline,'flow']=inverseflow-flow
						flows=flows.drop(line,axis=0)
				
				# compute with of arrows
				if flows['flow'].abs().max()>0: flows['width']=10.0*flows['flow']/flows['flow'].abs().max()
				else: flows['width']=flows['flow']
				
				flows['startpoint']=flows['start'].apply(lambda x: centers[numcenters[x]])
				flows['endpoint']=flows['end'].apply(lambda x: centers[numcenters[x]])
				flows['line']=flows.apply(lambda x: LineString([x['startpoint'], x['endpoint']]),axis=1)
				geoflows=gpd.GeoDataFrame(flows,geometry=flows['line'])

				# reduce length of arrows
				scaledgeometry=geoflows.scale(xfact=0.7,yfact=0.7,zfact=1.0,origin='center')
				geoflows.geometry=scaledgeometry

				# plot lines
				geoflows.plot(ax=axflows,column='flow',linewidth=geoflows.width,cmap='Reds',vmin=-100,vmax=500)

				# plot arrows
				for line in geoflows.index:
					if geoflows['flow'][line]!=0:
						plt.arrow(list(geoflows['geometry'][line].coords)[0][0],list(geoflows['geometry'][line].coords)[0][1],
							list(geoflows['geometry'][line].coords)[1][0]-list(geoflows['geometry'][line].coords)[0][0],
							list(geoflows['geometry'][line].coords)[1][1]-list(geoflows['geometry'][line].coords)[0][1],
							head_width=1,head_length=0.5,color='black',linewidth=0,zorder=2)

				axflows.set_title("Import/Exports (MWh)",fontsize=10)
				namefigpng='Scenario_'+str(NumScen)+'_ArrowImpExp.jpeg'
				plt.savefig(cfg['dirIMG']+namefigpng)	
				plt.close()
				
			if cfg['PostTreat']['SpecificPeriods']['latex']:
				writelatex=True
				if len(list_regions)>1:
					namefigpng='Scenario_'+str(NumScen)+'_ArrowImpExp.jpeg'
					bodylatex_Det=bodylatex_Det+figure(cfg['dirIMGLatex']+namefigpng,'Map of Interconnection uses over the whole period (MWh)',namefigpng)

		# draw chloromap of continent
			if cfg['map']:
				tablebilansNRJ=tablebilansNRJ.fillna(0)
				mytable=tablebilansNRJ
				mytable['name']=mytable.index
				mytable=mytable.reset_index()
				
				
				# draw chloromap for C02
				mycontinent['CO2']=tablebilansNRJ['CO2']
				mycontinent=mycontinent.fillna(0)
				figchloroeurope, axchloroeurope = plt.subplots(1,1)
				mycontinent.plot(column='CO2',ax=axchloroeurope,cmap='Greys',legend=True,legend_kwds={'label':"CO2 Emissions (t)",'orientation':"horizontal"})
				axchloroeurope.set_title('Emissions (t)',fontsize=10)
				namefigpng='Scenario_ChloroMap-CO2.jpeg'
				plt.savefig(cfg['dirIMG']+namefigpng)
				plt.close()
			
			if cfg['PostTreat']['SpecificPeriods']['latex']:		
				writelatex=True
				bodylatex_Det=bodylatex_Det+"\\subsection{Energy generation}\n"
				tablebilansNRJ=pd.read_csv(cfg['dirOUTScen']+'Balance.csv',index_col=0)	
				#bodylatex_Det=bodylatex_Det+tablelatex(tablebilansNRJ,'Total Energy Balance (MWh)')
				
				namefigpng='Scenario_ChloroMap-MeanEnergy.jpeg'
				bodylatex_Det=bodylatex_Det+figure(cfg['dirIMGLatex']+namefigpng,'Map of energy generated for Scenario '+str(NumScen)+' (MWh)',namefigpng)
				
			###########################################################
			#continent pies
			#############################################################
			#if cfg['map']: namefigpng=EuropePieMap(tablebilansNRJ.drop(['CO2','name','Non Served'],axis=1),'EnergyPieEurope-Scenario.jpeg',TechnoAggrColors)
			if cfg['PostTreat']['SpecificPeriods']['latex']:
				writelatex=True
				namefigpng='Scenario_EnergyPieEurope.jpeg'
				bodylatex_Det=bodylatex_Det+figure(cfg['dirIMGLatex']+namefigpng,'Map of countries Energy Generated',namefigpng)
				namefigpng='Scenario_ChloroMap-CO2.jpeg'
				bodylatex_Det=bodylatex_Det+"\\subsection{CO2 Emissions}\n"+figure(cfg['dirIMGLatex']+namefigpng,'Map of Carbon Emissions (tons)',namefigpng)
		indexScen=indexScen+1			

	# creation of the reports
	if isLatex:
		if writelatex:
			containerlatex=startlatex+bodylatex_IC+bodylatex_Sto+bodylatex_Det+endlatex
			filelatex=cfg['dirOUT']+cfg['namereport']+'.tex'
			if os.path.exists(filelatex): os.remove(filelatex)
			myfile=open(filelatex,"x")
			myfile.write(containerlatex)
			myfile.close()

log_and_exit(0, cfg['path'])

