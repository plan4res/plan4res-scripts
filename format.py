#!/usr/bin/env python
# -*- coding: utf-8 -*-

## Import packages
import pandas as pd
import numpy as np
import os
import yaml
from itertools import product
from netCDF4 import Dataset
import calendar
import sys

from p4r_python_utils import *

path = os.environ.get("PLAN4RESROOT")
logger.info('path='+path)

def abspath_to_relpath(path, basepath):
	return os.path.relpath(path, basepath) if os.path.abspath(path) else path

nbargs=len(sys.argv)
if nbargs>1: 
	settings_format=sys.argv[1]
else:
	settings_format="settings_format.yml"
settings_format = abspath_to_relpath(settings_format, path)
    
    
# read config file
cfg1={}
with open(os.path.join(path, settings_format),"r") as mysettings:
	cfg1=yaml.load(mysettings,Loader=yaml.FullLoader)
settings_create = None
if nbargs>2:
	settings_create=abspath_to_relpath(sys.argv[2], path)
	cfg2={}
	with open(os.path.join(path, settings_create),"r") as mysettings:
		cfg2=yaml.load(mysettings,Loader=yaml.FullLoader)
	cfg = {**cfg1, **cfg2}
else:
	cfg=cfg1
	
# replace name of current dataset by name given as input
if nbargs>3:
	namedataset=sys.argv[3]
	if 'path' in cfg:
		cfg['path']=cfg['path'].replace(cfg['path'].split('/')[len(cfg['path'].split('/'))-2],namedataset)
	else:
		cfg['path']='/data/local/'+namedataset+'/'
cfg['inputpath']=os.path.join(cfg['path'], cfg['inputDir'])
cfg['outputpath']=os.path.join(cfg['path'], cfg['outputDir'])
if 'timeseriespath' not in cfg: cfg['timeseriespath']=os.path.join(cfg['path'],'TimeSeries/')

if cfg['USEPLAN4RESROOT']:
	path = os.environ.get("PLAN4RESROOT")
	cfg['path']=path+cfg['path']
	cfg['outputpath']=path+cfg['outputpath']
	cfg['inputpath']=path+cfg['inputpath']
	cfg['timeseriespath']=path+cfg['timeseriespath']
logger.info('path: '+cfg['inputpath'])

cfg['treat']=cfg['csvfiles']
cfg['treat']['ZP']=cfg['csvfiles']['ZP_ZonePartition']
cfg['treat']['IN']=cfg['csvfiles']['IN_Interconnections']
cfg['treat']['ZV']=cfg['csvfiles']['ZV_ZoneValues']
cfg['treat']['TU']=cfg['csvfiles']['TU_ThermalUnits']
cfg['treat']['SS']=cfg['csvfiles']['SS_SeasonalStorage']
cfg['treat']['STS']=cfg['csvfiles']['STS_ShortTermStorage']
cfg['treat']['RES']=cfg['csvfiles']['RES_RenewableUnits']
if 'SYN_SynchCond' in cfg['csvfiles']:
	cfg['treat']['SYN']=cfg['csvfiles']['SYN_SynchCond']


if not os.path.isdir(cfg['outputpath']):os.mkdir(cfg['outputpath'])

format=cfg['inputformat']
if format=='excel':
	p4r_excel=cfg['inputpath']+cfg['excelfile']
	excelfile=pd.read_excel(p4r_excel,None)
	sheets=list(excelfile.keys())
else:
	sheets=cfg['csvfiles']

def read_input_timeseries(cfg, ts_name, which_path='timeseriespath', **kwargs):
	file = os.path.join(cfg[which_path], ts_name)
	if not os.path.isfile(file):
		logger.error('File '+file+' does not exist. Use key '+which_path+' in configuration file '+settings_format+(' or configuration file '+settings_create if settings_create is not None else '')+' to specify input directory.')
		log_and_exit(2, cfg['path'])
	logger.info('Read '+file)
	return pd.read_csv(file, **kwargs)

#################################################################################################################################################
#																																				#
#											Read main data file(s)																						#
#																																				#
#																																				#
#################################################################################################################################################

# some of the sheets may have an additional row on top (where usually units can be written)
if 'additionalRowInSheets' in cfg and cfg['additionalRowInSheets']:
	skiprow=True
	skip=1
else:
	skiprow=False
	skip=0


# Read Parameters
#################################################################################################################################################

# create table of coupling constraints
ListPossibleTypesCoupling=['ActivePowerDemand','PrimaryDemand','SecondaryDemand','InertiaDemand','PollutantBudget']

ListCouplingConstraints=[elem for elem in cfg['CouplingConstraints'] if 'Pollutant' not in elem]
ListPollutants=[]
if 'PollutantBudget' in cfg['CouplingConstraints']:
	for pollutant in cfg['CouplingConstraints']['PollutantBudget']:
		ListPollutants=ListPollutants+[pollutant]
Coupling=pd.DataFrame(index=ListCouplingConstraints+ListPollutants,dtype=object,columns=['Partition','Sum'])

NumberPollutants=0
for coupling_constraint in cfg['CouplingConstraints']:
	if coupling_constraint not in ListPossibleTypesCoupling:
		logger.error('Constraint '+coupling_constraint+'is not possible')
		logger.error('Possible constraints are: '+', '.join(ListPossibleTypesCoupling))
		log_and_exit(1, cfg['path'])
	if coupling_constraint=="PollutantBudget":
		NumberPollutants=len(cfg['CouplingConstraints']['PollutantBudget'])
		for pollutant in ListPollutants:
			Coupling.loc[pollutant]['Sum']=[cfg['CouplingConstraints']['PollutantBudget'][pollutant]['SumOf']]
			Coupling.loc[pollutant]['Partition']=cfg['CouplingConstraints']['PollutantBudget'][pollutant]['Partition']
	else:
		Coupling.loc[coupling_constraint,'Sum']=cfg['CouplingConstraints'][coupling_constraint]['SumOf']
		Coupling.loc[coupling_constraint,'Partition']=cfg['CouplingConstraints'][coupling_constraint]['Partition']
logger.info('Coupling constraints:'+', '.join(ListCouplingConstraints))
logger.info('Emissions constraints:'+', '.join(ListPollutants))

# get dates
dates=pd.Series()
beginTS=pd.to_datetime(cfg['Calendar']['BeginTimeSeries'],dayfirst=cfg['Calendar']['dayfirst'])
dates['UCBeginData']=pd.Timestamp(year=beginTS.year,month=beginTS.month,day=beginTS.day,hour=beginTS.hour,minute=beginTS.minute)
beginDS=pd.to_datetime(cfg['Calendar']['BeginDataset'],dayfirst=cfg['Calendar']['dayfirst'])
dates['UCBegin']=pd.Timestamp(year=beginDS.year,month=beginDS.month,day=beginDS.day,hour=beginDS.hour,minute=beginDS.minute)
endTS=pd.to_datetime(cfg['Calendar']['EndTimeSeries'],dayfirst=cfg['Calendar']['dayfirst'])
dates['UCEndData']=pd.Timestamp(year=endTS.year,month=endTS.month,day=endTS.day,hour=endTS.hour,minute=endTS.minute)
endDS=pd.to_datetime(cfg['Calendar']['EndDataset'],dayfirst=cfg['Calendar']['dayfirst'])
dates['UCEnd']=pd.Timestamp(year=endDS.year,month=endDS.month,day=endDS.day,hour=endDS.hour,minute=endDS.minute)
logger.info('dates: timeseries start: '+str(dates['UCBeginData'])+' end: '+str(dates['UCEndData']))
logger.info('plan4res dataset start : '+str(dates['UCBegin'])+' end: '+str(dates['UCEnd']))

dates['UCBeginDataYearP4R']=pd.Timestamp(year=dates['UCBegin'].year,month=dates['UCBeginData'].month,day=dates['UCBeginData'].day,hour=dates['UCBeginData'].hour,minute=dates['UCBeginData'].minute)
DurationTimeSeries=pd.Timedelta(dates['UCEndData']-dates['UCBeginData'])
if dates['UCBegin']<dates['UCBeginDataYearP4R']:
	dates['UCBeginDataYearP4R']=pd.Timestamp(year=dates['UCBegin'].year-1,month=dates['UCBeginData'].month,day=dates['UCBeginData'].day,hour=dates['UCBeginData'].hour,minute=dates['UCBeginData'].minute)
dates['UCEndDataYearP4R']=dates['UCBeginDataYearP4R']+DurationTimeSeries

dates['UCBeginExtendedData']=dates['UCBeginDataYearP4R']
dates['UCEndExtendedData']=dates['UCEndDataYearP4R']

SSVTimeStep=cfg['Calendar']['SSVTimeStep']['Duration']
UnitSSV=cfg['Calendar']['SSVTimeStep']['Unit']
if UnitSSV=='days': SSVTimeStep=SSVTimeStep*24
if UnitSSV=='weeks': SSVTimeStep=SSVTimeStep*168
UCTimeStep=cfg['Calendar']['TimeStep']['Duration']
UnitUC=cfg['Calendar']['TimeStep']['Unit']
if UnitUC=='days': UCTimeStep=UCTimeStep*24
if UnitUC=='weeks': UCTimeStep=UCTimeStep*168

# get scenarios and convert list to string
ListScenarios=[str(elem) for elem in cfg['ParametersFormat']['Scenarios'] ]
ScenarisedData=cfg['ParametersFormat']['ScenarisedData']

# other parameters
if 'ThermalMaxPowerTimeSpan' in cfg['ParametersFormat']:
	ThermalMaxPowerTimeSpan=cfg['ParametersFormat']['ThermalMaxPowerTimeSpan']['Duration']
	UnitTMPTS=cfg['ParametersFormat']['ThermalMaxPowerTimeSpan']['Unit']
	if UnitTMPTS=='days': ThermalMaxPowerTimeSpan=ThermalMaxPowerTimeSpan*24
	if UnitTMPTS=='weeks': ThermalMaxPowerTimeSpan=ThermalMaxPowerTimeSpan*168
else: 
	ThermalMaxPowerTimeSpan=SSVTimeStep
CoeffSpillage=cfg['ParametersFormat']['CoeffSpillage']

# dates treatments
############################
durationData=dates.UCEndData-dates.UCBeginData+pd.Timedelta(hours=1)
durationInstance=dates.UCEnd-dates.UCBegin+pd.Timedelta(hours=1)
NumberUCTimeSteps=int((durationInstance.days*24+durationInstance.seconds/3600)/UCTimeStep)
durationUCTimeStep=pd.Timedelta(str(UCTimeStep)+' hours')
durationSSVTimeStep=pd.Timedelta(str(SSVTimeStep)+' hours')
logger.info('Duration instance:'+str(durationInstance))

if not ((durationInstance.total_seconds()/3600)/SSVTimeStep)-int((durationInstance.total_seconds()/3600)/SSVTimeStep)==0:
	NewNumberUCTimeSteps=int((durationInstance.total_seconds()/3600)/SSVTimeStep)*SSVTimeStep/UCTimeStep
	NumberUCTimeStepsToDelete=NumberUCTimeSteps-NewNumberUCTimeSteps
	NumberUCTimeSteps=int(NewNumberUCTimeSteps)
	durationToDelete=pd.Timedelta(str(NumberUCTimeStepsToDelete*UCTimeStep)+' hours')
	dates['UCEnd']=dates['UCEnd']-durationToDelete
	durationInstance=dates.UCEnd-dates.UCBegin+pd.Timedelta(hours=1)
logger.info('Number of time steps:'+str(NumberUCTimeSteps)+' of duration:'+str(durationUCTimeStep))
	
TimeHorizonUC=int(SSVTimeStep/UCTimeStep)
NumberIntervals=TimeHorizonUC
NumberSSVTimeSteps=int(NumberUCTimeSteps*UCTimeStep/SSVTimeStep)
TimeHorizonSSV=NumberSSVTimeSteps
logger.info('Number of SSV time steps:'+str(NumberSSVTimeSteps)+' of duration:'+str(TimeHorizonUC)+' timesteps = '+str(durationSSVTimeStep))

# create dataframe with start and end dates of all SSV timesteps
datesSSV=pd.DataFrame(index=list(range(NumberSSVTimeSteps)),columns=['start','end'])
start=dates['UCBegin']
for i in range(NumberSSVTimeSteps):
	datesSSV.loc[i]=[start,start+durationSSVTimeStep-pd.Timedelta('1 hours')]
	start=start+durationSSVTimeStep
	
# create dataframe with start and end dates of all UC timesteps
datesUC=pd.DataFrame(index=list(range(NumberUCTimeSteps)),columns=['start','end'])
start=dates['UCBegin']
for i in range(NumberUCTimeSteps):
	datesUC.loc[i]=[start,start+durationUCTimeStep]
	start=start+durationUCTimeStep
	
# create dataframe with start and end dates of all UC timesteps
datesData=pd.DataFrame(columns=['start','end'])
start=dates['UCBeginData']
i=0
while start<=dates['UCEndData']:
	datesData.loc[i]=[start,start+durationUCTimeStep]
	start=start+durationUCTimeStep
	i=i+1

# check if csv files have to be modified to include results of investment_solver
CreateDataPostInvest=False # True if csv files are to be re-created former computation of investments
if 'RecomputeCSV' in cfg['ParametersFormat']:
	if cfg['ParametersFormat']['RecomputeCSV'] and os.path.isfile(cfg['path']+'results_invest/Solution_OUT.csv'):
		CreateDataPostInvest=True
		solInvest=check_and_read_csv(cfg, cfg['path']+'results_invest/Solution_OUT.csv',header=None)
		indexSolInvest=0
			
# Read sheet ZP_ZonePartition
#################################################################################################################################################
if 'ZP_ZonePartition' in sheets:
	if format=='excel':
		ZP=pd.read_excel(p4r_excel,sheet_name='ZP_ZonePartition',skiprows=0,index_col=None)
	else:
		ZP=read_input_csv(cfg, 'ZP_ZonePartition', skiprows=0, index_col=None)
	ZP=ZP.drop_duplicates()					
	Nodes=ZP[Coupling['Partition'].loc['ActivePowerDemand']]
	NumberNodes=len(Nodes)
	NumberSlackUnits=len(Nodes)  # there is one slack unit per node

	# create partitions
	# Partition['Leveli'] is a dictionnary whose keys are the regions at level i and values are the list of nodes in each region
	Partition=pd.Series(dtype=object)
	for level in ZP.columns:
		if level==Coupling['Partition'].loc['ActivePowerDemand']: 
			Partition.loc[level]=dict(zip(list(Nodes),[ [node] for node in list(Nodes)]))
		else:
			keys=list(set(list(ZP[level])))
			values= []
			for key in keys:
				values.append(list(set(list(ZP.loc[ZP[level]==key][Coupling['Partition'].loc['ActivePowerDemand']]))))
			Partition.loc[level]=dict(zip(keys,values))
			
	if 'PrimaryDemand' in Coupling.index:
		NumberPrimaryZones=len(Partition.loc[Coupling['Partition']['PrimaryDemand']])
	else: NumberPrimaryZones=0
	if 'SecondaryDemand' in Coupling.index:
		NumberSecondaryZones=len(Partition.loc[Coupling['Partition']['SecondaryDemand']])
	else: NumberSecondaryZones=0
	if 'InertiaDemand' in Coupling.index:
		NumberInertiaZones=len(Partition.loc[Coupling['Partition']['InertiaDemand']])
	else: NumberInertiaZones=0
	TotalNumberPollutantZones=sum(len(Partition.loc[Coupling['Partition'][elem]]) for elem in ListPollutants)
else:
	logger.error('ZP_Partition missing')
	log_and_exit(2, cfg['path'])

# Read sheet ZV_ZoneValues
#################################################################################################################################################
if 'ZV_ZoneValues' in sheets:
	if format=='excel':
		ZV=pd.read_excel(p4r_excel,sheet_name='ZV_ZoneValues',skiprows=0,index_col=['Type', 'Zone'])
	else:
		ZV=read_input_csv(cfg, 'ZV_ZoneValues',skiprows=0,index_col=['Type', 'Zone'])
	ZV=ZV.drop_duplicates()					
	if 'Profile_Timeserie' in ZV.columns:
		ZV['Profile_Timeserie']=ZV['Profile_Timeserie'].fillna('')
	else:
		ZV['Profile_Timeserie']=''
	ZV['Profile_Timeserie']=ZV['Profile_Timeserie'].fillna('')														   
	ZV=ZV.fillna(0)
else: 
	logger.error('ZV_ZoneValues missing')
	log_and_exit(1, cfg['path'])

InstalledCapacity=pd.DataFrame(index=Nodes)

# Read sheet SS_SeasonalStorage
#################################################################################################################################################if 'TU_ThermalUnits' in sheets:
if 'SS_SeasonalStorage' in sheets:
	if format=='excel':
		SS=pd.read_excel(p4r_excel,sheet_name='SS_SeasonalStorage',skiprows=skip,index_col=['Name','Zone'])
	else:
		SS=read_input_csv(cfg, 'SS_SeasonalStorage',skiprows=skip,index_col=['Name','Zone'])
	if not SS.empty:
		SS=SS.drop_duplicates()				 
		SS=SS.drop( SS[ SS['NumberUnits']==0 ].index )
		SS=SS.drop( SS[ SS['MaxPower']==0.0 ].index )
		SS['InflowsProfile']=SS['InflowsProfile'].fillna('')
		if 'WaterValues' in SS:
			SS['WaterValues']=SS['WaterValues'].fillna('')
		TotalNumberHydroUnits=SS['NumberUnits'].sum()
		
		if 'HydroSystem' in SS.columns:
			NumberHydroSystems=int(SS['HydroSystem'].max()+1)
		else:
			NumberHydroSystems=1
		SS['NumberReservoirs']=1
		SS['NumberArcs']=1
		for reservoir in SS.index:
			if 'MinPower' in SS.columns and SS['MinPower'][reservoir]<0:
				SS.loc[reservoir]['NumberReservoirs']=2
				SS.loc[reservoir]['NumberArcs']=2
		TotalNumberReservoirs=SS['NumberReservoirs'].sum()
		HSSS=pd.Series(dtype=object)
		for hs in range(NumberHydroSystems):
			HSSS.loc[hs]=SS[ SS['HydroSystem']==hs ]
	else: 
		logger.warning('No seasonal storage mix in this dataset')
		NumberHydroSystems=0	
else: 
	logger.warning('No seasonal storage mix in this dataset')
	NumberHydroSystems=0

  
  
  
# Read sheet TU_ThermalUnits
#################################################################################################################################################
if 'TU_ThermalUnits' in sheets:
	if format=='excel':
		TU=pd.read_excel(p4r_excel,sheet_name='TU_ThermalUnits',skiprows=skip,index_col=['Name','Zone'])
	else:
		TU=read_input_csv(cfg, 'TU_ThermalUnits', skiprows=skip, index_col=['Name','Zone'])
	TU=TU[TU['NumberUnits'] != 0]
	if CreateDataPostInvest:
		save_input_csv(cfg, 'TU_ThermalUnits',TU)
		for row in TU.index:
			if (TU.loc[row,'MaxAddedCapacity']>0)+(TU.loc[row,'MaxRetCapacity']>0):
				if solInvest[0].loc[indexSolInvest]>1:
					TU.loc[row,'MaxAddedCapacity']=TU.loc[row,'MaxAddedCapacity']-(solInvest[0].loc[indexSolInvest]-1)*TU.loc[row,'MaxPower']
					if TU.loc[row,'MaxAddedCapacity']<=cfg['ParametersCreate']['zerocapacity']:
						TU.loc[row,'MaxAddedCapacity']=0
				if solInvest[0].loc[indexSolInvest]<1:
					TU.loc[row,'MaxRetCapacity']=TU.loc[row,'MaxRetCapacity']-(1-solInvest[0].loc[indexSolInvest])*TU.loc[row,'MaxPower']
					if TU.loc[row,'MaxRetCapacity']<=cfg['ParametersCreate']['zerocapacity']:
						TU.loc[row,'MaxRetCapacity']=0
				logger.info('Added Capacity to TU '+str(row)+' :'+str(TU.loc[row,'MaxPower']*solInvest[0].loc[indexSolInvest]-TU.loc[row,'MaxPower']))
				for c in ['MaxPower', 'MinPower', 'Capacity']:
					if c in TU.columns:
						TU.loc[row,c] = np.round(TU.loc[row,c]*solInvest[0].loc[indexSolInvest], decimals=cfg['ParametersFormat']['RoundDecimals'])
				indexSolInvest=indexSolInvest+1
		write_input_csv(cfg, 'TU_ThermalUnits',TU)
	TU=TU.drop( TU[ TU['NumberUnits']==0 ].index )
	if ('MaxAddedCapacity' not in TU.columns and 'MaxRetCapacity' not in TU.columns):
		TU=TU.drop( TU[ TU['MaxPower']<=cfg['ParametersCreate']['zerocapacity'] ].index )
	else:
		TU=TU.drop( TU[ ( TU['MaxPower']        <=cfg['ParametersCreate']['zerocapacity'] )    \
			&           ( TU['MaxAddedCapacity']<=cfg['ParametersCreate']['zerocapacity'] )    \
			&           ( TU['MaxRetCapacity']  <=cfg['ParametersCreate']['zerocapacity'] ) ].index )

	TU=TU.drop_duplicates()
	if 'MaxPowerProfile' in TU.columns:
		TU['MaxPowerProfile']=TU['MaxPowerProfile'].fillna('')
	NumberThermalUnits=TU['NumberUnits'].sum()
	if ('MaxAddedCapacity' in TU.columns and 'MaxRetCapacity' in TU.columns):
		NumberInvestedThermalUnits=len(TU[ (TU['MaxRetCapacity']>0) | (TU['MaxAddedCapacity']>0) ])
	elif 'MaxAddedCapacity' in TU.columns:
		NumberInvestedThermalUnits=len(TU[ TU['MaxAddedCapacity']>0 ])
	elif 'MaxRetCapacity' in TU.columns:
		NumberInvestedThermalUnits=len(TU[ TU['MaxRetCapacity']>0 ])

	else:
		NumberInvestedThermalUnits=0
else: 
	logger.warning('No thermal mix in this dataset')
	NumberThermalUnits=0

# Read sheet RES_RenewableUnits
#################################################################################################################################################
if 'RES_RenewableUnits' in sheets:
	if format=='excel':
		RES=pd.read_excel(p4r_excel,sheet_name='RES_RenewableUnits',skiprows=skip,index_col=['Name','Zone'])
	else:
		RES=read_input_csv(cfg, 'RES_RenewableUnits',skiprows=skip,index_col=['Name','Zone'])
	if CreateDataPostInvest:
		save_input_csv(cfg, 'RES_RenewableUnits',RES)
		for row in RES.index:
			if (RES.loc[row,'MaxAddedCapacity']>0)+(RES.loc[row,'MaxRetCapacity']>0):
				if solInvest[0].loc[indexSolInvest]>1:
					RES.loc[row,'MaxAddedCapacity']=RES.loc[row,'MaxAddedCapacity']-(solInvest[0].loc[indexSolInvest]-1)*RES.loc[row,'MaxPower']
					if RES.loc[row,'MaxAddedCapacity']<=cfg['ParametersCreate']['zerocapacity']:
						RES.loc[row,'MaxAddedCapacity']=0
				if solInvest[0].loc[indexSolInvest]<1:
					RES.loc[row,'MaxRetCapacity']=RES.loc[row,'MaxRetCapacity']-(1-solInvest[0].loc[indexSolInvest])*RES.loc[row,'MaxPower']
					if RES.loc[row,'MaxRetCapacity']<=cfg['ParametersCreate']['zerocapacity']:
						RES.loc[row,'MaxRetCapacity']=0
				logger.info('Added Capacity to RES '+str(row)+' :'+str(RES.loc[row,'MaxPower']*solInvest[0].loc[indexSolInvest]-RES.loc[row,'MaxPower']))

				for c in ['MaxPower', 'MinPower', 'Capacity']:
					if c in RES.columns:
						RES.loc[row,c] = np.round(RES.loc[row,c]*solInvest[0].loc[indexSolInvest], decimals=cfg['ParametersFormat']['RoundDecimals'])
				indexSolInvest=indexSolInvest+1
		write_input_csv(cfg, 'RES_RenewableUnits',RES)

	RES=RES.drop( RES[ RES['NumberUnits']==0 ].index )
	if ('MaxAddedCapacity' not in RES.columns and 'MaxRetCapacity' not in RES.columns):
		RES=RES.drop( RES[ RES['MaxPower']<=cfg['ParametersCreate']['zerocapacity'] ].index )
	else:
		RES=RES.drop( RES[ (RES['MaxPower']<=cfg['ParametersCreate']['zerocapacity']) & ( RES['MaxAddedCapacity']<=cfg['ParametersCreate']['zerocapacity']) & (RES['MaxRetCapacity']<=cfg['ParametersCreate']['zerocapacity'] ) ].index )

	RES=RES.drop_duplicates()
	if 'Energy_Timeserie' in RES.columns and 'Energy' in RES.columns:
		RES['EnergyMaxPower']=RES.apply(lambda x: x.Energy if x.Name=="Hydro|Run of River" else x.Energy_Timeserie*x.MaxPower,axis=1)
	RES['MaxPowerProfile']=RES['MaxPowerProfile'].fillna('')
	NumberIntermittentUnits=RES['NumberUnits'].sum()
	if ('MaxAddedCapacity' in RES.columns and 'MaxRetCapacity' in RES.columns):
		NumberInvestedIntermittentUnits=len(RES[ (RES['MaxRetCapacity']>0) | (RES['MaxAddedCapacity']>0) ])
	elif 'MaxAddedCapacity' in RES.columns:
		NumberInvestedIntermittentUnits=len(RES[ RES['MaxAddedCapacity']>0 ])
	elif 'MaxRetCapacity' in RES.columns:
		NumberInvestedIntermittentUnits=len(RES[ RES['MaxRetCapacity']>0 ])
	else:
		NumberInvestedIntermittentUnits=0

else: 
	logger.warning('No intermittent generation mix in this dataset')
	NumberIntermittentUnits=0

# Read sheet STS_ShortTermStorage
#################################################################################################################################################
if 'STS_ShortTermStorage' in sheets:
	if format=='excel':
		STS=pd.read_excel(p4r_excel,sheet_name='STS_ShortTermStorage',skiprows=skip,index_col=None)
	else:
		STS=read_input_csv(cfg, 'STS_ShortTermStorage',skiprows=skip,index_col=None)
	if CreateDataPostInvest:
		save_input_csv(cfg, 'STS_ShortTermStorage',STS)
		for row in STS.index:
			if (STS.loc[row,'MaxAddedCapacity']>0)+(STS.loc[row,'MaxRetCapacity']>0):
				if solInvest[0].loc[indexSolInvest]>1:
					STS.loc[row,'MaxAddedCapacity']=STS.loc[row,'MaxAddedCapacity']-(solInvest[0].loc[indexSolInvest]-1)*STS.loc[row,'MaxPower']
					if STS.loc[row,'MaxAddedCapacity']<=cfg['ParametersCreate']['zerocapacity']:
						STS.loc[row,'MaxAddedCapacity']=0
				if solInvest[0].loc[indexSolInvest]<1:
					STS.loc[row,'MaxRetCapacity']=STS.loc[row,'MaxRetCapacity']-(1-solInvest[0].loc[indexSolInvest])*STS.loc[row,'MaxPower']
					if STS.loc[row,'MaxRetCapacity']<=cfg['ParametersCreate']['zerocapacity']:
						STS.loc[row,'MaxRetCapacity']=0
				logger.info('Added Capacity to STS '+str(row)+' :'+str(STS.loc[row,'MaxPower']*solInvest[0].loc[indexSolInvest]-STS.loc[row,'MaxPower']))

				for c in ['MaxPower', 'MinPower', 'Capacity','MaxVolume','MinVolume']:
					if c in STS.columns:
						STS.loc[row,c] = np.round(STS.loc[row,c]*solInvest[0].loc[indexSolInvest], decimals=cfg['ParametersFormat']['RoundDecimals'])
				indexSolInvest=indexSolInvest+1
		write_input_csv(cfg, 'STS_ShortTermStorage',STS)

	STS=STS.drop( STS[ STS['NumberUnits']==0 ].index )
	if ('MaxAddedCapacity' not in STS.columns and 'MaxRetCapacity' not in STS.columns):
		STS=STS.drop( STS[ STS['MaxPower']<=cfg['ParametersCreate']['zerocapacity'] ].index )
	else:
		STS=STS.drop( STS[ (STS['MaxPower']<=cfg['ParametersCreate']['zerocapacity']) & (STS['MaxAddedCapacity']<=cfg['ParametersCreate']['zerocapacity']) & (STS['MaxRetCapacity']<=cfg['ParametersCreate']['zerocapacity'] ) ].index )

	STS=STS.drop_duplicates()
	STS=STS.set_index(['Name','Zone'])
	NumberBatteryUnits=STS['NumberUnits'].sum()
	if ('MaxAddedCapacity' in STS.columns and 'MaxRetCapacity' in STS.columns):
		NumberInvestedBatteryUnits=len(STS[ (STS['MaxRetCapacity']>0) | (STS['MaxAddedCapacity']>0) ])
	elif 'MaxAddedCapacity' in STS.columns:
		NumberInvestedBatteryUnits=len(STS[ STS['MaxAddedCapacity']>0 ])
	elif 'MaxRetCapacity' in STS.columns:
		NumberInvestedBatteryUnits=len(STS[ STS['MaxRetCapacity']>0 ])
	else:
		NumberInvestedBatteryUnits=0
else: 
	logger.warning('No short term storage mix in this dataset')
	NumberBatteryUnits=0

# Read sheet SYN SynchCond
#################################################################################################################################################
if 'SYN_SynchCond' in sheets:
	if format=='excel':
		SYN=pd.read_excel(p4r_excel,sheet_name='SYN_SynchCond',skiprows=skip,index_col=None)
	else:
		SYN=pd.read_excel(cfg['inputpath']+cfg['csvfiles']['SYN_SynchCond'],skiprows=skip,index_col=None)
	SYN=SYN.drop( SYN[ SYN['NumberUnits']==0 ].index )
	SYN=SYN.drop( SYN[ SYN['MaxRotatingConsumption']==0.0 ].index )
	SYN=SYN.drop_duplicates()
	SYN=SYN.set_index(['Name','Zone'])
	NumberSyncUnits=SYN['NumberUnits'].sum()
else: 
	logger.warning('No synchronous condensers mix in this dataset')
	NumberSyncUnits=0

# Read sheet IN_Interconnections
#################################################################################################################################################
if 'IN_Interconnections' in sheets:
	if format=='excel':
		IN=pd.read_excel(p4r_excel,sheet_name='IN_Interconnections',skiprows=skip,index_col=0)
	else:
		IN=read_input_csv(cfg, 'IN_Interconnections', skiprows=skip,index_col=0)
	if CreateDataPostInvest:
		save_input_csv(cfg, 'IN_Interconnections',IN,index=True)
		for row in IN.index:
			if (IN.loc[row,'MaxAddedCapacity']>0)+(IN.loc[row,'MaxRetCapacity']>0):
				if solInvest[0].loc[indexSolInvest]>1:
					IN.loc[row,'MaxAddedCapacity']=IN.loc[row,'MaxAddedCapacity']-(solInvest[0].loc[indexSolInvest]-1)*IN.loc[row,'MaxPowerFlow']
					if IN.loc[row,'MaxAddedCapacity']<=cfg['ParametersCreate']['zerocapacity']:
						IN.loc[row,'MaxAddedCapacity']=0
				if solInvest[0].loc[indexSolInvest]<1:
					IN.loc[row,'MaxRetCapacity']=IN.loc[row,'MaxRetCapacity']-(1-solInvest[0].loc[indexSolInvest])*IN.loc[row,'MaxPowerFlow']
					if IN.loc[row,'MaxRetCapacity']<=cfg['ParametersCreate']['zerocapacity']:
						IN.loc[row,'MaxRetCapacity']=0
				logger.info('Added Capacity to IN '+str(row)+' :'+str(IN.loc[row,'MaxPowerFlow']*solInvest[0].loc[indexSolInvest]-IN.loc[row,'MaxPowerFlow']))

				for c in ['MaxPowerFlow', 'MinPowerFlow']:
					if c in IN.columns:
						IN.loc[row,c] = np.round(IN.loc[row,c]*solInvest[0].loc[indexSolInvest], decimals=cfg['ParametersFormat']['RoundDecimals'])
				indexSolInvest=indexSolInvest+1
		write_input_csv(cfg, 'IN_Interconnections',IN,index=True)
	
	if ('MaxAddedCapacity' not in IN.columns and 'MaxRetCapacity' not in IN.columns):
		IN=IN.drop( IN[ (IN['MaxPowerFlow']<=cfg['ParametersCreate']['zerocapacity']) & (IN['MinPowerFlow']>=(-1)*cfg['ParametersCreate']['zerocapacity'])  ].index )
	else:
		IN=IN.drop( IN[ (IN['MaxPowerFlow']<=cfg['ParametersCreate']['zerocapacity']) & (IN['MinPowerFlow']>=(-1)*cfg['ParametersCreate']['zerocapacity']) & (IN['MaxAddedCapacity']<=cfg['ParametersCreate']['zerocapacity']) & (IN['MaxRetCapacity']<=cfg['ParametersCreate']['zerocapacity'])  ].index )

	IN=IN.drop_duplicates()
	NumberLines=len(IN.index)
	if ('MaxAddedCapacity' in IN.columns and 'MaxRetCapacity' in IN.columns):
		NumberInvestedLines=len(IN[ (IN['MaxRetCapacity']>0) | (IN['MaxAddedCapacity']>0) ])
	elif 'MaxAddedCapacity' in IN.columns:
		NumberInvestedLines=len(IN[ IN['MaxAddedCapacity']>0 ])
	elif 'MaxRetCapacity' in IN.columns:
		NumberInvestedLines=len(IN[ IN['MaxRetCapacity']>0 ])
	else:
		NumberInvestedLines=0
else: 
	logger.warning('No interconnections in this dataset')
	NumberLines=0
	
#################################################################################################################################################
#																																				#
#											Compute data																						#
#																																				#
#																																				#
#################################################################################################################################################
NumberUnits=NumberHydroSystems+NumberThermalUnits+NumberBatteryUnits+NumberIntermittentUnits+NumberSyncUnits+NumberSlackUnits
if NumberHydroSystems>0: 
	NumberArcs=SS['NumberArcs'].sum()
	NumberHydroUnits=int(SS['NumberUnits'].sum())
else: 
	NumberArcs=0
	NumberHydroUnits=0
NumberElectricalGenerators=NumberArcs+NumberThermalUnits+NumberBatteryUnits+NumberIntermittentUnits+NumberSyncUnits+NumberSlackUnits
logger.info(str(NumberHydroSystems)+' Hydrosystems')
logger.info(str(NumberHydroUnits)+' Hydro Units with '+str(NumberArcs)+' generators')
logger.info(str(NumberThermalUnits)+' Thermal Units')
logger.info(str(NumberBatteryUnits)+' Short term Storage Units')
logger.info(str(NumberIntermittentUnits)+' Intermittent Units')
logger.info(str(NumberSyncUnits)+' Synchronous condensers')
logger.info(str(NumberSlackUnits)+' Slack Units')
logger.info(str(NumberUnits)+' units, '+str(NumberElectricalGenerators)+' generators')

#################################################################################################################################################
#																																				#
#											Read timeseries																						#
#																																				#
#																																				#
#################################################################################################################################################
def ExtendAndResample(name,TS,isEnergy=True):
	# change timeindex to adapt to dataset calendarb
	beginSerie=TS.index[0]
	if beginSerie>dates['UCBeginDataYearP4R']:
		datesDelta=pd.Timedelta(beginSerie-dates['UCBeginDataYearP4R'])
		TS.index=TS.index-datesDelta
	else:
		datesDelta=pd.Timedelta(dates['UCBeginDataYearP4R']-beginSerie)
		TS.index=TS.index+datesDelta
	
	# Extension is a copy of TS on the extended dates UCBeginExtendedData and UCEndExtendedData 
	Extension=TS[ TS.index>= dates['UCBeginExtendedData'] ]
	Extension=Extension[ Extension.index<= dates['UCEndExtendedData'] ]
	
	# resample
	newfreq=str(UCTimeStep)+'h'
	TS_freq=pd.infer_freq(TS.index)

	# upsample=True means that timeseries is given at a frequency
	# bigger than hour (eg 2 hours, daily, weekly)
	upsample=False
	
	# calcul de la frequence de la série en nombre d'heures
	if 'D' in TS_freq or 'W' in TS_freq or 'M' in TS_freq:
		upsample=True
		if 'D' in TS_freq:
			if len(TS_freq)>1:
				Hours_freq=int(TS_freq[:-1])*24
			else:
				Hours_freq=24
		if 'W' in TS_freq:
			if '-' in TS_freq: TS_freq=TS_freq.split('-')[0]
			if len(TS_freq)>1:
				Hours_freq=int(TS_freq[:-1])*168
			else:
				Hours_freq=168
		if 'M' in TS_freq:
			if len(TS_freq)>1:
				Hours_freq=int(TS_freq[:-1])*728
			else:
				Hours_freq=728
				
	if 'H' in TS_freq or 'h' in TS_freq:
		if len(TS_freq)>1:
			Hours_freq=int(TS_freq[:-1])
		else:
			Hours_freq=1
		if Hours_freq>int(UCTimeStep): upsample=True
	

	duration_TS_timestep=pd.Timedelta(str(Hours_freq)+' hours')
	
	if Hours_freq==1: 
		Extension.index=Extension.index+durationData
	else:
		Extension.index=Extension.index+pd.Timedelta(TS.index[-1]-TS.index[0])+pd.Timedelta(str(Hours_freq)+' hours')
	TS=pd.concat([TS,Extension])

	# case where timeserie is given at frequency bigger than hour
	# resample to hour frequecy before resampling to the required frequency
	if upsample:
		TS=TS.resample('h').ffill()		# convert to hourly frequency

		# extend with missing dates: duplicate last dates
		if TS.index[-1]< dates['UCEnd']:
			dur_missing=dates['UCEnd']-TS.index[-1] # compute duration of missing data
			Extension=TS[ TS.index> (TS.index[-1]-dur_missing) ] # take last period of TS of this duration
			Extension.index=Extension.index+dur_missing # shift over time
			TS=pd.concat([TS,Extension]) # add at end of serie

		TS=TS.resample(newfreq).sum()
	else:
		TS=TS.resample(newfreq).sum()

	# keep only period of dataset
	TS=TS[ TS.index>= dates['UCBegin'] ]
	TS=TS[ TS.index<= dates['UCEnd'] ]
	
	# case where there is only one value in TS_freq
	if len(TS.index)==1:
		TS2=TS[ TS.index>= dates['UCBegin'] ]
		TS2.index=TS2.index+pd.Timedelta(str(Hours_freq)+' hours')
		TS=pd.concat([TS,TS2])
	return TS

def read_deterministic_timeseries(IsDT):
	if IsDT:
		DeterministicTS=pd.read_csv(cfg['inputpath']+cfg['DeterministicTimeSeries'],index_col=0)
		DeterministicTS.index=pd.to_datetime(DeterministicTS.index,dayfirst=cfg['Calendar']['dayfirst'])
	
		DeterministicTS=ExtendAndResample('DET',DeterministicTS)

		# add constant serie
		DeterministicTS['One']=1.0
		DeterministicTS['Zero']=0.0
	else:
		DeterministicTS=pd.DataFrame(index=datesData['start'])
		DeterministicTS.index=pd.to_datetime(DeterministicTS.index)
		DeterministicTS['One']=1.0
		DeterministicTS['Zero']=0.0
		DeterministicTS=ExtendAndResample('DET',DeterministicTS)
	return DeterministicTS
	
def create_demand_scenarios():
	DemandScenarios=pd.Series(dtype=object)
	isEnergy=(cfg['ParametersFormat']['ScenarisedData']['ActivePowerDemand']['MultiplyTimeSerieBy']=='Energy')
	for node in Nodes:
		firstPart=True
		for component in Coupling.loc['ActivePowerDemand']['Sum']:
			nameTS=ZV.loc[component,node]['Profile_Timeserie']
			valTS=ZV.loc[component,node]['value']
			# case without a profile
			if nameTS=='':
				if firstPart:
					DemandScenarios.loc[node]=pd.DataFrame(columns=ListScenarios)
					for col in ListScenarios: DemandScenarios.loc[node][col]=DeterministicTimeSeries['One']
					firstPart=False
				else:
					for col in ListScenarios: DemandScenarios.loc[node][col]=DemandScenarios.loc[node][col]+DeterministicTimeSeries['One']
			# read timeserie if deterministic
			elif '.csv' in nameTS: # stochastic OR deterministic series
				isDeterministic=False
				TS=read_input_timeseries(cfg, nameTS, skiprows=0,index_col=0)
				if len(TS.columns)==1: isDeterministic=True # the serie is deterministic
					
				TS.index=pd.to_datetime(TS.index,dayfirst=cfg['Calendar']['dayfirst'])
				TS=ExtendAndResample(nameTS,TS,isEnergy)
					
				if firstPart:
					DemandScenarios[node]=pd.DataFrame(index=TS.index,columns=ListScenarios)
					if isDeterministic: 
						for col in ListScenarios: DemandScenarios.loc[node][col]=valTS*TS[TS.columns.tolist()[0]]
					else: DemandScenarios[node]=valTS*TS
					firstPart=False
				else:
					if isDeterministic:
						for col in ListScenarios: DemandScenarios.loc[node][col]=DemandScenarios.loc[node][col]+valTS*TS[TS.columns.tolist()[0]]
					else: DemandScenarios[node]=DemandScenarios.loc[node]+valTS*TS
					
			else:
				# deterministic serie in DeterministicTimeSeries dataframe
				DTS=valTS*DeterministicTimeSeries[nameTS]
				if firstPart:
					DemandScenarios.loc[node]=pd.DataFrame(columns=ListScenarios)
					for col in ListScenarios: DemandScenarios.loc[node][col]=DTS
					firstPart=False
				else:
					for col in ListScenarios: DemandScenarios.loc[node][col]=DemandScenarios.loc[node][col]+DTS
				
		
	return DemandScenarios
	
def create_inflows_scenarios():
	isEnergy=(cfg['ParametersFormat']['ScenarisedData']['Hydro:Inflows']['MultiplyTimeSerieBy']['reservoir']=='Energy')
	InflowsScenarios=pd.Series(dtype=object,index=SS.index)
	
	for reservoir in SS.index:
		nameTS=SS.loc[reservoir]['InflowsProfile']
		valTS=SS.loc[reservoir]['Inflows']
		# read timeserie 
		# case without a profile
		if nameTS=='':
			InflowsScenarios[reservoir]=pd.DataFrame(columns=ListScenarios)
			for col in ListScenarios: InflowsScenarios.loc[reservoir][col]=(valTS/cfg['ParametersFormat']['NumberHoursInYear'])*DeterministicTimeSeries['One'] # valTS is an energy per year
		elif '.csv' in nameTS:  # stochastic series
			TS=read_input_timeseries(cfg, nameTS, index_col=0)
			TS.index=pd.to_datetime(TS.index,dayfirst=cfg['Calendar']['dayfirst'])
			TS=ExtendAndResample(nameTS,TS,isEnergy)
			if len(TS.columns) > 1: # stochastic serie
				InflowsScenarios[reservoir]=pd.DataFrame(index=TS.index,columns=TS.columns)
				InflowsScenarios[reservoir]=valTS*TS
			else: # deterministic serie
				InflowsScenarios[reservoir]=pd.DataFrame(index=TS.index,columns=ListScenarios)
				for col in ListScenarios: InflowsScenarios.loc[reservoir][col]=valTS*TS[TS.columns.tolist()[0]]
		else:  # deterministic series
			DTS=valTS*DeterministicTimeSeries[nameTS]
			InflowsScenarios.loc[reservoir]=pd.DataFrame(columns=ListScenarios)
			for col in ListScenarios: InflowsScenarios.loc[reservoir][col]=DTS
					
	return InflowsScenarios
	
def create_res_scenarios():
	ResScenarios=pd.Series(dtype=object,index=RES.index)
	newIndex=[]
	for res in RES.index:
		nameTS=RES.loc[res]['MaxPowerProfile']
		valTS=RES.loc[res]['MaxPower']
		nameAsset=res[0]
		technoAsset = None
		for namekind in cfg['technos']:
			for nametechno in cfg['technos'][namekind]:
				if nametechno in nameAsset:
					technoAsset=namekind
		if technoAsset is None:
			logger.error('The techno for '+str(res)+' described in file RES_RenewableUnits.csv could not be identified. Check the data listed under the "technos" key in the configuration files.')
			log_and_exit(2, cfg['path'])
		
		isEnergy=(cfg['ParametersFormat']['ScenarisedData']['Renewable:MaxPowerProfile']['MultiplyTimeSerieBy'][technoAsset]=='Energy')
					
		# read timeserie 
		if '.csv' in nameTS:  # stochastic series
			TS=read_input_timeseries(cfg, nameTS, skiprows=0,index_col=0)
			TS.index=pd.to_datetime(TS.index,dayfirst=cfg['Calendar']['dayfirst'])
			TS=ExtendAndResample(nameTS,TS,isEnergy)
			if len(TS.columns) > 1: # stochastic serie
				ResScenarios[res]=pd.DataFrame(index=TS.index,columns=TS.columns)
				ResScenarios[res]=valTS*TS
				
			else: #deterministic
				ResScenarios[res]=pd.DataFrame(index=TS.index,columns=ListScenarios)
				for col in ListScenarios: ResScenarios[res][col]=valTS*TS[TS.columns.tolist()[0]]
				
			newIndex.append(res)
#
	
	ResScenarios=ResScenarios[ newIndex ]
	return ResScenarios
	
def create_thermal_scenarios():
	# define wether there is one profile per techno/region, or one profile per unit
	newIndex=[]
	isEnergy=(cfg['ParametersFormat']['ScenarisedData']['Thermal:MaxPowerProfile']['MultiplyTimeSerieBy']=='Energy')

	if True in np.column_stack(TU['MaxPowerProfile'].str.contains(r",", na=False)): # there exist multiple profiles for the same row
		multipleTHSeries=True
		listIndexes=[list(TU.index[0]),list(TU.index[1]),list(range(TU['NumberUnits'].max()))]
		ThermalScenarios=pd.Series(dtype=object,index=pd.MultiIndex.from_tuples((list(product(*listIndexes)))))
		for th in TU.index:
			for unit in range(TU.loc[th]['NumberUnits']):
				if ',' in str(TU.loc[th]['MaxPowerProfile']):
					nameTS=TU.loc[th]['MaxPowerProfile'].split(',')[unit]
				else:
					nameTS=TU.loc[th]['MaxPowerProfile']
				valTS=TU.loc[th]['MaxPower']
				
				# read timeserie 
				if '.csv' in nameTS:  # stochastic series
					TS=read_input_timeseries(cfg, nameTS,skiprows=0,index_col=0)
					TS.index=pd.to_datetime(TS.index,dayfirst=cfg['Calendar']['dayfirst'])
					TS=ExtendAndResample(nameTS,TS,isEnergy)
					
					if len(TS.columns) > 1: # stochastic serie
						ThermalScenarios[th[0],th[1],unit]=pd.DataFrame(index=TS.index,columns=TS.columns)
						ThermalScenarios[th[0],th[1],unit]=valTS*TS
					else:
						ThermalScenarios[th[0],th[1],unit]=pd.DataFrame(index=TS.index,columns=ListScenarios)
						for col in ListScenarios: ThermalScenarios[th[0],th[1],unit]=valTS*TS[TS.columns.tolist()[0]]
					newIndex.append( (th[0],th[1],unit) )
	else:
		multipleTHSeries=False
		ThermalScenarios=pd.Series(dtype=object,index=TU.index)
		for th in TU.index:
			nameTS=TU.loc[th]['MaxPowerProfile']
			valTS=TU.loc[th]['MaxPower']
			
			# read timeserie 
			if '.csv' in nameTS:  # stochastic series
				TS=read_input_timeseries(cfg,nameTS,skiprows=0,index_col=0)
				TS.index=pd.to_datetime(TS.index,dayfirst=cfg['Calendar']['dayfirst'])
				TS=ExtendAndResample(nameTS,TS)
				
				if len(TS.columns) > 1: # stochastic serie
					ThermalScenarios[th]=pd.DataFrame(index=TS.index,columns=TS.columns)
					ThermalScenarios[th]=valTS*TS
				else:
					ThermalScenarios[th]=pd.DataFrame(index=TS.index,columns=ListScenarios)
					for col in ListScenarios: ThermalScenarios[th]=valTS*TS[TS.columns.tolist()[0]]
			newIndex.append(th)
	
	ThermalScenarios=ThermalScenarios[newIndex]
	return ThermalScenarios
	
	
#################################################################################################################################################
#																																				#
#											Create netcDF files																					#
#																																				#
#																																				#
#################################################################################################################################################
# create the HydroUnitBlocks
def addHydroUnitBlocks(Block,indexUnitBlock,scenario,start,end,id):
	for hydrosystem in range(NumberHydroSystems):
		# create all hydrosystems
		HSBlock=Block.createGroup('UnitBlock_'+str(indexUnitBlock))
		HSBlock.type="HydroSystemUnitBlock"
		NumberHydroUnitsInHydroSystem=len(HSSS.loc[hydrosystem].index)
		HSBlock.createDimension("NumberHydroUnits",NumberHydroUnitsInHydroSystem)
		indexHU=0
		NumberReservoirsInHydroSystem=0
		for hydrounit in HSSS.loc[hydrosystem].index:
			# create all hydro units
			for index_subunit in range(HSSS.loc[hydrosystem]['NumberUnits'][hydrounit]):
				HBlock=HSBlock.createGroup('HydroUnitBlock_'+str(indexHU))
				HBlock.type="HydroUnitBlock"
				HBlock.setncattr("name",hydrounit[0]+'_'+hydrounit[1]+'_'+str(index_subunit))
				NumberReservoirs=HSSS.loc[hydrosystem]['NumberReservoirs'][hydrounit]
				NumberReservoirsInHydroSystem=NumberReservoirsInHydroSystem+NumberReservoirs
				HBlock.createDimension("NumberReservoirs",NumberReservoirs)
				NumberArcs=HSSS.loc[hydrosystem]['NumberArcs'][hydrounit]
				HBlock.createDimension("NumberArcs",NumberArcs)
				HBlock.createDimension("NumberIntervals",NumberIntervals)
				
				# create arcs
				StartArc=HBlock.createVariable("StartArc","u4",("NumberArcs"))
				EndArc=HBlock.createVariable("EndArc","u4",("NumberArcs"))
				if NumberArcs ==1:
					StartArc[:]=[0]
					EndArc[:]=[1]
				elif NumberArcs ==2:
					StartArc[:]=[0,1]
					EndArc[:]=[1,0]
				
				# create min and max volume
				MaxVolumeData=HSSS.loc[hydrosystem]['MaxVolume'][hydrounit]
				if type(MaxVolumeData)==str:
					MaxVolumetric=HBlock.createVariable("MaxVolumetric",np.double,("NumberReservoirs","NumberIntervals"))
					vmax=DeterministicTimeSeries[MaxVolumeData][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]
					if NumberReservoirs==1:
						MaxVolumetric[0,:]=vmax
					else:
						MaxVolumetric[0,:]=vmax
						MaxVolumetric[1,:]=cfg['DownReservoirVolumeMultFctor']*vmax
				else:
					MaxVolumetric=HBlock.createVariable("MaxVolumetric",np.double,("NumberReservoirs"))
					if NumberReservoirs==1:
						MaxVolumetric[0]=[MaxVolumeData]
					elif NumberReservoirs==2:
						MaxVolumetric[0]=[MaxVolumeData]
						MaxVolumetric[1]=[cfg['DownReservoirVolumeMultFctor']*MaxVolumeData]
					vmax=MaxVolumeData*DeterministicTimeSeries['One'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]
				if 'MinVolume' in SS.columns:
					MinVolumeData=HSSS.loc[hydrosystem]['MinVolume'][hydrounit]
					if type(MinVolumeData)==str:
						MinVolumetric=HBlock.createVariable("MinVolumetric",np.double,("NumberReservoirs","NumberIntervals"))
						if NumberReservoirs==1:
							MinVolumetric[0,:]=np.minimum(vmax,DeterministicTimeSeries[MinVolumeData][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
						elif NumberReservoirs==2:
							MinVolumetric[0,:]=np.minimum(vmax,DeterministicTimeSeries[MinVolumeData][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
							MinVolumetric[1,:]=np.array(DeterministicTimeSeries['Zero'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
					else:
						
						if type(MaxVolumeData)==str and  ( (MinVolumeData>vmax).isin([True]).sum()>0 ):
							MinVolumetric=HBlock.createVariable("MinVolumetric",np.double,("NumberReservoirs","NumberIntervals"))
							if NumberReservoirs==1:
								MinVolumetric[0,:]=np.minimum(vmax,MinVolumeData*DeterministicTimeSeries['One'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
							elif NumberReservoirs==2:
								MinVolumetric[0,:]=np.minimum(vmax,MinVolumeData*DeterministicTimeSeries['One'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
								MinVolumetric[1,:]=np.array(DeterministicTimeSeries['Zero'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
						else:
							MinVolumetric=HBlock.createVariable("MinVolumetric",np.double,("NumberReservoirs"))
							if NumberReservoirs==1:
								MinVolumetric[0]=[MinVolumeData]
							elif NumberReservoirs==2:
								MinVolumetric[0]=[MinVolumeData]
								MinVolumetric[1]=[0.0]

				# create inflows
				if ('Inflows' in SS.columns) and ('InflowsProfile' in SS.columns):
					Inflows=HBlock.createVariable("Inflows",np.double,("NumberReservoirs","NumberIntervals"))
					InflowsData=HSSS.loc[hydrosystem]['Inflows'][hydrounit]
					InflowsDataProfile=HSSS.loc[hydrosystem]['InflowsProfile'][hydrounit]
					if cfg['IncludeScenarisedData'] and 'Hydro:Inflows' in ScenarisedData:
						inflow=np.array(InflowsScenarios.loc[hydrounit][scenario][ ( InflowsScenarios.loc[hydrounit].index >= start ) & ( InflowsScenarios.loc[hydrounit].index <= end ) ])
					else:
						inflow=np.array(DeterministicTimeSeries['Zero'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
					if NumberReservoirs==1:
						Inflows[0,:]=inflow
					elif NumberReservoirs==2:
						Inflows[0,:]=inflow
						Inflows[1,:]=np.array(DeterministicTimeSeries['Zero'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
				elif ('Inflows' in SS.columns) and ('InflowsProfile' not in SS.columns):
					InflowsData=HSSS.loc[hydrosystem]['Inflows'][hydrounit]
					Inflows=HBlock.createVariable("Inflows",np.double,("NumberReservoirs"))
					if NumberReservoirs==1:
						Inflows[0]=[InflowsData]
					elif NumberReservoirs==2:
						Inflows[0]=[InflowsData]
						Inflows[1]=[0.0]
				elif ('Inflows' not in SS.columns) and ('InflowsProfile' in SS.columns):
					Inflows=HBlock.createVariable("Inflows",np.double,("NumberReservoirs","NumberIntervals"))
					InflowsDataProfile=HSSS.loc[hydrosystem]['InflowsProfile'][hydrounit]
					if cfg['IncludeScenarisedData'] and 'Hydro:Inflows' in ScenarisedData:
						inflow=np.array(InflowsScenarios.loc[hydrounit][scenario][ ( InflowsScenarios.loc[hydrounit].index >= start ) & ( InflowsScenarios.loc[hydrounit].index <= end ) ])
					else:
						inflow=np.array(DeterministicTimeSeries['Zero'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
					if NumberReservoirs==1:
						Inflows[0,:]=inflow
					elif NumberReservoirs==2:
						Inflows[0,:]=inflow
						Inflows[1,:]=np.array(DeterministicTimeSeries['Zero'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])

				# create min and max power and flow
				MaxPowerData=HSSS.loc[hydrosystem]['MaxPower'][hydrounit]
				if 'TurbineEfficiency' in SS.columns: TurbineEfficiency=HSSS.loc[hydrosystem]['TurbineEfficiency'][hydrounit]
				else: TurbineEfficiency=1
				if 'PumpingEfficiency' in SS.columns: PumpingEfficiency=HSSS.loc[hydrosystem]['PumpingEfficiency'][hydrounit]
				else: PumpingEfficiency=cfg['PumpingEfficiency']['reservoir']['Reservoir']
				if type(MaxPowerData)==str:
					MaxFlow=HBlock.createVariable("MaxFlow",np.double,("NumberIntervals","NumberArcs"))
					MaxPower=HBlock.createVariable("MaxPower",np.double,("NumberIntervals","NumberArcs"))
					pmax=DeterministicTimeSeries[MaxPowerData][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]
					Pmax=pd.DataFrame(pmax)
					if NumberArcs==1: 
						Pmax['Spillage']=CoeffSpillage*pmax
						Pmax=Pmax.transpose()
						Flow=(1+CoeffSpillage)*(1/TurbineEfficiency)*Pmax
						Pmax['Spillage']=DeterministicTimeSeries['Zero'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]
						for t in range(NumberIntervals): 
							MaxPower[t,:]=np.array(Pmax[t])
							MaxFlow[t,:]=np.array(Flow[t])
					elif NumberArcs==2:
						Pmax['zero']=0.0
						Pmax['Spillage']=CoeffSpillage*pmax
						Pmax=Pmax.transpose()
						Flow=(1+CoeffSpillage)*(1/TurbineEfficiency)*Pmax
						Pmax['Spillage']=DeterministicTimeSeries['Zero'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]
						for t in range(NumberIntervals): 
							MaxPower[t,:]=np.array(Pmax[t])
							MaxFlow[t,:]=np.array(Flow[t])
				else:
					MaxFlow=HBlock.createVariable("MaxFlow",np.double,("NumberArcs"))
					MaxPower=HBlock.createVariable("MaxPower",np.double,("NumberArcs"))
					if NumberArcs==1: 
						MaxPower[:]=[MaxPowerData*UCTimeStep]
						MaxFlow[:]=[(1+CoeffSpillage)*(1/TurbineEfficiency)*MaxPowerData*UCTimeStep]
					elif NumberArcs==2:
						MaxPower[:]=[MaxPowerData*UCTimeStep,0]
						MaxFlow[:]=[(1+CoeffSpillage)*(1/TurbineEfficiency)*MaxPowerData*UCTimeStep,(1/TurbineEfficiency)*CoeffSpillage*MaxPowerData*UCTimeStep]
					pmax=MaxPowerData*UCTimeStep*DeterministicTimeSeries['One'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]
				if 'MinPower' in SS.columns:
					MinPowerData=HSSS.loc[hydrosystem]['MinPower'][hydrounit]
					if type(MinPowerData)==str:
						MinFlow=HBlock.createVariable("MinFlow",np.double,("NumberIntervals","NumberArcs"))
						MinPower=HBlock.createVariable("MinPower",np.double,("NumberIntervals","NumberArcs"))
						pmin=np.maximum(DeterministicTimeSeries['Zero'],np.minimum(pmax,DeterministicTimeSeries[MinPowerData][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]))
						Pmin=pd.DataFrame(pmin)
						if NumberArcs==1: 
							Pmin['Zero']=DeterministicTimeSeries['Zero'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]
							Pmin=Pmin.transpose()
							for t in range(NumberIntervals): 
								MinPower[t,:]=np.array(Pmin[t])
								MinFlow[t,:]=np.array((1/PumpingEfficiency)*Pmin[t])
						elif NumberArcs==2:
							Pmin['Min']=np.minimum(DeterministicTimeSeries['Zero'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ],DeterministicTimeSeries[MinPowerData][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
							Pmin['Zero']=DeterministicTimeSeries['Zero'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]
							Pmin=Pmin.transpose()
							for t in range(NumberIntervals): 
								MinPower[t,:]=np.array(Pmin[t])
								MinFlow[t,:]=np.array((1/PumpingEfficiency)*Pmin[t])
					else:
						if type(MaxPowerData)==str and  ( (MinPowerData>pmax).isin([True]).sum()>0 ):
							MinFlow=HBlock.createVariable("MinFlow",np.double,("NumberIntervals","NumberArcs"))
							MinPower=HBlock.createVariable("MinPower",np.double,("NumberIntervals","NumberArcs"))
							# case where MaxPower may go under MinPower
							pmin=np.minimum(pmax,MinPowerData*DeterministicTimeSeries['One'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
							Pmin=pd.DataFrame(pmin)
							if NumberArcs==1: 
								Pmin['Zero']=DeterministicTimeSeries['Zero'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]
								for t in range(NumberIntervals): 
									MinPower[t,:]=np.array(Pmin[t])
									MinFlow[t,:]=np.array((1/PumpingEfficiency)*Pmin[t])
							elif NumberArcs==2:
								Pmin['Min']=MinPowerData*DeterministicTimeSeries['One'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]
								Pmin['Zero']=DeterministicTimeSeries['Zero'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]
								Pmin=Pmin.transpose()
								for t in range(NumberIntervals): 
									MinPower[t,:]=np.array(Pmin[t])
									MinFlow[t,:]=np.array((1/PumpingEfficiency)*Pmin[t])
						else:
							MinFlow=HBlock.createVariable("MinFlow",np.double,("NumberArcs"))
							MinPower=HBlock.createVariable("MinPower",np.double,("NumberArcs"))
							if NumberArcs==1: 
								MinPower[:]=[max(0,MinPowerData*UCTimeStep)]
								MinFlow[:]=[(1/TurbineEfficiency)*max(0,MinPowerData*UCTimeStep)]
							elif NumberArcs==2:
								MinPower[:]=[max(0,MinPowerData*UCTimeStep),(1/PumpingEfficiency)*MinPowerData*UCTimeStep]
								MinFlow[:]=[max(0,(1/TurbineEfficiency)*MinPowerData*UCTimeStep),(1/PumpingEfficiency)*MinPowerData*UCTimeStep]
								
				# create ramping constraints
				if 'DeltaRampUp' in SS.columns:
					DeltaRampUpData=HSSS.loc[hydrosystem]['DeltaRampUp'][hydrounit]
					if type(DeltaRampUpData)==str:
						DeltaRampUp=HBlock.createVariable("DeltaRampUp",np.double,("NumberIntervals","NumberArcs"))
						ramp=DeterministicTimeSeries[DeltaRampUpData][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]
						Ramp=pd.DataFrame(ramp)
						Ramp['zero']=0.0
						if NumberReservoirs==1:
							Ramp=Ramp.transpose()
							for t in range(NumberIntervals):
								DeltaRampUp[t,:]=[Ramp[t]]
						elif NumberReservoirs==2:
							Ramp['zero2']=0.0
							Ramp=Ramp.transpose()
							for t in range(NumberIntervals):
								DeltaRampUp[t,:]=[Ramp[t],cfg['ParametersFormat']['DownDeltaRampUpMultFactor']*Ramp[t]]
					else:
						DeltaRampUp=HBlock.createVariable("DeltaRampUp",np.double,("NumberArcs"))
						if NumberReservoirs==1:
							DeltaRampUp[:]=[DeltaRampUpData*UCTimeStep]
						elif NumberReservoirs==2:
							DeltaRampUp[:]=[DeltaRampUpData*UCTimeStep,cfg['ParametersFormat']['DownDeltaRampUpMultFactor']*DeltaRampUpData*UCTimeStep]
				if 'DeltaRampDown' in SS.columns:
					DeltaRampDownData=HSSS.loc[hydrosystem]['DeltaRampUp'][hydrounit]
					if type(DeltaRampDownData)==str:
						DeltaRampDown=HBlock.createVariable("DeltaRampDown",np.double,("NumberIntervals","NumberArcs"))
						ramp=DeterministicTimeSeries[DeltaRampDownData][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]
						Ramp=pd.DataFrame(ramp)
						Ramp['zero']=0.0
						if NumberReservoirs==1:
							Ramp=Ramp.transpose()
							for t in range(NumberIntervals):
								DeltaRampDown[t,:]=[Ramp[t]]
						elif NumberReservoirs==2:
							Ramp['zero2']=0.0
							Ramp=Ramp.transpose()
							for t in range(NumberIntervals):
								DeltaRampDown[t,:]=[Ramp[t],cfg['ParametersFormat']['DownDeltaRampDownMultFactor']*Ramp[t]]
					else:
						DeltaRampDown=HBlock.createVariable("DeltaRampDown",np.double,("NumberArcs"))
						if NumberReservoirs==1:
							DeltaRampDown[:]=[DeltaRampDownData*UCTimeStep]
						elif NumberReservoirs==2:
							DeltaRampDown[:]=[DeltaRampDownData*UCTimeStep,cfg['ParametersFormat']['DownDeltaRampDownMultFactor']*DeltaRampDownData*UCTimeStep]					
				
				# create primary and secondary rho
				if 'PrimaryRho' in SS.columns:
					PrimaryRhoData=HSSS.loc[hydrosystem]['PrimaryRho'][hydrounit]
					if type(PrimaryRhoData)==str:
						PrimaryRho=HBlock.createVariable("PrimaryRho",np.double,("NumberIntervals","NumberArcs"))
						ramp=DeterministicTimeSeries[PrimaryRhoData][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]
						Ramp=pd.DataFrame(ramp)
						Ramp['zero']=0.0
						Ramp=Ramp.transpose()
						if NumberReservoirs==1:
							for t in range(NumberIntervals):
								PrimaryRho[t,:]=[Ramp[t]]
						elif NumberReservoirs==2:
							for t in range(NumberIntervals):
								PrimaryRho[t,:]=[Ramp[t],0]
					else:
						PrimaryRho=HBlock.createVariable("PrimaryRho",np.double,("NumberArcs"))
						if NumberReservoirs==1:
							PrimaryRho[:]=[PrimaryRhoData]
						elif NumberReservoirs==2:
							PrimaryRho[:]=[PrimaryRhoData,0]		
				if 'SecondaryRho' in SS.columns:
					SecondaryRhoData=HSSS.loc[hydrosystem]['SecondaryRho'][hydrounit]
					if type(SecondaryRhoData)==str:
						SecondaryRho=HBlock.createVariable("SecondaryRho",np.double,("NumberIntervals","NumberArcs"))
						ramp=DeterministicTimeSeries[SecondaryRhoData][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]
						Ramp=pd.DataFrame(ramp)
						Ramp['zero']=0.0
						Ramp=Ramp.transpose()
						if NumberReservoirs==1:
							for t in range(NumberIntervals):
								SecondaryRho[t,:]=[Ramp[t]]
						elif NumberReservoirs==2:
							for t in range(NumberIntervals):
								SecondaryRho[t,:]=[Ramp[t],0]
					else:
						SecondaryRho=HBlock.createVariable("SecondaryRho",np.double,("NumberArcs"))
						if NumberReservoirs==1:
							SecondaryRho[:]=[SecondaryRhoData]
						elif NumberReservoirs==2:
							SecondaryRho[:]=[SecondaryRhoData,0]
				
				# create efficiency data : piecewise linear function (with only 1 piece)
				NumberPieces=HBlock.createVariable("NumberPieces","u4",("NumberArcs"))
				if NumberReservoirs==1:
					NumberPieces[:]=[1]
					TotalNumberPieces=1
				elif NumberReservoirs==2:
					NumberPieces[:]=[1,1]
					TotalNumberPieces=2
				LinearTerm=HBlock.createVariable("LinearTerm",np.double,("NumberArcs"))
				if NumberReservoirs==1:
					LinearTerm[:]=[TurbineEfficiency]
				elif NumberReservoirs==2:
					LinearTerm[:]=[TurbineEfficiency,PumpingEfficiency]
				ConstantTerm=HBlock.createVariable("ConstantTerm",np.double,("NumberArcs"))
				if NumberReservoirs==1:
					ConstantTerm[:]=[0]
				elif NumberReservoirs==2:
					ConstantTerm[:]=[0,0]
					
				# create inertia data
				if 'Inertia' in SS.columns:
					InertiaData=HSSS.loc[hydrosystem]['Inertia'][hydrounit]
					if type(InertiaData)==str:
						InertiaPower=HBlock.createVariable("InertiaPower",np.double,("NumberIntervals","NumberArcs"))
						inertia=cfg['ParametersFormat']['InertiaMultFactor']*DeterministicTimeSeries[InertiaData][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]
						if NumberReservoirs==1:
							for t in range(NumberIntervals):
								InertiaPower[t,:]=[Inertia[t]]
						elif NumberReservoirs==2:
							for t in range(NumberIntervals):
								InertiaPower[t,:]=[Inertia[t],0]
					else:
						InertiaPower=HBlock.createVariable("InertiaPower",np.double,("NumberArcs"))
						if NumberReservoirs==1:
							InertiaPower[:]=[cfg['ParametersFormat']['InertiaMultFactor']*InertiaData]
						elif NumberReservoirs==2:
							InertiaPower[:]=[cfg['ParametersFormat']['InertiaMultFactor']*InertiaData,0]
					
				# create initial conditions
				InitialVolumetric=HBlock.createVariable("InitialVolumetric",np.double,("NumberReservoirs"))
				if 'InitialVolume' in SS.columns:
					initialvolume=HSSS.loc[hydrosystem]['InitialVolume'][hydrounit]
				if NumberReservoirs==1:
					InitialVolumetric[:]=[initialvolume]
				elif NumberReservoirs==2:
					InitialVolumetric[:]=[initialvolume,cfg['ParametersFormat']['DownReservoirVolumeMultFactor']*initialvolume]
				
				InitialFlowRate=HBlock.createVariable("InitialFlowRate",np.double,("NumberArcs"))
				if 'InitialFlowRate' in SS.columns:
					InitialFlowRateData=HSSS.loc[hydrosystem]['InitialFlowRate'][hydrounit]
				else:
					InitialFlowRateData=0.0
				if NumberReservoirs==1:
					InitialFlowRate[:]=[InitialFlowRateData*UCTimeStep]
				elif NumberReservoirs==2:
					InitialFlowRate[:]=[InitialFlowRateData*UCTimeStep,0]
									
				indexHU=indexHU+1
		# add polyhedral function for water values
		PolyhedralFunctionBlock=HSBlock.createGroup('PolyhedralFunctionBlock')
		PolyhedralFunctionBlock.type="PolyhedralFunctionBlock"
		
		# compute size of function which is equal to number of reservoirs
		NumVar=HSSS.loc[hydrosystem]['NumberReservoirs'].sum()
		PolyhedralFunctionBlock.createDimension("PolyFunction_NumVar", NumVar)
		
		# create polyhedral function only when requested
		if ('WaterValues' in SS.columns) and (cfg['IncludeVU']!='None' and cfg['FormatVU']=='PerReservoir') and ((cfg['IncludeVU']=='Last' and id==NumberSSVTimeSteps-1) or cfg['IncludeVU']=='All' or (cfg['IncludeVU']=='Last' and cfg['FormatMode']=='SingleUC')):

			indexReservoir=0;
			NumberReservoirsInHydroSystem=HSSS.loc[hydrosystem]['NumberReservoirs'].sum()
			# water values are given per each reservoir
			for hydrounit in HSSS.loc[hydrosystem].index:
				for index_subunit in range(HSSS.loc[hydrosystem]['NumberUnits'][hydrounit]):
					
					BVfile=HSSS.loc[hydrosystem]['WaterValues'][hydrounit]
					if BVfile!='':
						logger.info('Add Bellman values from file')
						BVdata=read_input_timeseries(cfg,BVfile,'inputpath',index_col=0,skiprows=skip)
						BVdata.index=pd.to_datetime(BVdata.index,dayfirst=cfg['Calendar']['dayfirst'])
						# keep only data included in the period of the block
						BVdata=BVdata[ (BVdata.index >=start) & (BVdata.index <=end)   ]
						
						# keep only the last date within this pediod
						BVdata=BVdata[ BVdata.index== BVdata.index.max() ]
						BVsize=len(BVdata.index)-1
						
						if BVsize>0:
							# compute coefficients a
							if indexReservoir==0:
								PolyhedralFunctionBlock.createDimension("PolyFunction_NumRow", BVsize)
								Data_A=np.zeros(shape=(BVsize,NumberReservoirsInHydroSystem))
								Data_B=np.zeros(shape=(BVsize))
							# order BVData per growing volumes
							BVdata=BVdata.sort_values(by=['Volume'], ascending=True)
							
							for i in range(BVsize):
								Data_A[i][indexReservoir]=(BVdata.iloc[i+1]['Value']-BVdata.iloc[i]['Value'])/(BVdata.iloc[i+1]['Volume']-BVdata.iloc[i]['Volume'])
								Data_B[i]=Data_B[i]+BVdata.iloc[i]['Value']-Data_A[i][indexReservoir]*BVdata.iloc[i]['Volume']
						else:
							if indexReservoir==0:
								PolyhedralFunctionBlock.createDimension("PolyFunction_NumRow", 1)
								Data_A=np.zeros(shape=(BVsize,NumberReservoirsInHydroSystem))
								Data_B=np.zeros(shape=(BVsize))
					else:
						if indexReservoir==0:
							PolyhedralFunctionBlock.createDimension("PolyFunction_NumRow", 1)
							Data_A=np.zeros(shape=(1,NumberReservoirsInHydroSystem))
							Data_B=np.zeros(shape=1)
							
					indexReservoir=indexReservoir+HSSS.loc[hydrosystem]['NumberReservoirs'][hydrounit]
			
			PolyFunction_A=PolyhedralFunctionBlock.createVariable("PolyFunction_A",np.double,("PolyFunction_NumRow","PolyFunction_NumVar"))
			PolyFunction_b=PolyhedralFunctionBlock.createVariable("PolyFunction_b",np.double,("PolyFunction_NumRow"))
			for row in range(BVsize):
				PolyFunction_A[row,:]=Data_A[row]
				PolyFunction_b[row]=Data_B[row]
					
		# cases polyhedral water values
		###############################
		
		# case with a cuts file
		elif (cfg['FormatVU']=='Polyhedral') and ('WaterValues' in SS.columns) and ( (cfg['IncludeVU']=='Last' and id==NumberSSVTimeSteps-1) or cfg['IncludeVU']=='All' ):
			# water values are given once for the hydrosystem
			
			# find file
			ListWVfile=[elem for elem in list(HSSS.loc[hydrosystem]['WaterValues']) if len(elem)>0 ]					
			if len(ListWVfile)>0:
				WVfile=ListWVfile[0]
				logger.info('Add Bellman values from file')
				WVdata=read_input_timeseries(cfg,WVfile,'inputpath', index_col=0,skiprows=0,dayfirst=cfg['Calendar']['dayfirst'])
				# keep only data included in the period of the block
				WVdata= WVdata[WVdata.index < NumberSSVTimeSteps]

				# replace index column by start dates of SSV stages
				WVdata['new_index'] = datesSSV['start'].iloc[WVdata.index].values
				WVdata.set_index('new_index', inplace=True)
								
				WVdata.index=pd.to_datetime(WVdata.index,dayfirst=cfg['Calendar']['dayfirst'])
				# keep only data included in the period of the block
				WVdata=WVdata[ (WVdata.index >=start) & (WVdata.index <=end)   ]
				# keep only the last date within this pediod
				WVdata=WVdata[ WVdata.index == WVdata.index.max() ]
				WVsize=len(WVdata.index)
				
				if WVsize>0: 
					NumRowPolyFunction=WVsize
				else:
					NumRowPolyFunction=1
				PolyhedralFunctionBlock.createDimension("PolyFunction_NumRow", NumRowPolyFunction)
				PolyFunction_A=PolyhedralFunctionBlock.createVariable("PolyFunction_A",np.double,("PolyFunction_NumRow","PolyFunction_NumVar"))
				PolyFunction_b=PolyhedralFunctionBlock.createVariable("PolyFunction_b",np.double,("PolyFunction_NumRow"))
				
				if WVsize>0: 
					Data_A=np.array( WVdata[ [elem for elem in WVdata.columns if 'a_' in elem  ] ] )
					Data_B=np.array( WVdata[ [elem for elem in WVdata.columns if elem=='b' ] ] )
				elif ('LastStepPolyFunctionB' in cfg['ParametersFormat'] and 'LastStepPolyFunctionA' in cfg['ParametersFormat']):
					Data_A=np.full(shape=(1,NumberReservoirsInHydroSystem),fill_value=cfg['ParametersFormat']['LastStepPolyFunctionA'])
					Data_B=np.full(shape=1,fill_value =cfg['ParametersFormat']['LastStepPolyFunctionB'])
				else:
					Data_A=np.zeros(shape=(1,NumberReservoirsInHydroSystem))
					Data_B=np.zeros(shape=1)
					
				for row in range(NumRowPolyFunction):
					PolyFunction_A[row,:]=Data_A[row]
					PolyFunction_b[row]=Data_B[row]
		
		# case without a cuts file
		elif (cfg['FormatVU']=='Polyhedral') and ('WaterValues' not in SS.columns) and ( cfg['IncludeVU']=='Last' and id==NumberSSVTimeSteps-1 ):
			PolyhedralFunctionBlock.createDimension("PolyFunction_NumRow", 1)
			PolyFunction_A=PolyhedralFunctionBlock.createVariable("PolyFunction_A",np.double,("PolyFunction_NumRow","PolyFunction_NumVar"))
			PolyFunction_b=PolyhedralFunctionBlock.createVariable("PolyFunction_b",np.double,("PolyFunction_NumRow"))
			Data_A=np.full(shape=(1,NumberReservoirsInHydroSystem),fill_value=cfg['ParametersFormat']['LastStepPolyFunctionA'])
			Data_B=np.full(shape=1,fill_value =cfg['ParametersFormat']['LastStepPolyFunctionB'])
			PolyFunction_A[0,:]=Data_A[0]
			PolyFunction_b[0]=Data_B[0]

		indexUnitBlock=indexUnitBlock+1
	return indexUnitBlock
	
def addThermalUnitBlocks(Block,indexUnitBlock,scenario,start,end):
	for tu in TU.index:
		# create all thermal units
		for index_subunit in range(TU['NumberUnits'][tu]):
			TBlock=Block.createGroup('UnitBlock_'+str(indexUnitBlock))
			TBlock.type="ThermalUnitBlock"
			TBlock.setncattr("name",tu[0]+'_'+tu[1]+'_'+str(index_subunit))
			TBlock.createDimension("NumberIntervals",NumberIntervals)

			# add minpower and maxpower
			MaxPowerData=TU['MaxPower'][tu]
			if 'MaxPowerProfile' in TU.columns and len(TU['MaxPowerProfile'][tu])>0:
				MaxPowerProfile=TU['MaxPowerProfile'][tu]
				MaxPower=TBlock.createVariable("MaxPower",np.double,("NumberIntervals"))
				
				if '.csv' in MaxPowerProfile:
					if 'Thermal:MaxPowerProfile' in ScenarisedData:
						if cfg['IncludeScenarisedData']:
							MaxPower[:]=np.array(ThermalScenarios.loc[tu][scenario][ ( ThermalScenarios.loc[tu].index >= start ) & ( ThermalScenarios.loc[tu].index <= end ) ])
							pmax=ThermalScenarios.loc[tu][scenario][ ( ThermalScenarios.loc[tu].index >= start ) & ( ThermalScenarios.loc[tu].index <= end ) ]
						else:
							MaxPower[:]=np.array(DeterministicTimeSeries['Zero'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
							pmax=MaxPowerData*DeterministicTimeSeries['Zero'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]
					else:
						logger.warning('deterministic timeseries for maxpower not implemented')
						# not implemented
						# read det time serie
						# extend and resample
				else:
					MaxPower[:]=np.array(MaxPowerData*DeterministicTimeSeries[MaxPowerProfile][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
					pmax=MaxPowerData*DeterministicTimeSeries[MaxPowerProfile][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]
			else:
				pmax=MaxPowerData*DeterministicTimeSeries['One'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]
				MaxPowerProfile=''
				MaxPower=TBlock.createVariable("MaxPower",np.double,())
				MaxPower[:]=MaxPowerData*UCTimeStep
			if 'MinPower' in TU.columns:
				MinPowerData=TU['MinPower'][tu]
				if type(MinPowerData)==str :
					MinPower=TBlock.createVariable("MinPower",np.double,("NumberIntervals"))
					pmin=np.minimum(DeterministicTimeSeries[MinPowerData][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ],pmax)
					MinPower[:]=pmin
				else:
					if len(MaxPowerProfile)>0 and  ( (MinPowerData>pmax).isin([True]).sum()>0 ):
						pmin=np.minimum(MinPowerData*DeterministicTimeSeries['One'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ],pmax)
						MinPower=TBlock.createVariable("MinPower",np.double,("NumberIntervals"))
						MinPower[:]=pmin
					else:
						MinPower=TBlock.createVariable("MinPower",np.double,())
						MinPower[:]=MinPowerData*UCTimeStep

			# create cost variables
			if 'QuadTerm' in TU.columns:
				QuadTermData=TU['QuadTerm'][tu]
				if type(QuadTermData)==str:
					QuadTerm=TBlock.createVariable("QuadTerm",np.double,("NumberIntervals"))
					QuadTerm[:]=np.array(DeterministicTimeSeries[QuadTermData][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
				elif QuadTermData>0 or QuadTermData<0 :
					QuadTerm=TBlock.createVariable("QuadTerm",np.double,())
					QuadTerm[:]=[QuadTermData]
					
			if 'StartUpCost' in TU.columns:
				StartUpCostData=TU['StartUpCost'][tu]
				if type(StartUpCostData)==str:
					StartUpCost=TBlock.createVariable("StartUpCost",np.double,("NumberIntervals"))
					StartUpCost[:]=np.array(DeterministicTimeSeries[StartUpCostData][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
				elif StartUpCostData>0 or StartUpCostData<0:
					StartUpCost=TBlock.createVariable("StartUpCost",np.double,())
					StartUpCost[:]=[StartUpCostData]
					
			if 'VariableCost' in TU.columns:
				VariableCostData=TU['VariableCost'][tu]
				if type(VariableCostData)==str:
					LinearTerm=TBlock.createVariable("LinearTerm",np.double,("NumberIntervals"))
					LinearTerm[:]=np.array(DeterministicTimeSeries[VariableCostData][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
				else :
					LinearTerm=TBlock.createVariable("LinearTerm",np.double,())
					LinearTerm[:]=[VariableCostData]
			else:
				LinearTerm=TBlock.createVariable("LinearTerm",np.double,())
				LinearTerm[:]=[0.0]
					
			if 'FixedCost' in TU.columns:
				FixedCostData=TU['FixedCost'][tu]*int(UCTimeStep)/8760 # FixedCost is in €/MW/year
				if type(FixedCostData)==str:
					ConstTerm=TBlock.createVariable("ConstTerm",np.double,("NumberIntervals"))
					ConstTerm[:]=UCTimeStep*np.array(DeterministicTimeSeries[FixedCostData][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
				elif FixedCostData>0 or FixedCostData<0:
					ConstTerm=TBlock.createVariable("ConstTerm",np.double,())
					ConstTerm[:]=[FixedCostData]
						
			# create ramping constraints
			if 'DeltaRampUp' in TU.columns:
				DeltaRampUpData=TU['DeltaRampUp'][tu]
				if type(DeltaRampUpData)==str:
					DeltaRampUp=TBlock.createVariable("DeltaRampUp",np.double,("NumberIntervals"))
					DeltaRampUp[:]=np.array(DeterministicTimeSeries[DeltaRampUpData][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
				else:
					DeltaRampUp=HBlock.createVariable("DeltaRampUp",np.double,())
					DeltaRampUp[:]=[DeltaRampUpData*UCTimeStep]
			if 'DeltaRampDown' in TU.columns:
				DeltaRampDownData=TU['DeltaRampDown'][tu]
				if type(DeltaRampDownData)==str:
					DeltaRampDown=TBlock.createVariable("DeltaRampDown",np.double,("NumberIntervals"))
					DeltaRampDown[:]=np.array(DeterministicTimeSeries[DeltaRampDownData][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
				else:
					DeltaRampDown=TBlock.createVariable("DeltaRampDown",np.double,())
					DeltaRampDown[:]=[DeltaRampDownData*UCTimeStep]			
			
			# create primary and secondary rho
			if 'PrimaryRho' in TU.columns:
				PrimaryRhoData=TU['PrimaryRho'][tu]
				if type(PrimaryRhoData)==str:
					PrimaryRho=TBlock.createVariable("PrimaryRho",np.double,("NumberIntervals"))
					PrimaryRho[:]=np.array(DeterministicTimeSeries[PrimaryRhoData][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
				else:
					PrimaryRho=TBlock.createVariable("PrimaryRho",np.double,())
					PrimaryRho[:]=[PrimaryRhoData]

			if 'SecondaryRho' in TU.columns:
				SecondaryRhoData=TU['SecondaryRho'][tu]
				if type(SecondaryRhoData)==str:
					SecondaryRho=TBlock.createVariable("SecondaryRho",np.double,("NumberIntervals"))
					SecondaryRho[:]=np.array(DeterministicTimeSeries[SecondaryRhoData][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
				else:
					SecondaryRho=TBlock.createVariable("SecondaryRho",np.double,())
					SecondaryRho[:]=[SecondaryRhoData]

			# create duration constraints
			if 'MinUpTime' in TU.columns:
				MinUpTimeTimeData=TU['MinUpTime'][tu]
				if MinUpTimeTimeData>0:
					MinUpTime=TBlock.createVariable("MinUpTime","u4",())
					MinUpTime[:]=int(MinUpTimeData/UCTimeStep)
				
			if 'MinDownTime' in TU.columns:
				MinDownTimeTimeData=TU['MinDownTime'][tu]
				if MinDownTimeTimeData>0:
					MinDownTime=TBlock.createVariable("MinDownTime","u4",())
					MinDownTime[:]=int(MinDownTimeData/UCTimeStep)
			
			# create initial conditions
			if 'InitialPower' in TU.columns:
				InitialPowerData=TU['InitialPower'][tu]*UCTimeStep
			else: InitialPowerData=MaxPowerData*UCTimeStep
			InitialPower=TBlock.createVariable("InitialPower",np.double,())
			InitialPower[:]=InitialPowerData
			
			if 'InitUpDownTime' in TU.columns:
				InitUpDownTimeData=int(TU['InitUpDownTime'][tu]/UCTimeStep)
			elif 'MinUpTime' in TU.columns and 'MinDownTime' in TU.columns:
				InitUpDownTimeData=max(int(MinUpTimeTimeData/UCTimeStep),int(MinDownTimeTimeData/UCTimeStep))
			elif 'MinUpTime' in TU.columns and 'MinDownTime' not in TU.columns:
				InitUpDownTimeData=int(MinUpTimeTimeData/UCTimeStep)
			elif 'MinUpTime' not in TU.columns and 'MinDownTime' in TU.columns:
				InitUpDownTimeData=int(MinDownTimeTimeData/UCTimeStep)
			else: InitUpDownTimeData=1
			InitUpDownTime=TBlock.createVariable("InitUpDownTime","u4",())
			InitUpDownTime[:]=InitUpDownTimeData
				
			# create consumption when started but producing 0
			if 'Pauxiliary' in TU.columns:
				FixedConsumptionData=TU['Pauxiliary'][tu]*UCTimeStep
				if type(FixedConsumptionData)==str:
					FixedConsumption=TBlock.createVariable("FixedConsumption",np.double,("NumberIntervals"))
					FixedConsumption[:]=UCTimeStep*np.array(DeterministicTimeSeries[FixedConsumptionData][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
				else:
					FixedConsumption=TBlock.createVariable("FixedConsumption",np.double,())
					FixedConsumption[:]=[FixedConsumptionData]
					
			# create consumption when started but producing 0
			if 'Inertia' in TU.columns:
				InertiaData=TU['Inertia'][tu]
				if type(InertiaData)==str:
					InertiaCommitment=TBlock.createVariable("InertiaCommitment",np.double,("NumberIntervals"))
					InertiaCommitment[:]=np.array(cfg['ParametersFormat']['InertiaMultFactor']*DeterministicTimeSeries[InertiaData][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
				else:
					InertiaCommitment=TBlock.createVariable("InertiaCommitment",np.double,())
					InertiaCommitment[:]=[cfg['ParametersFormat']['InertiaMultFactor']*InertiaData]
		
			indexUnitBlock=indexUnitBlock+1
	return indexUnitBlock
			
def addIntermittentUnitBlocks(Block,indexUnitBlock,scenario,start,end):
	for tu in RES.index:
		# create all thermal units
		for index_subunit in range(RES['NumberUnits'][tu]):
			IBlock=Block.createGroup('UnitBlock_'+str(indexUnitBlock))
			IBlock.type="IntermittentUnitBlock"
			IBlock.setncattr("name",tu[0]+'_'+tu[1]+'_'+str(index_subunit))
			IBlock.createDimension("NumberIntervals",NumberIntervals)

			# add minpower and maxpower
			MaxPowerData=RES['MaxPower'][tu]
			if 'MaxPowerProfile' in RES.columns and len(RES['MaxPowerProfile'][tu])>0:
				MaxPowerProfile=RES['MaxPowerProfile'][tu]
				MaxPower=IBlock.createVariable("MaxPower",np.double,("NumberIntervals"))
				
				if '.csv' in MaxPowerProfile:
					if 'Renewable:MaxPowerProfile' in ScenarisedData:
						if cfg['IncludeScenarisedData']:
							MaxPower[:]=np.array(RESScenarios.loc[tu][scenario][ ( RESScenarios.loc[tu].index >= start ) & ( RESScenarios.loc[tu].index <= end ) ])
							pmax=RESScenarios.loc[tu][scenario][ ( RESScenarios.loc[tu].index >= start ) & ( RESScenarios.loc[tu].index <= end ) ]
					else:
						MaxPower[:]=np.array(DeterministicTimeSeries['Zero'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
						pmax=MaxPowerData*DeterministicTimeSeries['Zero'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]
				else:
					MaxPower[:]=np.array(MaxPowerData*DeterministicTimeSeries[MaxPowerProfile][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
					pmax=MaxPowerData*DeterministicTimeSeries[MaxPowerProfile][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]
			else:
				pmax=MaxPowerData*DeterministicTimeSeries['One'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]
				MaxPowerProfile=''
				MaxPower=IBlock.createVariable("MaxPower",np.double,())
				MaxPower[:]=[MaxPowerData*UCTimeStep]
			
			if 'MinPower' in RES.columns:
				MinPowerData=RES['MinPower'][tu]
				if type(MinPowerData)==str:
					MinPower=IBlock.createVariable("MinPower",np.double,("NumberIntervals"))
					pmin=np.minimum(DeterministicTimeSeries[MinPowerData][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ],pmax)
					MinPower[:]=pmin
				else:
					if len(MaxPowerProfile)>0 and  ( (MinPowerData>pmax).isin([True]).sum()>0 ):
						pmin=np.minimum(MinPowerData*DeterministicTimeSeries['One'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ],pmax)
						MinPower=IBlock.createVariable("MinPower",np.double,("NumberIntervals"))
						MinPower[:]=pmin
					else:
						MinPower=IBlock.createVariable("MinPower",np.double,())
						MinPower[:]=MinPowerData*UCTimeStep
						
			# create gamma and kappa
			if 'Gamma' in RES.columns:
				GammaData=RES['Gamma'][tu]
			else:
				GammaData=0.0
			Gamma=IBlock.createVariable("Gamma",np.double,())
			Gamma[:]=[GammaData]
			
			if 'Kappa' in RES.columns:
				KappaData=RES['Kappa'][tu]
				Kappa=IBlock.createVariable("Kappa",np.double,())
				Kappa[:]=[KappaData]

			# create inertia
			if 'Inertia' in RES.columns:
				InertiaData=RES['Inertia'][tu]
				if type(InertiaData)==str:
					InertiaPower=IBlock.createVariable("InertiaPower",np.double,("NumberIntervals"))
					InertiaPower[:]=np.array(cfg['ParametersFormat']['InertiaMultFactor']*DeterministicTimeSeries[InertiaData][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
				else:
					InertiaPower=IBlock.createVariable("InertiaPower",np.double,())
					InertiaPower[:]=[cfg['ParametersFormat']['InertiaMultFactor']*InertiaData]
		
			indexUnitBlock=indexUnitBlock+1
	return indexUnitBlock

def addSynchConsUnitBlocks(Block,indexUnitBlock,start,end):
	for tu in SYN.index:
		# create all thermal units
		for index_subunit in range(SYN['NumberUnits'][tu]):
			SBlock=Block.createGroup('UnitBlock_'+str(indexUnitBlock))
			SBlock.type="ThermalUnitBlock"
			SBlock.setncattr("name",tu[0]+'_'+tu[1]+'_'+str(index_subunit))
			SBlock.createDimension("NumberIntervals",NumberIntervals)

			# add minpower and maxpower
			MaxConsoData=SYN['MaxRotatingConsumption'][tu]
			if 'StartUpCost' in SYN.columns:
				StartUpCostData=SYN['StartUpCost'][tu]
			else:
				StartUpCostData=0
			if 'FixedCost' in SYN.columns:
				FixedCostData=SYN['FixedCost'][tu]
			else:
				FixedCostData=0
			MaxConsoData=SYN['MaxRotatingConsumption'][tu]
			InertiaData=SYN['Inertia'][tu]
			MaxPower=SBlock.createVariable("MaxPower",np.double,())
			MinPower=SBlock.createVariable("MaxPower",np.double,())
			FixedConsumption=SBlock.createVariable("FixedConsumption",np.double,())
			StartUpCost=SBlock.createVariable("StartUpCost",np.double,())
			FixedCost=SBlock.createVariable("FixedCost",np.double,())
			InertiaCommitment.createVariable("InertiaCommitment",np.double,())
			MaxConso[:]=[MaxConsoData]
			Inertia[:]=[cfg['ParametersFormat']['InertiaMultFactor']*InertiaData]
			MaxPower[:]=[0.0]
			MinPower[:]=[0.0]
			FixedConsumption[:]=[MaxConsoData]
			StartUpCost[:]=[StartUpCostData]
			FixedCost[:]=[FixedCostData]

			indexUnitBlock=indexUnitBlock+1
	return indexUnitBlock
			
def addBatteryUnitBlocks(Block,indexUnitBlock,start,end):
	for tu in STS.index:
		# create all thermal units
		for index_subunit in range(STS['NumberUnits'][tu]):
			TBlock=Block.createGroup('UnitBlock_'+str(indexUnitBlock))
			TBlock.type="BatteryUnitBlock"
			TBlock.setncattr("name",tu[0]+'_'+tu[1]+'_'+str(index_subunit))
			TBlock.createDimension("NumberIntervals",NumberIntervals)

			nameAsset=tu[0]
			for nametechno in cfg['technos']['hydrostorage']+cfg['technos']['battery']:
				if nametechno in nameAsset:
					technoAsset=nametechno

			# add minpower and maxpower
			MaxPowerData=STS['MaxPower'][tu]
			if 'MaxPowerProfile' in STS.columns and len(STS['MaxPowerProfile'][tu])>0:
				MaxPowerProfile=STS['MaxPowerProfile'][tu]
				MaxPower=TBlock.createVariable("MaxPower",np.double,("NumberIntervals"))
				MaxPower[:]=np.array(MaxPowerData*DeterministicTimeSeries[MaxPowerProfile][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
				pmax=MaxPowerData*DeterministicTimeSeries[MaxPowerProfile][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]
			else:
				pmax=MaxPowerData*DeterministicTimeSeries['One'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]
				MaxPowerProfile=''
				MaxPower=TBlock.createVariable("MaxPower",np.double,())
				MaxPower[:]=[MaxPowerData*UCTimeStep]
			
			if 'MinPower' in STS.columns:
				MinPowerData=STS['MinPower'][tu]
				if type(MinPowerData)==str:
					MinPower=TBlock.createVariable("MinPower",np.double,("NumberIntervals"))
					pmin=np.minimum(DeterministicTimeSeries[MinPowerData][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ],pmax)
					MinPower[:]=pmin
				else:
					if len(MaxPowerProfile)>0 and  ( (MinPowerData*UCTimeStep>pmax).isin([True]).sum()>0 ):
						pmin=np.minimum(MinPowerData*DeterministicTimeSeries['One'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ],pmax)
						MinPower=TBlock.createVariable("MinPower",np.double,("NumberIntervals"))
						MinPower[:]=pmin
					else:
						MinPower=TBlock.createVariable("MinPower",np.double,())
						MinPower[:]=[MinPowerData*UCTimeStep]
						
			# add minstorage and maxstorage
			MaxStorageData=STS['MaxVolume'][tu]
			if 'MaxStorageProfile' in STS.columns and len(STS['MaxStorageProfile'][tu])>0:
				MaxStorageProfile=STS['MaxStorageProfile'][tu]
				MaxStorage=TBlock.createVariable("MaxStorage",np.double,("NumberIntervals"))
				MaxStorage[:]=np.array(MaxStorageData*DeterministicTimeSeries[MaxStorageProfile][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
				vmax=MaxStorageData*DeterministicTimeSeries[MaxStorageProfile][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]
			else:
				vmax=MaxStorageData*DeterministicTimeSeries['One'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]
				MaxStorageProfile=''
				MaxStorage=TBlock.createVariable("MaxStorage",np.double,())
				MaxStorage[:]=MaxStorageData
			
			if 'MinVolume' in STS.columns:
				MinStorageData=STS['MinVolume'][tu]
				if type(MinStorageData)==str:
					MinStorage=TBlock.createVariable("MinStorage",np.double,("NumberIntervals"))
					vmin=np.minimum(DeterministicTimeSeries[MinStorageData][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ],vmax)
					MinStorage[:]=vmin
				else:
					if 'VolumeLevelTarget' in STS.columns:
						vmin=np.minimum(MinStorageData*DeterministicTimeSeries['One'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ],vmax)
						vmin.loc[vmin.tail(1).index.item()]=STS['VolumeLevelTarget'][tu]
						vmin.loc[vmin.head(1).index.item()]=STS['VolumeLevelTarget'][tu]
						MinStorage=TBlock.createVariable("MinStorage",np.double,("NumberIntervals"))
						MinStorage[:]=vmin
					elif len(MaxStorageProfile)>0:
						vmin=np.minimum(MinStorageData*DeterministicTimeSeries['One'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ],vmax)
						MinStorage=TBlock.createVariable("MinStorage",np.double,("NumberIntervals"))
						MinStorage[:]=vmin
					else:
						MinStorage=TBlock.createVariable("MinStorage",np.double,())
						MinStorage[:]=MinStorageData
						
			# create Initial Power
			if 'InitialPower' in STS.columns:
				InitialPowerData=STS['InitialPower'][tu]*UCTimeStep
			else:
				InitialPowerData=0.0
			InitialPower=TBlock.createVariable("InitialPower",np.double,())
			InitialPower[:]=[InitialPowerData]
			
			# create Initial Storage
			if 'InitialStorage' in STS.columns:
				InitialStorageData=STS['InitialStorage'][tu]
			elif 'VolumeLevelTarget' in STS.columns:
				InitialStorageData=STS['VolumeLevelTarget'][tu]
			else:
				InitialStorageData=0.0
			InitialStorage=TBlock.createVariable("InitialStorage",np.double,())
			InitialStorage[:]=[InitialStorageData]
			
			# create max primary and secondary power
			if 'MaxPrimaryPower' in STS.columns:
				MaxPrimaryPowerData=STS['MaxPrimaryPower'][tu]*UCTimeStep
				if type(MaxPrimaryPowerData)==str:
					MaxPrimaryPower=TBlock.createVariable("MaxPrimaryPower",np.double,("NumberIntervals"))
					MaxPrimaryPower[:]=np.array(DeterministicTimeSeries[MaxPrimaryPowerData][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
				else:
					MaxPrimaryPower=TBlock.createVariable("MaxPrimaryPower",np.double,())
					MaxPrimaryPower[:]=MaxPrimaryPowerData
					
			if 'MaxSecondaryPower' in STS.columns:
				MaxSecondaryPowerData=STS['MaxSecondaryPower'][tu]*UCTimeStep
				if type(MaxSecondaryPowerData)==str:
					MaxSecondaryPower=TBlock.createVariable("MaxSecondaryPower",np.double,("NumberIntervals"))
					MaxSecondaryPower[:]=np.array(DeterministicTimeSeries[MaxSecondaryPowerData][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
				else:
					MaxSecondaryPower=TBlock.createVariable("MaxSecondaryPower",np.double,())
					MaxSecondaryPower[:]=MaxSecondaryPowerData
						
			# create ramping constraints
			if 'DeltaRampUp' in STS.columns:
				DeltaRampUpData=STS['DeltaRampUp'][tu]*UCTimeStep
				if type(DeltaRampUpData)==str:
					DeltaRampUp=TBlock.createVariable("DeltaRampUp",np.double,("NumberIntervals"))
					DeltaRampUp[:]=np.array(DeterministicTimeSeries[DeltaRampUpData][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
				else:
					DeltaRampUp=HBlock.createVariable("DeltaRampUp",np.double,())
					DeltaRampUp[:]=[DeltaRampUpData]
			if 'DeltaRampDown' in STS.columns:
				DeltaRampDownData=STS['DeltaRampUp'][tu]*UCTimeStep
				if type(DeltaRampDownData)==str:
					DeltaRampDown=TBlock.createVariable("DeltaRampDown",np.double,("NumberIntervals"))
					DeltaRampDown[:]=np.array(DeterministicTimeSeries[DeltaRampDownData][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
				else:
					DeltaRampDown=TBlock.createVariable("DeltaRampDown",np.double,())
					DeltaRampDown[:]=[DeltaRampDownData]		
					
			# create rho
			if 'PumpingEfficiency' in STS.columns:
				PumpingEfficiency=STS['PumpingEfficiency'][tu]
				if type(PumpingEfficiency)==str:
					StoringBatteryRho=TBlock.createVariable("StoringBatteryRho",np.double,("NumberIntervals"))
					StoringBatteryRho[:]=np.array(DeterministicTimeSeries[PumpingEfficiency][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
				else:
					StoringBatteryRho=TBlock.createVariable("StoringBatteryRho",np.double,())
					StoringBatteryRho[:]=[PumpingEfficiency]
			else:
				StoringBatteryRho=TBlock.createVariable("StoringBatteryRho",np.double,())
				StoringBatteryRho[:]=[ cfg['PumpingEfficiency'][technoAsset][nameAsset] ]
				
			if 'TurbineEfficiency' in STS.columns:
				TurbineEfficiency=STS['TurbineEfficiency'][tu]
				if type(TurbineEfficiency)==str:
					ExtractingBatteryRho=TBlock.createVariable("ExtractingBatteryRho",np.double,("NumberIntervals"))
					ExtractingBatteryRho[:]=np.array(1/DeterministicTimeSeries[TurbineEfficiency][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
				else:
					ExtractingBatteryRho=TBlock.createVariable("ExtractingBatteryRho",np.double,())
					if TurbineEfficiency!=0:
						ExtractingBatteryRho[:]=[1/TurbineEfficiency]
					else:
						ExtractingBatteryRho[:]=[1]
			else:
				ExtractingBatteryRho=TBlock.createVariable("ExtractingBatteryRho",np.double,())
				ExtractingBatteryRho[:]=[1]
			
			# create cost
			if 'Cost' in STS.columns:
				CostData=STS['Cost'][tu]
				if type(Cost)==str:
					Cost=TBlock.createVariable("Cost",np.double,("NumberIntervals"))
					Cost[:]=np.array(DeterministicTimeSeries[CostData][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
				else:
					Cost=TBlock.createVariable("Cost",np.double,())
					Cost[:]=[CostData]
				
			# create cost
			if 'Inflows' in STS.columns:
				Inflow=STS['Inflows'][tu]
				if type(Inflow)==str:
					Demand=TBlock.createVariable("Demand",np.double,("NumberIntervals"))
					Demand[:]=np.array((-1)*DeterministicTimeSeries[Inflow][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
				else:
					Demand=TBlock.createVariable("Demand",np.double,("NumberIntervals"))					
					Demand[:]=np.array((-1)*Inflow*DeterministicTimeSeries['One'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])

			if 'Kappa' in STS.columns:
				KappaData=STS['Kappa'][tu]
				Kappa=TBlock.createVariable("Kappa",np.double,())
				Kappa[:]=[KappaData]

			# create inertia
			if 'Inertia' in STS.columns:
				InertiaData=STS['Inertia'][tu]
				if type(InertiaData)==str:
					InertiaPower=TBlock.createVariable("InertiaPower",np.double,("NumberIntervals"))
					InertiaPower[:]=np.array(cfg['ParametersFormat']['InertiaMultFactor']*DeterministicTimeSeries[InertiaData][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
				else:
					InertiaPower=TBlock.createVariable("InertiaPower",np.double,())
					InertiaPower[:]=[cfg['ParametersFormat']['InertiaMultFactor']*InertiaData]
		
			indexUnitBlock=indexUnitBlock+1
	return indexUnitBlock

def addSlackUnitBlocks(Block,indexUnitBlock,start,end):
	for node in Nodes:
		# create all thermal units
		SBlock=Block.createGroup('UnitBlock_'+str(indexUnitBlock))
		SBlock.type="SlackUnitBlock"
		SBlock.setncattr("name",'SlackUnit_'+str(node))
		SBlock.createDimension("NumberIntervals",NumberIntervals)
		
		# add  maxpower for demand constraint
		if ('MaxActivePowerDemand',node) in ZV.index:
			MaxPowerData=ZV.loc['MaxActivePowerDemand',node]['value']
			MaxPowerProfile=ZV.loc['MaxActivePowerDemand',node]['Profile_Timeserie']
		elif 'MaxPower' in cfg['CouplingConstraints']['ActivePowerDemand']:
			MaxPowerData=cfg['CouplingConstraints']['ActivePowerDemand']['MaxPower']
			MaxPowerProfile=''
		else:
			MaxPowerData=100000
			MaxPowerProfile=''
		if len(MaxPowerProfile)>1:
			MaxPower=SBlock.createVariable("MaxPower",np.double,("NumberIntervals"))
			MaxPower[:]=np.array(MaxPowerData*DeterministicTimeSeries[MaxPowerProfile][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
		else:
			MaxPower=SBlock.createVariable("MaxPower",np.double,())
			MaxPower[:]=MaxPowerData*UCTimeStep
		
		# add cost for demand constraint
		if ('CostActivePowerDemand',node) in ZV.index:
			CostPowerData=ZV.loc['CostActivePowerDemand',node]['value']
			CostPowerProfile=ZV.loc['CostActivePowerDemand',node]['Profile_Timeserie']
		elif 'Cost' in cfg['CouplingConstraints']['ActivePowerDemand']:
			CostPowerData=cfg['CouplingConstraints']['ActivePowerDemand']['Cost']
			CostPowerProfile=''
		else:
			CostPowerData=0.0
			CostPowerProfile=''
		if len(CostPowerProfile)>1:
			ActivePowerCost=SBlock.createVariable("ActivePowerCost",np.double,("NumberIntervals"))
			ActivePowerCost[:]=np.array(CostPowerData*DeterministicTimeSeries[CostPowerProfile][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
		else:
			ActivePowerCost=SBlock.createVariable("ActivePowerCost",np.double,())
			ActivePowerCost[:]=CostPowerData
		
		isPrimary=('PrimaryDemand' in ListCouplingConstraints)
		if isPrimary:
			if ('MaxPrimaryDemand',node) in ZV.index:
				MaxPrimaryPowerData=ZV.loc['MaxPrimaryDemand',node]['value']
				MaxPrimaryPowerProfile=ZV.loc['MaxPrimaryDemand',node]['Profile_Timeserie']
			elif 'PrimaryDemand' in cfg['CouplingConstraints']:
				if 'MaxPower' in cfg['CouplingConstraints']['PrimaryDemand']:
					MaxPrimaryPowerData=cfg['CouplingConstraints']['PrimaryDemand']['MaxPower']
					MaxPrimaryPowerProfile=''
				else:
					MaxPrimaryPowerData=10000
					MaxPrimaryPowerProfile=''
				
			if len(MaxPrimaryPowerProfile)>1:
				MaxPrimaryPower=SBlock.createVariable("MaxPrimaryPower",np.double,("NumberIntervals"))
				MaxPrimaryPower[:]=np.array(MaxPrimaryPowerData*DeterministicTimeSeries[MaxPrimaryPowerProfile][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
			else:
				MaxPrimaryPower=SBlock.createVariable("MaxPrimaryPower",np.double,())
				MaxPrimaryPower[:]=MaxPrimaryPowerData*UCTimeStep
			if ('CostPrimaryDemand',node) in ZV.index:
				CostPrimaryPowerData=ZV.loc['CostPrimaryDemand',node]['value']
				CostPrimaryPowerProfile=ZV.loc['CostPrimaryDemand',node]['Profile_Timeserie']
			elif 'PrimaryDemand' in cfg['CouplingConstraints']:
				if 'Cost' in cfg['CouplingConstraints']['PrimaryDemand']:
					PrimaryCost=cfg['CouplingConstraints']['PrimaryDemand']['Cost']
					MaxPrimaryPowerProfile=''
				else:
					MaxPrimaryPowerData=10000
					MaxPrimaryPowerProfile=''
			if len(CostPrimaryPowerProfile)>1:
				PrimaryCost=SBlock.createVariable("PrimaryCost",np.double,("NumberIntervals"))
				PrimaryCost[:]=np.array(CostPrimaryPowerData*DeterministicTimeSeries[CostPrimaryPowerProfile][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
			else:
				PrimaryCost=SBlock.createVariable("PrimaryCost",np.double,())
				PrimaryCost[:]=CostPrimaryPowerData
				
		isSecondary=('SecondaryDemand' in ListCouplingConstraints)
		if isSecondary:
			if ('MaxSecondaryDemand',node) in ZV.index:
				MaxSecondaryPowerData=ZV.loc['MaxSecondaryDemand',node]['value']
				MaxSecondaryPowerProfile=ZV.loc['MaxSecondaryDemand',node]['Profile_Timeserie']
			elif 'SecondaryDemand' in cfg['CouplingConstraints']:
				if 'MaxPower' in cfg['CouplingConstraints']['SecondaryDemand']:
					MaxSecondaryPowerData=cfg['CouplingConstraints']['SecondaryDemand']['MaxPower']
					MaxSecondaryPowerProfile=''
				else:
					MaxSecondaryPowerData=10000
					MaxSecondaryPowerProfile=''
				
			if len(MaxSecondaryPowerProfile)>1:
				MaxSecondaryPower=SBlock.createVariable("MaxSecondaryPower",np.double,("NumberIntervals"))
				MaxSecondaryPower[:]=np.array(MaxSecondaryPowerData*DeterministicTimeSeries[MaxSecondaryPowerProfile][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
			else:
				MaxSecondaryPower=SBlock.createVariable("MaxSecondaryPower",np.double,())
				MaxSecondaryPower[:]=MaxSecondaryPowerData*UCTimeStep
			if ('CostSecondaryDemand',node) in ZV.index:
				CostSecondaryPowerData=ZV.loc['CostSecondaryDemand',node]['value']
				CostSecondaryPowerProfile=ZV.loc['CostSecondaryDemand',node]['Profile_Timeserie']
			elif 'SecondaryDemand' in cfg['CouplingConstraints']:
				if 'Cost' in cfg['CouplingConstraints']['SecondaryDemand']:
					SecondaryCost=cfg['CouplingConstraints']['SecondaryDemand']['Cost']
					MaxSecondaryPowerProfile=''
				else:
					MaxSecondaryPowerData=10000
					MaxSecondaryPowerProfile=''
			if len(CostSecondaryPowerProfile)>1:
				SecondaryCost=SBlock.createVariable("SecondaryCost",np.double,("NumberIntervals"))
				SecondaryCost[:]=np.array(CostSecondaryPowerData*DeterministicTimeSeries[CostSecondaryPowerProfile][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
			else:
				SecondaryCost=SBlock.createVariable("SecondaryCost",np.double,())
				SecondaryCost[:]=CostSecondaryPowerData
				
		isInertia=('InertiaDemand' in ListCouplingConstraints)
		if isInertia:
			if ('InertiaDemand',node) in ZV.index:
				MaxInertiaData=ZV.loc['MaxInertiaDemand',node]['value']
				MaxInertiaProfile=ZV.loc['MaxInertiaDemand',node]['Profile_Timeserie']
			elif 'InertiaDemand' in cfg['CouplingConstraints']:
				if 'MaxPower' in cfg['CouplingConstraints']['InertiaDemand']:
					MaxInertiaData=cfg['CouplingConstraints']['InertiaDemand']['MaxPower']
					MaxInertiaProfile=''
				else:
					MaxInertiaData=10000
					MaxInertiaProfile=''
				
			if len(MaxInertiaProfile)>1:
				MaxInertia=SBlock.createVariable("MaxInertia",np.double,("NumberIntervals"))
				MaxInertia[:]=np.array(MaxInertiaData*DeterministicTimeSeries[MaxInertiaProfile][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
			else:
				MaxInertia=SBlock.createVariable("MaxInertia",np.double,())
				MaxInertia[:]=MaxInertiaData
			if ('CostInertiaDemand',node) in ZV.index:
				CostInertiaData=ZV.loc['CostInertiaDemand',node]['value']
				CostInertiaProfile=ZV.loc['CostInertiaDemand',node]['Profile_Timeserie']
			elif 'InertiaDemand' in cfg['CouplingConstraints']:
				if 'Cost' in cfg['CouplingConstraints']['InertiaDemand']:
					InertiaCost=cfg['CouplingConstraints']['InertiaDemand']['Cost']
					MaxInertiaProfile=''
				else:
					MaxInertiaData=10000
					MaxInertiaProfile=''
			if len(CostInertiaProfile)>1:
				InertiaCost=SBlock.createVariable("InertiaCost",np.double,("NumberIntervals"))
				InertiaCost[:]=np.array(CostInertiaData*DeterministicTimeSeries[CostInertiaProfile][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ])
			else:
				InertiaCost=SBlock.createVariable("InertiaCost",np.double,())
				InertiaCost[:]=CostInertiaData
			
		indexUnitBlock=indexUnitBlock+1
	return indexUnitBlock


def createUCBlock(filename,id,scenario,start,end):
	UCBlock=Dataset(filename,'w',format='NETCDF4')
		
	# general attribute of the file
	UCBlock.setncattr('SMS++_file_type',np.int64(1))
	
	# add the unique group Block of type UCBlock 
	Block=UCBlock.createGroup("Block_0")
	Block.id=str(id)
	Block.type="UCBlock"
	
	# add dimensions
	Block.createDimension("TimeHorizon",TimeHorizonUC)
	Block.createDimension("NumberUnits",NumberUnits)
	Block.createDimension("NumberElectricalGenerators",NumberElectricalGenerators)
	Block.createDimension("NumberNodes",NumberNodes)
	if NumberLines>0:
		Block.createDimension("NumberLines",NumberLines)
	if NumberPrimaryZones>0:
		Block.createDimension("NumberPrimaryZones",NumberPrimaryZones)
	if NumberSecondaryZones >0:
		Block.createDimension("NumberSecondaryZones",NumberSecondaryZones)
	if NumberInertiaZones >0:
		Block.createDimension("NumberInertiaZones",NumberInertiaZones)
	if NumberPollutants >0:
		Block.createDimension("NumberPollutants",NumberPollutants)
	
	# create generatornode
	if NumberNodes>1:
		GeneratorNode=Block.createVariable("GeneratorNode","u4" ,("NumberElectricalGenerators"))
		indexGen=0
		GeneratorNodeData=np.array([0]*NumberElectricalGenerators)
		for data in listData:
			if data[1]=='HSSS':
				for unit in data[0].index:
					NbUnits=data[0]['NumberUnits'][unit]
					unitnode=unit[1]
					Node_index=Nodes[ Nodes==unitnode ].index[0]
					for i in range(NbUnits):
						# case with pumping: 3 generators
						if ('MinPower' in data[0].columns)  and (data[0]['MinPower'][unit]<0): 
							GeneratorNodeData[indexGen]=int(Node_index)
							GeneratorNodeData[indexGen+1]=int(Node_index)
							GeneratorNodeData[indexGen+2]=int(Node_index)
							indexGen=indexGen+2
						else:
							GeneratorNodeData[indexGen]=int(Node_index)
							GeneratorNodeData[indexGen+1]=int(Node_index)
							indexGen=indexGen+1
			else:
				for unit in data[0].index:
					NbUnits=data[0]['NumberUnits'][unit]
					unitnode=unit[1]
					Node_index=Nodes[ Nodes==unitnode ].index[0]
					for i in range(NbUnits):
						GeneratorNodeData[indexGen]=int(Node_index)
						indexGen=indexGen+1
		# include slack units
		for node in range(NumberNodes):
			GeneratorNodeData[indexGen]=node
			indexGen=indexGen+1
		GeneratorNode[:]=GeneratorNodeData
	
	# create nodename
	NodeName=Block.createVariable("NodeName",str,("NumberNodes",))
	NodeName[:]=np.array(Nodes)
	
	# create and fill variables
	ActivePowerDemand=Block.createVariable("ActivePowerDemand",np.double ,("NumberNodes","TimeHorizon"))
	if cfg['IncludeScenarisedData'] and 'ActivePowerDemand' in ScenarisedData:
		for node in range(NumberNodes):
			ActivePowerDemand[node,:]=np.array(DemandScenarios.loc[Nodes[node]][scenario][ ( DemandScenarios.loc[Nodes[node]].index >= start ) & ( DemandScenarios.loc[Nodes[node]].index <= end ) ])
	else:
		for node in range(NumberNodes):
			ActivePowerDemand[node,:]=np.zeros(int(TimeHorizonUC))
	
	# create other coupling constraints
	if NumberPrimaryZones>0:
		PrimaryZones=Block.createVariable("PrimaryZones",str ,("NumberNodes"))
		PrimaryZones[:]=np.array(list(Partition[Coupling.loc['PrimaryDemand']['Partition']]))
		PrimaryDemand=Block.createVariable("PrimaryDemand",np.double ,("NumberPrimaryZones","TimeHorizon"))
		parts=Coupling.loc['PrimaryDemand']['Sum']
		indexZone=0
		for zone in list(Partition[Coupling.loc['PrimaryDemand']['Partition']]):
			firstpart=True
			for part in parts:
				value=ZV.loc[part,zone]['value']
				profile=ZV.loc[part,zone]['Profile_Timeserie']
				if len(profile)>0:
					part_constraint=value*DeterministicTimeSeries[profile][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]
				else:
					part_constraint=value*DeterministicTimeSeries['One'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]
				if firstpart:
					constraint=part_constraint
				else:
					constraint=constraint+part_constraint
				firstpart=False
			PrimaryDemand[indexZone,:]=constraint
			indexZone=indexZone+1
	if NumberSecondaryZones>0:
		SecondaryZones=Block.createVariable("SecondaryZones",str ,("NumberNodes"))
		SecondaryZones[:]=np.array(list(Partition[Coupling.loc['SecondaryZones']['Partition']]))
		SecondaryDemand=Block.createVariable("SecondaryDemand",np.double,("NumberSecondaryZones","TimeHorizon"))
		parts=Coupling.loc['SecondaryDemand']['Sum']
		indexZone=0
		for zone in list(Partition[Coupling.loc['SecondaryDemand']['Partition']]):
			firstpart=True
			for part in parts:
				value=ZV.loc[part,zone]['value']
				profile=ZV.loc[part,zone]['Profile_Timeserie']
				if len(profile)>0:
					part_constraint=value*DeterministicTimeSeries[profile][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]
				else:
					part_constraint=value*DeterministicTimeSeries['One'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]
				if firstpart:
					constraint=part_constraint
				else:
					constraint=constraint+part_constraint
				firstpart=False
			SecondaryDemand[indexZone,:]=constraint
			indexZone=indexZone+1
	if NumberInertiaZones>0:
		InertiaZones=Block.createVariable("InertiaZones",str ,("NumberNodes"))
		InertiaZones[:]=np.array(list(Partition[Coupling.loc['InertiaZones']['Partition']]))
		InertiaDemand=Block.createVariable("InertiaDemand",np.double ,("NumberInertiaZones","TimeHorizon"))
		parts=Coupling.loc['InertiaDemand']['Sum']
		indexZone=0
		for zone in list(Partition[Coupling.loc['InertiaDemand']['Partition']]):
			firstpart=True
			for part in parts:
				value=ZV.loc[part,zone]['value']
				profile=ZV.loc[part,zone]['Profile_Timeserie']
				if len(profile)>0:
					part_constraint=value*DeterministicTimeSeries[profile][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]
				else:
					part_constraint=value*DeterministicTimeSeries['One'][ ( DeterministicTimeSeries.index >= start ) & ( DeterministicTimeSeries.index <= end ) ]
				if firstpart:
					constraint=part_constraint
				else:
					constraint=constraint+part_constraint
				firstpart=False
			InertiaDemand[indexZone,:]=constraint
			indexZone=indexZone+1
	if NumberPollutants>0:
		NumberPollutantZones=Block.createVariable("NumberPollutantZones","u4" ,("NumberPollutants"))
		NumberPollutantZonesData=np.array([0]*NumberPollutants)
		PollutantBudgetData=np.array([])
		PollutantRhoData=np.array([])
		PollutantZonesData=np.array(['']*NumberNodes)
		PollutantZones=Block.createVariable("PollutantZones",str ,("NumberPollutants","NumberNodes"))
		Block.createDimension("TotalNumberPollutantZones",TotalNumberPollutantZones)
		PollutantBudget=Block.createVariable("PollutantBudget",np.double ,("TotalNumberPollutantZones"))
		PollutantRho=Block.createVariable("PollutantRho",np.double ,("NumberPollutants","NumberElectricalGenerators"))
		indexPollutant=0
		indexP=0
		for p in ListPollutants:
			NumberPollutantZonesData[indexP]=len(list(Partition[Coupling.loc[p]['Partition']]))
			PollutantZonesData=np.array(list(Partition[Coupling.loc[p]['Partition']]))
			PollutantZones[indexPollutant,:]=PollutantZonesData
			indexPollutant=indexPollutant+1
			for z in list(Partition[Coupling.loc[p]['Partition']]):
				PollutantBudgetData=np.append(PollutantBudgetData,ZV.loc[p,z]['value'])
			
			for data in listData:
				for unit in data[0].index:
					NbUnits=data[0]['NumberUnits'][unit]
					if data[1]=='HSSS':
						for i in range(NbUnits):
							PollutantRhoData=np.append(PollutantRhoData,0.0)
							PollutantRhoData=np.append(PollutantRhoData,0.0)
					elif p in data[0].columns:
						pollutantrho=data[0][p][unit]
						for i in range(NbUnits):
							PollutantRhoData=np.append(PollutantRhoData,pollutantrho)
					else: 
						for i in range(NbUnits):
							PollutantRhoData=np.append(PollutantRhoData,0.0)
					
			# slack units
			for n in range(NumberNodes):PollutantRhoData=np.append(PollutantRhoData,0.0)
			indexP=indexP+1
		NumberPollutantZones[:]=NumberPollutantZonesData
		PollutantBudget[:]=PollutantBudgetData
		if NumberThermalUnits>0:
			PollutantRho[:]=PollutantRhoData
	
	# create grid	
	if NumberLines>0:
		StartLine=Block.createVariable("StartLine","u4",("NumberLines"))
		EndLine=Block.createVariable("EndLine","u4" ,("NumberLines"))
		StartData=np.array([])
		EndData=np.array([])
		MaxData=np.array(IN['MaxPowerFlow']*UCTimeStep)
		MinData=np.array(IN['MinPowerFlow']*UCTimeStep)
		NameData=np.array(IN.index)
		CostData=np.array([0.0]*NumberLines)
		SusceptanceData=np.array([0.0]*NumberLines)
		if 'Impedance' in IN.columns:
			IN['Susceptance']=IN['Impedance'].apply(lambda x: -1/x if x>0 else 0)
			SusceptanceData=np.array(IN['Susceptance'])
		if 'Cost' in IN.columns:
			CostData=np.array(IN['Cost'])
		for line in IN.index:
			StartData=np.append(StartData,Nodes[ Nodes==IN['StartLine'][line] ].index[0])
			EndData=np.append(EndData,Nodes[ Nodes==IN['EndLine'][line] ].index[0])
		StartLine[:]=StartData
		EndLine[:]=EndData
		
		MinPowerFlow=Block.createVariable("MinPowerFlow",np.double ,("NumberLines"))
		MinPowerFlow[:]=MinData
		MaxPowerFlow=Block.createVariable("MaxPowerFlow",np.double ,("NumberLines"))
		MaxPowerFlow[:]=MaxData
		Susceptance=Block.createVariable("Susceptance",np.double,("NumberLines"))
		Susceptance[:]=SusceptanceData
		NetworkCost=Block.createVariable("NetworkCost",np.double ,("NumberLines"))
		NetworkCost[:]=CostData
		LineName=Block.createVariable("LineName",str,("NumberLines"))
		LineName[:]=NameData
		
	indexUnitBlock=0
	if NumberHydroSystems>0:
		indexUnitBlock=addHydroUnitBlocks(Block,indexUnitBlock,scenario,start,end,id)
	if NumberThermalUnits>0:
		indexUnitBlock=addThermalUnitBlocks(Block,indexUnitBlock,scenario,start,end)
	if NumberIntermittentUnits>0:
		indexUnitBlock=addIntermittentUnitBlocks(Block,indexUnitBlock,scenario,start,end)
	if NumberSyncUnits>0:
		indexUnitBlock=addSynchCondUnitBlocks(Block,indexUnitBlock,start,end)
	if NumberBatteryUnits>0:
		indexUnitBlock=addBatteryUnitBlocks(Block,indexUnitBlock,start,end)
	if NumberSlackUnits>0:
		indexUnitBlock=addSlackUnitBlocks(Block,indexUnitBlock,start,end)

	UCBlock.close()


def createSDDPBlock(filename,id):
	SDDPBlock=Dataset(filename,'w',format='NETCDF4')
	# general attribute of the file
	SDDPBlock.setncattr('SMS++_file_type',1)
	Block=SDDPBlock.createGroup("Block_"+str(id))
	Block.type="SDDPBlock"
	Block.createDimension("NumPolyhedralFunctionsPerSubBlock",NumberHydroSystems)
	Block.createDimension("TimeHorizon",TimeHorizonSSV)
	Block.createDimension("NumberScenarios",len(ListScenarios))
	SubScenarioSize = ThermalMaxPowerSize * Nb_TPP + (SSVTimeStep/UCTimeStep) * (Nb_APD + Nb_SS + Nb_RGP)
	ScenarioSize=NumberSSVTimeSteps*SubScenarioSize
	Block.createDimension("SubScenarioSize",SubScenarioSize)
	Block.createDimension("ScenarioSize",ScenarioSize)
	NumberRandomDataGroups=0
	SizeRandomDataGroupsData=[]
	if Nb_APD>0: 
		NumberRandomDataGroups=NumberRandomDataGroups+1
		SizeRandomDataGroupsData=SizeRandomDataGroupsData+[Nb_APD*SSVTimeStep/UCTimeStep]
	if Nb_SS>0: 
		NumberRandomDataGroups=NumberRandomDataGroups+1
		SizeRandomDataGroupsData=SizeRandomDataGroupsData+[Nb_SS*SSVTimeStep/UCTimeStep]
	if Nb_TPP>0: 
		NumberRandomDataGroups=NumberRandomDataGroups+1
		SizeRandomDataGroupsData=SizeRandomDataGroupsData+[Nb_TPP*ThermalMaxPowerSize]
	if Nb_RGP>0: 
		NumberRandomDataGroups=NumberRandomDataGroups+1
		SizeRandomDataGroupsData=SizeRandomDataGroupsData+[Nb_RGP*SSVTimeStep/UCTimeStep]
	Block.createDimension("NumberRandomDataGroups",NumberRandomDataGroups)
	AdmissibleStateSize=Block.createDimension("AdmissibleStateSize",NumberSSVTimeSteps*TotalNumberReservoirs)
	InitialStateSize=Block.createDimension("InitialStateSize",TotalNumberReservoirs)

	SizeRandomDataGroups=Block.createVariable("SizeRandomDataGroups",'u4',("NumberRandomDataGroups"))
	SizeRandomDataGroups[:]=SizeRandomDataGroupsData
	
	# create scenarios
	logger.info('fill scenarios')
	Scenarios=Block.createVariable("Scenarios",np.double,("NumberScenarios","ScenarioSize"))
	indexScenario=0
	for scenario in ListScenarios:
		logger.info('scenario '+scenario)
		ScenarioData=pd.Series(index=range( NumberSSVTimeSteps ),dtype=object)
		for t in range(NumberSSVTimeSteps):
			start=datesSSV.loc[t]['start']
			if Nb_APD>0:
				ScenarioData.loc[t]=np.concatenate([ DemandScenarios.loc[node][scenario][ (DemandScenarios.loc[node][scenario].index >=datesSSV.loc[t]['start'] ) & (DemandScenarios.loc[node][scenario].index <=datesSSV.loc[t]['end'] ) ] for node in Nodes ]  )
				if Nb_SS>0:
					ScenarioData.loc[t]=np.concatenate([ScenarioData.loc[t],np.concatenate([ InflowsScenarios.loc[reservoir][scenario][ (InflowsScenarios.loc[reservoir][scenario].index >=datesSSV.loc[t]['start'] ) & (InflowsScenarios.loc[reservoir][scenario].index <=datesSSV.loc[t]['end'] ) ] for reservoir in InflowsScenarios.index  ]  )])
				if Nb_TPP>0:
					ScenarioData.loc[t]=np.concatenate([ScenarioData.loc[t],np.concatenate([ ThermalScenarios.loc[th][scenario][ (ThermalScenarios.loc[th][scenario].index >=datesSSV.loc[t]['start'] ) & (ThermalScenarios.loc[th][scenario].index <=datesSSV.loc[t]['end'] ) ] for th in ThermalScenarios.index ]  )])
				if Nb_RGP>0:
					ScenarioData.loc[t]=np.concatenate([ScenarioData.loc[t],np.concatenate([ RESScenarios.loc[res][scenario][ (RESScenarios.loc[res][scenario].index >=datesSSV.loc[t]['start'] ) & (RESScenarios.loc[res][scenario].index <=datesSSV.loc[t]['end'] ) ] for res in RESScenarios.index ]  )])
			elif Nb_SS>0:
				ScenarioData.loc[t]=np.concatenate([ InflowsScenarios.loc[reservoir][scenario][ (InflowsScenarios.loc[reservoir][scenario].index >=datesSSV.loc[t]['start'] ) & (InflowsScenarios.loc[reservoir][scenario].index <=datesSSV.loc[t]['end'] ) ] for reservoir in InflowsScenarios.index  ]  )
				if Nb_TPP>0:
					ScenarioData.loc[t]=np.concatenate([ScenarioData.loc[t],np.concatenate([ ThermalScenarios.loc[th][scenario][ (ThermalScenarios.loc[th][scenario].index >=datesSSV.loc[t]['start'] ) & (ThermalScenarios.loc[th][scenario].index <=datesSSV.loc[t]['end'] ) ] for th in ThermalScenarios.index ]  )])
				if Nb_RGP>0:
					ScenarioData.loc[t]=np.concatenate([ScenarioData.loc[t],np.concatenate([ RESScenarios.loc[res][scenario][ (RESScenarios.loc[res][scenario].index >=datesSSV.loc[t]['start'] ) & (RESScenarios.loc[res][scenario].index <=datesSSV.loc[t]['end'] ) ] for res in RESScenarios.index ]  )])
			elif Nb_TPP>0:
				ScenarioData.loc[t]=np.concatenate([ ThermalScenarios.loc[th][scenario][ (ThermalScenarios.loc[th][scenario].index >=datesSSV.loc[t]['start'] ) & (ThermalScenarios.loc[th][scenario].index <=datesSSV.loc[t]['end'] ) ] for th in ThermalScenarios.index ]  )
				if Nb_RGP>0:
					ScenarioData.loc[t]=np.concatenate([ScenarioData.loc[t],np.concatenate([ RESScenarios.loc[res][scenario][ (RESScenarios.loc[res][scenario].index >=datesSSV.loc[t]['start'] ) & (RESScenarios.loc[res][scenario].index <=datesSSV.loc[t]['end'] ) ] for res in RESScenarios.index ]  )])
			elif Nb_RGP>0:
				ScenarioData.loc[t]=np.concatenate([ RESScenarios.loc[res][scenario][ (RESScenarios.loc[res][scenario].index >=datesSSV.loc[t]['start'] ) & (RESScenarios.loc[res][scenario].index <=datesSSV.loc[t]['end'] ) ] for res in RESScenarios.index ]  )
			
		datascenario=np.concatenate([ScenarioData.loc[t] for t in range(NumberSSVTimeSteps)   ])

		Scenarios[indexScenario,:]=datascenario
		indexScenario=indexScenario+1
	
	StateSize=Block.createVariable("StateSize",'u4',())
	StateSize[:]=TotalNumberReservoirs
	AdmissibleState=Block.createVariable("AdmissibleState",np.double,("AdmissibleStateSize"))
	InitialState=Block.createVariable("InitialState",np.double,("InitialStateSize"))
	AdmissibleStateData=np.zeros(shape=(int(SS['NumberUnits'].sum()),NumberSSVTimeSteps))
	AdmissibleStateDownData=np.zeros(shape=(int(SS['NumberUnits'].sum()),NumberSSVTimeSteps))
	indexHydroUnit=0
	indexInitialState=0
	for hs in HSSS.index:
		for hu in HSSS.loc[hs].index:
			for t in range( NumberSSVTimeSteps ):
				if 'InitialVolume' in SS.columns:
					if type(HSSS.loc[hs].loc[hu]['NumberReservoirs'])==str:
						AdmissibleStateData[indexHydroUnit][t]=DeterministicTimeSeries[HSSS.loc[hs].loc[hu]['InitialVolume']].loc[datesSSV.loc[t]['start'] ]
					else:
						AdmissibleStateData[indexHydroUnit][t]=HSSS.loc[hs].loc[hu]['InitialVolume']
					AdmissibleStateDownData[indexHydroUnit][t]=HSSS.loc[hs].loc[hu]['MaxVolume']
				if t==0:
					if 'InitialVolume' in SS.columns:
						if type(HSSS.loc[hs].loc[hu]['NumberReservoirs'])==str:
							if HSSS.loc[hs].loc[hu]['NumberReservoirs']==1:
								InitialState[indexInitialState]=DeterministicTimeSeries[HSSS.loc[hs].loc[hu]['InitialVolume']].loc[datesSSV.loc[0]['start'] ]
								indexInitialState=indexInitialState+1
							elif HSSS.loc[hs].loc[hu]['NumberReservoirs']==2:
								InitialState[indexInitialState]=DeterministicTimeSeries[HSSS.loc[hs].loc[hu]['InitialVolume']].loc[datesSSV.loc[0]['start'] ]
								InitialState[indexInitialState+1]=HSSS.loc[hs].loc[hu]['MaxVolume']
								indexInitialState=indexInitialState+2
						else:
							if HSSS.loc[hs].loc[hu]['NumberReservoirs']==1:
								InitialState[indexInitialState]=HSSS.loc[hs].loc[hu]['InitialVolume']
								indexInitialState=indexInitialState+1
							elif HSSS.loc[hs].loc[hu]['NumberReservoirs']==2:
								InitialState[indexInitialState]=HSSS.loc[hs].loc[hu]['InitialVolume']
								InitialState[indexInitialState+1]=HSSS.loc[hs].loc[hu]['MaxVolume']
								indexInitialState=indexInitialState+1+2
					else:
						if HSSS.loc[hs].loc[hu]['NumberReservoirs']==1:
							InitialState[indexInitialState]=0
							indexInitialState=indexInitialState+1
						elif HSSS.loc[hs].loc[hu]['NumberReservoirs']==2:
							InitialState[indexInitialState]=0
							InitialState[indexInitialState+1]=HSSS.loc[hs].loc[hu]['MaxVolume']
							indexInitialState=indexInitialState+1+2
			indexHydroUnit=indexHydroUnit+1
	indexHydroUnit=0
	indexAdmissibleState=0
	for t in range( NumberSSVTimeSteps ):
		indexHydroUnit=0
		for hs in HSSS.index:
			for hu in HSSS.loc[hs].index:
				if HSSS.loc[hs].loc[hu]['NumberReservoirs']==1:
					AdmissibleState[indexAdmissibleState]=AdmissibleStateData[indexHydroUnit][t]
					indexAdmissibleState=indexAdmissibleState+1
				elif HSSS.loc[hs].loc[hu]['NumberReservoirs']==2:
					AdmissibleState[indexAdmissibleState]=AdmissibleStateData[indexHydroUnit][t]
					AdmissibleState[indexAdmissibleState+1]=AdmissibleStateDownData[indexHydroUnit][t]
					indexAdmissibleState=indexAdmissibleState+2
				indexHydroUnit=indexHydroUnit+1
	
	##################################################################################
	# create stochastic blocks and benders blocks
	##################################################################################
	logger.info('Create StochasticBlocks')
	for indexSSV in range(NumberSSVTimeSteps):
		# create blocks
		StochasticBlocks=Block.createGroup("StochasticBlock_"+str(indexSSV))
		StochasticBlocks.type="StochasticBlock"
		StochasticBlocks.createDimension("NumberDataMappings",NumberDataMappings)
		StochasticBlocks.createDimension("SetSizeSize",2*NumberDataMappings)
		SetElementsSize=4 * (Nb_SS + Nb_TPP + Nb_RGP)+2 * Nb_SS 
		if Nb_SS == 0: SetElementsSize=SetElementsSize+4
		if Nb_APDTo>0: SetElementsSize=SetElementsSize+NumberNodes * (Nb_APDTo + 2)
		StochasticBlocks.createDimension("SetElementsSize",SetElementsSize)

		# fill stochastic blocks
		DataType=StochasticBlocks.createVariable("DataType",'S1' ,("NumberDataMappings"))
		DataType[:]=NumberDataMappings*'D'

		SetSize=StochasticBlocks.createVariable("SetSize","u4",("SetSizeSize"))
		SetSizeData=NumberNodes*[Nb_APDTo,0]+2 * (Nb_SS + Nb_TPP + Nb_RGP)*[0]+2*[Nb_SS]
		SetSize[:]=np.array(SetSizeData)
		
		SetElements=StochasticBlocks.createVariable("SetElements","u4",("SetElementsSize"))
		SetElementsData=[]
		if Nb_APDTo > 0:
			for i in range(NumberNodes):
				for t in range(Nb_APDTo):
					SetElementsData=SetElementsData+[ t+ i * SSVTimeStep/UCTimeStep ]
				SetElementsData=SetElementsData+[ i * Nb_APDTo, (i+1)*Nb_APDTo]
		begin = NumberNodes * Nb_APDTo
		inflow_begin=begin
		for i in range(Nb_SS):
			SetElementsData=SetElementsData+ [ begin + i*SSVTimeStep/UCTimeStep, begin + (i+1)*SSVTimeStep/UCTimeStep,0,SSVTimeStep/UCTimeStep ]
		begin=begin+Nb_SS*SSVTimeStep/UCTimeStep
		for i in range(Nb_TPP):
			SetElementsData=SetElementsData+[begin + i * ThermalMaxPowerSize,begin + (i + 1) * ThermalMaxPowerSize,0,SSVTimeStep/UCTimeStep  ]
		begin=begin+Nb_TPP * ThermalMaxPowerSize
		for i in range(Nb_RGP):
			SetElementsData=SetElementsData+[ begin + i*SSVTimeStep/UCTimeStep, begin + (i+1)*SSVTimeStep/UCTimeStep,0,SSVTimeStep/UCTimeStep]
		if Nb_SS==0: SetElementsData=SetElementsData+[0,0,0,0]
		else:
			for i in range(Nb_SS):
				SetElementsData=SetElementsData+[inflow_begin + i * SSVTimeStep/UCTimeStep]
			reservoir_index = 0
			for reservoir in SS.index:
				SetElementsData=SetElementsData+[reservoir_index]
				reservoir_index=reservoir_index+SS.loc[reservoir]['NumberReservoirs']
		SetElements[:]=np.array(SetElementsData)
		
		FunctionName=StochasticBlocks.createVariable("FunctionName",str,("NumberDataMappings"))
		FunctionNameData=[]
		if Nb_APDTo > 0: FunctionNameData=FunctionNameData+NumberNodes*["UCBlock::set_active_power_demand"]
		FunctionNameData=FunctionNameData+Nb_SS*["HydroUnitBlock::set_inflow"]+Nb_TPP*["ThermalUnitBlock::set_maximum_power"]+Nb_RGP*["IntermittentUnitBlock::set_maximum_power"]+["BendersBFunction::modify_constants"]
		FunctionName[:]=np.array(FunctionNameData)

		Caller=StochasticBlocks.createVariable("Caller",'S1',("NumberDataMappings"))
		CallerData=(NumberDataMappings-1)*"B"+"F"
		Caller[:]=np.array(list(CallerData))

		# create abstract path of stochastic blocks
		APSB=StochasticBlocks.createGroup("AbstractPath")
		APSB.createDimension("PathDim",NumberDataMappings)
		totalLength = (4 * Nb_SS) + (3 * Nb_TPP) + (3 * Nb_RGP) + 1
		demandPathTotalLength = 0
		if Nb_APDTo>0: 
			demandPathTotalLength=2*NumberNodes
			totalLength=totalLength+demandPathTotalLength
		APSB.createDimension("TotalLength",totalLength)
		
		PathStart=APSB.createVariable("PathStart",'u4',("PathDim"))
		PathStartData=[]
		if Nb_APDTo>0: 
			for i in range(NumberNodes): 
				PathStartData=PathStartData+[2*i]
		for s in range(Nb_SS): 
			PathStartData=PathStartData+[demandPathTotalLength + 4*s]
		for s in range(Nb_TPP+Nb_RGP): 
			PathStartData=PathStartData+[demandPathTotalLength+ 4*Nb_SS  +3*s]
		PathStartData=PathStartData+[demandPathTotalLength+ 4*Nb_SS +3*Nb_TPP+3*Nb_RGP]
		PathStart[:]=PathStartData		
		
		PathNodeTypes=APSB.createVariable("PathNodeTypes",'S1',("TotalLength"))
		PathNodeTypesData=""
		if Nb_APDTo>0: 
			PathNodeTypesData=NumberNodes*"OB"
		PathNodeTypesData=PathNodeTypesData+Nb_SS*"OBBB"+(Nb_TPP+Nb_RGP)*"OBB"+"O"
		PathNodeTypes[:]=np.array(list(PathNodeTypesData))
			
		PathGroupIndices=APSB.createVariable("PathGroupIndices",'u4',("TotalLength"))
		PathGroupIndicesIndex=0
		for n in range(NumberNodes):
			PathGroupIndices[PathGroupIndicesIndex+1]=0
			PathGroupIndicesIndex=PathGroupIndicesIndex+2
		num_previous_units=0
		for h in range(NumberHydroSystems): # loop on hydrosystems
			indexUnitInHS=0
			for u in range(int(HSSS.loc[h]['NumberUnits'].sum())):  # loop on hydrounits of hydrosystem h
				PathGroupIndices[PathGroupIndicesIndex+1]=0
				PathGroupIndices[PathGroupIndicesIndex+2]=h
				PathGroupIndices[PathGroupIndicesIndex+3]=u+num_previous_units
				indexUnitInHS=indexUnitInHS+1
				PathGroupIndicesIndex=PathGroupIndicesIndex+4
			num_previous_units=HSSS.loc[h]['NumberUnits'].sum()-1
		for i in range(Nb_TPP):
			PathGroupIndices[PathGroupIndicesIndex+1]=0
			PathGroupIndices[PathGroupIndicesIndex+2]=NumberHydroSystems+i
			PathGroupIndicesIndex=PathGroupIndicesIndex+3
		for i in range(Nb_RGP):
			PathGroupIndices[PathGroupIndicesIndex+1]=0
			PathGroupIndices[PathGroupIndicesIndex+2]=NumberHydroSystems+NumberThermalUnits+i
			PathGroupIndicesIndex=PathGroupIndicesIndex+3
		
		BendersBlocks = StochasticBlocks.createGroup("Block")
		BendersBlocks.type="BendersBlock"
		
		# fill benders blocks
		BendersBlocks.createDimension("NumVar",TotalNumberReservoirs)
		
		# add BendersBFunction to BendersBlocks
		BBF=BendersBlocks.createGroup("BendersBFunction")
		BBF.createDimension("NumVar",TotalNumberReservoirs)
		BBF.createDimension("NumRow",TotalNumberReservoirs)
		BBF.createDimension("NumNonzero",TotalNumberReservoirs)
		
		# include UCBlocks
		UCB=BBF.createGroup("Block")
		UCB.id=str(indexSSV)
		UCB.filename="Block_"+str(indexSSV)+".nc4"
		
		# create abstract path of bendersbfunction
		AP=BBF.createGroup("AbstractPath")
		PathDim=TotalNumberReservoirs
		TotalLength=3*TotalNumberReservoirs
		AP.createDimension("PathDim",PathDim)
		AP.createDimension("TotalLength",TotalLength)
		
		PathNodeTypes=AP.createVariable("PathNodeTypes",'S1',("TotalLength"))
		PathNodeTypesData=TotalNumberReservoirs*"BBC"
		PathNodeTypes[:]=np.array(list(PathNodeTypesData))
		
		PathGroupIndices=AP.createVariable("PathGroupIndices",'u4',("TotalLength"))
		PathGroupIndicesData=[]
		
		for h in range(NumberHydroSystems): # loop on hydrosystems
			indexUnitInHS=0
			for u in HSSS.loc[h].index:  # loop on hydrounits of hydrosystem h
				for r in range(SS.loc[u]['NumberReservoirs']):
					PathGroupIndicesData=PathGroupIndicesData+[h,indexUnitInHS,0]
				indexUnitInHS=indexUnitInHS+1
		PathGroupIndices[:]=PathGroupIndicesData
		
		PathElementIndices=AP.createVariable("PathElementIndices",'u4',("TotalLength"))
		PathElementIndicesIndex=2
		indexUnitInHS=0
		for h in range(NumberHydroSystems): # loop on hydrosystems
			for u in HSSS.loc[h].index:  # loop on hydrounits of hydrosystem h
				for r in range(SS.loc[u]['NumberReservoirs']):
					PathElementIndices[PathElementIndicesIndex]=r
					PathElementIndicesIndex=PathElementIndicesIndex+3
				indexUnitInHS=indexUnitInHS+1
		
		
		PathStart=AP.createVariable("PathStart",'u4',("PathDim"))
		PathStartData=[]
		for r in range(TotalNumberReservoirs):
			PathStartData=PathStartData+[3*r]
		PathStart[:]=PathStartData
		
	# create SDDPBlock abstract path
	# create abstract path of bendersbfunction
	AP=Block.createGroup("AbstractPath")
	if NumberHydroSystems>0: PathDim=NumberHydroSystems
	else: PathDim=1
	TotalLength=4*PathDim
	AP.createDimension("PathDim",PathDim)
	AP.createDimension("TotalLength",TotalLength)
	
	PathNodeTypes=AP.createVariable("PathNodeTypes",'S1',("TotalLength"))
	PathNodeTypesData=PathDim*"OBBB"
	PathNodeTypes[:]=np.array(list(PathNodeTypesData))

	PathGroupIndices=AP.createVariable("PathGroupIndices",'u4',("TotalLength"))
	PathGroupIndicesIndex=1
	for h in range(NumberHydroSystems): # loop on hydrosystems
		PathGroupIndices[PathGroupIndicesIndex]=0
		PathGroupIndices[PathGroupIndicesIndex+1]=h
		PathGroupIndices[PathGroupIndicesIndex+2]=HSSS.loc[h]['NumberUnits'].sum()
		PathGroupIndicesIndex=PathGroupIndicesIndex+4
		
	PathStart=AP.createVariable("PathStart",'u4',("PathDim"))
	PathStartData=[]
	start=0
	for r in range(PathDim):
		PathStartData=PathStartData+[start]
		start=start+4
	PathStart[:]=PathStartData
	
	SDDPBlock.close()

def createInvestmentBlock(filename):
	InvestBlock=Dataset(filename,'w',format='NETCDF4')
	# general attribute of the file
	InvestBlock.setncattr('SMS++_file_type',1)
	Block=InvestBlock.createGroup("InvestmentBlock")
	Block.type="InvestmentBlock"
	NumAssets=NumberInvestedLines+NumberInvestedThermalUnits+NumberInvestedBatteryUnits+NumberInvestedIntermittentUnits
	Block.createDimension("NumAssets",NumAssets)

	if NumberInvestedThermalUnits>0:
		TUDisagr=TU.loc[TU.index.repeat(TU.NumberUnits)]
		TUDisagr['IndexUnit']=TUDisagr.reset_index().index+1
		if ('MaxAddedCapacity' in TUDisagr.columns and 'MaxRetCapacity' in TUDisagr.columns):
			TUInvested=pd.DataFrame(TUDisagr[ (TUDisagr['MaxRetCapacity']>0) | (TUDisagr['MaxAddedCapacity']>0) ])
			
		elif 'MaxAddedCapacity' in TUDisagr.columns:
			TUInvested=pd.DataFrame(TUDisagr[ TUDisagr['MaxAddedCapacity']>0 ])
			TUInvested['MaxRetCapacity']=0.0
		elif 'MaxRetCapacity' in TUDisagr.columns:
			TUInvested=pd.DataFrame(TUDisagr[ TUDisagr['MaxRetCapacity']>0 ])
			TUInvested['MaxAddedCapacity']=0.0
		TUNotInvested=pd.DataFrame(TUDisagr[ (TUDisagr['MaxRetCapacity']==0) & (TUDisagr['MaxAddedCapacity']==0) ])
		TUInvested['UpperBound']=TUInvested.apply(lambda x: (x.MaxPower+x.MaxAddedCapacity)/x.MaxPower if x.MaxAddedCapacity>0 else 1,axis=1)
		TUInvested['LowerBound']=TUInvested.apply(lambda x: max(x.MaxPower-x.MaxRetCapacity,0)/x.MaxPower if x.MaxRetCapacity>0 else 1,axis=1)
		TUInvested['Cost']=TUInvested.apply(lambda x: (x.InvestmentCost)*x.MaxPower,axis=1)
		TUInvested['AssetType']=0
		TUInvested=TUInvested.reset_index()
		TUNotInvested=TUNotInvested.reset_index()
		
	if NumberInvestedIntermittentUnits>0:
		RESDisagr=RES.loc[RES.index.repeat(RES.NumberUnits)]
		RESDisagr['IndexUnit']=RESDisagr.reset_index().index+NumberThermalUnits+1
		if ('MaxAddedCapacity' in RES.columns and 'MaxRetCapacity' in RES.columns):
			RESInvested=pd.DataFrame(RESDisagr[ (RESDisagr['MaxRetCapacity']>0) | (RESDisagr['MaxAddedCapacity']>0) ])
		elif 'MaxAddedCapacity' in RES.columns:
			RESInvested=pd.DataFrame(RESDisagr[ RESDisagr['MaxAddedCapacity']>0 ])
			RESInvested['MaxRetCapacity']=0.0
		elif 'MaxRetCapacity' in RES.columns:
			RESInvested=pd.DataFrame(RESDisagr[ RESDisagr['MaxRetCapacity']>0 ])
			RESInvested['MaxAddedCapacity']=0.0
		RESNotInvested=pd.DataFrame(RESDisagr[ (RESDisagr['MaxRetCapacity']==0) & (RESDisagr['MaxAddedCapacity']==0) ])
		RESInvested['UpperBound']=RESInvested.apply(lambda x: (x.MaxPower+x.MaxAddedCapacity)/x.MaxPower if x.MaxAddedCapacity>0 else 1,axis=1)
		RESInvested['LowerBound']=RESInvested.apply(lambda x: max(x.MaxPower-x.MaxRetCapacity,0)/x.MaxPower if x.MaxRetCapacity>0 else 1,axis=1)
		RESInvested['Cost']=RESInvested.apply(lambda x: x.InvestmentCost*x.MaxPower,axis=1)
		RESInvested['AssetType']=0
		RESInvested=RESInvested.reset_index()
		RESNotInvested=RESNotInvested.reset_index()

	if NumberInvestedBatteryUnits>0:
		ISts=True
		STSDisagr=STS.loc[STS.index.repeat(STS.NumberUnits)]
		STSDisagr['IndexUnit']=STSDisagr.reset_index().index+NumberThermalUnits+NumberIntermittentUnits+1
		if ('MaxAddedCapacity' in STS.columns and 'MaxRetCapacity' in STS.columns):
			STSInvested=pd.DataFrame(STSDisagr[ (STSDisagr['MaxRetCapacity']>0) | (STSDisagr['MaxAddedCapacity']>0) ])
		elif 'MaxAddedCapacity' in STS.columns:
			STSInvested=pd.DataFrame(STSDisagr[ STSDisagr['MaxAddedCapacity']>0 ])
			STSInvested['MaxRetCapacity']=0.0
		elif 'MaxRetCapacity' in STS.columns:
			STSInvested=pd.DataFrame(STSDisagr[ STSDisagr['MaxRetCapacity']>0 ])
			STSInvested['MaxAddedCapacity']=0.0
		STSInvested['UpperBound']=STSInvested.apply(lambda x: (x.MaxPower+x.MaxAddedCapacity)/x.MaxPower if x.MaxAddedCapacity>0 else 1,axis=1)
		STSInvested['LowerBound']=STSInvested.apply(lambda x: max(x.MaxPower-x.MaxRetCapacity,0)/x.MaxPower if x.MaxRetCapacity>0 else 1,axis=1)
		STSInvested['Cost']=STSInvested.apply(lambda x: x.InvestmentCost*x.MaxPower,axis=1)
		STSInvested['AssetType']=0
		STSInvested=STSInvested.reset_index()


	if NumberInvestedLines>0:
		Iint=True
		IN['IndexUnit']=IN.reset_index().index
		if ('MaxAddedCapacity' in IN.columns and 'MaxRetCapacity' in IN.columns):
			INInvested=pd.DataFrame(IN[ (IN['MaxRetCapacity']>0) | (IN['MaxAddedCapacity']>0) ])
		elif 'MaxAddedCapacity' in IN.columns:
			INInvested=pd.DataFrame(IN[ IN['MaxAddedCapacity']>0 ])
			INInvested['MaxRetCapacity']=0.0
		elif 'MaxRetCapacity' in IN.columns:
			INInvested=pd.DataFrame(IN[ IN['MaxRetCapacity']>0 ])
			INInvested['MaxAddedCapacity']=0.0
		INInvested['UpperBound']=INInvested.apply(lambda x: (x.MaxPowerFlow+x.MaxAddedCapacity)/x.MaxPowerFlow if x.MaxAddedCapacity>0 else 1,axis=1)
		INInvested['LowerBound']=INInvested.apply(lambda x: max(x.MaxPowerFlow-x.MaxRetCapacity,0)/x.MaxPowerFlow if x.MaxRetCapacity>0 else 1,axis=1)
		INInvested['Cost']=INInvested.apply(lambda x: x.InvestmentCost*x.MaxPowerFlow,axis=1)
		INInvested['AssetType']=1
		INInvested=INInvested.reset_index()

	Assets=Block.createVariable("Assets",'u4',("NumAssets"))
	AssetType=Block.createVariable("AssetType",'u4',("NumAssets"))
	UpperBound=Block.createVariable("UpperBound",np.double,("NumAssets"))
	LowerBound=Block.createVariable("LowerBound",np.double,("NumAssets"))
	Cost=Block.createVariable("Cost",np.double,("NumAssets"))
	
	listAssets=[]
	listAssetType=[]
	listUpperBound=[]
	listLowerBound=[]
	listCost=[]
	if NumberInvestedThermalUnits>0: 
		listAssets.append(np.array(TUInvested['IndexUnit']))
		listAssetType.append(np.array(TUInvested['AssetType']))
		listUpperBound.append(np.array(TUInvested['UpperBound']))
		listLowerBound.append(np.array(TUInvested['LowerBound']))
		listCost.append(np.array(TUInvested['Cost']))
	if NumberInvestedIntermittentUnits>0: 
		listAssets.append(np.array(RESInvested['IndexUnit']))
		listAssetType.append(np.array(RESInvested['AssetType']))
		listUpperBound.append(np.array(RESInvested['UpperBound']))
		listLowerBound.append(np.array(RESInvested['LowerBound']))
		listCost.append(np.array(RESInvested['Cost']))
	if NumberInvestedBatteryUnits>0: 
		listAssets.append(np.array(STSInvested['IndexUnit']))
		listAssetType.append(np.array(STSInvested['AssetType']))
		listUpperBound.append(np.array(STSInvested['UpperBound']))
		listLowerBound.append(np.array(STSInvested['LowerBound']))
		listCost.append(np.array(STSInvested['Cost']))
	if NumberInvestedLines>0:
		listAssets.append(np.array(INInvested['IndexUnit']))
		listAssetType.append(np.array(INInvested['AssetType']))
		listUpperBound.append(np.array(INInvested['UpperBound']))
		listLowerBound.append(np.array(INInvested['LowerBound']))
		listCost.append(np.array(INInvested['Cost']))
	
	if NumberInvestedThermalUnits+NumberInvestedIntermittentUnits+NumberInvestedBatteryUnits+NumberInvestedLines>0:
		Assets[:]=np.concatenate(listAssets)
		AssetType[:]=np.concatenate(listAssetType)
		UpperBound[:]=np.concatenate(listUpperBound)
		LowerBound[:]=np.concatenate(listLowerBound)
		Cost[:]=np.concatenate(listCost)
	else:
		Assets[:]=np.array([1])
		AssetType[:]=np.array([0])
		UpperBound[:]=np.array([1])
		LowerBound[:]=np.array([0])
		Cost[:]=np.array([0])
			
	if cfg['Invest']=='NRJ':
		Block.createDimension("NumConstraints",NumberNodes)
		Constraints_A=Block.createVariable("Constraints_A",np.double,("NumConstraints","NumAssets"))
		Constraints_LowerBound=Block.createVariable("Constraints_LowerBound",np.double,("NumConstraints"))
		# each region has enough energy
		indexNode=0
		for node in Nodes:
			Adata=np.array(TUInvested.apply(lambda x: cfg['ParametersFormat']['NumberHoursInYear']*x['MaxPower'] if x['Zone']==node else 0,axis=1))
			Adata=np.concatenate([Adata,np.array(RESInvested.apply(lambda x: x['Energy_Timeserie']*x['MaxPower'] if x['Zone']==node else 0,axis=1))])
			Adata=np.concatenate([Adata,np.zeros(shape=NumberInvestedBatteryUnits+NumberInvestedLines)])
			if len(TUNotInvested[ TUNotInvested['Zone']==node ].index)>0:
				ThermalEnergyNotInvested=cfg['ParametersFormat']['NumberHoursInYear']*TUNotInvested[ TUNotInvested['Zone']==node ]['MaxPower'].sum()
			else:
				ThermalEnergyNotInvested=0
			if len(RESNotInvested[ RESNotInvested['Zone']==node ].index) >0:
				RESEnergyNotInvested=RESNotInvested[ RESNotInvested['Zone']==node ]['EnergyMaxPower'].sum()
			else:
				RESEnergyNotInvested=0
			if ('Reservoir',node) in SS.index:
				SSEnergy=SS.loc['Reservoir',node]['Inflows'].sum()
			else:
				SSEnergy=0
			if ('Pumped Storage',node) in STS.index:
				STSEnergy=STS.loc['Pumped Storage',node]['Energy'].sum()
			else:
				STSEnergy=0
			LowerBound=ZV.loc['Total',node]['value'] - ThermalEnergyNotInvested - RESEnergyNotInvested - SSEnergy -  STSEnergy

			Constraints_A[indexNode,:]=Adata
			Constraints_LowerBound[indexNode]=LowerBound
			indexNode=indexNode+1
			
	if cfg['Invest']=='TargetRES':
		Block.createDimension("NumConstraints",NumberNodes)
		Constraints_A=Block.createVariable("Constraints_A",np.double,("NumConstraints","NumAssets"))
		Constraints_LowerBound=Block.createVariable("Constraints_LowerBound",np.double,("NumConstraints"))
		# each region has enough energy
		indexNode=0
		for node in Nodes:
			Adata=np.zeros(shape= NumberInvestedThermalUnits)
			Adata=np.concatenate([Adata,np.array(RESInvested.apply(lambda x: x['Energy_Timeserie']*x['MaxPower'] if x['Zone']==node else 0,axis=1))])
			Adata=np.concatenate([Adata,np.zeros(shape=NumberInvestedBatteryUnits+NumberInvestedLines)])
			
			if len(RESNotInvested[ RESNotInvested['Zone']==node ].index) >0:
				RESEnergyNotInvested=RESNotInvested[ RESNotInvested['Zone']==node ]['EnergyMaxPower'].sum()
			else:
				RESEnergyNotInvested=0
			if ('Reservoir',node) in SS.index:
				SSEnergy=SS.loc['Reservoir',node]['Inflows'].sum()
			else:
				SSEnergy=0
			if ('Pumped Storage',node) in STS.index:
				STSEnergy=STS.loc['Pumped Storage',node]['Energy'].sum()
			else:
				STSEnergy=0
			LowerBound=ZP[ ZP['Level1']== node ].iloc[0]['RhoTargetRES']*ZV.loc['Total',node]['value'] - RESEnergyNotInvested - SSEnergy - STSEnergy
			
			Constraints_A[indexNode,:]=Adata
			Constraints_LowerBound[indexNode]=LowerBound
			indexNode=indexNode+1

	# add group SDDP
	SDDPBlock = Block.createGroup("SDDPBlock")
	SDDPBlock.id=str(0)
	SDDPBlock.filename="SDDPBlock.nc4"

# read all timeseries data
if 'DeterministicTimeSeries' in cfg:
	logger.info('read deterministic timeseries')
	DeterministicTimeSeries=read_deterministic_timeseries(True)
else:
	DeterministicTimeSeries=read_deterministic_timeseries(False)

if cfg['IncludeScenarisedData'] or 'SDDP' in cfg['FormatMode'] or 'INVEST' in cfg['FormatMode']:
	logger.info('read demand timeseries')
	Nb_APD=Nb_APDTo=Nb_SS=Nb_TPP=Nb_RGP=0
	if 'ActivePowerDemand' in ScenarisedData:
		DemandScenarios=create_demand_scenarios()
		Nb_APDTo=TimeHorizonUC
		Nb_APD=len(DemandScenarios.index)
	if NumberIntermittentUnits>0 and 'Renewable:MaxPowerProfile' in ScenarisedData: 
		logger.info('read renewable generation timeseries')
		RESScenarios=create_res_scenarios()
		Nb_RGP=len(RESScenarios.index)
	if NumberHydroSystems>0 and 'Hydro:Inflows' in ScenarisedData: 
		logger.info('read inflows timeseries')
		InflowsScenarios=create_inflows_scenarios()
		Nb_SS=len(InflowsScenarios.index)
	if NumberThermalUnits>=0 and 'Thermal:MaxPowerProfile' in ScenarisedData: 
		logger.info('read thermal maxpower timeseries')
		ThermalScenarios=create_thermal_scenarios()
		Nb_TPP=len(ThermalScenarios.index)
	ThermalMaxPowerSize=SSVTimeStep/ThermalMaxPowerTimeSpan
	NumberDataMappings = Nb_SS + Nb_TPP + Nb_RGP + 1
	if Nb_APDTo > 0:
		NumberDataMappings = NumberDataMappings+NumberNodes;
	
listData=[]
if NumberHydroSystems>0:
	#listData.append((SS,'SS'))
	for hs in range(NumberHydroSystems): listData.append((HSSS.loc[hs],'HSSS'))
if NumberThermalUnits>0:
	listData.append((TU,'TU'))
if NumberIntermittentUnits>0:
	listData.append((RES,'RES'))
if NumberSyncUnits>0:
	listData.append((SYN,'SYN'))
if NumberBatteryUnits>0:
	listData.append((STS,'STS'))
if cfg['FormatMode']=='SingleUC':
	logger.info('create single UCBlock on the whole period =>'+cfg['outputpath']+'UCBlock.nc4')
	createUCBlock(cfg['outputpath']+'UCBlock.nc4',0,ListScenarios[0],datesSSV.loc[0]['start'],datesSSV.loc[0]['end'])

elif cfg['FormatMode']=='UC':
	logger.info('create One UCBlock per SSV timestep =>'+cfg['outputpath']+'UCBlock_*.nc4')
	for i in range(NumberSSVTimeSteps):
		logger.info('Create UCBlock '+str(i)+' from '+str(datesSSV.loc[i]['start'])+' to '+str(datesSSV.loc[i]['end'])+' => '+cfg['outputpath']+'Block_'+str(i)+'.nc4')
		createUCBlock(cfg['outputpath']+'Block_'+str(i)+'.nc4',i,ListScenarios[0],datesSSV.loc[i]['start'],datesSSV.loc[i]['end'])

elif cfg['FormatMode']=='SDDP':
	logger.info('create One SDDPBlock =>'+cfg['outputpath']+'SDDPBlock.nc4')
	createSDDPBlock(cfg['outputpath']+'SDDPBlock.nc4',0)

elif cfg['FormatMode']=='SDDPandUC':
	logger.info('create One SDDPBlock and One UCBlock per SSV timestep =>'+cfg['outputpath']+'SDDPBlock.nc4')
	createSDDPBlock(cfg['outputpath']+'SDDPBlock.nc4',0)
	for i in range(NumberSSVTimeSteps):
		logger.info('Create UCBlock '+str(i)+' from '+str(datesSSV.loc[i]['start'])+' to '+str(datesSSV.loc[i]['end'])+' => '+cfg['outputpath']+'Block_'+str(i)+'.nc4')
		createUCBlock(cfg['outputpath']+'Block_'+str(i)+'.nc4',i,ListScenarios[0],datesSSV.loc[i]['start'],datesSSV.loc[i]['end'])

elif cfg['FormatMode']=='INVEST':
	logger.info('create One InvestmentBlock, => '+cfg['outputpath']+'InvestmentBlock.nc4')
	createInvestmentBlock(cfg['outputpath']+'InvestmentBlock.nc4')

elif cfg['FormatMode']=='INVESTandSDDP':
	logger.info('create One SDDPBlock and one InvestmentBlock =>'+cfg['outputpath']+'SDDPBlock.nc4  and '+cfg['outputpath']+'InvestmentBlock.nc4'  )
	createSDDPBlock(cfg['outputpath']+'SDDPBlock.nc4',0)
	createInvestmentBlock(cfg['outputpath']+'InvestmentBlock.nc4')

elif cfg['FormatMode']=='INVESTandSDDPandUC':
	logger.info('create One InvestmentBlock, One SDDPBlock and One UCBlock per SSV timestep => '+cfg['outputpath']+'SDDPBlock.nc4 , '+cfg['outputpath']+'InvestmentBlock.nc4' )
	createSDDPBlock(cfg['outputpath']+'SDDPBlock.nc4',0)
	createInvestmentBlock(cfg['outputpath']+'InvestmentBlock.nc4')
	for i in range(NumberSSVTimeSteps):
		logger.info('Create UCBlock '+str(i)+' from '+str(datesSSV.loc[i]['start'])+' to '+str(datesSSV.loc[i]['end'])+' => '+cfg['outputpath']+'Block_'+str(i)+'.nc4')
		createUCBlock(cfg['outputpath']+'Block_'+str(i)+'.nc4',i,ListScenarios[0],datesSSV.loc[i]['start'],datesSSV.loc[i]['end'])

log_and_exit(0, cfg['path'])
