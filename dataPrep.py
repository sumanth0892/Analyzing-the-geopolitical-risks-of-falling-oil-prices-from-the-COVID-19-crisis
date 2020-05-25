import os 
import numpy as np
import pandas as pd 

#There are three files to be processed here and create three time series dataframes 

def getEPUdata():
	epuFile = 'USEPU.csv'
	usepu = pd.read_csv(epuFile,parse_dates = True,index_col = 'DATE')
	return usepu 

def getGPRdata():
	gprFile = 'GPRdaily.csv'
	gprData = pd.read_csv(gprFile,parse_dates = True,index_col = 'DATE')
	gprData = gprData.drop(['DAY_NUMBERS','EVENT',
		'*value of the index on the day the newspaper was published'],axis = 1)
	return gprData 

def getCrudeOilData():
	crudeFile = 'CrudeOilFutures.csv'
	crudeFutures = pd.read_csv(crudeFile,parse_dates = True,
					index_col = 'Date')
	crudeFutures = crudeFutures.sort_index(ascending = True)
	return crudeFutures 
