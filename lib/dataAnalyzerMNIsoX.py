##!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 
#Last Updated: March 15, 2024

import copy
import itertools
import pandas as pd
import numpy as np
import os
import organizeData
import re 
import numpy as np
import readCSVAndSimulate as sim
import dataScreenIsoX

def readIsoX(filePath):
    '''
    Read in the 'combined' output from isoX; rename & drop columns as desired. 

    Inputs:
        fileName: A string; the name of the file ('combined.isox')

    Outputs:
        IsoXDf: A dataframe with the data from the .csv output. 
    '''
    IsoXDf = pd.read_csv(filePath, sep = '\t')

    IsoXDf.drop(columns=['ions.incremental','peakResolution','basePeakIntensity','rawOvFtT','intensCompFactor','agc','analyzerTemperature','numberLockmassesFound'], inplace = True)
    IsoXDf.rename(columns={'scan.no':'scanNumber','time.min':'retTime','it.ms':'integTime','mzMeasured':'mass'},inplace = True)
    #IsoX does not return TIC*IT by default. We can calculate it this way (divide by 1000 to convert from ms). The values calculated differ from the FTStat TIC*IT by <0.02%, in the files I checked. 
    IsoXDf['TIC*IT'] = IsoXDf['tic'] * IsoXDf['integTime'] / 1000

    return IsoXDf

def calculate_Counts_And_ShotNoise(peakDf,resolution=120000,CN=4.4,z=1):
    '''
    Calculate counts of each scan peak
    
    Inputs: 
        peakDf: An individual dataframe consisting of a single peak extracted by FTStatistic.
        CN: A factor from the 2017 paper to convert intensities into counts
        resolution: A reference resolution, for the same purpose (do not change!)
        z: The charge of the ion, used to convert intensities into counts
        Microscans: The number of scans a cycle is divided into, typically 1.
        
    Outputs: 
        The inputDF, with a column for 'counts' added. 
    '''
    peakDf['counts'] = (peakDf['intensity'] /
                  peakDf['peakNoise']) * (CN/z) *(resolution/peakDf['resolution'])**(0.5) * peakDf['microscans']**(0.5)
    return peakDf

def findMostAbundantSub(mergedDf, subNameList):
    counts = []
    for sub in subNameList:
        thisCounts = mergedDf['counts'+sub].sum()
        counts.append(thisCounts)

    max_index = np.argmax(counts)
    max_sub = subNameList[max_index]

    return max_sub

def calc_Append_Ratios(mergedDf, subNameList, mostAbundant = True):
    '''
    Calculates the ratios for each combination of substitutions. Calculates all ratios in the order such that they are < 1.

    Inputs:
        mergedDf: A dataframe with all information for a single file. 
        subNameList: A list of substitution names, e.g. ['13C','18O','D']

    Outputs: 
        mergedDf: The same dataframe with ratios added. 
    '''
    max_sub = findMostAbundantSub(mergedDf, subNameList)

    for sub1, sub2 in itertools.combinations(subNameList,2):
        if ((mostAbundant) and (sub1 != max_sub) and (sub2 != max_sub)): 
            continue
        if mergedDf['counts' + sub1].sum() <= mergedDf['counts' + sub2].sum():
            mergedDf[sub1 + '/' + sub2] = mergedDf['counts' + sub1] / mergedDf['counts' + sub2]
        else:
            mergedDf[sub2 + '/' + sub1] = mergedDf['counts' + sub2] / mergedDf['counts' + sub1]
                            
    return mergedDf
    
def calc_MN_Rel_Abundance(mergedDf, subNameList):
    mergedDf['total Counts'] = 0
    for sub in subNameList:
        mergedDf['total Counts'] += mergedDf['counts' + sub]
        
    for sub in subNameList:
        mergedDf['MN Relative Abundance ' + sub] = mergedDf['counts' + sub] / mergedDf['total Counts']
        
    return mergedDf
    
def cull_By_Time(mergedDf, timeBounds, scanNumber = False):
    if scanNumber:
        mergedDf = mergedDf[mergedDf['scanNumber'].between(timeBounds[0], timeBounds[1], inclusive='both')]
    else:
        mergedDf = mergedDf[mergedDf['retTime'].between(timeBounds[0], timeBounds[1], inclusive='both')]
    return mergedDf
    
def combine_Substituted_Peaks(splitDf, cullOn = None, cullAmt = 3, scanNumber = False, timeBounds = None,MNRelativeAbundance = False):
    '''
    splitDf: A Pandas groupBy object; can be iterated through; iterations give tuples of isotopolog strings ('D','M0', etc.) and dataframes corresponding to the data for that substitution. 
    cullByTime: If True, only include data within certain set of timepoints (by scan # or retTime)
    scanNumber: If True & cullByTime is True, then cull based on scan #. OTherwise, cull by retTime.
    timeBounds: A tuple; the retTime or scanNumber limits to use. 
    MNRelativeAbundance: Calculate MN Relative Abundances rather than ratios. 
    '''
    combinedData = {}
    subIdx = 0 
    subNameList = []
    for subName, subData in splitDf:
        if subName == 'M0':
            subName = 'Unsub'
        subNameList.append(subName)
        subData = calculate_Counts_And_ShotNoise(subData)
        subData.rename(columns={'counts':'counts' + subName,'mass':'mass'+subName,'intensity':'intensity'+subName,'peakNoise':'peakNoise'+subName},inplace = True)

        if subIdx == 0:
            baseDf = copy.deepcopy(subData)
        else:
            subData.drop(columns=['filename','retTime','compound','isotopolog','tic','TIC*IT','integTime','resolution','agcTarget','microscans'],axis = 1, inplace=True)

            baseDf = pd.merge_ordered(baseDf, subData,on='scanNumber',suffixes =(False,False))

        subIdx += 1

    #If there are no entries for a given scan, IsoX will not include that scan in the output. 
    #This fills in 0s for all entries between the minimum and maximum scan
    maxScan = baseDf['scanNumber'].max()
    minScan = baseDf['scanNumber'].min()
    baseDf = baseDf.set_index('scanNumber').reindex(range(minScan,maxScan)).fillna(0).reset_index()

    #Fill in additional 0s
    for subName in subNameList:
        baseDf.loc[baseDf['mass' + subName].isnull(), 'mass' + subName] = 0
        baseDf.loc[baseDf['intensity' + subName].isnull(), 'intensity' + subName] = 0
        baseDf.loc[baseDf['peakNoise' + subName].isnull(), 'peakNoise' + subName] = 0
        baseDf.loc[baseDf['counts' + subName].isnull(), 'counts' + subName] = 0 

    if MNRelativeAbundance: 
        baseDf = calc_MN_Rel_Abundance(baseDf, subNameList)
    else:
        baseDf = calc_Append_Ratios(baseDf, subNameList)

    if timeBounds: 
        baseDf = cull_By_Time(baseDf, timeBounds, scanNumber = scanNumber)

    if cullOn != None:
        upperLim = baseDf[cullOn].mean() + baseDf[cullOn].std() * cullAmt
        lowerLim = baseDf[cullOn].mean() - baseDf[cullOn].std() * cullAmt
        baseDf = baseDf[(baseDf[cullOn] > lowerLim) & (baseDf[cullOn] < upperLim)]
        
    combinedData['subNameList'] = subNameList
    combinedData['mergedDf'] = baseDf
    return combinedData

def processIsoXDf(IsoXDf, cullOn = None, cullAmt = 3, scanNumber = False, timeBounds = None,MNRelativeAbundance = False, subsToExclude = []):
    '''
    Takes in the IsoXDataframe; splits it based on filename. For each filename, combines the data from all substitutions into a single dataframe. Adds these dataframes to an output list. 

    Inputs:
        IsoXDf: A dataframe with the data from the .csv output.

    Outputs:
        mergedList: mergedList, is a list of dataframes. Each dataframe corresponds to one file & has all information about each scan on a single line. 
    '''
    thisFileData = IsoXDf
    
    #Check to see if there is any issue with multiple entries being found for the same substitution
    splitDf = thisFileData.groupby('isotopolog')
    
    for subName, subData in splitDf:
        #IsoX will return multiple values if multiple peaks are observed with closely similar masses. This prints a warning if this is the case and where. Not as useful as it could be, because low abundance (e.g. start of scan) periods will often have warnings as well. 
        g = splitDf.get_group(subName)['scanNumber'].value_counts()
        if len(g[g>1]):
            a=2
            #print("Multiple Entries Found for file " + thisFileName + " sub " + subName)
            #print(g[g>1])

     #For each isotopolog & each scan, select only the observation with highest intensity
    selectTopScan = thisFileData.loc[thisFileData.groupby(['isotopolog','scanNumber'])["intensity"].idxmax()]

    #Filter out substitutions based on a user-provided list
    filtered_df = selectTopScan[~selectTopScan['isotopolog'].isin(subsToExclude)]

    #sort by isotopolog to prepare for combine_substituted
    topScanDf = filtered_df.groupby('isotopolog')
    thisCombined = combine_Substituted_Peaks(topScanDf, cullOn = cullOn, cullAmt = cullAmt, scanNumber = scanNumber, timeBounds = timeBounds,MNRelativeAbundance = MNRelativeAbundance)

    return thisCombined

def output_Raw_File_MN_Rel_Abundance(mergedDf, subNameList, massStr = None):
    #Initialize output dictionary 
    rtnDict = {}
      
    if massStr == None:
        #Try to set massStr programatically; first look for unsub
        try:
            massStr = str(round(mergedDf['massUnsub'].median(),1))
        #If this fails, find the first substitution and use that instead. 
        except:
            mergedEntries = list(mergedDf.keys())
            subNameList = [x[6:] for x in mergedEntries if 'counts' in x]
            massStr = str(round(mergedDf['mass' + subNameList[0]].median(),1))
        
    rtnDict[massStr] = {}
    
    for sub in subNameList:
        rtnDict[massStr][sub] = {}
        rtnDict[massStr][sub]['MN Relative Abundance'] = np.mean(mergedDf['MN Relative Abundance ' + sub])
        rtnDict[massStr][sub]['StDev'] = np.std(mergedDf['MN Relative Abundance ' + sub])
        rtnDict[massStr][sub]['StError'] = rtnDict[massStr][sub]['StDev'] / np.power(len(mergedDf),0.5)
        rtnDict[massStr][sub]['RelStError'] = rtnDict[massStr][sub]['StError'] / rtnDict[massStr][sub]['MN Relative Abundance']
        rtnDict[massStr][sub]['tic'] = mergedDf['tic'].mean()
        rtnDict[massStr][sub]['TICVar'] = 0
        rtnDict[massStr][sub]['TIC*ITVar'] = 0
        rtnDict[massStr][sub]['TIC*ITMean'] = 0
                        
        a = mergedDf['counts' + sub].sum()
        b = mergedDf['total Counts'].sum()
        shotNoiseByQuad = np.power((1./a + 1./b), 0.5)
        rtnDict[massStr][sub]['ShotNoiseLimit'] = shotNoiseByQuad
        
    return rtnDict        

def output_Raw_File_Ratios(mergedDf, subNameList, mostAbundant = True, massStr = None, omitRatios = [], debug = True, MNRelativeAbundance = False):
    '''
    For each ratio of interest, calculates mean, stdev, SErr, RSE, and ShotNoise based on counts. 
    Outputs these in a dictionary which organizes by fragment (i.e different entries for fragments at 119 and 109).
    
    Inputs:
        df: A list of merged data frames from the _combineSubstituted function. Each dataframe constitutes one fragment.
        isotopeList: A list of isotopes corresponding to the peaks extracted by FTStat, in the order they were extracted. 
                    This must be the same for each fragment. This is used to determine all ratios of interest, i.e. 13C/UnSub, and label them in the proper order. 
        omitRatios: A list of ratios to ignore. I.e. by default, the script will report 13C/15N ratios, which one may not care about. 
        weightByNLHeight: Used to say whether to weight by NL height or not. Tim: 20210330: Added this as the weightByNLHeight does not yet work with M+N routines. 
        debug: Set false to suppress print output
        percentAbundance: Optionally calculate the percent abundance within each fragment, for M+N experiments. 
        
         
    Outputs: 
        A dictionary giving mean, stdev, StandardError, relative standard error, and shot noise limit for all peaks.  
    '''
    #Initialize output dictionary 
    rtnDict = {}
      
    #Adds the peak mass to the output dictionary
    if massStr == None:
        #Try to set massStr programatically; first look for unsub
        try:
            massStr = str(round(mergedDf['massUnsub'].median(),1))
        #If this fails, find the first substitution and use that instead. 
        except:
            mergedEntries = list(mergedDf.keys())
            subNameList = [x[6:] for x in mergedEntries if 'counts' in x]
            massStr = str(round(mergedDf['mass' + subNameList[0]].median(),1))
        
    rtnDict[massStr] = {}
        
    maxSub = findMostAbundantSub(mergedDf, subNameList)

    for sub1, sub2 in itertools.combinations(subNameList,2):
        if ((mostAbundant) and (sub1 != maxSub) and (sub2 != maxSub)): 
            continue
        if sub1 + '/' + sub2 in mergedDf:
            header = sub1 + '/' + sub2
        else:
            header = sub2 + '/' + sub1
          
        #perform calculations and add them to the dictionary     
        rtnDict[massStr][header] = {}
        rtnDict[massStr][header]['Ratio'] = np.mean(mergedDf[header])
        rtnDict[massStr][header]['StDev'] = np.std(mergedDf[header])
        rtnDict[massStr][header]['StError'] = rtnDict[massStr][header]['StDev'] / np.power(len(mergedDf),0.5)
        rtnDict[massStr][header]['RelStError'] = rtnDict[massStr][header]['StError'] / rtnDict[massStr][header]['Ratio']
        
        a = mergedDf['counts' + sub1].sum()
        b = mergedDf['counts' + sub2].sum()
        shotNoiseByQuad = np.power((1./a + 1./b), 0.5)
        rtnDict[massStr][header]['ShotNoiseLimit'] = shotNoiseByQuad

        averageTIC = np.mean(mergedDf['tic'])
        valuesTIC = mergedDf['tic']
        rtnDict[massStr][header]['tic'] = averageTIC
        rtnDict[massStr][header]['TICVar'] = np.sqrt(
            np.mean((valuesTIC-averageTIC)**2))/np.mean(valuesTIC)

        averageTICIT = np.mean(mergedDf['TIC*IT'])
        valuesTICIT = mergedDf['TIC*IT']
        rtnDict[massStr][header]['TIC*ITMean'] = averageTICIT
        rtnDict[massStr][header]['TIC*ITVar'] = np.sqrt(
            np.mean((valuesTICIT-averageTICIT)**2))/np.mean(valuesTICIT)
                        
    return rtnDict

def read_IsoX_Folder(isoXFileNames, cullOn = None, cullAmt = 3, scanNumber = False, timeBounds = None, MNRelativeAbundance = False, zeroCountsScreen = True, zeroCountsThreshold = 0, peakDriftScreen = True, peakDriftThreshold = 2, exclusionDict = {}):
    mergedDictByFrags = {}

    for isoXFileName in isoXFileNames:
        print("Reading in " + str(isoXFileName))
        thisIsoX = readIsoX(isoXFileName)
        fragKeys = set(thisIsoX['compound'])

        for fragKey in fragKeys:
            if str(fragKey) not in mergedDictByFrags:
                mergedDictByFrags[str(fragKey)] = {}

            if str(fragKey) in exclusionDict:
                thisSubsToExclude = exclusionDict[str(fragKey)]
            else:
                thisSubsToExclude = []

            thisFragDf = thisIsoX[thisIsoX['compound'] == fragKey]
            thisMergedDict = processIsoXDf(thisFragDf, cullOn = cullOn, cullAmt = cullAmt, scanNumber = scanNumber, timeBounds = timeBounds, MNRelativeAbundance = MNRelativeAbundance, subsToExclude = thisSubsToExclude)

            if zeroCountsScreen:
                dataScreenIsoX.zeroCountsScreen(thisMergedDict, isoXFileName, fragKey, threshold = zeroCountsThreshold)

            if peakDriftScreen:
                dataScreenIsoX.peakDriftScreen(thisMergedDict, isoXFileName, fragKey, threshold = peakDriftThreshold)

            mergedDictByFrags[str(fragKey)][isoXFileName] = thisMergedDict

    return mergedDictByFrags

def calc_Folder_Output(isoXFileNames, smpStdOrdering = None, cullOn = None, cullAmt = 3, debug = False, scanNumber = False, timeBounds = None, MNRelativeAbundance = False, RSESNScreen = True, zeroCountsScreen = True, zeroCountsThreshold = 0, peakDriftScreen = True, peakDriftThreshold = 2, exclusionDict = {}):
    '''
    Calculates the output for many isoX files (NOT combined.isox. Files should be processed individually). 

    Inputs:
        isoXFilePaths: A list of .isox files, including both standard and sample files, to read in. 
        smpStdOrdering: A user specified list which explicitly indicates whether each file is a sample or standard. 
        cullOn: 
        cullAmt: 
        debug: 
        cullByTime: 
        scanNumber:
        timeBounds:
        MNRelativeAbundance:

    Outputs: 
        Two different versions of the output. These are:
        rtnAllFilesDF: A dataframe containing the processed data means, stdevs, serrs, rses, and SN. 
        mergedDict: A dictionary, containing 'subNameList' (metadata used in computations) and 'mergedDf', a dataframe containing the scan-by-scan data for that file. 
    '''    
    mergedDictByFrags = read_IsoX_Folder(isoXFileNames, cullOn = cullOn, cullAmt = cullAmt, scanNumber = scanNumber, timeBounds = timeBounds, MNRelativeAbundance = MNRelativeAbundance, zeroCountsScreen=zeroCountsScreen, zeroCountsThreshold=zeroCountsThreshold, peakDriftScreen=peakDriftScreen, peakDriftThreshold=peakDriftThreshold, exclusionDict=exclusionDict)
    
    rtnAllFilesDF = []
    for fragKey, fragData in mergedDictByFrags.items():
        for thisFileIdx, (thisFileName, thisFileData) in enumerate(fragData.items()):
            if smpStdOrdering == None:
                thisFileSmpStd = 'N/A'
            else:
                thisFileSmpStd = smpStdOrdering[thisFileIdx]

            thisMergedDf = thisFileData['mergedDf']
            thisSubNameList = thisFileData['subNameList']

            if MNRelativeAbundance:
                header = ["FileName", "Fragment", "MN Relative Abundance", "Average", "StdDev", "StdError", "RelStdError",'ShotNoise','Tic','TicVar', 'File Type']
                thisFileOutput = output_Raw_File_MN_Rel_Abundance(thisMergedDf, thisSubNameList, massStr=fragKey)

            else:
                header = ["FileName", "Fragment", "IsotopeRatio", "Average", "StdDev", "StdError", "RelStdError",'ShotNoise','Tic','TicVar', 'File Type']
                thisFileOutput = output_Raw_File_Ratios(thisMergedDf, thisSubNameList, massStr=fragKey)

            currentOutput = thisFileOutput[fragKey]
            for subKey, subData in currentOutput.items():
            #add subkey to each separate df for isotope specific 
                if MNRelativeAbundance:
                    thisRVal = subData["MN Relative Abundance"]
                else: 
                    thisRVal = subData["Ratio"]

                thisStdDev = subData["StDev"]
                thisStError = subData["StError"] 
                thisRelStError = subData["RelStError"]
                thisShotNoise = subData["ShotNoiseLimit"]
                thisTic = subData['tic']
                thisTicVar = subData['TICVar']
                thisRow = [thisFileName, fragKey, subKey, thisRVal, thisStdDev, thisStError,thisRelStError, thisShotNoise, thisTic, thisTicVar, thisFileSmpStd] 
                rtnAllFilesDF.append(thisRow)

    rtnAllFilesDF = pd.DataFrame(rtnAllFilesDF)

    # set the header row as the df header
    rtnAllFilesDF.columns = header 

    if MNRelativeAbundance:
        rtnAllFilesDF = rtnAllFilesDF.sort_values(by=['Fragment', 'MN Relative Abundance'], axis=0, ascending=True)

    else:
        #sort by fragment and isotope ratio, output to csv
        rtnAllFilesDF = rtnAllFilesDF.sort_values(by=['Fragment', 'IsotopeRatio'], axis=0, ascending=True)

    if RSESNScreen:
        dataScreenIsoX.RSESNScreen(rtnAllFilesDF, MNRelativeAbundance = MNRelativeAbundance)

    return rtnAllFilesDF, mergedDictByFrags

def folderOutputToDict(rtnAllFilesDF,MNRelativeAbundance = False):
    '''
    takes the output dataframe and processes to a dictionary, for .json output
    '''

    sampleOutputDict = {}
    fragmentList = []
    for i, info in rtnAllFilesDF.iterrows():
        fragment = info['Fragment']
        file = info['FileName']

        if MNRelativeAbundance:
            ratio = info['MN Relative Abundance']
        else:
            ratio = info['IsotopeRatio']
        avg = info['Average']
        std = info['StdDev']
        stderr = info['StdError']
        rse = info['RelStdError']
        SN = info['ShotNoise']

        if file not in sampleOutputDict:
            sampleOutputDict[file] = {}

        if fragment not in sampleOutputDict[file]:
            sampleOutputDict[file][fragment] = {}

        sampleOutputDict[file][fragment][ratio] = {'Average':avg,'StdDev':std,'StdError':stderr,'RelStdError':rse,'ShotNoise':SN}
        
    return sampleOutputDict

def processIndividualAndAverageIsotopeRatios(fragmentFolderPaths, cwd, outputToCSV=False, csvOutputPath = 'output.csv', file_extension = '.isox', processed_data_subfolder='Processed Data', time_bounds = (0,0), RSESNScreen = True, zeroCountsScreen = True, zeroCountsThreshold = 0, peakDriftScreen = True, peakDriftThreshold = 2, exclusionDict = {}):
    '''
    Process statistics on isox files and output processed data results. Prepare data to run M+1 model. If you have multiple input files, it will take their average relative standard error and divide this by the square root of the number of files to use for future computations. 

    Inputs:
        fragmentFolderPaths: A list of strings representing the paths of subfolders within the main folder. Output of organizeData.get_subfolder_paths()
        cwd: The current working directory. 
        outputToCSV: If True, output the processed data as a .csv. 
        csvOutputPath: Specifies the output path for that .csv. 
        file_extension: The extension of the input data files. Should be .isox. 
        processed_data_subfolder: Directory name for processed data. 
        acquisition_length: Cull by time, within these bounds. 
    
    Outputs:
        rtnMeans: A dataframe containing information about the mean sample and standard values for each isotope of each fragment. 
        rtnMergedDict: A dictionary, where keys are fileNames, and values are dictionaries. The inner dictionaries contain 'subNameList' (metadata used in computations) and 'mergedDf', a dataframe containing the scan-by-scan data for that file. 
    '''
    #Initialize outputs
    allDataReturnedDFList = []
    MN_RELATIVE_ABUNDANCE = False
    allSortedMeanIsotopeRatios = []
    allMergedDict = []
    
    #Iterate through each folder
    for thisFolder in fragmentFolderPaths:
        #Get the filenames for this folder
        isoXFileNames, smpStdOrdering = organizeData.get_file_paths_in_subfolders(thisFolder, file_extension)
        thisFolderName = os.path.basename(thisFolder) 

        #Determine whether we are processing M+1 data or molecular average data
        if thisFolderName == 'full_molecular_average':
            MN_RELATIVE_ABUNDANCE = False
        else:
            MN_RELATIVE_ABUNDANCE = True

        #Compute output and append to lists
        rtnAllFilesDF, mergedDict = calc_Folder_Output(isoXFileNames, smpStdOrdering = smpStdOrdering, cullOn = None, cullAmt = 3, debug = False, scanNumber = False, timeBounds = time_bounds, MNRelativeAbundance = MN_RELATIVE_ABUNDANCE, RSESNScreen = RSESNScreen, zeroCountsScreen = zeroCountsScreen, zeroCountsThreshold = zeroCountsThreshold, peakDriftScreen = peakDriftScreen, peakDriftThreshold = peakDriftThreshold, exclusionDict = exclusionDict)
        allDataReturnedDFList.append(rtnAllFilesDF)
        allMergedDict.append(mergedDict)
    
        #Compute means of sample and stanadrd
        if  MN_RELATIVE_ABUNDANCE == True:
            thisSortedAverageDF = rtnAllFilesDF.sort_values(by=['Fragment','MN Relative Abundance', 'File Type'])
            nFiles = thisSortedAverageDF.groupby(['Fragment','MN Relative Abundance', 'File Type'])['RelStdError'].transform('size')[0]
            means = thisSortedAverageDF.groupby(['Fragment','MN Relative Abundance', 'File Type']).mean(numeric_only = True).reset_index()
            means['RelStdError'] /= np.sqrt(nFiles)
            means['StdError'] /= np.sqrt(nFiles)

        else:
            thisSortedAverageDF = rtnAllFilesDF.sort_values(by=['IsotopeRatio', 'File Type'])
            nFiles = thisSortedAverageDF.groupby(['IsotopeRatio', 'File Type'])['RelStdError'].transform('size')[0]
            means = thisSortedAverageDF.groupby(['IsotopeRatio', 'File Type']).mean(numeric_only = True).reset_index()
            means['RelStdError'] /= np.sqrt(nFiles)
            means['StdError'] /= np.sqrt(nFiles)
            means['Fragment'] = thisFolderName
        
        #add to output list
        allSortedMeanIsotopeRatios.append(means)
    
    #Store scan-by-scan data for future use.
    rtnMergedDict = {}
    for d in allMergedDict:
        rtnMergedDict.update(d)

    #Add the means from each folder to a combined dataframe
    rtnData = pd.concat(allDataReturnedDFList, ignore_index=True)
    rtnMeans = pd.concat(allSortedMeanIsotopeRatios, ignore_index=True)

    #Output data frame to csv
    if outputToCSV:
        rtnData = rtnData.sort_values(by=['Fragment'], axis=0, ascending=True)
        rtnData.to_csv(str(cwd) + '/Processed Data/'+ csvOutputPath, index = False, header=True)
    
    return rtnMeans, rtnMergedDict

def prepareDataForM1(rtnMeans):
    '''
    Processes the data from rtnMeans into a format that can be read into the M+1 algorithm. 

    Inputs:
        rtnMeans: A dataframe containing information about the mean sample and standard values for each isotope of each fragment. The output of processIndividualAndAverageIsotopeRatios.

    Outputs:
        forM1Algo: A dictionary, structured like {'Std':{'44':'13C':{'Average':avg,      StdDev':std,'StdError':stderr,'RelStdError':rse,'ShotNoise':SN
        }}}

        The parent dictionary contains 'Std' and 'Smp'. The final dictionary contains data about that acquisition. 
    '''
    forM1Algo = {'Std':{},'Smp':{}}

    # Group by columns 'Fragment' and 'MN Relative Abundance' and iterate over groups
    organizeBy = ['File Type', 'Fragment','MN Relative Abundance']
    for (smpStd, fragment, isotope), group_df in rtnMeans.groupby(organizeBy):
        if fragment == 'full_molecular_average':
            continue
        # Extract values from the group
        thisIsotopeData = group_df.drop(columns=organizeBy).to_dict(orient='records')[0]
        
        # Assign values to the nested dictionary
        forM1Algo[smpStd].setdefault(fragment, {})[isotope] = thisIsotopeData

    return forM1Algo
