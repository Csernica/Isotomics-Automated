import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import dataAnalyzerMNIsoX as dA
from tqdm import tqdm
    
SUB_MASS_DICT = {'d':1.00627674587,'2h':1.00627674587,'15n':0.997034886,'13c':1.003354835,
                   'unsub':0,'33s':0.999387735,'34s':1.995795825,'36s':3.995009525,
                  '18o':2.0042449924,'17o':1.0042171364}

def RSESNScreen(rtnAllFilesDF, MNRelativeAbundance = False, threshold = 2):
    '''
    Screen all peaks and print any which have RSE/SN > threshold

    Inputs:
        rtnAllFilesDf: A dataframe containing all of the output ratios, from dataAnalyzerMN.calc_Folder_Output

    Outputs:
        None. Prints flags for peaks that exceed the threshold. 
    '''
    MN_Relative_Abundance = True
    threshold = 2
    rtnAllFilesDF['RSE over SN'] = rtnAllFilesDF['RelStdError'] / rtnAllFilesDF['ShotNoise']
    subColumn = 'MN Relative Abundance' if MNRelativeAbundance else 'IsotopeRatio'
    failSN = rtnAllFilesDF.loc[rtnAllFilesDF['RSE over SN'] > threshold]

    for index, row in failSN.iterrows():
        fileName = row['FileName']
        fragment = row['Fragment']
        sub = row[subColumn]
        thisRSEOverSN = row['RSE over SN']
        print(fileName + ' ' + fragment + ' ' + sub + " fails RSE/SN test with a value of " + f'{thisRSEOverSN:.3f}')

def zeroCountsScreen(mergedDict, fileName, fragKey, threshold = 0):
    '''
    Iterates through all peaks and prints those with zero counts higher than a certain relative threshold.

    Inputs: 
        folderPath: The directory containing the .txt or .csv files from FTStatistic.
        fragmentDict: A dictionary containing information about the fragments and substitutions present in the FTStat output file. 
        mergedList: A list of lists; each outer list corresponds to a file, each inner list to a fragment; elements of this inner list are dataframes giving the scans and data for that fragment of that file. 
        fileExt: The file extension of the FTStat output file, either '.txt' or '.csv'
        threshold: The relative number of zero scans to look for

    Outputs:
        None. Prints the name of peaks with more than the threshold number of zero scans. 
    '''
    subNameList = mergedDict['subNameList']
    thisMergedDf = mergedDict['mergedDf']

    for iso in subNameList:
        thisIsoFileZeros = len(thisMergedDf['counts' + iso].values) - np.count_nonzero(thisMergedDf['counts' + iso].values)

        thisIsoFileZerosFraction = thisIsoFileZeros / len(thisMergedDf['counts' + iso])
        if thisIsoFileZerosFraction > threshold:
            print('Zero Counts: ' + fileName + ' ' + str(fragKey) + ' ' + iso + ' has ' + str(thisIsoFileZeros) + ' zero scans, out of ' + str(len(thisMergedDf['counts' + iso])) + ' scans (' + str(thisIsoFileZerosFraction) + ')') 

def getThisSubMass(subKey):
    '''
    Given a subKey, like '13C', (or 'D-13C', split by '-', for multiple substitutions), get the mass change relative to no substitutions. 
    '''
    computedMass = 0
    thisSubs = subKey.split('-')

    #Find the increase in mass due to substitutions
    for sub in thisSubs:
        try:
            computedMass += SUB_MASS_DICT[sub.lower()]
        except:
            print("Could not look up substitution " + sub + " correctly.")
            computedMass += 0

    return computedMass
    
def peakDriftScreen(mergedDict, fileName, fragKey, threshold = 2):
    '''
    Check each peak to determine if it has drifted relative to the most abundant isotope. Threshold is the amount of drift (in ppm) required to raise a warning.  
    '''
    subNameList = mergedDict['subNameList']
    mergedDf = mergedDict['mergedDf']
  
    mostAbundantIso = dA.findMostAbundantSub(mergedDf, subNameList)
    idxMostAbundant = subNameList.index(mostAbundantIso)

    observedMasses = []
    for iso in subNameList:
        observedMassIso = mergedDf[mergedDf['mass' + iso]!=0]['mass' + iso].mean()
        observedMasses.append(observedMassIso)

    mostAbundantMassObs = observedMasses[idxMostAbundant]
    mostAbundantMassTheory = getThisSubMass(mostAbundantIso)
    
    for thisIsoIdx, thisIso in enumerate(subNameList):
        if thisIso != mostAbundantIso:
            thisMassTheory = getThisSubMass(thisIso)
            thisMassObs = observedMasses[thisIsoIdx]

            #compute observed and theoretical mass differences
            massDiffTheory = thisMassTheory - mostAbundantMassTheory
            massDiffObserve = thisMassObs - mostAbundantMassObs

            peakDrift = np.abs(massDiffObserve - massDiffTheory)

            peakDriftppm = peakDrift / observedMassIso * 10**6

            if peakDriftppm > threshold:
                print("Peak Drift Observed for " + fileName + " " + str(fragKey) + ' ' + thisIso + " with size " + str(peakDriftppm))

def subsequenceOutlierDetection(timeSeries, priorSubsequenceLength = 1000, testSubsequenceLength = 1000):
    '''
    Calculates the anomaly score for a timeseries and a given subsequence length. 

    Inputs:
        timeSeries: A univariate time series, i.e. a pandas series.
        subsequenceLength: The length of the subsequence to use..

    Outputs:
        allDev: The euclidian distance between the subsequence of interest and the mean of the previous subsequenceLength observations.
    '''
    allDev = []
    for i in tqdm(range(priorSubsequenceLength,len(timeSeries)-testSubsequenceLength)):
        thisSubsequence = timeSeries[i:i+testSubsequenceLength]
        thisPrediction = timeSeries[i-priorSubsequenceLength:i].mean()
        meanZScore = np.abs(((thisSubsequence.values - thisPrediction) / thisSubsequence.std()).mean())
        allDev.append(meanZScore)

    return np.array(allDev)

def internalStabilityScreenSubsequence(mergedDict, MNRelativeAbundance = False, priorSubsequenceLength = 1000, testSubsequenceLength = 1000, thresholdConstant = 0.2):
    '''
    Screens all peaks for subsequence outlier detection; prints those that fail.

    Inputs: 
        folderPath: The directory containing the .txt or .csv files from FTStatistic.
        fragmentDict: A dictionary containing information about the fragments and substitutions present in the FTStat output file. 
        fragmentMostAbundant: A list giving the most abundant peak in each fragment.
        mergedList: A list of lists; each outer list corresponds to a file, each inner list to a fragment; elements of this inner list are dataframes giving the scans and data for that fragment of that file. 
        MNRelativeAbundance: Whether to calculate results as MN Relative Abundances or ratios.
        fileExt: The file extension of the FTStat output file, either '.txt' or '.csv'
        subsequenceLength: The length of the subsequence to use. 
        thresholdConstant: A constant to multiply by to determine the appropriate threshold for screening. 

    Outputs:
        No output; prints the identity of any file failing the test. 
    '''
    firstKey = list(mergedDict.keys())[0]
    subNameList = mergedDict[firstKey]['subNameList']
  
    for iso in subNameList:
        print(iso)
        for fileIdx, (fileName, fileData) in enumerate(mergedDict.items()):
            thisMergedDf = fileData['mergedDf']

            mostAbundantSub = dA.findMostAbundantSub(thisMergedDf, subNameList)
            if MNRelativeAbundance:
                series = thisMergedDf['MN Relative Abundance ' + iso]
            elif iso != mostAbundantSub:
                series = thisMergedDf[iso + '/' + mostAbundantSub]
            else:
                continue

            allDev = subsequenceOutlierDetection(series, priorSubsequenceLength = priorSubsequenceLength, testSubsequenceLength = testSubsequenceLength)
            if max(allDev) > thresholdConstant:
                print("Failed Subsequence Detection " +  fileName + " " + iso + " with a value of " + "{:.2f}".format(max(allDev)))


def visualizeTICVersusTime(fileName, mergedDict, targetRatio='13C/Unsub',scan_averaging_number = 50):
    '''
    Visualize the TIC versus time for a single file in the merged dict
    '''
    fig, axes = plt.subplots(nrows = 1, ncols = 2, figsize = (9*0.8,3*0.8), dpi = 300, gridspec_kw={'width_ratios': [2, 1]})
    thisDf = mergedDict[fileName]['mergedDf']
    series = thisDf[targetRatio]
    l = len(series)
    serrs = []

    cAx = axes[0]
    movingAvg = []
    for i in tqdm(range(scan_averaging_number,len(series)-scan_averaging_number)):
        current = series[i-scan_averaging_number:i+scan_averaging_number]
        mean = current.mean()
        movingAvg.append(mean)

    cAx.scatter(range(len(series)),series, s = 8, color = 'darkgreen', alpha = 0.3)
    cAx.plot(range(50,len(movingAvg)+50),movingAvg, label = "Moving Average of" + str(scan_averaging_number*2)+ "Scans", color = 'tab:orange')
    cAx.set_ylabel(targetRatio, fontsize = 16)
    cAx.set_xlabel("Scan Number")
    cAx.legend(loc = 'upper left')

    cAx = axes[1]
    cAx.hist(thisDf[thisDf['massUnsub']!=0]['massUnsub'], 30, density=False, facecolor='w',edgecolor = 'k', alpha=1, orientation=u'horizontal')

    cAx.yaxis.set_label_position("right")
    cAx.set_xlabel("Number of Scans")

    cAx.legend()

    cAx.spines['right'].set_visible(False)
    cAx.spines['top'].set_visible(False)

    cAx.spines['top'].set_visible(False)
    cAx.spines['right'].set_visible(False)
    
    plt.tight_layout()