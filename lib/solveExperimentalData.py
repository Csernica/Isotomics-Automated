import os
import re

import numpy as np

import readInput as ri
import solveSystem as ss
import readCSVAndSimulate as sim
import basicDeltaOperations as op
import fragmentAndSimulate as fas
import dataAnalyzerMNIsoX

def defineProcessFragKeys(fragmentationDictionary, fragSubset = None):
    '''
    processFragKeys is a dictionary used to help read in experimental data. The idea is, your experimental data is read in with some fragment keys (e.g., '44', 'full_relative_abundance') which may not correspond to those used for simulated data, and processFragKeys converts between the two. Based on the way we set this up in Isotomics-Automated they should match except for the 'full_relative_abundance'. This fills them in. 

    Inputs:
        fragmentationDictionary: Basic information about the fragmentation, from the input CSV. 

    Outputs:
        processFragKeys: A dictionary, where fragment keys for experimental data are keyed to fragment keys for simulated data. 
    '''
    processFragKeys = {'full_relative_abundance':'full'}
    for fragKey in fragmentationDictionary:
        if fragKey != 'full':
            if fragSubset == None:
                processFragKeys[fragKey] = fragKey
            elif fragKey in fragSubset:
                processFragKeys[fragKey] = fragKey
                

    return processFragKeys

def experimentalDataM1(rtnMeans, cwd, MOLECULE_INPUT_PATH, std_deltas, UValue = '13C/Unsub', mAObs = None, mARelErr = None, perturbTheoryOAmt = 0.001, MonteCarloN = 1000, outputPrecision = 3, resultsFileName = 'M1Output.csv', plot = True, fragSubset = None):
    '''
    Parent function to process experimental M+1 Data and return results. 

    Inputs:
        rtnMeans: A dataframe, returned via dataAnalyzerMNIsoX.processIndividualAndAverageIsotopeRatios, which gives the mean values of samples and standards. 
        cwd: The current working directory. 
        MOLECULE_INPUT_PATH: A path to the .csv file used to initialize data about your molecule. 
        std_deltas: Known or approximated delta values for the standard, used for forward model standardization. 
        UValue, mAObs, mARelErr: Used to set the molecular average U Value of the sample used to convert M+1 relative abundancees to site-specific deltas. See getUVal for details. 
        perturbTheoryOAmt: Used for the 'observed abundance correction' to low abundance peaks. See the appendix to Csernica and Eiler 2023 for details about this correction. A typical value of 0.001 (1 per mil) is a good estimate. Test different values with simulated data to see if your choice is appropriate. 
        MonteCarloN: The number of iterations used for the M+1 monte carlo solver. 
        outputPrecision: The number of decimals to include in the output .csv. 
        resultsFileName: Filename to export results to. 
        plot: If True, return a plot. 

    Outputs: 
        cleanExperimentalOutput: A dataframe containing basic information about the molecule and the results of the M+1 algorithm. 
        Also returns a plot (if plot) and an output .csv (of cleanExperimentalOutput)
    '''
    #GET FORWARD MODEL STANDARD
    initializedMolecule = sim.moleculeFromCsv(os.path.join(cwd, MOLECULE_INPUT_PATH), deltas = std_deltas, fragSubset = fragSubset)
    processFragKeys = defineProcessFragKeys(initializedMolecule['fragmentationDictionary'])
    mDf = initializedMolecule['molecularDataFrame']
    predictedMeasurement, MNDict, fractionationFactors = sim.simulateMeasurement(initializedMolecule, massThreshold = 5)

    #GET U VALUE
    UValuesSmp = getUVal(rtnMeans, initializedMolecule, UValue = UValue, mAObs = mAObs, mARelErr = mARelErr)

    #GET M+1 DATA
    preparedData = dataAnalyzerMNIsoX.prepareDataForM1(rtnMeans)
    replicateData = ri.readObservedData(preparedData, theory = predictedMeasurement,
                                                    standard = [True, False],
                                                    processFragKeys = processFragKeys)
    
    #Generate observed abundance ('O') correction factors
    OValueCorrection = ss.OValueCorrectTheoretical(predictedMeasurement,
                                                replicateData['Smp'],
                                                massThreshold = 1)

    isotopologuesDict = fas.isotopologueDataFrame(MNDict, mDf)

    rare_sub = UValue.split('/')[0]
    #Run the M+1 algorithm and process the results
    M1Results = ss.M1MonteCarlo(replicateData['Std'], replicateData['Smp'], 
                                OValueCorrection, 
                                isotopologuesDict, 
                                initializedMolecule['fragmentationDictionary'], 
                                N = MonteCarloN,
                                perturbTheoryOAmt = perturbTheoryOAmt, 
                                disableProgress = True)
    
    processedResults = ss.processM1MCResults(M1Results, UValuesSmp,
                                             isotopologuesDict, 
                                             mDf, 
                                             UMNSub = [rare_sub],
                                             disableProgress = True)

    mDf = ss.updateSiteSpecificDfM1MC(processedResults, mDf)
    #END M1 ALGORITHM

    #output nice csv
    cleanExperimentalOutput = mDf.drop(['deltas','M1 M+N Relative Abundance', 'M1 M+N Relative Abundance Error', 'UM1','UM1 Error','Calc U Values','Calc U Values Error'], axis=1, inplace=False)
    #Round off decimals
    toRound = ['VPDB etc. Deltas', 'VPDB etc. Deltas Error', 'Relative Deltas','Relative Deltas Error']
    cleanExperimentalOutput[toRound] = cleanExperimentalOutput[toRound].round(decimals=outputPrecision) 
    cleanExperimentalOutput.to_csv(resultsFileName, index=False)

    if plot:
        sim.plotOutput(mDf)

    return cleanExperimentalOutput

def getUVal(rtnMeans, initializedMolecule, UValue = '13C/Unsub', mAObs = None, mARelErr = None):
    '''
    Computes the molecular Average U Value for the sample used to convert M+1 relative abundances to site-specific deltas. Either calculates this from Orbitrap data or takes it as an input (e.g., if it is constrained externally via EA). 

    Inputs:
        rtnMeans: A dataframe, returned via dataAnalyzerMNIsoX.processIndividualAndAverageIsotopeRatios, which gives the mean values of samples and standards. 
        initializedMolecule: A read in molecular csv containing basic details about the molecule. 
        UValue, mAObs, mARelErr: Information about the SAMPLE's ratio used to convert from M+1 relative abundances to site-specific deltas. mAObs is the delta value. mARelErr is the error. For example, if I have an EA measurement OF MY SAMPLE of 13C/Unsub that gives -30.0 +/- 0.1 per mil vs VPDB, then I set: UValue = '13C/Unsub', mAObs = -30.0, mARelErr = 0.1. If I do not have this information for the sample, I can use Orbitrap data. In this case, I must include appropriate Orbitrap files in 'full_molecular_average', specify the substitution I want to use, and include 'None' for mAObs and mARelErr. 

    Outputs: 
        UValuesSmp: A dictionary. The first key is the rare substitution used (e.g., '13C'.) Values are a dictionary. In that dictionary, you have 'Observed' and 'Error'. 'Observed' gives the U Value as a ratio (e.g., 0.01118 for 13C/Unsub for a molecule with 1 carbon, or 2*0.01118 if there are two carbons, etc.) and 'Error' gives the error on that ratio (e.g., 0.00001118 if there is 1 per mil error). 
    '''
    #Various representations of U_VAL; rare_sub gives '13C', U_ValID gives 'C'
    rare_sub = UValue.split('/')[0]
    U_ValID = re.sub(r'\d+', '', rare_sub)
    mDf = initializedMolecule['molecularDataFrame']

    #Optionally pass a tuple, in which case return it directly
    if mAObs:
        nAtomsThisU = mDf[mDf['IDS'] == U_ValID]['Number'].sum()
        UOrbi = op.concentrationToM1Ratio(op.deltaToConcentration(U_ValID,mAObs)) * nAtomsThisU

        if mARelErr == None:
            print("No error given for molecular average measurement, defaulting to 0")
            return {rare_sub:{'Observed': UOrbi, 'Error': UOrbi * 0}}
        else:
            return {rare_sub:{'Observed': UOrbi, 'Error': UOrbi * mARelErr / 1000}}

    #Get actual standard enrichment vs reference frame
    stdRatioVSRef = getAverageRatio(mDf, U_ValID)
    nAtomsThisU = mDf[mDf['IDS'] == U_ValID]['Number'].sum()
    
    #Filter by the U Value of interest and perform a sample standard comparison.
    MAData = rtnMeans[rtnMeans['Fragment'] == 'full_molecular_average']
    filtered_data = MAData[MAData['IsotopeRatio'] == UValue]
    grouped = filtered_data.groupby('File Type')[['Average', 'RelStdError']].first()
    thisUVal = stdRatioVSRef * grouped.loc['Smp', 'Average'] / grouped.loc['Std', 'Average']
    thisUVal *= nAtomsThisU
    thisURelativeErr = np.sqrt((grouped.loc['Smp', 'RelStdError'])**2 + (grouped.loc['Std', 'RelStdError'])**2)

    #Convert the delta value to concentration, multiply by 3 to go from R to U value
    UValuesSmp = {rare_sub:{'Observed': thisUVal, 'Error': thisUVal * thisURelativeErr}}

    return UValuesSmp

def getAverageRatio(molecularDataFrame, element):
    #pull out this element only
    thisElementDf = molecularDataFrame[molecularDataFrame['IDS'] == element]
    numbers = thisElementDf['Number']
    deltas = thisElementDf['deltas']

    #Loop through sites. Get the concentration of each and add to total
    totalConc = [0,0,0,0]
    for siteIdx, siteDelta in enumerate(deltas):
        thisConc = np.array(op.deltaToConcentration(element,siteDelta))
        totalConc += thisConc * numbers[siteIdx]

    #convert total concentration to ratio and delta
    avgRatio = op.concentrationToM1Ratio(totalConc)

    return avgRatio