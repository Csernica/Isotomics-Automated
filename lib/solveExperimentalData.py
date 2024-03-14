import os
import re

import numpy as np

import readInput as ri
import solveSystem as ss
import readCSVAndSimulate as sim
import basicDeltaOperations as op
import fragmentAndSimulate as fas
import organizeData


def experimentalDataM1(rtnMeans, cwd, MOLECULE_INPUT_PATH, UValue = '13C/Unsub', UOrbiAndErr = False, processFragKeys = {'44':'44','full_relative_abundance':'full'}, perturbTheoryOAmt = 0.001, MonteCarloN = 100, outputPrecision = 3, resultsFileName = 'M1Output.csv'):
    #GET FORWARD MODEL STANDARD
    initializedMolecule = sim.moleculeFromCsv(os.path.join(cwd, MOLECULE_INPUT_PATH), deltas = [0] * 6)
    mDf = initializedMolecule['molecularDataFrame']
    predictedMeasurement, MNDict, fractionationFactors = sim.simulateMeasurement(initializedMolecule, massThreshold = 5)

    #GET U VALUE
    UValuesSmp = getUVal(rtnMeans, initializedMolecule, UValue = UValue, UOrbiAndErr = UOrbiAndErr)

    #GET M+1 DATA
    preparedData = organizeData.prepareDataForM1(rtnMeans)
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

    return cleanExperimentalOutput

def getUVal(rtnMeans, initializedMolecule, UValue = '13C/Unsub', UOrbiAndErr = False):
    '''
    
    '''
    #Various representations of U_VAL; rare_sub gives '13C', U_ValID gives 'C'
    rare_sub = UValue.split('/')[0]
    U_ValID = re.sub(r'\d+', '', rare_sub)

    #Optionally pass a tuple, in which case return it directly
    if UOrbiAndErr:
        return {rare_sub:{'Observed': UOrbiAndErr[0], 'Error': UOrbiAndErr[0] * UOrbiAndErr[1]}}

    #Filter by the U Value of interest and perform a sample standard comparison.
    MAData = rtnMeans[rtnMeans['Fragment'] == 'full_molecular_average']
    filtered_data = MAData[MAData['IsotopeRatio'] == UValue]
    grouped = filtered_data.groupby('File Type')[['Average', 'RelStdError']].first()
    thisUVal = 1000 * (grouped.loc['Smp', 'Average'] / grouped.loc['Std', 'Average'] - 1)
    thisURelativeErr = np.sqrt((grouped.loc['Smp', 'RelStdError'])**2 + (grouped.loc['Std', 'RelStdError'])**2)

    #Convert the delta value to concentration, multiply by 3 to go from R to U value
    mDf = initializedMolecule['molecularDataFrame']
    nAtomsThisU = mDf[mDf['IDS'] == U_ValID]['Number'].sum()
    UOrbi = op.concentrationToM1Ratio(op.deltaToConcentration(U_ValID,thisUVal)) * nAtomsThisU
    UValuesSmp = {rare_sub:{'Observed': UOrbi, 'Error': UOrbi * thisURelativeErr}}

    return UValuesSmp
