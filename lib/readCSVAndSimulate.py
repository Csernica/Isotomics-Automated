import pandas as pd
import fragmentAndSimulate as fas
import readInput as ri
import solveSystem as ss
import calcIsotopologues as ci
import matplotlib.pyplot as plt
import basicDeltaOperations as op

import json
import copy

def process_fragment_column(col, fragSubset, molecularDf, fragmentationDictionary, renameCols):
    '''
    Read and process a single fragment column. Helper function for moleculeFromCsv.

    Inputs: 
        col: A string, format like 'Fragment 44' or 'Fragment 44_01'.
        fragSubset: A list of fragment names, like ['44','full']
        molecularDf: The molecular DataFrame
        fragmentationDictionary: A dictionary specifying the fragments and their geometries
        renameCols: A dictionary, where keys are the column names in the input csv and values are new desired column names for the molecular Dataframe

    Outputs:
        fragmentationDictionary, molecularDf, renameCols

        All as above, but either filled with more information for frags in fragSubset (fragmentationDictionary and renameCols) or with information removed (molecularDf, if there are fragments not in 'fragSubset').     
    '''
    #Get just the fragment name
    newName = col.split(' ')[1]
    #fragments could be given as 44_01, 44_02, or just as 44; split into cases
    if '_' not in newName:
        #If the fragment is desired, fill in the dictionary
        if newName in fragSubset:
            fragmentationDictionary[newName] = {'01':{'subgeometry':list(molecularDf[col]),'relCont':1}}
        #otherwise, drop it
        else:
            molecularDf.drop(columns=col, inplace=True)

    else:
        #In this format, split into parent frag name and number
        parentFrag, geometryNumber = newName.split('_')
        #If the fragment is desired, fill in the dictionary
        if parentFrag in fragSubset:
            renameCols[col] = newName
            fragmentationDictionary[parentFrag] = {geometryNumber:{'subgeometry':list(molecularDf[col]), 'relCont':1}}
        #otherwise, drop it
        else:
            molecularDf.drop(columns=col, inplace=True)

    return fragmentationDictionary, molecularDf, renameCols

def otherFragmentReps(fragmentationDictionary, molecularDf):
    """
    Computes alternative representations of the fragments for easier processing.

    Inputs:
        fragmentationDictionary: A dictionary specifying the fragments and their geometries

    Outputs:
        expandedFrags: An ATOM depiction of each fragment, where an ATOM depiction has one entry for each atom (rather than for each site). See fragmentAndSimulate for details.
        fragSubgeometryKeys: A list of strings, e.g. 44_01, 44_02, corresponding to each subgeometry of each fragment. A fragment will have multiple subgeometries if there are multiple fragmentation pathways to form it.
    """
    condensedFrags = []
    fragSubgeometryKeys = []
    for fragKey, subFragDict in fragmentationDictionary.items():
        for subFragNum, subFragInfo in subFragDict.items():
            #for each subgeometry, get which sites it samples and its key. Add to list.
            condensedFrags.append(subFragInfo['subgeometry'])
            fragSubgeometryKeys.append(fragKey + '_' + subFragNum)

    #Expand to ATOM depictions of the fragments.
    expandedFrags = [fas.expandFrag(x, molecularDf['Number']) for x in condensedFrags]

    return expandedFrags, fragSubgeometryKeys

def moleculeFromCsv(path, deltas = False, fragSubset = None):
    '''
    Reads in a .csv file to get basic information about that molecule.

    Inputs:
        path: a string locating the .csv, e.g. 'Example Input.csv'
        fragSubset: A list of fragments to include. Defaults to None, in which case it is filled programatically to contain all fragments. 

    Outputs:
        initializedMolecule: A dictionary, containing:

            molecularDataFrame: A dataframe containing basic information about the molecule. 
            expandedFrags: An ATOM depiction of each fragment, where an ATOM depiction has one entry for each atom (rather than for each site). See fragmentAndSimulate for details.
            fragSubgeometryKeys: A list of strings, e.g. 44_01, 44_02, corresponding to each subgeometry of each fragment. A fragment will have multiple subgeometries if there are multiple fragmentation pathways to form it.
            fragmentationDictionary: A dictionary like the allFragments variable, but only including the subset of fragments selected by fragSubset.
    '''
    #Read csv
    molecularDf = pd.read_csv(path, index_col='Site Names')
    molecularDf.index.name = None

    #Columns must be renamed from the input csv to the format the code expects. We add to this in 'process_fragment_column'.
    renameCols = {'Element': 'IDS', 'Number Atoms': 'Number'}
    molecularDf.rename(columns=renameCols, inplace=True)

    #Fill in programatically in the default case
    if fragSubset == None:
        fragSubset = [col.split(' ')[1] for col in molecularDf.columns if 'Fragment' in col]

    #Fill a dictionary with the fragment information and which sites each fragment samples
    fragmentationDictionary = {}
    for col in molecularDf.columns:
        #get fragment columns
        if 'Fragment' in col:
            molecularDf[col] = molecularDf[col].replace(0, 'x')
            fragmentationDictionary, molecularDf, renameCols = process_fragment_column(col, fragSubset, molecularDf, fragmentationDictionary, renameCols)
        
        #calculate some other objects relating to fragments that are used in computations
        expandedFrags, fragSubgeometryKeys = otherFragmentReps(fragmentationDictionary, molecularDf)

    #rename the columns
    molecularDf.rename(columns=renameCols, inplace=True)

    molecularDf['deltas'] = deltas
    
    #Construct output
    initializedMolecule = {'molecularDataFrame':molecularDf,
                           'expandedFrags':expandedFrags,
                           'fragSubgeometryKeys':fragSubgeometryKeys,
                           'fragmentationDictionary':fragmentationDictionary}
    
    return initializedMolecule

def extractExperimentalErrors(errorPath, UValueError = 0.001, MNError = 0.001):
    '''
    Read in error bars from a user provided csv. If no CSV is provided, instead return a single float giving error for the molecular average measurement (U Value) and for the M+N Experiment. A value of 0.001 corresponds to a 1 per mil error. 

    Inputs:
        errorPath: Path to the CSV. 
        UValueError: A float giving the error on the U Value. Only used if no errorPath provided.
        MNError: A float giving the error on the M+N measurements. Only used if no errorPath provided.

    Outputs:
        UValueErrors, MNErrors: Dictionaries. UValueErrors looks like: {'13C':0.0001}, while MNErrors looks like: {'M1':{'44':{'15N':0.002,'13C':0.001,'Unsub':0.001}}}

        Alternatively:
        UValueError, MNError: The same as the inputs. 
    '''
    if errorPath:
        errorDf = pd.read_csv(errorPath, header = None)

        targetHeaders = ['Molecular Average']
        fragmentPrefix = 'Fragment'

        errorDict = {}

        for rowIdx, row in errorDf.iterrows():
            headerCell = row[0]
            if any(header in headerCell for header in targetHeaders) or headerCell.startswith(fragmentPrefix):
                #Clean up the header to not inlude 'Fragment' prefix
                keyName = headerCell[len(fragmentPrefix):].strip() if headerCell.startswith(fragmentPrefix) else headerCell

                #Initialize Dictionary
                errorDict[keyName] = {}

                # Loop through subsequent rows to gather data until another target header or the end of the DataFrame
                for j in range(rowIdx+1, len(errorDf)):
                    nextRow = errorDf.iloc[j]
                    nextHeaderCell = nextRow[0]
                    # Check if the next row is another header
                    if any(header in nextHeaderCell for header in targetHeaders) or nextHeaderCell.startswith(fragmentPrefix):
                        break  # Found another header, stop processing this section
                    else:
                        # This row is data to be included for the current header
                        errorDict[keyName][nextHeaderCell] = float(nextRow[1]) / 1000

        #Input file has '13C/Unsub' for U Values. Later code expects '13C'. 
        UValueErrors = {subKey.split('/')[0]: subValue for subKey, subValue in errorDict.get('Molecular Average', {}).items()}
        errorDict.pop('Molecular Average', None)

        #Fill in MNErrors in expected format
        MNErrors = {'M1':errorDict}

        return UValueErrors, MNErrors

    else:
        return UValueError, MNError
    

def addClumps(byAtom, molecularDataFrame, clumpD):
    '''
    A function which adjusts the abundances of target isotopologues to include clumps. 

    Inputs:
        byAtom: A dictionary where keys are isotoplogues and values are dictionaries containing details about those isotopologues. 
        molecularDataFrame: A dataframe containing basic information about the molecule. 
        clumpD: Specifies information about clumps to add; otherwise the isotome follows the stochastic assumption. Currently works only for mass 1 substitutions (e.g. 1717, 1317, etc.) See ci.introduceClump for details.

    Outputs:
        byAtom: The same dictionary, with concentrations modified to include clumps. 
    '''
    stochD = copy.deepcopy(byAtom)
    
    for clumpNumber, clumpInfo in clumpD.items():
        byAtom = ci.introduceClump(byAtom, clumpInfo['Sites'], clumpInfo['Amount'], molecularDataFrame)
        
    for clumpNumber, clumpInfo in clumpD.items():
        ci.checkClumpDelta(clumpInfo['Sites'], molecularDataFrame, byAtom, stochD)

    return byAtom

def simulateMeasurement(initializedMolecule, abundanceThreshold = 0, UValueList = [],
                        massThreshold = 1, clumpD = {}, outputPath = None, disableProgress = True, calcFF = False, fractionationFactors = {}, omitMeasurements = {}, ffstd = 0.05, unresolvedDict = {}, outputFull = False):
    '''
    Simulates M+N measurements of an alanine molecule with input deltas specified by the input dataframe molecularDataFrame. 

    Inputs:
        molecularDataFrame: A dataframe containing basic information about the molecule. 
        expandedFrags: An ATOM depiction of each fragment, where an ATOM depiction has one entry for each atom (rather than for each site). See fragmentAndSimulate for details.
        fragSubgeometryKeys: A list of strings, e.g. 44_01, 44_02, corresponding to each subgeometry of each fragment. A fragment will have multiple subgeometries if there are multiple fragmentation pathways to form it.
        fragmentationDictionary: A dictionary like the allFragments variable, but only including the subset of fragments selected by fragSubset.
        abundanceThreshold: A float; Does not include measurements below this M+N relative abundance, i.e. assuming they will not be  measured due to low abundance. 
        UValueList: A list giving specific substitutions to calculate molecular average U values for ('13C', '15N', etc.)
        massThreshold: An integer; will calculate M+N relative abundances for N <= massThreshold
        clumpD: Specifies information about clumps to add; otherwise the isotome follows the stochastic assumption. Currently works only for mass 1 substitutions (e.g. 1717, 1317, etc.) See ci.introduceClump for details.
        outputPath: A string, e.g. 'output', or None. If it is a string, outputs the simulated spectrum as a json. 
        disableProgress: Disables tqdm progress bars when True.
        calcFF: When True, computes a new set of fractionation factors for this measurement.
        fractionationFactors: A dictionary, specifying a fractionation factor to apply to each ion beam. This is used to apply fractionation factors calculated previously to this predicted measurement (e.g. for a sample/standard comparison with the same experimental fractionation)
        omitMeasurements: omitMeasurements: A dictionary, {}, specifying measurements which I will not observed. For example, omitMeasurements = {'M1':{'61':'D'}} would mean I do not observe the D ion beam of the 61 fragment of the M+1 experiment, regardless of its abundance. 
        ffstd: A float; if new fractionation factors are calculated, they are pulled from a normal distribution centered around 1, with this standard deviation.
        unresolvedDict: A dictionary, specifying which unresolved ion beams add to each other.
        outputFull: A boolean. Typically False, in which case beams that are not observed are culled from the dictionary. If True, includes this information; this should only be used for debugging, and will likely break the solver routine. 
        
    Outputs:
        predictedMeasurement: A dictionary giving information from the M+N measurements. 
        MN: A dictionary where keys are mass selections ("M1", "M2") and values are dictionaries giving information about the isotopologues of each mass selection.
        fractionationFactors: The calculated fractionation factors for this measurement (empty unless calcFF == True)
    '''
    #Unpack initialized molecule dictionary
    molecularDataFrame = initializedMolecule['molecularDataFrame']
    expandedFrags = initializedMolecule['expandedFrags']
    fragSubgeometryKeys = initializedMolecule['fragSubgeometryKeys']
    fragmentationDictionary = initializedMolecule['fragmentationDictionary']

    byAtom = ci.inputToAtomDict(molecularDataFrame, disable = disableProgress)
    if clumpD != {}:
        byAtom = addClumps(byAtom, molecularDataFrame, clumpD)

    #bySub is an representation of data, a dictionary where keys are substitutions (e.g., '13C'), and values are their abundances. 
    bySub = ci.calcSubDictionary(byAtom, molecularDataFrame, atomInput = True)
    
    #Initialize Measurement output
    allMeasurementInfo = {}
    allMeasurementInfo = fas.UValueMeasurement(bySub, allMeasurementInfo, massThreshold = massThreshold,subList = UValueList)

    MN = ci.massSelections(byAtom, massThreshold = massThreshold)
    MN = fas.trackMNFragments(MN, expandedFrags, fragSubgeometryKeys, molecularDataFrame, unresolvedDict = unresolvedDict)
        
    predictedMeasurement, fractionationFactors = fas.predictMNFragmentExpt(allMeasurementInfo, MN, expandedFrags, fragSubgeometryKeys, molecularDataFrame, 
                                                 fragmentationDictionary, calcFF = calcFF, ffstd = ffstd,
                                                 abundanceThreshold = abundanceThreshold, fractionationFactors = fractionationFactors, omitMeasurements = omitMeasurements, unresolvedDict = unresolvedDict, outputFull = outputFull)
    
    if outputPath != None:
        output = json.dumps(predictedMeasurement)

        f = open(outputPath + ".json","w")
        f.write(output)
        f.close()
        
    return predictedMeasurement, MN, fractionationFactors



def plotOutput(simulationOutput, relativeDeltasActual, ylim = False):
    '''
    Simple function to plot the output of a simulation.

    Inputs:
        simulationOutput: A dataframe containing the results of a M+1 simulation. returned by simulateSmpStd. 
        relativeDeltasActual: A list, giving the actual delta values for the sample standard comparison. The ordering of the list should be the same as the indices of simulationOutput. 

    Outputs:
        None. Constructs a plot. 
    '''
    fig, ax = plt.subplots(figsize = (8,4))

    plt.scatter(range(len(relativeDeltasActual)),relativeDeltasActual, label = 'Actual', marker = 's', color = 'r')
    plt.errorbar(range(len(relativeDeltasActual)),simulationOutput['Relative Deltas'],simulationOutput['Relative Deltas Error'], fmt = 'o',c='k',capsize = 5,
                label = 'Computed')
    plt.legend()
    ax.set_xticks(range(len(relativeDeltasActual)))
    ax.set_xticklabels(simulationOutput.index, rotation = 45)
    ax.set_ylabel("$\delta_{STD}$")
    if ylim:
        ax.set_ylim(*ylim)

def simulateSmpStd(path, deltasStd, deltasSmp, deltasStdAppx, abundanceThreshold = 0, UValueList = [], massThreshold = 1,  disableProgress = True, calcFF = False, omitMeasurements = {}, ffstd = 0.05, plot = True, MonteCarloN = 100, perturbTheoryOAmt = 0, errorPath = False, MNError = 0, UValueError = 0, resultsFileName = 'output.csv', outputPrecision = 3, UMNSub = '13C'):
    '''
    Parent function which constructs and runs a full sample standard comparison. 

    Inputs:
        path: The path to a CSV file containing data about a molecule. 
        deltaStd: Delta values for the standard, used in the simulation. 
        deltaSmp: Delta values for the sample used in the simulation. 
        deltasStdAppx: Delta values for the approximation of the standard used for the forward model of the standard. 
        abundanceThreshold, UValueList, massThreshold, disableProgress, calcFF, omitMeasurements, ffstd: Inputs to simulateMeasurement. See that function for documentation. 
        plot: Construct a plot showing the comparison. 
        MonteCarloN, perturbTheoryOAmt: Used to run the Monte Carlo routine; see ss.M1MonteCarlo.
        errorPath: The path to a CSV containing information about the experimental errors. 
        MNError: IF no errorPath, then use this value for the error on every observed fragment. 
        UValueError: IF no errorPath, use this value for the error on the molecular average measurement. 
        resultsFileName: Name the output CSV. 
        outputPrecision: Adjust precision on the output CSV. 
        UMNSub: Choose which isotope to anchor to when calculating site-specific data. 

    Outputs:
        simulationOutput: A dataframe containing the results of the simulation. Also produces a cleaned CSV containing the most relevant pieces of information. Also plots these data (if plot). 
    '''
    simArgs = {
    'abundanceThreshold': abundanceThreshold,
    'UValueList': UValueList,
    'massThreshold': massThreshold,
    'disableProgress': disableProgress,
    'calcFF': calcFF,
    'omitMeasurements': omitMeasurements,
    'ffstd': ffstd}

    #Simulate observations of the standard
    thisStd = moleculeFromCsv(path, deltas = deltasStd)
    stdMeasurement, stdMNDict, stdFF = simulateMeasurement(thisStd, **simArgs)

    #Simulate observations of the sample
    thisSmp = moleculeFromCsv(path, deltas = deltasSmp)
    smpMeasurement, smpMNDict, smpFF = simulateMeasurement(thisSmp, **simArgs)

    #Simulate a forward model of the standard. 
    forwardModel = moleculeFromCsv(path, deltas = deltasStdAppx)
    forwardModelPredictions, MNDict, FF = simulateMeasurement(forwardModel,abundanceThreshold = 0,
                                                                        massThreshold = 1,
                                                                            unresolvedDict = {})

    #Get error bars used for the simulation from an input CSV or use the same specified value for each.
    UValueError, MNError = extractExperimentalErrors(errorPath, UValueError = UValueError, MNError = MNError)

    #'Read in' the 'experimental' data (from smpMeasurement, stdMeasurement)
    processStandard = ri.readComputedData(stdMeasurement, error = MNError, theory = forwardModelPredictions)
    processSample = ri.readComputedData(smpMeasurement, error = MNError)
    UValuesSmp = ri.readComputedUValues(smpMeasurement, error = UValueError, UMNSub = UMNSub)
    
    #Get estimates for the OValueCorrection
    OCorrection = ss.OValueCorrectTheoretical(forwardModelPredictions, processSample, massThreshold = 1)
    isotopologuesDict = fas.isotopologueDataFrame(MNDict, forwardModel['molecularDataFrame'])

    #Solve the system, update the dataframe
    M1Results = ss.M1MonteCarlo(processStandard, processSample, OCorrection, isotopologuesDict,
                                forwardModel['fragmentationDictionary'], 
                                N = MonteCarloN, perturbTheoryOAmt = perturbTheoryOAmt, disableProgress = disableProgress)

    processedResults = ss.processM1MCResults(M1Results, UValuesSmp, isotopologuesDict, forwardModel['molecularDataFrame'], UMNSub = [UMNSub], disableProgress = disableProgress)

    simulationOutput = ss.updateSiteSpecificDfM1MC(processedResults, forwardModel['molecularDataFrame'])

    #output nice csv
    cleanSimulationOutput = simulationOutput.drop(['deltas','M1 M+N Relative Abundance', 'M1 M+N Relative Abundance Error', 'UM1','UM1 Error','Calc U Values','Calc U Values Error'], axis=1, inplace=False)
    #Round off decimals
    toRound = ['VPDB etc. Deltas', 'VPDB etc. Deltas Error', 'Relative Deltas','Relative Deltas Error']
    cleanSimulationOutput[toRound] = cleanSimulationOutput[toRound].round(decimals=outputPrecision) 
    cleanSimulationOutput.to_csv(resultsFileName, index=False)

    #Construct Plot
    if plot:
        relativeDeltasActual = [op.compareRelDelta(atomID, delta1, delta2) for atomID, delta1, delta2 in zip(forwardModel['molecularDataFrame']['IDS'], deltasStd, deltasSmp)]

        plotOutput(simulationOutput, relativeDeltasActual)

    return simulationOutput