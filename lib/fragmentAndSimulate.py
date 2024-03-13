import numpy as np
import pandas as pd

import basicDeltaOperations as op
import calcIsotopologues as ci

'''
This code extracts the concentrations of isotopologues of interest from the dictionary of all isotopologues   
in order to predict the outcomes of meaurements. It also allows one to fragment the isotopologues to compute the     
outcome of fragment measurements.                                                                                                            

It assumes one has access to a dictionary with information about the isotopologues. See calcIsotopologues.py. 
'''

#Gives an easy way to recover an isotope from an element and its cardinal mass representation. 
subDict = {'C':{'0':'','1':'13C'},
           'N':{'0':'','1':'15N'},
           'H':{'0':'','1':'D'},
           'O':{'0':'','1':'17O','2':'18O'},
           'S':{'0':'','1':'33S','2':'34S','4':'36S'}}

#An easy way to recover the mass of an isotope from an element and its cardinal mass representation
massDict = {'C':{'0':12,'1':13.00335484},
            'N':{'0':14.003074,'1':15.00010889},
            'H':{'0':1.007825032,'1':2.014101778},
            'O':{'0':15.99491462,'1':16.99913175,'2':17.9991596},
            'S':{'0':31.97207117,'1':32.9714589,'2':33.96786701,'4':35.9670807}}

def UValueMeasurement(bySub, allMeasurementInfo, massThreshold = 3, subList = []):
    '''
    Simulates measurements with no fragmentation. Extracts the concentration of all isotopologues with mass below some threshold for easy reference. 
    
    Inputs:
        bySub: A dictionary with information about all isotopologues of a molecule, sorted by substitution. 
        allMeasurementInfo: A dictionary containing information from many types of measurements. 
        massThreshold: A mass cutoff; isotopologues with cardinal mass change above this will not be included unless indicated in subList. 
        subList: A list giving specific substitutions to calculate U values for ('13C', '15N', etc.). If substitutions are given, calculates U values only for these substitutions. Otherwise, calculates all U values below the mass threshold. 
        
    Outputs:
        allMeasurementInfo: A dictionary, updated to include information from direct measurements.
    '''
    unsubConc = bySub['']['Conc']
    
    if 'Full Molecule' not in allMeasurementInfo:
        allMeasurementInfo['Full Molecule'] = {}
        
    for sub, info in bySub.items():
        if info['Mass'][0] <= massThreshold and subList == []:
            allMeasurementInfo['Full Molecule'][sub] = info['Conc'] / unsubConc
        elif sub in subList:
            allMeasurementInfo['Full Molecule'][sub] = info['Conc'] / unsubConc
            
    return allMeasurementInfo
    

def fragMult(z, y):
    '''
    Fragments an individual site of an isotopologue. z should be 1 or 'x'. 
    
    Inputs:
        z: specifies whether a site is retained (1) or lost ('x')
        y: The mass of a substitution at that site
        
    Outputs:
        'x', specifying that the site is lost, or y, specifying that the site remains. 
    '''
    if z not in [1,'x']:
        raise Exception("Cannot fragment successfully, each site must be lost ('x') or retained (1)")
    if z == 'x' or y == 'x':
        return 'x'
    else:
        return y
    
def expandFrag(siteDepict, numberAtSite):
    '''
    Creates an ATOM depiction of a fragment from a SITE depiction of a fragment. Expands the fragmentation vector according to the number of atoms at each site. For example, if I fragment [0,(0,1)] with fragmentation vector [0,1], I do so by applying the fragmentation vector [011] to the isotopologue [001], expanding the tuple. This function expands the fragmentation vector. 
    
    Inputs:
        siteDepict: SITE depiction of fragmentation vector.
        
    Outputs:
        atomDepict: expanded depiction of fragmentation vector
    '''
    atomDepict = []
    for i, v in enumerate(siteDepict):
        atomDepict += [v] * numberAtSite[i]
    
    return atomDepict

def fragmentOneIsotopologue(atomFrag, isotopologue):
    '''
    Applies the ATOM fragmentation vector to the ATOM depiction of an isotopologue. Raises a warning if they are not the same length. Returns the ATOM depiction of the isotopologue with "x" in positions that are lost.
    
    Inputs:
        atomFrag: The ATOM depiction of the fragmentation vector
        isotopologue: The ATOM depiction of the isotopologue, a string
        
    Outputs:
        A string giving the ATOM depiction of a fragmented isotopologue. 
    '''
    #important to raise this--otherwise one may inadvertantly fragment incorrectly. 
    if len(atomFrag) != len(isotopologue):
           raise Exception("Cannot fragment successfully, as the fragment and the isotopologue you want to fragment have different lengths")
            
    a = [fragMult(x,y) for x, y in zip(atomFrag, isotopologue)]
    
    if len(a) != len(isotopologue):
        raise Exception("Cannot fragment successfully, the resulting fragment has a different length than the input isotopologue.")
    
    return ''.join(a)

def fragmentIsotopologueDict(atomIsotopologueDict, atomFrag, relContribution = 1):
    '''
    Applies the same fragmentation vector to all isotopologues of an input isotopologue dict and stores the results. This operation corresponds to the "fragmentation" operation from the M+N paper. Combines isotopologues which fragment to yield the same product. For the version which does track, see "fragmentAndTrackIsotopologues"
    
    Inputs:
        atomIsotopologueDict: A dictionary containing some set of isotopologues, often a M1, M2, ... set, keyed by their ATOM depiction. 
        atomFrag: An ATOM depiction of a fragment
        relContribution: A float between 0 and 1, giving the relative contribution of this fragmentation geometry to the observed ion beam at that mass
        
    Outputs: 
        fragmentedDict: A dictionary where the keys are the ATOM isotopologues after fragmentation (i.e. "0000x") and the values are the concentrations of those isotopologues. Note that this may combine isotopologues from the input dictionary which fragment in the same way; i.e. 001 and 002 both fragment to yield "00x". 
    '''
    fragmentedDict = {}
    for isotopologue, value in atomIsotopologueDict.items():
        newIsotopologue = fragmentOneIsotopologue(atomFrag, isotopologue)
        if newIsotopologue not in fragmentedDict:
            fragmentedDict[newIsotopologue] = 0
        fragmentedDict[newIsotopologue] += (value['Conc'] * relContribution)
        
    return fragmentedDict
    
def computeSubs(isotopologue, IDs):
    '''
    Given an ATOM depiction of an isotopologue, computes which substitutions are present. 
    
    Inputs:
        isotopologue: The ATOM string depiction of an isotopologue
        IDs: The string of site elements, i.e. the output of strSiteElements
        
    Outputs:
        A string giving substitutions present in that isotopologue, separated by "-". I.e. "17O-17O"
    '''
    subs = []
    for i in range(len(isotopologue)):
        if isotopologue[i] != 'x':
            element = IDs[i]
            if subDict[element][str(isotopologue[i])] != '':
                subs.append(subDict[element][str(isotopologue[i])])
                
    if subs == []:
        return "Unsub"
        
    return '-'.join(subs)

def computeMass(isotopologue, IDs):
    '''
    Used to predict and generate spectra with exact masses. 
    
    Inputs:
        isotopologue: A string, the ATOM depiction of an isotopologue.
        IDs: A string, the ATOM depiction of element IDs.
        
    Outputs:
        mass: A float, giving the exact mass of the isotopologue. 
    '''
    mass = 0
    for i in range(len(isotopologue)):
        if isotopologue[i] != 'x':
            element = IDs[i]
            mass += massDict[element][str(isotopologue[i])]
        
    return mass

def predictMNFragmentExpt(allMeasurementInfo, MNDict, atomFragList, fragSubgeometryKeys, molecularDataFrame, fragmentationDictionary, abundanceThreshold = 0, omitMeasurements = {}, fractionationFactors = {}, calcFF = False, ffstd = 0.05, randomseed = 25, unresolvedDict = {}, outputFull = False):
    '''
    Predicts the results of several M+N experiements across a range of mass selected populations and fragments. It incorporates the preceding functions into a whole, so you can just call this and get results.
    
    Inputs:
        allMeasurementInfo: A dictionary containing information from many types of measurements. 
        MNDict: A dictionary where the keys are "M0", "M1", etc. and the values are dictionaries containing all isotopologues from the ATOM dictionary with a specified cardinal mass difference. See massSelections function in calcIsotopologues.py
        atomFragList: A list of expanded fragments, one for each subgeometry, e.g. [[1, 1, 1, 1, 'x'], ['x', 1, 1, 1, 'x']]. See expandFrags function.
        fragSubgeometryKeys: A list of strings, indicating the identity of each fragment subgeometry. I.e. ['54_01','42_01']
        molecularDataFrame: A dataFrame containing information about the molecule.
        fragmentationDictionary: A dictionary, e.g. {'full': {'01': {'subgeometry': [1, 1, 1, 1, 1, 1], 'relCont': 1}},
                                                     '44': {'01': {'subgeometry': [1, 'x', 'x', 1, 1, 'x'], 'relCont': 1}}}
                                 which gives information about the fragments, their subgeometries and relative contributions.
        abundanceThreshold: Does not include measurements below a certain relative abundance, i.e. assuming they will not be  measured due to low abundance. 
        omitMeasurements: A dictionary, {}, specifying measurements which I will not observed. For example, omitMeasurements = {'M1':{'61':'D'}} would mean I do not observe the D ion beam of the 61 fragment of the M+1 experiment, regardless of its abundance. 
        fractionationFactors: A dictionary, specifying a fractionation factor to apply to each ion beam. This is used to apply fractionation factors calculated previously to this predicted measurement (e.g. for a sample/standard comparison with the same experimental fractionation). 
        calcFF: A boolean, specifying whether new fractionation factors should be calculated via this function. If True, fractionFactors should be left empty. 
        ffstd: A float. If new fractionation factors are calculated, they are generated from a normal distribution with mean 1 and standard deviation of ffstd. 
        randomseed: An integer. If new fractionation factors are calculated, we initialize this random seed; this allows us to generate the same factors if we run multiple times. 
        unresolvedDict: A dictionary, specifying which unresolved ion beams add to each other. 
        outputFull: A boolean. Typically False, in which case beams that are not observed are culled from the dictionary. If True, includes this information; this should only be used for debugging, and will likely break the solver routine. 
        
    Outputs: 
        allMeasurementInfo: A dictionary, updated to include information from the M+N measurements. 
        calculatedFF: The calculated fractionation factors for this measurement (empty unless calcFF == True)
    '''
    calculatedFF = {}
    siteElements = ci.strSiteElements(molecularDataFrame)
    np.random.seed(randomseed)
    #For each population (M1, M2, M3) that we mass select
    for massSelection, MN in MNDict.items():
        #add a key to output dictionary
        if massSelection not in allMeasurementInfo:
            allMeasurementInfo[massSelection] = {}
            
        if calcFF == True:
            calculatedFF[massSelection] = {}

        #For each fragment we will observe
        for j, fragment in enumerate(atomFragList):

            #add a key to output dictionary
            if fragSubgeometryKeys[j] not in allMeasurementInfo[massSelection]:
                allMeasurementInfo[massSelection][fragSubgeometryKeys[j]] = {}
                
            if calcFF == True:
                calculatedFF[massSelection][fragSubgeometryKeys[j]] = {}
 
            #fragment the mass selection accordingly 
            fragKey, fragNum = fragSubgeometryKeys[j].split('_')
            relContribution = fragmentationDictionary[fragKey][fragNum]['relCont']
            fragmentedIsotopologues = fragmentIsotopologueDict(MN, fragment, relContribution = relContribution)

            #compute the absolute abundance of each substitution
            predictSpectrum = {}

            for key, item in fragmentedIsotopologues.items():
                sub = computeSubs(key, siteElements)
                    
                if sub not in predictSpectrum:
                    predictSpectrum[sub] = {'Abs. Abundance':0}
                predictSpectrum[sub]['Abs. Abundance'] += item
            
            #Fractionate
            if calcFF == True:
                for sub in predictSpectrum.keys():
                    ff = np.random.normal(1,ffstd)
                    calculatedFF[massSelection][fragSubgeometryKeys[j]][sub] = ff
                    predictSpectrum[sub]['Abs. Abundance'] *= ff
            
            elif fractionationFactors != {}:
                for sub in predictSpectrum.keys():
                    predictSpectrum[sub]['Abs. Abundance'] *= fractionationFactors[massSelection][fragSubgeometryKeys[j]][sub]
                    
            allMeasurementInfo[massSelection][fragSubgeometryKeys[j]] = predictSpectrum
    
    allMeasurementInfo = combineFragmentSubgeometries(allMeasurementInfo, fragmentationDictionary)
    
    allMeasurementInfo = computeMNRelAbundances(allMeasurementInfo, omitMeasurements = omitMeasurements, abundanceThreshold = abundanceThreshold, unresolvedDict = unresolvedDict, outputFull = outputFull)
                             
    return allMeasurementInfo, calculatedFF

def combineFragmentSubgeometries(allMeasurementInfo, fragmentationDictionary):
    '''
    Takes fragments with multiple subgeometries and combines their measurements. For example, if frag 82 is made via 82_01 (relCont = 0.4) and 82_02 (relCont = 0.6) this function adds the values of these subfragments to give the actual measurement. 
    
    Inputs:
        allMeasurementInfo: A dictionary containing information about the measurement including fragment subgeometries.
        fragmentationDictionary: A dictionary giving information about the fragments and their subgeometries. 
        
    Outputs:
        combinedAllMeasurementInfo: A dictionary containing information about the measurement including only full fragments. 
    '''
    combinedAllMeasurementInfo = {}
    for massSelection, fragmentData in allMeasurementInfo.items():
        #only take MN experiments
        if massSelection[0] != 'M':
            combinedAllMeasurementInfo[massSelection] = fragmentData
        else:
            combinedAllMeasurementInfo[massSelection] = {}
            for fullFragKey, isotopicData in fragmentData.items():
                fragKey, fragNum = fullFragKey.split('_')
               
                if fragKey not in combinedAllMeasurementInfo[massSelection]:
                    combinedAllMeasurementInfo[massSelection][fragKey] = {}

                for isotopicSub, subData in isotopicData.items():
                    if isotopicSub not in combinedAllMeasurementInfo[massSelection][fragKey]:
                        combinedAllMeasurementInfo[massSelection][fragKey][isotopicSub] = {'Abs. Abundance':0}

                    combinedAllMeasurementInfo[massSelection][fragKey][isotopicSub]['Abs. Abundance'] += subData['Abs. Abundance']
                    
    return combinedAllMeasurementInfo
        
def computeMNRelAbundances(allMeasurementInfo, omitMeasurements = {}, abundanceThreshold = 0, unresolvedDict = {}, outputFull = False):
    '''
    Compute relative abundances from a MN experiment.
    
    Inputs:
        allMeasurementInfo: A dictionary containing information about the absolute abundance of peaks observed in the measurement. 
        omitMeasurements: Allows a user to manually specify ion beams to not measure. For example, omitMeasurements = {'M1':{'61':'D'}} would mean I do not observe the D ion beam of the 61 fragment of the M+1 experiment, regardless of its abundance. 
        abundanceThreshold: gives a relative abundance threshold (e.g. 0.01) below which peaks will not be observed. If a simulated ion beam has relative abundance below this threshold, it is culled from the predicted measurement. 
        unresolvedDict: A dictionary, specifying which unresolved ion beams add to each other. 
        outputFull: False by default. Can be set True to include information about all ion beams, not only the observed ones. This is useful for debugging. outputFull: False by default. Can be set True to include information about all ion beams, not only the observed ones. This is useful for debugging. 
        
    Outputs:
        allMeasurementInfo: A dictionary, containing information about the relative abundances of peaks observed in the measurement. 
    '''
    
    for massSelection, fragmentData in allMeasurementInfo.items():
        #only take MN experiments
        if massSelection[0] == 'M':
            #By fragment
            for fragKey, isotopicData in fragmentData.items():
                #compute relative abundance of each substitution
                totalAbundance = 0
                
                #Get abundance of each sub
                for isotopicSub, subData in isotopicData.items():
                    totalAbundance += subData['Abs. Abundance']
            
                #compute relative abundances
                for isotopicSub, subData in isotopicData.items():
                    subData['Rel. Abundance'] = subData['Abs. Abundance'] / totalAbundance
                
                #Coalescing peaks--if we are moving abundance from one substitution to another
                for isotopicSub, subData in isotopicData.items():
                    #check to see if we have to 
                    try:
                        #if we do, set the coalesced relative abundance of the old sub to 0
                        newSub = unresolvedDict[massSelection][fragKey][isotopicSub]
                        subData['Combined Rel. Abundance'] = 0
                    except:
                        newSub = isotopicSub
                        
                    #Then find the new substitution
                    newSubData = allMeasurementInfo[massSelection][fragKey][newSub]
                    
                    #Add the old subs relative abundance to the new sub
                    if 'Combined Rel. Abundance' not in newSubData:
                        newSubData['Combined Rel. Abundance'] = subData['Rel. Abundance']
                    else:
                        newSubData['Combined Rel. Abundance'] += subData['Rel. Abundance']
                            
                #Calculate adjusted relative abundance, which does not include contributions from peaks below some
                #threshold
                shortSpectrum = {}
                totalRelAbund = 0
                try:
                    forbiddenPeaks = omitMeasurements[massSelection][fragKey]
                except:
                    forbiddenPeaks = []

                for isotopicSub, subData in isotopicData.items():
                    #If the peak is observed, count it
                    if subData['Combined Rel. Abundance'] > abundanceThreshold and isotopicSub not in forbiddenPeaks:
                        shortSpectrum[isotopicSub] = subData
                        totalRelAbund += subData['Combined Rel. Abundance']
                    #Otherwise, either 1) set Adj. Rel. Abundance to 0, keeping it in the spectrum or
                    #                  2) cull it from the spectrum 
                    else:
                        if outputFull:
                            shortSpectrum[isotopicSub] = subData
                            shortSpectrum[isotopicSub]['Adj. Rel. Abundance'] = 0

                #calculate adj. rel. abundance for the qualifying peaks
                for isotopicSub, subData in shortSpectrum.items():
                    #If we added adj. rel. abundance = 0 the previous step, we don't want to repeat that calculation
                    if 'Adj. Rel. Abundance' not in subData:
                        subData['Adj. Rel. Abundance'] = subData['Combined Rel. Abundance'] / totalRelAbund
                    
                allMeasurementInfo[massSelection][fragKey] = shortSpectrum
                
    return allMeasurementInfo

def trackMNFragments(MN, expandedFrags, fragSubgeometryKeys, molecularDataFrame, unresolvedDict = {}):
    '''
    Fragments and tracks isotopologues across a range of mass selections and fragments. 
    
    Inputs:
        MN: A dictionary, where the key is an MNKey and the values give information about all isotopologues associated with that mass selection. 
        expandedFrags: A list of the expanded fragments
        fragSubgeometryKeys: A list of the fragment subgeometry keys 
        molecularDataFrame: The initial dataframe with information about the molecule
        unresolvedDict: A dictionary, specifying which unresolved ion beams add to each other. 
        
    Outputs:
        MN: The same dictionary, with information about fragmentation added. 
    '''
    unsubString = list(MN['M0'].keys())[0]
    UnsubConc = MN['M0'][unsubString]['Conc']

    for key in list(MN.keys()):
        massSelection = MN[key]
        for i, fragment in enumerate(expandedFrags):
            fragmentAndTrackIsotopologues(massSelection, fragment, fragSubgeometryKeys[i], UnsubConc, molecularDataFrame, unresolvedDict = unresolvedDict)
            
    return MN

def fragmentAndTrackIsotopologues(massSelection, atomFrag, fragmentKey, unsubConc, molecularDataFrame, unresolvedDict = {}):
    '''
    Fragments isotopologues and tracks which parent isotopologues end up in which product. For the version that combines isotopologues, for simulating measurements, see fragmentIsotopologueDict (that is, if 001 and 002 both form 00x on fragmentation, this function tracks 001 and 002 explictly; fragmentIsotopologueDict only reports 00x). This function fills in a dictionary with the isotopologues introduced to be fragmented by identifying the product and substitutions of each. 
    
    inputs:
        massSelection: A subset of isotopologues indexed using the ATOM depiction. 
        atomFrag: An ATOM depiction fragment
        fragmentKey: A string giving the identity of the fragment. 
        unsubConc: The concentration of the unsubstituted isotopologue. 
        molecularDataFrame: A dataFrame containing information about the molecule.
        unresolvedDict: {'133':{'17O':'13C'}}
        
    outputs:
        massSelection: The same dictionary, updated to include information about fragmentation. 
    '''
    siteElements = ci.strSiteElements(molecularDataFrame)
    
    fragmentedDict = {}
    for isotopologue, value in massSelection.items():
        value['Stochastic U'] = value['Conc'] / unsubConc
        frag = [fragMult(x,y) for x, y in zip(atomFrag, isotopologue)]
        newIsotopologue = ''.join(frag)
        massSelection[isotopologue][fragmentKey + ' Identity'] = newIsotopologue
        
        sub = computeSubs(newIsotopologue, siteElements)
            
        #If unresolved peaks are a problem
        if fragmentKey in unresolvedDict:
            if sub in unresolvedDict[fragmentKey]:
                sub = unresolvedDict[fragmentKey][sub]
            
        massSelection[isotopologue][fragmentKey + ' Subs'] = sub
        
    return massSelection

def isotopologueDataFrame(MNDictionary, molecularDataFrame):
    '''
    Given a dictionary containing different mass selections, iterates through each mass selection. Extracts the isotopologues from each and puts them into a dataframe, identifying their concentration, substitution, as well as a long string giving a "precise identity", i.e. including explicit labels. Returns these as a dictionary with keys "M0", "M1", etc. where the values are dataFrames of the isotopologues. 
    
    Inputs:
        MNDictionary: A dictionary containing different mass selections, i.e. the output of fragmentAndTrackIsotopologues
        molecularDataFrame: A dataFrame containing information about the molecule.
        
    Outputs:
        isotopologuesDict: A dictionary where the keys are "M0", "M1", etc. and the values are dataFrames giving the isotopologues with those substitutions. 
    '''
    
    isotopologuesDict = {}
    siteElements = ci.strSiteElements(molecularDataFrame)
    
    for key in list(MNDictionary.keys()):
        massSelection = MNDictionary[key]
    
        Isotopologues = pd.DataFrame.from_dict(massSelection).T
        Isotopologues.rename(columns={'Conc':'Stochastic',"Subs": "Composition"},inplace = True)
        
        preciseStrings = []
        
        expandedIndices = []
        for i, n in enumerate(molecularDataFrame.Number):
            expandedIndices += n * [molecularDataFrame.index[i]]
        
        for i, v in Isotopologues.iterrows():
            Subs = [ci.uEl(element, int(number)) for element, number in zip(siteElements, i)]
           
            Precise = [x + " " + y for x, y in zip(Subs, expandedIndices) if x != '']
            output = '   |   '.join(Precise)
            preciseStrings.append(output)
        Isotopologues['Precise Identity'] = preciseStrings
        Isotopologues.sort_values('Composition',inplace = True)
        
        isotopologuesDict[key] = Isotopologues
        
    return isotopologuesDict