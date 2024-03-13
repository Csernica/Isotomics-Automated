import itertools
import copy
import math

import numpy as np
import pandas as pd
from tqdm import tqdm

import basicDeltaOperations as op

'''                                                 
This code calculates a dictionary giving all possible isotopologues of a molecule and their concentrations, based on input information about the sites and their isotopic composition.                                                 
                                                                                                                      
The theory for this section is developed in the working M+N paper. Contact Tim for details (tcsernic@caltech.edu).

It assumes one has access to a dataframe specifying details about a molecule. See the tutorial.
'''

#The possible substitutions, by cardinal mass, for each element. 
setsOfElementIsotopes = {'H':(0,1),'N':(0,1),'C':(0,1),'O':(0,1,2),'S':(0,1,2,4)}

def calculateSetsOfSiteIsotopes(molecularDataFrame):
    '''
    Every site has some set of possible isotopes. For single-atomic sites, this is equal to the set of element isotopes value for the relevant element. For multiatomic sites, it is given by a multinomial expansion of the set of element isotopes. For example, a nitrogen site with 2 atoms can have (00), (01), or (11) as possible isotopes. The number of ways to make these combinations are 1, 2, and 1 respectively. This function calculates the possible substitutions and multinomial coefficients. 
    
    Inputs:
        molecularDataFrame: A dataFrame containing information about the molecule.
        
    Outputs: 
        setsOfSiteIsotopes: A list of tuples, where tuple i gives the possible combinations of substitutions at site i. 
        multinomialCoefficients: A list of tuples, where tuple i gives the multinomial coefficients of substitutions at site i. 
    '''
    elIDs = molecularDataFrame['IDS'].values
    numberAtSite = molecularDataFrame['Number'].values

    siteList = [(x,y) for x,y in zip(elIDs, numberAtSite)]

    #Determine the set of Site Isotopes for every site
    setsOfSiteIsotopes = []
    multinomialCoefficients = []

    for site in siteList:
        el = site[0]
        n = site[1]

        if n == 1:
            setsOfSiteIsotopes.append(setsOfElementIsotopes[el])
            multinomialCoefficients.append([1] * len(setsOfElementIsotopes[el]))

        else:
            siteIsotopes = tuple(itertools.combinations_with_replacement(setsOfElementIsotopes[el], n))
            setsOfSiteIsotopes.append(siteIsotopes)

            #Determining the multinomial coefficients takes a bit of work. There must be a more elegant way to do so.
            #First we generate all possible isotopic structures
            siteIsotopicStructures = itertools.product(setsOfElementIsotopes[el],repeat = n)

            #Then we sort each structure, i.e. so that [0,1] and [1,0] both become [0,1]
            sort = [sorted(list(x)) for x in siteIsotopicStructures]

            #Then we count how many times each structure occurs; i.e. [0,1] may appear twice
            counts = []
            for siteComposition in siteIsotopes:
                c = 0
                for isotopeStructure in sort:
                    if list(siteComposition) == isotopeStructure:
                        c += 1
                counts.append(c)

            #Finally, we expand the counts. Suppose we have a set of site Isotopes with [0,0], [0,1], [1,1] with 
            #multinomial Coefficients of [1,2,1], respectively. We want to output [[1,1],[2,2],[1,1]] rather than 
            #[1,2,1], because doing so will allow us to take advantage of the optimized itertools.product function
            #to calculate the symmetry number of isotopologues with many multiatomic sites.

            #One can check that by doing so, "multinomialCoefficients" and "setsOfSiteIsotopes", the output variables
            #from this section, have the same form. 
            processedCounts = [[x] * n for x in counts]

            multinomialCoefficients.append(processedCounts)
            
    return setsOfSiteIsotopes, multinomialCoefficients

def calcAllIsotopologues(setsOfSiteIsotopes, multinomialCoefficients, M1Only = False):
    '''
    Compute all isotopologues of a molecule. For much larger molecules (>1 million isotopologues), we will want to avoid this step and instead just calculate the MN populations we are most interested in. 
    
    Inputs:
        setsOfSiteIsotopes: A list of tuples, where tuple i gives the possible combinations of substitutions at site i. 
        multinomialCoefficients: A list of tuples, where tuple i gives the multinomial coefficients of substitutions at site i. 
        
    Outputs: 
        setOfAllIsotopologues: A list of tuples, where each tuple is an isotopologue of a molecule.
        symmetryNumbers: A list of ints, where int i gives the number of ways to construct isotopologue i. Follows same indexing as setOfAllIsotopologues. 
    '''
    if M1Only:
        setOfM1Isotopologues, symmetryNumbers = calcThroughM1Isotopologues(setsOfSiteIsotopes)
        return setOfM1Isotopologues, symmetryNumbers
    
    i = 0
    setOfAllIsotopologues = []
    symmetryNumbers = []
    for isotopologue in itertools.product(*setsOfSiteIsotopes):
        setOfAllIsotopologues.append(isotopologue)

    #As setsOfSiteIsotopes and multinomialCoefficients are in the same form, we can use the optimized itertools.product
    #again to efficiently calculate the symmetry numbers
    for isotopologue in itertools.product(*multinomialCoefficients):
        flat = [x[0] if type(x) == list else x for x in isotopologue]
        n = np.array(flat).prod()

        symmetryNumbers.append(n)
                             
    return setOfAllIsotopologues, symmetryNumbers

def calcThroughM1Isotopologues(setsOfSiteIsotopes):
    '''
    A workaround to compute only the M1 population of isotopologues (and the unsubstituted isotopologue). This will speed calculation for M+1 experiments. Variants could be written for M+2, M+3, etc.; this is a future goal. 
    
    Inputs:
        setsOfSiteIsotopes: A list of tuples, where tuple i gives the possible combinations of substitutions at site i. 
        
    Outputs: 
        setOfAllIsotopologues: A list of tuples, where each tuple is an isotopologue of a molecule.
        symmetryNumbers: A list of ints, where int i gives the number of ways to construct isotopologue i. Follows same indexing as setOfAllIsotopologues. 
    '''
    symmetryNumbers = []
    setOfM1Isotopologues = []
    unsubstitutedIsotopologue = []

    for siteIsotope in setsOfSiteIsotopes:
        #siteIsotope will contain integers (for single sites) or tuples (for multiatomic sites). So here we are
        #checking only single sites
        if 0 in siteIsotope:
            unsubstitutedIsotopologue.append(0)

        else:
            n = len(siteIsotope[0])
            zeroTup = tuple([0]*n)
            unsubstitutedIsotopologue.append(zeroTup)

    setOfM1Isotopologues.append(unsubstitutedIsotopologue)
    symmetryNumbers.append(1)

    for index, siteIsotope in enumerate(setsOfSiteIsotopes):
        isotopologue = copy.copy(unsubstitutedIsotopologue)
        
        if 1 in siteIsotope:
            #siteIsotope will contain integers (for single sites) or tuples (for multiatomic sites). So here we are
            #checking only single sites
            isotopologue[index] = 1
            symmetryNumbers.append(1)
            setOfM1Isotopologues.append(tuple(isotopologue))

        else:
            n = len(siteIsotope[0])

            singleM1Sub = [0] * n
            singleM1Sub[-1] = 1
            singleM1SubTup = tuple(singleM1Sub)

            if singleM1SubTup in siteIsotope:
                isotopologue[index] = singleM1SubTup
                setOfM1Isotopologues.append(tuple(isotopologue))
                symmetryNumbers.append(n)
            
    return setOfM1Isotopologues, symmetryNumbers

def siteSpecificConcentrations(molecularDataFrame):
    '''
    Calculates all site-specific concentrations and puts them in an array for easy access. Note that at present, it only works for C,N,O,S,H. If we add new elements, we may need to play with the structure of this function. 
    
    The basic structure of the array is: array[i][j] gives the concentration of an isotope with cardinal mass difference i at position j. 
    
    Inputs:
        molecularDataFrame: A dataFrame containing information about the molecule.
        
    Outputs:
        concentrationArray: A numpy array giving the concentration of each isotope at each site. 
    '''
    elIDs = molecularDataFrame['IDS'].values
    numberAtSite = molecularDataFrame['Number'].values
    deltas = molecularDataFrame['deltas'].values
    
    concentrationList = []
    for index in range(len(elIDs)):
        element = elIDs[index]
        delta = deltas[index]
        concentration = op.deltaToConcentration(element, delta)
        concentrationList.append(concentration)

    #put site-specific concentrations into a workable form
    unsub = []
    M1 = []
    M2 = []
    M3 = []
    M4 = []
    for concentration in concentrationList:
        unsub.append(concentration[0])
        M1.append(concentration[1])
        M2.append(concentration[2])
        M3.append(0)
        M4.append(concentration[3])

    concentrationArray = np.array(unsub), np.array(M1), np.array(M2), np.array(M3), np.array(M4)
    
    return concentrationArray

def calculateIsotopologueConcentrations(setOfAllIsotopologues, symmetryNumbers, concentrationArray, disable = False):
    '''
    Puts information about the isotopologues of a molecule, their symmetry numbers, and concentrations of individual isotopes together in order to calculate the concentration of each isotopologue. Does so under the stochastic assumption, i.e. assuming that isotopes are distributed stochastically across all isotopologues.
    
    This is a computationally expensive step--ways to improve would be welcome. For molecules where it is too expensive, it would be expedient to avoid calculating all isotopologues and only calculate the M1, M2, etc populations of interest. 
    
    Inputs:
        setOfAllIsotopologues: A list of tuples, where each tuple is an isotopologue of a molecule.
        symmetryNumbers: A list of ints, where int i gives the number of ways to construct isotopologue i. Follows same indexing as setOfAllIsotopologues. 
        concentrationArray: A numpy array giving the concentration of each isotope at each site. 
        disable: Disables the tqdm progress bar if True. 
        
    Outputs:
        d: A dictionary where the keys are string representations of each isotopologue and the values are dictionaries. For example, a string could be '00100', where there is an M1 substitution at position 2 and M0 isotopes at all other sites. The value dictionaries include "Conc", or concentration, and "num", giving the number of isotopologues of that form. The sum of all concentrations should be 1.  
        
        The keys can be "expanded" strings, i.e. including multiple atomic sites in parentheses. For example, N1/N2 and O3 would appear as (0,1)0. 
    '''
    d = {}
    for i, isotopologue in enumerate(tqdm(setOfAllIsotopologues, disable = disable)):
        number = symmetryNumbers[i]
        isotopeConcList = []
        for index, value in enumerate(isotopologue):
            if type(value) == tuple:
                isotopeConc = [concentrationArray[sub][index] for sub in value]
                isotopeConcList += isotopeConc
            else:
                isotopeConc = concentrationArray[value][index]
                isotopeConcList.append(isotopeConc)      

        isotopologueConc = np.array(isotopeConcList).prod()

        string = ''.join(map(str, isotopologue))

        d[string] = {'Conc':isotopologueConc * number,'num':number}
        
    return d

def condenseStr(text):
    '''
    Takes the "expanded" string depictions, i.e. "(0,1)0" for multiatomic sites and transforms them into "ATOM" depictions, i.e. "010". This makes it easy to pick out the element for a particular substitution by finding the index of the ATOM depiction and looking at that same index in strSiteElements. 
    
    Inputs:
        text: A string, the "expanded" string depiction. 
        
    Outputs:
        text: A string, the "ATOM" string depiction. 
    '''
    text = text.replace('(', '')
    text = text.replace(')', '')
    text = text.replace(',', '')
    text = text.replace(' ', '')
    
    return text
 
def uEl(el, n):
    '''
    Returns the type of substitution, given a chemical element and cardinal mass of isotope.
    
    Inputs:
        el: A string, giving the element of interest
        n: An int, giving the cardinal mass of the isotope
        
    Returns: 
        A string identifying the isotope substitution. 
    '''
    if n == 0:
        return ''
    if n == 'x':
        return ''
    if el == 'C':
        if n == 1:
            return '13C'
    if el == 'H':
        if n == 1:
            return 'D'
    if el == 'O':
        if n == 1:
            return '17O'
        if n == 2:
            return '18O'
    if el == 'N':
        if n == 1:
            return '15N'
    if el == 'S':
        if n == 1:
            return '33S'
        if n == 2:
            return '34S'
        if n == 4:
            return '36S'

def strSiteElements(molecularDataFrame):
    '''
    Our dataframe may include multiatomic sites--for example, we may define site N1/N2 to include two nitrogens and site O3 to have one oxygen. It is useful to have a string where we can index in by position--i.e. "NNO"--to determine the chemical element at a given position. This function defines that string. 
    
    Inputs:
        molecularDataFrame: A dataFrame containing information about the molecule.
        
    Outputs: 
        siteElements: A string giving the chemical element by position, expanding multiatomic sites. 
    '''
    elIDs = molecularDataFrame['IDS'].values
    numberAtSite = molecularDataFrame['Number'].values

    siteList = [(x,y) for x,y in zip(elIDs, numberAtSite)]
    siteElementsList = [site[0] * site[1] for site in siteList]
    siteElements = ''.join(siteElementsList)
    
    return siteElements

def calcAtomDictionary(isotopologueConcentrationDict, molecularDataFrame, disable = False):
    '''
    Given the dictionary from calculateIsotopologueConcentrations, calculates another dictionary with more complete information. Takes the "expanded" string depictions i.e. "(0,1)0" to "ATOM" depictions i.e. "010" and makes these the keys. Stores the expanded depictions, number, and concentration for each isotopologue, then additionally calculates their mass and relevant substitutions. 
    
    Computationally expensive; 20 seconds for methionine.
    
    An example entry from methionine for the unsubstituted isotopologue is shown below:  
    
    '000000000000000000000': {'Number': 1,
      'full': '00(0, 0)00000(0, 0, 0)(0, 0)(0, 0, 0)(0, 0)00',
      'Conc': 0.8906400439358315,
      'Mass': 0,
      'Subs': ''}
      
    Inputs: 
        isotopologueConcentrationDict: The output from calculateIsotopologueConcentrations. 
        molecularDataFrame: A dataFrame containing information about the molecule.
        disable: Disables the tqdm progress bar if True
        
    Outputs: 
        byAtom: A new dictionary containing more complete information about the isotopologues. 
    '''
    siteElements = strSiteElements(molecularDataFrame)
    
    byAtom = {}
    for i, v in tqdm(isotopologueConcentrationDict.items(), disable = disable):
        ATOM = condenseStr(i)
        byAtom[ATOM] = {}
        byAtom[ATOM]['Number'] = v['num']
        byAtom[ATOM]['Full'] = i
        byAtom[ATOM]['Conc'] = v['Conc']
        byAtom[ATOM]['Mass'] = np.array(list(map(int,ATOM))).sum()
        byAtom[ATOM]['Subs'] = '-'.join(filter(None, [uEl(element, int(number)) for element, number in zip(siteElements, ATOM)]))
    
    return byAtom

def calcSubDictionary(isotopologueConcentrationDict, molecularDataFrame, atomInput = False):
    '''
    Similar to the "byAtom" dictionary, a more complete depiction of all isotopologues of a molecule. In this case, rather than index in by ATOM string, index in by substitution--i.e., the key '17O' gives information for all isotopologues with the substituion '17O'. This is a better way to index into this information for certain mass spectrometry experiments, e.g. a molecular average measurement of the ratio between two substitutions. 
    
    This "bySub" dictionary can be calculated directly from the isotopologueConcentrationDict, or from a precalculated ATOM dictionary. 
    
    Computationally expensive; 20 seconds for methionine. 
    
    An example entry from methionine is shown below.
    
    'D': {'Number': 12,
      'Full': ['00(0, 0)00000(0, 0, 0)(0, 0)(0, 0, 0)(0, 0)01',
       '00(0, 0)00000(0, 0, 0)(0, 0)(0, 0, 0)(0, 0)10',
       '00(0, 0)00000(0, 0, 0)(0, 0)(0, 0, 0)(0, 1)00',
       '00(0, 0)00000(0, 0, 0)(0, 0)(0, 0, 1)(0, 0)00',
       '00(0, 0)00000(0, 0, 0)(0, 1)(0, 0, 0)(0, 0)00',
       '00(0, 0)00000(0, 0, 1)(0, 0)(0, 0, 0)(0, 0)00'],
      'Conc': 0.0015953500722996194,
      'Mass': [1, 1, 1, 1, 1, 1],
      'ATOM': ['000000000000000000001',
       '000000000000000000010',
       '000000000000000000100',
       '000000000000000010000',
       '000000000000010000000',
       '000000000001000000000']},
       
    'ATOM' gives the list of all ATOM strings with this substitution. 
      
    Inputs: 
        isotopologueConcentrationDict: The output from calculateIsotopologueConcentrations or a "byAtom" dictionary.
        molecularDataFrame: A dataFrame containing information about the molecule.
        atomInput: Specifies whether the input dictionary is the output of calculateIsotopologueConcentrations or a "byAtom" dictionary
        
    Outputs: 
        bySub: A new dictionary containing giving the same information as the ATOM dictionary but indexed via substitution. 
    '''
    siteElements = strSiteElements(molecularDataFrame)
    if atomInput == False:
        bySub = {}
        for i, v in isotopologueConcentrationDict.items():
            ATOM = condenseStr(i)
            Subs = ''.join([uEl(element, int(number)) for element, number in zip(siteElements, ATOM)])
            if Subs not in bySub:
                bySub[Subs] = {'Number': 0, 'Full': [],'Conc': 0, 'Mass': [], 'ATOM': []}
            bySub[Subs]['Number'] += v['num']
            bySub[Subs]['Full'].append(i)
            bySub[Subs]['Conc'] += v['Conc']
            bySub[Subs]['Mass'].append(np.array(list(map(int,ATOM))).sum())
            bySub[Subs]['ATOM'].append(ATOM)
    
    else:
        bySub = {}
        for i, v in isotopologueConcentrationDict.items():
            Subs = v['Subs']
            if Subs not in bySub:
                bySub[Subs] = {'Number': 0, 'Full': [],'Conc': 0, 'Mass': [], 'ATOM': []}
            bySub[Subs]['Number'] += v['Number']
            bySub[Subs]['Full'].append(v['Full'])
            bySub[Subs]['Conc'] += v['Conc']
            bySub[Subs]['Mass'].append(v['Mass'])
            bySub[Subs]['ATOM'].append(i)
                
    return bySub

def inputToAtomDict(molecularDataFrame, disable = False, M1Only = False):
    '''
    A function wrapper to combine several of the basic tasks leading to construction of the isotopologue dictionary. If you are trying to understand how this works, run each of these functions individually. 
    
    Inputs:
        molecularDataFrame: A dataFrame containing information about the molecule.
        disable: If True, disables tqdm progress bars for the dictionary calculations (which can be time-intensive). 
        M1Only: If True, only calculates the M+1 population. 
    
    Outputs: 
        byAtom: A dictionary where keys are "ATOM strings" (i.e. '0000100010') corresponding to different isotopologues and values are dictionaries, listing information about the concentration, number, composition, etc. of those isotopologues. 
    '''
    if disable == False:
        print("Calculating Isotopologue Concentrations")
    siteElements = strSiteElements(molecularDataFrame)
    siteIsotopes, multinomialCoeff = calculateSetsOfSiteIsotopes(molecularDataFrame)
    bigA, SN = calcAllIsotopologues(siteIsotopes, multinomialCoeff, M1Only = M1Only)
    concentrationArray = siteSpecificConcentrations(molecularDataFrame)
    d = calculateIsotopologueConcentrations(bigA, SN, concentrationArray, disable = disable)

    if disable == False:
        print("Compiling Isotopologue Dictionary")
    byAtom = calcAtomDictionary(d, molecularDataFrame, disable = disable)
    
    return byAtom

def massSelections(atomDictionary, massThreshold = 4):
    '''
    Pulls out M0, M1, etc. populations from the ATOM dictionary, up to specified threshold. Packages them into a dictionary, where keys are "M0", "M1", etc. and values are dictionaries giving the isotopologues associated with that population. 
    
    Inputs:
        atomDictionary: A dictionary with information about all isotopologues, keyed by ATOM strings. The output of calcAtomDictionary.
        massThreshold: An int. Does not include populations with cardinal mass difference above this threshold. 
        
    Outputs:
        A dictionary where the keys are "M0", "M1", etc. and the values are dictionaries containing all isotopologues from the ATOM dictionary with a specified cardinal mass difference. 
    '''
    MNDict = {}
    
    for i in range(massThreshold+1):
        MNDict['M' + str(i)] = {}
        
    for i, v in atomDictionary.items():
        for j in range(massThreshold+1):
            if v['Mass'] == j:
                MNDict['M' + str(j)][i] = v
            
    return MNDict

def introduceClump(clumpD, siteList, clumpAmount, molecularDataFrame):
    '''
    Introduce a clump between any number of sites while keeping the site-specific concentrations the same. This is a complicated operation--we must add concentration to the clumped isotopologue, remove it from the singly-substituted isotopologues, and add it to the unsubstituted isotopologue.
    
    Note--this function only works for mass 1 substitutions. A more general function should be built for mass 2, 3, etc. 
    
    Inputs:
        clumpD: a "byAtom" dictionary including all isotopologues.
        siteList: A list of sites to introduce a clump at, e.g. ['Cmethyl','Cgamma']
        clumpAmount: The amount, in concentration space (not CAP Delta), of the clump to introduce. 
        molecularDataFrame: The initial molecular info dataFrame. 
        
    Outputs:
        clumpD: A byAtom dictionary with the clump added.
    '''
    siteNameList = []

    siteName = molecularDataFrame.index
    siteNumber = molecularDataFrame.Number

    for i, name in enumerate(siteName):
        siteNameList += [name] * siteNumber[i]
        
    unsub = "0" * len(siteNameList)

    clumpD[unsub]['Conc'] += clumpAmount
    clumpSubList = list(unsub)
    for site in siteList:
        #For multiatomic sites, we need to pick the last index, as this is what all variants of the singly-substituted
        #isotopologue are indexed under
        reverseIndex = siteNameList[::-1].index(site)
        index = len(siteNameList) - reverseIndex - 1

        singleSubList = list(unsub)
        singleSubList[index] = "1"
        clumpSubList[index] = "1"

        singleSubString = "".join(singleSubList)
        #Have to multiply by number of occurences of that isotopologue; each one sees the increase/decrease in Conc
        n = clumpD[singleSubString]['Number']
        clumpD[singleSubString]['Conc'] -= clumpAmount * n


    clumpSubString = "".join(clumpSubList)
    n = clumpD[clumpSubString]['Number']
    clumpD[clumpSubString]['Conc'] += clumpAmount * n
    
    return clumpD

def checkClumpDelta(siteList, molecularDataFrame, clumpD, stochD):
    '''
    Checks the CAP Delta value for substitutions at some set of sites and prints these. 
    
    Inputs:
        siteList: A list of sites to introduce a clump at, e.g. ['Cmethyl','Cgamma']
        molecularDataFrame: The initial molecular info dataFrame. 
        clumpD: a "byAtom" dictionary with the clumps added.
        stochD: A "byAtom" dictionary without the clumps present. 
        
    Outputs:
        None. Prints the size of the clump at the input sites. 
    '''
    siteNameList = []

    siteName = molecularDataFrame.index
    siteNumber = molecularDataFrame.Number

    for i, name in enumerate(siteName):
        siteNameList += [name] * siteNumber[i]
        
    unsubStr = "0" * len(siteNameList)
    clumpSubList = ['0'] * len(siteNameList)
    
    for site in siteList:
        #For multiatomic sites, we need to pick the last index, as this is what all variants of the singly-substituted
        #isotopologue are indexed under
        reverseIndex = siteNameList[::-1].index(site)
        index = len(siteNameList) - reverseIndex - 1

        clumpSubList[index] = "1"

    clumpSubString = "".join(clumpSubList)
    
    a = stochD[clumpSubString]['Conc'] / stochD[unsubStr]['Conc']
    b = clumpD[clumpSubString]['Conc'] / clumpD[unsubStr]['Conc']
    capDelta = 1000 * (b/a -1)
    print("Clump at " + ' '.join(siteList) + ' is ' + str(capDelta))