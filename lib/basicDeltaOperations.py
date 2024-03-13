import numpy as np
import pandas as pd

'''
This code takes care of basic manipulations between delta, ratio, and concentration space. Currently, it works for H, C, N, O, and S. 
'''
                                                                                                                       

#DEFINE STANDARDS
#Doubly isotopic atoms are given as standard ratios. The standards are: VPDB for carbon, AIR for nitrogen, VSMOW for H/O, and CDT for sulfur. 
STD_Rs = {"H": 0.00015576, "C": 0.011180, "N": 0.003676, "17O": 0.0003799, "18O": 0.0020052,
         "33S":0.007877,"34S":0.0441626,"36S":0.000105274}
    
def deltaToConcentration(atomIdentity,delta):
    '''
    Converts an input delta value for a given type of atom in some reference frame to a 4-tuple containing the concentration of the unsubstituted, M+1, M+2, and all other versions of the atom. For atoms with multiple rare isotopes, sets the additional isotopes following a mass-scaling law. 
    
    Inputs:
        atomIdentity: A string giving the isotope of interest. Note "S" gives 33S. 
        delta: The input delta value. See STD_Rs, above, for the standard ratios. 
        
    Outputs:
        A 4-tuple giving the concentration vector for this delta value. The first entry gives the unsubstituted; successive entries give higher mass substitutions. E.g. for S the entries are (32S, 33S, 34S, 36S). 
    '''
    if type(delta) == tuple:
        return twoDeltasToConcentration(atomIdentity, delta)
    
    if atomIdentity == 'C' or atomIdentity == '13C':
        ratio = (delta/1000+1)*STD_Rs['C']
        concentrationSub = ratio/(1+ratio)
        
        return (1-concentrationSub,concentrationSub,0,0)
    
    if atomIdentity == 'H' or atomIdentity == 'D':
        ratio = (delta/1000+1)*STD_Rs['H']
        concentrationSub = ratio/(1+ratio)
        
        return (1-concentrationSub,concentrationSub,0,0)
    
    if atomIdentity == 'N' or atomIdentity == '15N':
        ratio = (delta/1000+1)*STD_Rs['N']
        concentrationSub = ratio/(1+ratio)
        
        return (1-concentrationSub,concentrationSub,0,0)
    
    elif atomIdentity == '17O' or atomIdentity == 'O':
        r17 = (delta/1000+1)*STD_Rs['17O']
        delta18 = 1/0.52 * delta
        r18 = (delta18/1000+1)*STD_Rs['18O']
        
        o17 = r17/(1+r17+r18)
        o18 = r18/(1+r17+r18)
        o16 = 1-(o17+o18)
        
        return (o16,o17,o18,0)
    
    elif atomIdentity == '33S' or atomIdentity == 'S':
        r33 = (delta/1000+1)*STD_Rs['33S']
        delta34 = delta/0.515
        r34 = (delta34/1000+1)*STD_Rs['34S']
        delta36 = delta34*1.9
        r36 = (delta36/1000+1)*STD_Rs['36S']
        
        s33 = r33/(1+r33 + r34 + r36)
        s34 = r34/(1+r33+r34+r36)
        s36 = r36/(1+r33+r34+r36)

        s32 = 1-(s33+s34+s36)
        
        return (s32,s33,s34,s36)
        
    elif atomIdentity == '34S':
        r34 = (delta/1000+1)*STD_Rs['34S']
        delta33 = delta * 0.515
        r33 = (delta/1000+1)*STD_Rs['33S']
        delta36 = delta*1.9
        r36 = (delta36/1000+1)*STD_Rs['36S']
        
        s33 = r33/(1+r33 + r34 + r36)
        s34 = r34/(1+r33+r34+r36)
        s36 = r36/(1+r33+r34+r36)

        s32 = 1-(s33+s34+s36)
        
        return (s32,s33,s34,s36)
                
    else:
        raise Exception('Sorry, I do not know how to deal with ' + atomIdentity)

def twoDeltasToConcentration(atomIdentity, deltaTuple):
    '''
    A special version of the deltaToConcentration function, for 17/18O or 33/34S, when both are constrained via experiment. In this case, we do not set the additional isotopes via a mass scaling law, and instead calculate explicitly.
    
    Inputs:
        atomIdentity: "O" or "S", for oxygen or sulfur. 
        deltaTuple: A 2-tuple. The first entry is 17O or 33S, the second is 18O or 34S. 
        
    Outputs: 
        A 4-tuple giving the concentration vector for this delta value. The first entry gives the unsubstituted; successive entries give higher mass substitutions. E.g. for S the entries are (32S, 33S, 34S, 36S). 
    '''
    if atomIdentity == 'O':
        delta17 = deltaTuple[0]
        delta18 = deltaTuple[1]
        
        r17 = (delta17/1000+1)*STD_Rs['17O']
        r18 = (delta18/1000+1)*STD_Rs['18O']
        
        o17 = r17/(1+r17+r18)
        o18 = r18/(1+r17+r18)
        o16 = 1-(o17+o18)
        
        return (o16,o17,o18,0)
    
    elif atomIdentity == 'S':
        delta33 = deltaTuple[0]
        delta34 = deltaTuple[1]
        
        r33 = (delta33/1000+1)*STD_Rs['33S']
        r34 = (delta34/1000+1)*STD_Rs['34S']
        delta36 = delta34*1.9
        r36 = (delta36/1000+1)*STD_Rs['36S']
        
        s33 = r33/(1+r33 + r34 + r36)
        s34 = r34/(1+r33+r34+r36)
        s36 = r36/(1+r33+r34+r36)

        s32 = 1-(s33+s34+s36)
        
        return (s32,s33,s34,s36)
    
def concentrationToM1Ratio(concentrationTuple):
    '''
    Gives the ratio for the mass 1 rare isotope of an atom.
    
    Inputs:
        concentrationTuple: A 4-tuple giving the concentration vector for this delta value. The first entry gives the unsubstituted; successive entries give higher mass substitutions. E.g. for S the entries are (32S, 33S, 34S, 36S). 
        
    Outputs:
        A float, giving the ratio between the mass 1 isotope and the unsubstituted isotope. 
    '''
    return concentrationTuple[1]/concentrationTuple[0] 

def ratioToDelta(atomIdentity, ratio):
    '''
    Converts an input ratio for a given atom to a delta value.
    
    Inputs:
        atomIdentity: A string giving the isotope of interest
        ratio: The isotope ratio.
        
    outputs: 
        delta: The delta value for that isotope and ratio.
    '''
    delta = 0
    if atomIdentity == 'D':
        atomIdentity = 'H'
        
    if atomIdentity in 'HCN' or atomIdentity in ['13C','15N']:
        #in case atomIdentity is 2H, 13C, 15N, take last character only
        delta = (ratio/STD_Rs[atomIdentity[-1]]-1)*1000
        
    elif atomIdentity == 'O' or atomIdentity == '17O':
        delta = (ratio/STD_Rs['17O']-1)*1000
        
    elif atomIdentity == '18O':
        delta = (ratio/STD_Rs['18O']-1)*1000
        
    elif atomIdentity == 'S' or atomIdentity == '33S':
        delta = (ratio/STD_Rs['33S']-1)*1000
        
    elif atomIdentity == '34S':
        delta = (ratio/STD_Rs['34S']-1)*1000
        
    elif atomIdentity == '36S':
        delta = (ratio/STD_Rs['36S']-1)*1000
        
    else:
        raise Exception('Sorry, I do not know how to deal with ' + atomIdentity)
        
    return delta

def compareRelDelta(atomID, deltaStd, deltaSmp):
    '''
    Given at atom ID and two deltas in VPDB etc. space, finds their relative difference (not in VPDB etc. space). This is useful for making sample standard comparisons. 
    
    Inputs:
        atomID: A string, to be fed to deltaToConcentration
        deltaStd: The delta value of the "standard" (denominator)
        deltaSmp: The delta value of the "sample" (numerator)
        
    Outputs:
        relDelta: The delta value of the sample relative to the standard.
    '''
    rStd = concentrationToM1Ratio(deltaToConcentration(atomID,deltaStd))
    rSmp = concentrationToM1Ratio(deltaToConcentration(atomID,deltaSmp))

    relDelta = 1000*(rSmp/rStd-1)
    
    return relDelta