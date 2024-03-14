import json

def readJSON(filepath):
    '''
    A basic function wrapper to read in a JSON file.
    
    Inputs:
        filepath: A string, pointing to the .json file of interest.
        
    Returns:
        read: A dictionary of the .json file. 
    '''
    with open(filepath, 'r') as j:
        read = json.loads(j.read())
    return read

def readComputedData(JSONInput, error = 0, theory = {}):
    '''
    Takes a .json file saved from a computed dataset, extracts the M+N information, and processes it to get ready to solve the M+N measurement. Allows the user to specify some amount of relative experimental error (e.g. 1 per mil) for all observed beams. For standards, includes the predicted measurements in the abscence of resolution, abundance, fractionation, etc. errors. 
    
    Inputs:
        JSONInput: The .json dictionary to input. 
        error: The relative error to assign to each observed peak. Either an int/float, in which case it is the same for all observed peaks, or a dictionary, in which the relative error of each peak is set individually. If it is a dictionary, it is keyed via 'massSelection' then 'fragKey' then 'subKey', e.g. error = {'M1':{'44':{'13C':0.001,'2H':0.005,'15N':0.002}}}
        theory: A dictionary giving the predicted measurement in the abscence of any experimental errors. This is used to standardize in later steps. If the input is a Sample, should be left empty. 
 
    Outputs:
        process: A dictionary giving the processed .json file. Each fragment of each mass selection has a dictionary, specifying the substitutions, their observed abundances and errors, and for standards only, their predicted abundances. 
    '''
    process = {}
    for massSelection, massSelectionData in JSONInput.items():
        #skip non MN measurements
        if massSelection != 'Full Molecule' and massSelection != "M0":
            process[massSelection] = {}
            for fragKey, fragData in massSelectionData.items():
                process[massSelection][fragKey] = {'Subs':[],
                                                           'Predicted Abundance':[],
                                                           'Observed Abundance':[],
                                                           'Error':[]}

                for subKey, subData in fragData.items():
                    process[massSelection][fragKey]['Subs'].append(subKey)
                    process[massSelection][fragKey]['Observed Abundance'].append(subData['Adj. Rel. Abundance'])

                    if type(error) == float or type(error) == int:
                        thisError = error
                    else:
                        try:
                            thisError = error[massSelection][fragKey][subKey]
                        except: 
                            raise Exception("Expected to have error bar for " + massSelection + ' ' + fragKey + ' ' + subKey + ' but do not')

                    process[massSelection][fragKey]['Error'].append(thisError * subData['Adj. Rel. Abundance'])
                    
                    if theory != {}:
                        try:
                            t = theory[massSelection][fragKey][subKey]['Rel. Abundance']
                        except:
                            raise Exception('Measurement has value for ' + fragKey + ' ' + subKey + ' but theory does not')
                        process[massSelection][fragKey]['Predicted Abundance'].append(t)
                    
    return process

def readComputedUValues(JSONInput, error = 0, UMNSub = '13C'):
    '''
    Takes a .json file from a computed dataset, extracts the full molecule U Values, and packages them into a dictionary. Allows the user to specify the relative error for each U value. 
    
    Inputs:
        JSONInput: The .json dictionary to input. 
        error: The relative error to assign to each observed peak. 
        
    Outputs:
        U Values: A dictionary specifying the molecular average U values and their errors. 
        
    '''
    UValues = {}
    if isinstance(error, float):
        # If input_value is a float, multiply every value in k by input_value
        for subKey, subValue in JSONInput['Full Molecule'].items():
            if subKey == UMNSub:
                UValues[subKey] = {'Observed':subValue,'Error':subValue * error}

    elif isinstance(error, dict):
        # If input_value is a dictionary, iterate through k and multiply by the corresponding value
        for subKey, subValue in JSONInput['Full Molecule'].items():
            if subKey == UMNSub:
                if subKey in error:
                    UValues[subKey] = {'Observed':subValue,'Error':subValue * error[subKey]}
                else:
                    print(subKey + ' not found in input, defaulting to 1 per mil')
                    UValues[subKey] = {'Observed':subValue,'Error':subValue * 0.001}
    else:
        raise ValueError("input_value must be either a float or a dictionary.")
            
    return UValues

def checkSampleStandard(sample, standard):
    '''
    A useful feature for checking that sample and standard measurements include the same observations. This can be run when using actual data or using experimental data with fractionation factors which may cause certain beams to be unobserved. In the latter case, passing these peaks to the omitMeasurements option of simulateMeasurement will help. 
    
    Inputs:
        sample: Observed data for the sample
        standard: Observed data for the standard
        
    Outputs:
        None. Prints which peaks are observed in one dataset but not the other. 
    '''
    ####If they include different measurements, will raise errors. Check here, then go back and omit those 
    ####measurements from synthetic dataset. 

    ####If the fractionation factors change, the beams below the abundance cutoff may change and this may
    ####have to be modified. 
    for massSelection, fragData in standard.items():
        for fragKey, subData in fragData.items():
            for subKey in subData['Subs']:
                if subKey not in sample[massSelection][fragKey]['Subs']:
                    print(massSelection + " " + fragKey + " " + subKey + " in standard but not sample")

    for massSelection, fragData in sample.items():
        for fragKey, subData in fragData.items():
            for subKey in subData['Subs']:
                if subKey not in standard[massSelection][fragKey]['Subs']:
                    print(massSelection + " " + fragKey + " " + subKey + " in sample but not standard")
                    
def readObservedData(sampleOutputDict, MNKey = "M1", theory = {}, standard = [], processFragKeys = {}):
    '''
    Process an observed .json file, where keys are files, then fragments, then isotopes, into the "processStandard" or "processSample" input dictionaries. Standardization (by adding predicted values to the dictionary) may or may not be performed; we determine this by setting a list, standard, which is true or false for each file in order.
    
    Inputs:
        sampleOutputDict: A processed .json file containing experimental data.
        MNKey: "M1", "M2", etc. 
        theory: Predicted measurements for a standard. 
        standard: A list of booleans, determining whether each file (in order) is a standard or a sample. 
        processFragKeys: A dictionary. Keys are the fragment Keys assigned to the data in the experimental file; values are the same in the theory data. E.g. {'134.0':'133'}
        
    Outputs:
        process: A dictionary, containing the experimental data processed into the "processStandard" or "processSample" forms, ready for the MC algorithm. 
    '''
    process = {}
    fileIdx = 0
    for fileKey, fileInfo in sampleOutputDict.items():
        shortFileKey = fileKey.split('/')[-1]
        shortFileKey = shortFileKey.split('.')[0]
        process[shortFileKey] = {}
        process[shortFileKey][MNKey] = {}
        for fragKey, fragInfo in fileInfo.items():
            newFragKey = fragKey
            if fragKey in processFragKeys:
                newFragKey = processFragKeys[fragKey] 
                process[shortFileKey][MNKey][newFragKey] = {'Subs':[],
                                           'Predicted Abundance':[],
                                           'Observed Abundance':[],
                                           'Error':[]}

                for subKey, subInfo in fragInfo.items():
                    process[shortFileKey][MNKey][newFragKey]['Subs'].append(subKey)
                    process[shortFileKey][MNKey][newFragKey]['Observed Abundance'].append(subInfo['Average'])
                    process[shortFileKey][MNKey][newFragKey]['Error'].append(subInfo['StdError'])

                    if standard[fileIdx] == True:
                        if theory != {}:
                            try:
                                t = theory[MNKey][newFragKey][subKey]['Adj. Rel. Abundance']
                                process[shortFileKey][MNKey][newFragKey]['Predicted Abundance'].append(t)
                            except:
                                raise Exception('Measurement has value for ' + newFragKey + ' or ' + fragKey + ' ' + subKey + ' but theory does not')
        fileIdx += 1

    return process