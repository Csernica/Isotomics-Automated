import fragmentAndSimulate as fas
import pandas as pd
import matplotlib.pyplot as plt
import calcIsotopologues as ci

MASS_CHANGE = {'13C': 13.00335484 - 12,
            '15N': 15.00010889 - 14.003074,
            'D': 2.014101778 - 1.007825032,
            '17O': 16.99913175 - 15.99491462,
            '18O': 17.9991596 - 15.99491462,
            '33S': 32.9714589 - 31.97207117,
            '34S':33.96786701 - 31.97207117,
            '36S':35.96708707 - 31.97207117}

def massChangeVsUnsub(subKey):
    '''
    Given a subKey (e.g., 'D-D-13C'), returns the mass change caused by those isotope substitutions.

    Inputs:
        subKey: A string, like 'D-D-13C'.

    Outputs:
        mc: A float, the mass change associated with those rare isotope substitutions. 
    '''
    individualSubs = subKey.split('-')

    mc = 0
    for thisSub in individualSubs:
        mc += MASS_CHANGE[thisSub]

    return mc

def fullSpectrumVis(fullMolecule, molecularDataFrame, figsize = (8,4), massError = 0, lowAbundanceCutOff = 0, xlim = (), ylim = (), outputIsotopologs = False, tolerance = 0.1):
    '''
    Takes in some predicted spectral data for the full molecule and outputs a plot of that spectrum. Plots these as relative abundances. 

    (Note: The 'fullMolecule' input is typically from simulateMeasurement, which gives U values (e.g., 13C/Unsub) rather than abundances (e.g., 13C abundance). Because we renormalize to relative abundances and the Unsub in each denominator cancels, the spectrum is still correct). 

    Inputs:
        fullMolecule: A dictionary, where keys are substitutions (e.g., '' for Unsub, '13C', '13C-D'), and values are abundances. 
        molecularDataFrame: A dataframe containing basic information about the molecule. 
        figsize: Desired size of the output figure. 
        massError: Apply any mass error to the computed masses. 
        lowAbundanceCutOff: Do not show peaks below this abundance. 
        xlim: Override the xlimits of the plot
        ylim: Override the ylimits of the plot. 
        outputIsotopologs: If true, outputs a isotoplogs.tsv file with the peaks in the spectrum. 
        tolerance: The tolerance setting to use for an output isotopologs.tsv. 

    Outputs:
        None. Displays a plot.
    '''
    spectrumToPlot = {}

    #Calculate the mass of the unsubstituted isotopologue
    unsubStr = molecularDataFrame['Number'].sum() * '0'
    strSiteElements = ci.strSiteElements(molecularDataFrame)
    unsubMass = fas.computeMass(unsubStr, strSiteElements)

    #Go through the spectral data and pull out abundances. 
    for subKey, subData in fullMolecule.items():
        mass = unsubMass + massChangeVsUnsub(subKey) if subKey else unsubMass
        correctedMass = mass + massError

        # Initialize or update the dictionary for correctedMass
        spectrumToPlot.setdefault(correctedMass, {'U Value': 0, 'Sub': subKey})
        spectrumToPlot[correctedMass]['U Value'] = subData

    #Calculates the relative abundances, and places these, the masses, and substitutions into lists to plot. 
        totalAbundance = sum(item['U Value'] for item in spectrumToPlot.values())

        massPlot = []
        relAbundPlot = []
        subPlot = []
        for mass, data in spectrumToPlot.items():
            relAbundance = data['U Value'] / totalAbundance
            if relAbundance > lowAbundanceCutOff:
                massPlot.append(mass)
                relAbundPlot.append(relAbundance)
                subPlot.append(data['Sub'])
            
    #Constructs a figure; does not plot peaks below the relative abundance cut off. 
    fig, ax = plt.subplots(figsize = figsize)
    ax.vlines(massPlot, 0, relAbundPlot)

    # Set x-axis labels to mass and substitution info
    labels = [f"{round(mass, 5)}\n{sub}" for mass, sub in zip(massPlot, subPlot)]
    ax.set_xticks(massPlot)
    ax.set_xticklabels(labels, rotation=45)
    # Apply optional axis limits
    if xlim:
        ax.set_xlim(*xlim)
    if ylim:
        ax.set_ylim(*ylim)
    ax.set_ylabel("Relative Abundance")

    if outputIsotopologs:
        constructIsotopologsTSV('Molecular Average', subPlot, massPlot, tolerance = tolerance, file_path = 'isotopologsMolecularAverage.tsv')
            
    plt.show()

def MNSpectrumVis(molecularDataFrame, fragKey, predictedMeasurement, MNKey, MNDict, lowAbundanceCutOff = 0, massError = 0, xlim = (), ylim = (), outputIsotopologs = False, tolerance = 0.1):
    '''
    Visualizes the fragmented spectrum of an M+N experiment based on the abundances of fragment peaks.

    Inputs:
        molecularDataFrame: A dataframe containing basic information about the molecule. 
        fragKey: A string identifying the fragment, e.g. '133', '44'.
        predictedMeasurement: A dictionary giving information about the abundance of isotopic peaks in the fragment. See fragmentAndSimulate.predictMNFragmentExpt
        MNKey: A string identifying the mass selection to visualize; e.g. 'M1', 'M2'. This population of isotopologues is selected prior to fragmentation.
        MNDict: A dictionary; the keys are MNKeys ("M1", "M2") and the values are dictionaries containing the isotopologues present in each mass selection. See calcIsotopologues.massSelections
        lowAbundanceCutOff: Do not show peaks below this relative abundance.
        massError: In amu, shifts all observed peaks by this amount. 
        xlim: Set an xlim for the plot, as (xlow, xhigh)
        ylim: as xlim. 
        outputIsotopologs: If true, outputs a isotoplogs.tsv file with the peaks in the spectrum. 
        tolerance: The tolerance setting to use for an output isotopologs.tsv. 

    Outputs:
        None. Displays plot. 
    '''
    isotopologueData = predictedMeasurement[MNKey][fragKey]
    siteElements = ci.strSiteElements(molecularDataFrame)
    
    massPlot, relAbundPlot, subPlot = [], [], []
    for subKey, observation in isotopologueData.items():
        # Simplified isotopologue mass calculation
        isotopologues_df = pd.DataFrame.from_dict(MNDict[MNKey], orient='index')
        isotopologue_str = isotopologues_df.loc[isotopologues_df[fragKey + '_01 Subs'] == subKey,
                                                fragKey + '_01 Identity'].iloc[0]
        mass = fas.computeMass(isotopologue_str, siteElements) + massError

        if observation['Rel. Abundance'] > lowAbundanceCutOff:
            massPlot.append(mass)
            relAbundPlot.append(observation['Rel. Abundance'])
            subPlot.append(subKey)

    fig, ax = plt.subplots(figsize = (10,4))
    ax.vlines(massPlot, 0, relAbundPlot)
    ax.set_title(MNKey + ' ' + fragKey)

    ax.set_xticks(massPlot)
    labels = [str(round(x,5)) +'\n' + y for x,y in zip(massPlot,subPlot)]
    ax.set_xticklabels(labels,rotation = 45)
    
    if xlim:
        ax.set_xlim(*xlim)
    if ylim:
        ax.set_ylim(*ylim)
    ax.set_ylabel("Relative Abundance")

    if outputIsotopologs:
        constructIsotopologsTSV(fragKey, subPlot, massPlot, tolerance = tolerance, file_path = 'isotopologs' + fragKey + '.tsv')
    
    plt.show()

def constructIsotopologsTSV(parentName, subList, massList, tolerance = 0.1, file_path = 'isotopologs.tsv'):
    '''
    Construct an 'isotopologs.tsv' file that can be used for isoX. This function is housed in 'spectrumVis' because it generates the .tsv based on the predictions used in the spectrum. 

    Inputs:
        parentName: becomes 'compound' in isoX. Written as 'parent' instead of compound because it could be used to correspond to a fragment. 
        subList: A list of the substitutions; 'isotopolog' column in isoX. 
        massList: A list of masses; 'm/z' in isoX. 
        tolerance: Applies the same tolerance to each peak. 
        file_path: Desired output file path. 
    '''
    compoundList = [parentName] * len(subList)
    toleranceList = [tolerance] * len(subList)
    zList = [1] * len(subList)

    # Create a DataFrame
    isotopologsDf = pd.DataFrame({
        'Compound': compoundList,
        'Isotopolog': subList,
        'm/z': massList,
        'Tolerance [mmu]': toleranceList,
        'z': zList
    })

    isotopologsDf['m/z'] = isotopologsDf['m/z'].apply(lambda x: f"{x:.5f}")

    # Write to a .tsv file with a custom header row
    with open(file_path, 'w', newline='') as f:
        f.write('#Compound\tIsotopolog\tm/z\tTolerance [mmu]\tz\n')  # Write the header row
        isotopologsDf.to_csv(f, sep='\t', index=False, header = False)  # Append the DataFrame below the header row