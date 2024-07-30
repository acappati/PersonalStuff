#!/bin/env python3

## short script to read a ROOT file in the nanoAOD format (from CJSLT framework) and select some events

import pandas
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True



def take_event(samplename: str, filename: str):

    # define ROOT dataframe and access the tree 'Events'
    rdf = ROOT.RDataFrame('Events',filename)

    # GetEntries
    nEntries = rdf.Count()
    print(samplename)
    print(f'Events {nEntries.GetValue()}')


    # Check that the event contains a selected candidate, and that
    # passes the required triggers (which is necessary for samples
    # processed with TRIGPASSTHROUGH=True)
    rdf = rdf.Filter("bestCandIdx != -1 && HLT_passZZ4l == true")  # content of the filter must be C++ code

    # Define the best ZZ candidate quantities and add them as columns to the dataframe
    rdf = rdf.Define('ZZMass', 'ZZCand_mass[bestCandIdx]')
    rdf = rdf.Define('Z1flav', 'ZZCand_Z1flav[bestCandIdx]')
    rdf = rdf.Define('Z2flav', 'ZZCand_Z2flav[bestCandIdx]')
    rdf = rdf.Define('KD', 'ZZCand_KD[bestCandIdx]')
    rdf = rdf.Define('nRun', 'run')
    rdf = rdf.Define('nLumi', 'luminosityBlock')
    rdf = rdf.Define('nEvent', 'event')

    #leptons = getLeptons(rdf, 'event')

    # define pandas dataframe columns
    variables = ['ZZMass', 'Z1flav', 'Z2flav', 'KD', 'nRun', 'nLumi', 'nEvent']

    # define pandas dataframe
    pdf = pandas.DataFrame(rdf.AsNumpy(columns=variables))


    # --- select only events with ZZMass close to 125 GeV
    pdf = pdf.loc[(pdf['ZZMass'] > 120) & (pdf['ZZMass'] < 130)]

    # select only events with 2e2mu final state
    # pdf = pdf.loc[((pdf['Z1flav'].abs() == 121) & (pdf['Z2flav'].abs() == 169))]

    # --- select event with specific run:lumi:event numbers
    # pdf = pdf[(pdf['nRun'] == 359694) &
    #           (pdf['nLumi'] == 126) &
    #           (pdf['nEvent'] == 234858814)]
    # this event has ZZMass 199 GeV, not good for display

    # print the first 10 rows of the dataframe
    print(pdf.head(10))



if __name__ == '__main__':

    # MuonEG
    #filename = '/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII_byZ1Z2/2022EE/Data/MuonEG2022F/ZZ4lAnalysis.root'

    # EGamma
    #filename = '/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII_byZ1Z2/2022EE/Data/EGamma2022E/ZZ4lAnalysis.root'

    # Muon
    filename = '/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII_byZ1Z2/2022EE/Data/Muon2022E/ZZ4lAnalysis.root'

    # total data
    #filename = '/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII_byZ1Z2/2022EE/Data/ZZ4lAnalysis.root'

    take_event('Data', filename)
