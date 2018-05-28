#!/bin/bash

# *******************
# usage: 
#    - choose the period of your dataset
#    - ./CalcJSON.sh 
#
# structure:
#    - read file JSON
#    - convert JSON file in txt and compute luminosity
# ********************



# **********************
# Golden Json 2018 data 
# **********************

JSON_path=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_314472-316271_13TeV_PromptReco_Collisions18_JSON.txt


JSON_localName=GoldenJSON_2018data.txt
Output=LumiCalc_2018data.txt


cp $JSON_path $JSON_localName

# normtag file is not available yet
brilcalc lumi -b "STABLE BEAMS" -u /fb -i $JSON_path -o $Output

# norm tag now available 
# brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i $JSON_path -u /fb -o $Output




echo DONE
