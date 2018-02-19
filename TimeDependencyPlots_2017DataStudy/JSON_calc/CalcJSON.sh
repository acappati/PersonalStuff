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



# ****************
# Moriond2017 Json; full 2016 dataset; 05 Dec, 2016; 36.8/fb
# ****************
# JSON_path=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt 

# JSON_localName=GoldenJSON_2016data.txt
# Output=LumiCalc_2016data.txt

# NORMTAG_path=/afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json

# cp $JSON_path $JSON_localName
# brilcalc lumi -b "STABLE BEAMS" -i $JSON_localName --normtag $NORMTAG_path -u /fb -o $Output


# ****************
# Golden Json 2017 data 
# ****************
#JSON_path=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-305636_13TeV_PromptReco_Collisions17_JSON.txt # Nov 10, 2017; 35.88/fb
#JSON_path=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Final/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt # Dec 15, 2017; 41.86/fb
JSON_path=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt # Feb 2, 2018; 41.37/fb


JSON_localName=GoldenJSON_2017data.txt
Output=LumiCalc_2017data.txt


cp $JSON_path $JSON_localName
# normtag file is not available yet
brilcalc lumi -b "STABLE BEAMS" -u /fb -i $JSON_path -o $Output



# ****************************
# 2017 data different periods
# ****************************

# *** Run2017B  
# JSON_path=Cert_294927-301141_13TeV_PromptReco_Collisions17_JSON_B.txt 
# JSON_localName=GoldenJSON_2017data_B.txt
# Output=LumiCalc_2017data_B.txt

# *** Run2017C - PromptReco-v1 
# JSON_path=Cert_294927-301141_13TeV_PromptReco_Collisions17_JSON_Cv1.txt 
# JSON_localName=GoldenJSON_2017data_Cv1.txt
# Output=LumiCalc_2017data_Cv1.txt

# *** Run2017C - PromptReco-v2
#JSON_path=Cert_294927-301141_13TeV_PromptReco_Collisions17_JSON_Cv2.txt  
#JSON_localName=GoldenJSON_2017data_Cv2.txt
#Output=LumiCalc_2017data_Cv2.txt

# cp $JSON_path $JSON_localName
# # normtag file is not available yet
# brilcalc lumi -b "STABLE BEAMS" -u /fb -i $JSON_path -o $Output

echo DONE
