*******************************
 To perform Data vs Lumi study
*******************************

-  enter in JSONCalc/ directory and set 
       export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.1.7/bin:$PATH

-  then run the scripts:
   -   ./CalcJSON.sh                       to convert JSON file in txt
   -   python SplitPerLumi.py              to divide JSON into lumi block 

-  exit the JSONCalc/ directory and run the scripts:
   -   python ExtractHistos.py             to store data into histo for each lumi block (needs some time)
   -   python DataVsLumi_modifies.py       to plot Data vs Lumi


if you need to fit DATA or MC histos, not divided per lumisection:

-  python FitDATA.py
-  python FitMC.py

Then you can use:

-  python yellowPlots.py           to plot Data vs MC histos 

************************	
