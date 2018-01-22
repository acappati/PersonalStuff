##############################################################
#                                                            #
# Helper script that takes the output of CalcJSON.sh script  #
# sh CalcJSON.sh                                             #
# and splits it into blocks of approximately the same size   #
#                                                            #
##############################################################

LumiBlock = 0.5 #define wanted size in /fb
debug = False

# 2016 data
# input = "LumiCalc_2016data.txt"
# out   = "SplittedBlocks_2016data_0p5_new.txt"
# eraBoundaries = [275376, 276283, 276811, 277420, 278808, 280385] #LAST runs of 2016B, 2016C, 2016D, 2016E, 2016F, 2016G


# 2017 data
input = "LumiCalc_2017data.txt"
out   = "SplittedBlocks_2017data_0p5_new.txt"
eraBoundaries = [299329, 302029, 303434, 304797] #LAST runs of 2017B, 2017C, 2017D, 2017E 


print("---------------------------------------------------------------------------------------------------")
print("Splitting the JSON file in segments with approximate luminosity of " + str(LumiBlock) +" /fb ..." )
print("---------------------------------------------------------------------------------------------------")

with open(input, "r") as Inp:
    print("Opening " + input + ".")
    run = []
    fill = []
    time = []
    nls = []
    ncms = []
    delivered = []
    recorded = []
    
    for line in Inp:
        if(line[0] != "#"):
        
            cleared_line = line.split(",")
            run_fill = cleared_line[0].split(":")
            run.append(run_fill[0])
            fill.append(run_fill[1])
            time.append(cleared_line[1])
            nls.append(cleared_line[2])
            ncms.append(cleared_line[3])
            delivered.append((cleared_line[4]))
            recorded.append(float(cleared_line[5]))

sum = 0
run_B = []
run_E = []
fill_B = []
fill_E = []
LumiBlock_A = []
n_runs = 0
totalLumi = 0

currentEra=0


for i in range(0, len(recorded)):
    sum += recorded[i]
    n_runs += 1
    #Force end of block at era boundaries
    newEra=False
    if ((i < len(recorded)-1) and (currentEra < len(eraBoundaries)) and (int(run[i+1]) > eraBoundaries[currentEra])) :
        print "Era: ", currentEra, "Up to run:", eraBoundaries[currentEra], " up to lumi:", totalLumi+sum
        currentEra +=1
        newEra=True
        
    if ( (newEra) or (i == len(recorded)-1) or (sum > 0.95*LumiBlock and sum < 1.1*LumiBlock and sum + recorded[i+1] > 1.2*LumiBlock) or (sum + recorded[i+1] > 1.3*LumiBlock)):
        LumiBlock_A.append(sum)
        run_B.append(run[i+1-n_runs])
        run_E.append(run[i])
        fill_B.append(fill[i+1-n_runs])
        fill_E.append(fill[i])
        n_runs = 0
        totalLumi += sum
        sum = 0
print(min(LumiBlock_A))
print(max(LumiBlock_A))
if(debug):
    for i in range(0, len(LumiBlock_A)):
        print("{0}:{1} - {2}:{3} block with LUMI: {4:2f}".format(run_B[i],fill_B[i],run_E[i],fill_E[i],LumiBlock_A[i]))

output = open(out,"w")
for i in range(0, len(LumiBlock_A)):
    output.write(str(run_B[i]) + " " + str(fill_B[i]) + " " + str(run_E[i]) + " " + str(fill_E[i]) + " " + str(LumiBlock_A[i]) + "\n")

print("Done!\n\n" + str(totalLumi) + " /fb splitted in " + str(len(LumiBlock_A)) + " blocks.\nOutput written in " + out )

