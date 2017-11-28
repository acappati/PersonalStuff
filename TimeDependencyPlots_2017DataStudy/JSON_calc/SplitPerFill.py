##############################################################
#                                                            #
# Helper script that takes the output of CalcJSON.sh script  #
# sh CalcJSON.sh                                             #
# and splits it into fill-blocks                             #
#                                                            #
##############################################################

debug = False
input = "LumiCalc_Moriond17.txt"
out = "SplittedFills_Moriond17.txt"

print("--------------------------------------")
print("  Splitting the JSON file in fills ..." )
print("--------------------------------------")

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

cleared_fill = []
for i in fill:
    if i not in cleared_fill:
        cleared_fill.append(i)

sum = 0
run_B = []
run_E = []
fill_B = []
fill_E = []
LumiBlock_A = []
n_runs = 0
totalLumi = 0
j=0

for i in range(0, len(recorded)):
    sum += recorded[i]
    n_runs += 1
    if ( i == len(recorded)-1 or ((fill[i] == cleared_fill[j]) and fill[i+1] != cleared_fill[j])):
        LumiBlock_A.append(sum)
        run_B.append(run[i+1-n_runs])
        run_E.append(run[i])
        fill_B.append(fill[i+1-n_runs])
        fill_E.append(fill[i])
        n_runs = 0
        totalLumi += sum
        sum = 0
        if(j < len(cleared_fill) - 1):
            j+=1

print(min(LumiBlock_A))
print(max(LumiBlock_A))
if(debug):
    for i in range(0, len(LumiBlock_A)):
        print("{0}:{1} - {2}:{3} block with LUMI: {4:2f}".format(run_B[i],fill_B[i],run_E[i],fill_E[i],LumiBlock_A[i]))

output = open(out,"w")
for i in range(0, len(LumiBlock_A)):
    output.write(str(run_B[i]) + " " + str(fill_B[i]) + " " + str(run_E[i]) + " " + str(fill_E[i]) + " " + str(LumiBlock_A[i]) + "\n")

print("Done!\n\n" + str(totalLumi) + " /fb splitted in " + str(len(LumiBlock_A)) + " blocks.\nOutput written in " + out + ".")

