"""
Parameter sweep tool: Invoke Canvas across a smoke test, using various parameter combinations, to review overall accuracy.
"""
import os
import traceback

# Set up which parameters to vary, and which values to try:
CanvasCallerParameters = []
#CanvasCallerParameters.append(("-D", [1.25, 1.5, 1.75, 2.0, 2.25, 2.5]))
#CanvasCallerParameters.append(("-C", [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0]))
#CanvasCallerParameters.append(("-M", [500, 2000, 5000, 10000, 20000]))

CanvasCallerParameters.append(("-D", [1.25, 1.75, 2.5]))
CanvasCallerParameters.append(("-C", [0.25, 1.0, 2.0]))
CanvasCallerParameters.append(("-M", [500, 2000, 50000]))

CommandLines = []
Descriptions = []

# Main test driver: Consider all combinations
CurrentIndexes = [0] * len(CanvasCallerParameters)
while 1:
    #print CurrentIndexes
    # Build the command-line for this parameter combo:
    CommandLine = "-x \""
    Description = ""
    for Index in range(len(CanvasCallerParameters)):
        Value = CanvasCallerParameters[Index][1][CurrentIndexes[Index]]
        CommandLine += "%s %s "%(CanvasCallerParameters[Index][0], Value)
        Description += "%s\t"%(Value)
    CommandLine += "\""
    Descriptions.append(Description)
    CommandLines.append(CommandLine)
    # Iterate to the next combination:
    Index = len(CanvasCallerParameters) - 1
    while Index >= 0:
        CurrentIndexes[Index] += 1
        #print Index, CurrentIndexes[Index],
        if CurrentIndexes[Index] < len(CanvasCallerParameters[Index][1]):
            break
        CurrentIndexes[Index] = 0
        Index -= 1
        #print CurrentIndexes
    if Index < 0:
        break

SweepFolder = "ParamSweep"
if not os.path.exists(SweepFolder):
    os.makedirs(SweepFolder)

OutputFile = open(os.path.join(SweepFolder, "ResultsSummary.txt"), "wb")    
for SetIndex in range(len(CommandLines)):
    print CommandLines[SetIndex], Descriptions[SetIndex]
    # Invoke test suite and save results:
    TestOutputPath = os.path.join(SweepFolder, "Results%s.txt"%SetIndex)
    CommandLine = "TestCanvasSomaticCaller.py %s -o %s"%(CommandLines[SetIndex], TestOutputPath)
    print CommandLine
    OutputFile.write("%s\t%s"%(SetIndex, Descriptions[SetIndex]))
    try:
        os.system(CommandLine)
    except:
        traceback.print_exc()
        OutputFile.write("FAIL\n")
        continue
    File = open(TestOutputPath, "rb")
    Count = 0
    AccuracyMin = 9999
    AccuracyMean = 0
    EventMean = 0
    for FileLine in File.xreadlines():
        Bits = FileLine.split("\t")
        if len(Bits) < 2 or len(Bits[0]) == 0 or Bits[0][0] == "#":
            continue
        print Bits[0], Bits[1]
        try:
            Value = float(Bits[1])
        except:
            continue 
        AccuracyMin = min(AccuracyMin, Value)
        Count += 1
        AccuracyMean += Value
        EventMean += float(Bits[3])
    AccuracyMean /= max(1, Count) 
    EventMean /= max(1, Count)
    OutputFile.write("%s\t%s\t%s\t%s\t\n"%(Count, AccuracyMin, AccuracyMean, EventMean))
    print
    print ">>>%s\t%s\t%s\t%s\t%s\t%s\t"%(SetIndex, Descriptions[SetIndex], Count, AccuracyMin, AccuracyMean, EventMean)
    OutputFile.flush()
                         
