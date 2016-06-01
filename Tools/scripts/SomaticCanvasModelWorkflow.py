#!/usr/bin/env python
import os
import sys
import re
import json
import numpy
import operator
import shutil

ScriptDir=os.path.abspath(os.path.dirname(__file__))
pyFlowPath=os.path.join(ScriptDir,"redist","pyflow","src")
sys.path.append(pyFlowPath)
sys.path.append(ScriptDir)
from pyflow import WorkflowRunner

def isint(x):
    try:
        a = float(x)
        b = int(a)
    except ValueError:
        return False
    else:
        return a == b

def ensureDir(d):
    """
    make directory if it doesn't already exist, raise exception if
    something else is in the way:
    """
    if os.path.exists(d):
        if not os.path.isdir(d) :
            raise Exception("Can't create directory: %s" % (d))
    else :
        os.makedirs(d, mode=0777)


def readTestCorpus(params, testcorpus):
    """
    read in training samples data    
    """
    inFile = open(testcorpus, "rb")
    params.sampleNames = []
    params.sampleDataPath = []
    params.sampleReferenceGenome = []
    params.truthFiles  = []
    params.sampleFilterBed = []
    params.exlcudeRegions = []
    params.sampleSize = 0
    baseDir = os.path.dirname(testcorpus)


    headers = inFile.next().lstrip('#').rstrip('\r\n').split('\t')
    for line in inFile:
        columns = line.rstrip('\r\n').split("\t")
        if len(columns) < 6:
            continue
        if (columns[0][0] == "#"):
            continue
        row = dict(zip(headers, columns))
        params.sampleNames.append(row['sampleNames'])
        params.sampleDataPath.append(os.path.join(baseDir, row['sampleDataPath']))
        params.sampleReferenceGenome.append(os.path.join(baseDir, row['sampleReferenceGenome']))
        params.sampleFilterBed.append(os.path.join(baseDir, row['sampleFilterBed']))
        params.truthFiles.append(os.path.join(baseDir, row['sampleDataPath'], row['truthFile']))
        params.exlcudeRegions.append(os.path.join(baseDir, row['exlcudeRegions']))
        params.sampleSize += 1


def readModelParameters(model_parameters_file):
    with open(model_parameters_file) as parameterFile:
        config = json.load(parameterFile)
        return config


def writeModelParameters(outFilePath, modelParameters):
    with open(outFilePath, 'w') as outFile:
        json.dump(modelParameters, outFile, indent=4)


def updateModelParameters(parameterConfig, bestOptimizationParameters, bestOptimizationParameterValues):
    for i in range(len(bestOptimizationParameters)): 
        parameterConfig[bestOptimizationParameters[i]] = bestOptimizationParameterValues[i]


def mutateModelParameters(parameterConfig, parametername, parameterValueMin, parameterValueMax, parameterConfigFile):
    if  (parameterConfig.get(parametername)!= None):
        newValue = numpy.random.uniform(parameterValueMin, parameterValueMax)
        if isint(parameterValueMax) & isint(parameterValueMin):
            newValue = int(round(newValue))
        else:
            newValue = round(newValue, 3)
        parameterConfig[parametername] = newValue
        return newValue
    else:
        return sys.exit("Parameter %s does not exist in configuration file %s\n" % parametername, parameterConfigFile)


def parseEvaluateCNV(params, tmpPath, outFilePath):
    # Get results from the output file:
    outFile =  open(outFilePath, 'a')
    for sampleIndex in range(params.sampleSize):
        inputPath = os.path.join(tmpPath, params.sampleNames[sampleIndex], 'Results.txt')
        inFile = open(inputPath, 'r')
        Accuracy = "N/A"
        DirectionAccuracy = "N/A"
        Recall = "N/A"
        Precision = "N/A"
        for fileLine in inFile.xreadlines():
            line = fileLine.strip().split("\t")
            if line[0] == "Accuracy":
                Accuracy = line[1]
            if line[0] == "DirectionAccuracy":
                DirectionAccuracy = line[1]
            if line[0] == "Recall":
                Recall = line[1]
            if line[0] == "Precision":
                Precision = line[1]
        # Report results:
        ResultString = params.sampleNames[sampleIndex] + "\t" + str(Accuracy) + "\t" + str(DirectionAccuracy) + "\t" + str(Recall) + "\t" + str(Precision)
        outFile.write(ResultString + "\n")
    outFile.close()


def getBestparameters(resultsPath, nbestParams):
    """
    Best new parameters for a given iteration   
    """
    medianMeanAccuracy = []
    meanAccuracy = []
    optimizationStepParameter = []
    optimizationStepParameterValue = []
    inFile = open(resultsPath, "r")
    for fileLine in inFile.xreadlines():
        line = fileLine.strip().split("\t")
        medianMeanAccuracy.append((float(line[1]) + float(line[2]))/2.0)
        optimizationStepParameter.append(line[0].split("_")[0])
        if isint(line[0].split("_")[1]):
            optimizationStepParameterValue.append(int(line[0].split("_")[1]))
        else:
            optimizationStepParameterValue.append(float(line[0].split("_")[1]))
    sortedAccuracy = sorted(range(len(medianMeanAccuracy)), key=lambda k: medianMeanAccuracy[k],reverse=True)
    selectedOptimizationStepParameters = [optimizationStepParameter[i] for i in sortedAccuracy[0:nbestParams]]
    selectedOptimizationStepParameterValues = [optimizationStepParameterValue[i] for i in sortedAccuracy[0:nbestParams]]
    return selectedOptimizationStepParameters, selectedOptimizationStepParameterValues


def parseOptimizationStep(params, inFilePath, outFilePath):
    """
    Average over EvaluateCNV metrics for all samples for a given optimization step   
    """

    inFile =  open(inFilePath, 'r') 
    outFile =  open(outFilePath, 'a') 
    Accuracy = []
    DirectionAccuracy = []
    Recall = []
    Precision = []
    for fileLine in inFile.xreadlines():
        line = fileLine.strip().split("\t")
        Accuracy.append(float(line[1]))
        DirectionAccuracy.append(float(line[2]))
        Recall.append(float(line[3]))
        Precision.append(float(line[4]))

    medianAccuracy = numpy.median(Accuracy)
    meanAccuracy = numpy.mean(Accuracy)
    medianDirectionAccuracy = numpy.median(DirectionAccuracy)
    meanDirectionAccuracy = numpy.mean(DirectionAccuracy)
    medianRecall = numpy.median(Recall)
    meanRecall = numpy.mean(Recall)
    medianPrecision = numpy.median(Precision)
    meanPrecision = numpy.mean(Precision)
    ResultString = params.currentParameter + "_" + str(params.currentParameterValue) + "\t" + str(medianAccuracy) + "\t" + str(meanAccuracy) + "\t" + str(medianDirectionAccuracy) + "\t" + str(meanDirectionAccuracy) + "\t" + str(medianRecall) + "\t" + str(meanRecall) + "\t" + str(medianPrecision) + "\t" + str(meanPrecision)
    outFile.write(ResultString + "\n")
    outFile.close()


class RunSomaticCanvasWorkflow(WorkflowRunner):
    """Workflow for running CanvasSomaticCaller"""

    def __init__(self, params):
        self.params = params

    def workflow(self):
        parameterStep = "Parameter_" + self.params.currentParameter + "_" + str(self.params.currentParameterValue)
        otimizationTaskID = "Iteration_" + str(self.params.iteration) + parameterStep
        tmpPath = os.path.join(self.params.outputPath, "Iteration_" + str(self.params.iteration), parameterStep)
        configPath = os.path.join(tmpPath, "SomaticCallerParameters.json")
        ensureDir(tmpPath)
        canvasBinary = os.path.join(self.params.executablePath, "CanvasSomaticCaller.exe")
        print tmpPath
        for sampleIndex in range(self.params.sampleSize):
            canvasTaskID = "canvasTaskID_" + self.params.sampleNames[sampleIndex]
            outputPath = os.path.join(tmpPath, self.params.sampleNames[sampleIndex])
            ensureDir(outputPath)
            canvasTask = "/illumina/sync/software/unofficial/Isas/packages/mono-4.0.2/bin/mono %s" % canvasBinary
            canvasTask += " -v %s/VFResultsTumor.txt.gz" % self.params.sampleDataPath[sampleIndex]
            canvasTask += " -i %s/Tumor.partitioned" % self.params.sampleDataPath[sampleIndex]
            canvasTask += " -o %s" % os.path.join(outputPath, "CNV.vcf.gz")
            canvasTask += " -b %s" % self.params.sampleFilterBed[sampleIndex]
            canvasTask += " -n Tumor"
            canvasTask += " -r %s" % self.params.sampleReferenceGenome[sampleIndex]
            canvasTask += " -c %s" % configPath
            canvasTask += " -t %s" % self.params.truthFiles[sampleIndex]
            self.addTask(canvasTaskID, canvasTask, memMb = self.params.memoryCore, retryMax = 3, retryMode = "all")


class RunEvaluateCNV(WorkflowRunner):
    """Workflow for running RunEvaluateCNV"""

    def __init__(self, params):
        self.params = params

    def workflow(self):
        parameterStep = "Parameter_" + self.params.currentParameter + "_" + str(self.params.currentParameterValue)
        otimizationTaskID = "Iteration_" + str(self.params.iteration) + parameterStep
        tmpPath = os.path.join(self.params.outputPath, "Iteration_" + str(self.params.iteration), parameterStep)
        ensureDir(tmpPath)
        for sampleIndex in range(self.params.sampleSize):
            evaluateCNVID = "evaluateCNVTaskID_" + self.params.sampleNames[sampleIndex]
            outputPath = os.path.join(tmpPath, self.params.sampleNames[sampleIndex])
            ensureDir(outputPath)
            inputPath = os.path.join(outputPath, "CNV.vcf.gz")
            evaluateCNVtask = "/illumina/sync/software/unofficial/Isas/packages/mono-4.0.2/bin/mono %s" % self.params.evaluateCNVPath + " "
            evaluateCNVtask += self.params.truthFiles[sampleIndex] + " "
            evaluateCNVtask += inputPath + " "
            evaluateCNVtask += self.params.exlcudeRegions[sampleIndex] + " "
            evaluateCNVtask +=  os.path.join(outputPath, "Results.txt")
            self.addTask(evaluateCNVID, evaluateCNVtask, memMb = self.params.memoryCore, retryMax = 3, retryMode = "all")


class ParseEvaluateCNV(WorkflowRunner):
    """Workflow for parsing RunEvaluateCNV """

    def __init__(self, params):
        self.params = params

    def workflow(self):
        parameterStep = "Parameter_" + self.params.currentParameter + "_" + str(self.params.currentParameterValue)
        tmpPath = os.path.join(self.params.outputPath, "Iteration_" + str(self.params.iteration), parameterStep)
        outFilePathEvaluateCNV = os.path.join(tmpPath, 'Results.txt')
        outFilePathOptimizationStep = os.path.join(self.params.outputPath, "Iteration_" + str(self.params.iteration), 'Results.txt')
        parseEvaluateCNV(self.params, tmpPath, outFilePathEvaluateCNV)
        parseOptimizationStep(self.params, outFilePathEvaluateCNV, outFilePathOptimizationStep)


class FullWorkflow(WorkflowRunner):
    """Complete workflow for training Canvas somatic model"""

    def __init__(self, params):
        self.params = params

    def workflow(self):
        somaticCanvasWorkflow = self.addWorkflowTask("runSomaticCanvasWorkflow", RunSomaticCanvasWorkflow(self.params))
        runEvaluateCNVWorkflow = self.addWorkflowTask("runEvaluateCNV", RunEvaluateCNV(self.params), dependencies=somaticCanvasWorkflow)
        self.addWorkflowTask("parseEvaluateCNV", ParseEvaluateCNV(self.params), dependencies=runEvaluateCNVWorkflow)


class OptimizeSomaticCanvasWorkflow(WorkflowRunner):
    """
    Alter parameters at each iteration of optimization steps 
    """

    def __init__(self, params):
        self.params = params
        self.params.otimizationStep = 1


    def workflow(self):
        ensureDir(self.params.outputPath)
        if self.params.iteration > 0:
            configPath = os.path.join(self.params.outputPath, "Iteration_" + str(self.params.iteration - 1))
        else:
            configPath = self.params.configPath
        parameterConfigFile = os.path.join(configPath, "SomaticCallerParameters.json")
        for keyModelParameters, valueModelParameters in self.params.modelParametersSet.iteritems():
            parameterConfig = readModelParameters(parameterConfigFile)
            self.params.currentParameter  = keyModelParameters
            newModelParameter = mutateModelParameters(parameterConfig, keyModelParameters, valueModelParameters[0], valueModelParameters[1], parameterConfigFile) 
            self.params.currentParameterValue = newModelParameter
            parameterStep = "Parameter_" + self.params.currentParameter + "_" + str(self.params.currentParameterValue)
            parameterConfigPath = os.path.join(self.params.outputPath, "Iteration_" + str(self.params.iteration), parameterStep)
            ensureDir(parameterConfigPath)
            tmpParameterConfigFile = os.path.join(parameterConfigPath, "SomaticCallerParameters.json")
            writeModelParameters(tmpParameterConfigFile, parameterConfig)
            taskID = "WorkflowOptimizatioNStep_" + str(self.params.otimizationStep)
            self.addWorkflowTask(taskID, FullWorkflow(self.params))
            self.params.otimizationStep += 1


class ParseOptimizeSomaticCanvasWorkflow(WorkflowRunner):
    """
    Identify best parameter combination and update SomaticCallerParameters configuration file 
    """

    def __init__(self, params):
        self.params = params

    def workflow(self):
        resultsPath     = os.path.join(self.params.outputPath, "Iteration_" + str(self.params.iteration), "Results.txt")
        bestOptimizationParameters, bestOptimizationParameterValues = getBestparameters(resultsPath, self.params.nbestParams)
        currentIterationPath  = os.path.join(self.params.outputPath, "Iteration_" + str(self.params.iteration))
        if self.params.iteration > 0:
            previousIterationPath = os.path.join(self.params.outputPath, "Iteration_" + str(self.params.iteration - 1))
        else:
            previousIterationPath = self.params.configPath
        parameterConfigFile = os.path.join(previousIterationPath, "SomaticCallerParameters.json")
        parameterConfig = readModelParameters(parameterConfigFile)
        updateModelParameters(parameterConfig, bestOptimizationParameters, bestOptimizationParameterValues)
        parameterConfigFile = os.path.join(currentIterationPath, "SomaticCallerParameters.json")
        writeModelParameters(parameterConfigFile, parameterConfig)


class OptimizeSomaticCanvasFullWorkflow(WorkflowRunner):
    """
    Runs FullWorkflow across all training samples.
    Each optimization step uses different alterations of model parameters.
    """

    def __init__(self, params):
        self.params = params
        self.params.otimizationStep = 1

    def workflow(self):     
        totalIterations = 40
        currentIteration = 1
        optimizeWorkflows = list()
        parseWorkflows = list()
        taskID = "OptimizeSomaticCanvasWorkflow_0"
        self.params.iteration = 0
        print taskID
        optimizeWorkflows.append(self.addWorkflowTask(taskID, OptimizeSomaticCanvasWorkflow(self.params)))
        taskID = "ParseOptimizeSomaticCanvasWorkflow_0" 
        print taskID
        parseWorkflows.append(self.addWorkflowTask(taskID, ParseOptimizeSomaticCanvasWorkflow(self.params), optimizeWorkflows[currentIteration-1]))
        while (currentIteration < totalIterations):
            self.params.iteration = currentIteration
            taskID = "OptimizeSomaticCanvasWorkflow_" + str(currentIteration)
            optimizeWorkflows.append(self.addWorkflowTask(taskID, OptimizeSomaticCanvasWorkflow(self.params), dependencies = parseWorkflows[currentIteration-1]))
            taskID = "ParseOptimizeSomaticCanvasWorkflow_" + str(currentIteration)
            parseWorkflows.append(self.addWorkflowTask(taskID, ParseOptimizeSomaticCanvasWorkflow(self.params), dependencies = optimizeWorkflows[currentIteration]))
            currentIteration += 1