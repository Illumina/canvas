#!/usr/bin/env python

import os, sys, json, copy
import multiprocessing
from psutil import virtual_memory
from SomaticCanvasModelWorkflow import *
import argparse


class Params:
    """Object to store all the parameters of the SomaticCanvasModel workflow"""
    
    pass


def getParams():
    """Get SomaticCanvasModel parameters from command line and configuration file into a SimParams object"""

    params = Params()

    # Parameters from command line
    options = getCommandlineOptions()
    params.executablePath = options.executablePath
    params.outputPath = options.outputPath
    params.configPath = options.configPath
    params.trainingSamples = options.trainingSamples
    params.evaluateCNVPath = options.evaluateCNVPath
    with open(options.modelParametersSet) as modelParametersSetFile:
        params.modelParametersSet = json.load(modelParametersSetFile)
    params.mode = options.mode
    params.nbestParams = options.nbestParams
    params.crossValidation = options.crossValidation
    params.currentParameter = ""
    params.currentParameterValue = 0
    # read training samples 
    readTestCorpus(params, params.trainingSamples)
    return params

def getCommandlineOptions() :
    """
    get user options from the command-line and validate required arguments
    """

    parser = argparse.ArgumentParser(prog='optimizeSomaticCanvasModel 0.1', conflict_handler='resolve')

    # Full build required arguments
    parser.add_argument('-i', '--input', dest='trainingSamples', action='store', help='Isis executable path')
    parser.add_argument('-e', '--executablePath', dest='executablePath', action='store', help='Isis executable path')
    parser.add_argument('-o', '--outputPath', dest='outputPath', help='Output directory for Canvas builds')
    parser.add_argument('-c', '--configPath', dest='configPath', help='Canvas somatic model parameters config file')
    parser.add_argument('-v', '--evaluateCNV', dest='evaluateCNVPath', help='Canvas somatic model parameters configFile')
    parser.add_argument('-m', '--mode', dest='mode', default='sge', choices=['sge', 'local'], help="select run mode (local|sge, default=sge)")
    parser.add_argument('-p', '--modelParametersSet', dest='modelParametersSet', help="Path to .json file with Canvas model parameters")
    parser.add_argument('--email', dest='mailTo', action='store', help='e-mail to notify on job completion or error (may be specified more than once)')
    parser.add_argument('--nbestParams', dest='nbestParams', action='store', default = 2, help='Number of best best parameters to keep at each iteration')
    parser.add_argument('--crossValidationFraction', dest='crossValidation', action='store', default = 0.2, help='Proportion of samples to use in testing')


    options = parser.parse_args()

    # check for any errors or missing requirements from argument parsing:
    if (options.trainingSamples is None) :
        print "\nTraining samples file is not specified!\n\n"
        parser.print_help()
        sys.exit(2)

    if (options.executablePath is None):
        print "\nCanvas executable path is not specified!\n\n"
        parser.print_help()
        sys.exit(2)

    if (options.outputPath is None):
        print "\nOutput directory is not specified!\n\n"
        parser.print_help()
        sys.exit(2)

    if (options.configPath is None):
        print "\nCanvas somatic model parameters config file is not specified!\n\n"
        parser.print_help()
        sys.exit(2)

    if (options.evaluateCNVPath is None):
        print "\nEvaluateCNV path is not specified!\n\n"
        parser.print_help()
        sys.exit(2)

    if (options.modelParametersSet is None):
        print "\nModelParametersSet path is not specified!\n\n"
        parser.print_help()
        sys.exit(2)

    return options


def main():

    params = getParams()
    wflow = OptimizeSomaticCanvasFullWorkflow(params)
    if params.mode == "local":  
        nCores=multiprocessing.cpu_count()
        memoryTotal = virtual_memory().total >> 20
        params.memoryCore  = int(memoryTotal/nCores)
    else:
        nCores = 128
        memoryTotal = "unlimited"
        params.memoryCore = 3048
    wflow.run(mode = params.mode, dataDirRoot = params.outputPath, isContinue = "Auto", isForceContinue = True, isQuiet = True, nCores = nCores, memMb = memoryTotal, retryMode = "all")


if __name__=="__main__":
    main()
