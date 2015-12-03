#!/usr/bin/env python
"""
BuildIsisRNA_Tools builds IsisRNA_Tools for each Isis RNA workflow.
"""
import os
import os.path
import sys
import getopt
import getpass
import glob
import shutil
import traceback
import subprocess

class Workflows:
    "List the RNA workflow names once and for all, to avoid magic strings"
    RNAQuantWorker = "RNAQuantWorker"
    SmallRNAWorker = "SmallRNAWorker"
    TargetedRnaSeqWorker = "TargetedRnaSeqWorker"
    WholeGenomeRnaSeqWorker = "WholeGenomeRnaSeqWorker"
    ZodiacRNAWorker = "ZodiacRNAWorker"
    RNA = [TargetedRnaSeqWorker, WholeGenomeRnaSeqWorker, SmallRNAWorker, RNAQuantWorker, ZodiacRNAWorker]

class Preparer:
    def __init__(self):
        self.ReleaseFolder = None
        self.DirectoryName = None
        self.CleanOnly = False

    def SafeDelete(self, Path):
        "Delete a file/dir (or list of files/dirs), if it exists, swallowing any exceptions"
        if type(Path) == type([]):
            for FilePath in Path:
                self.SafeDelete(FilePath)
            return
        if not os.path.exists(Path):
            return
        try:
            if os.path.isdir(Path):
                shutil.rmtree(Path)
            else:
                os.remove(Path)
        except:
            print "Error deleting %s: %s" % (Path, sys.exc_info()[0])

    def SafeDeleteAll(self, folder, skip=[]):
        for filen in os.listdir(folder):
            if filen in skip:
                continue
            path = os.path.join(folder, filen)
            print "Removing", path
            self.SafeDelete(path)

    def PurgeUnusedWorkflows(self):
        "If the workflow was not included in the list above, remove it from this release."
        WorkflowsFolder = os.path.join(self.ReleaseFolder, "Workflows")
        for WorkflowName in os.listdir(WorkflowsFolder):
            if WorkflowName not in Workflows.RNA:
                print "Removing unused workflow", WorkflowName
                shutil.rmtree(os.path.join(WorkflowsFolder, WorkflowName))

    def MakeRNATools(self):
        print "Building IsisRNA_Tools for RNA workflows at: %s" % self.ReleaseFolder
        for Workflow in Workflows.RNA:
            workDir = os.path.join(self.ReleaseFolder, "Workflows", Workflow)
            buildScript = "{0}/build_release.sh".format(workDir)
            basespaceBuildScript = "{0}/build_basespace_release.sh".format(workDir)
            if os.path.exists(basespaceBuildScript):
                buildScript = basespaceBuildScript
            if os.path.exists(buildScript):
                print "Building RNA tools for ", Workflow
                exitCode = os.system("/bin/sh " + buildScript)
                if exitCode != 0:
                    print "Failed building RNA tools for", Workflow
                    sys.exit(-1)
                print "Done building RNA tools for ", Workflow
                # Remove everything but IsisRNA_Tools
                for filen in os.listdir(workDir):
                    if filen == 'IsisRNA_Tools':
                        continue
                    path = os.path.join(workDir, filen)
                    print 'Removing', path
                    self.SafeDelete(path)
            else:
                # Remove workDir
                print 'IsisRNA_Tools build script not found in %s. Removing the workflow.' % (workDir, )
                self.SafeDelete(workDir)
                pass
            
    def CleanupUnusedFiles(self):
        print "Cleaning up files not used to build IsisRNA_Tools at: %s" % self.ReleaseFolder
        skip = ['Workflows', 'BuildIsisRNA_Tools.py']
        self.SafeDeleteAll(self.ReleaseFolder, skip=skip)
        self.PurgeUnusedWorkflows()
        skip = ['build_release.sh', 'build_basespace_release.sh', 'IsisRNA_ToolsSrc']
        for Workflow in Workflows.RNA:
            workDir = os.path.join(self.ReleaseFolder, "Workflows", Workflow)
            self.SafeDeleteAll(workDir, skip=skip)

    def ChangeDirectoryName(self):
        if self.DirectoryName != None:
            newFolder = os.path.join(os.path.dirname(self.ReleaseFolder), self.DirectoryName)
            if self.ReleaseFolder != newFolder:
                print "Changing directory name from %s to %s" % (self.ReleaseFolder, newFolder)
                os.rename(self.ReleaseFolder, newFolder)
                self.ReleaseFolder = newFolder

    def Main(self, Arguments):
        "Main entry point!  Parse command-line arguments, prepare a Linux release folder, and do other goodies like RPM generation"
        self.OriginalFolder = os.getcwd()
        self.ParseCommandLine(Arguments) 
        if (self.ReleaseFolder == None):
            self.PrintUsageInfo()
            return
        if not os.path.exists(self.ReleaseFolder):
            print "Error - no release folder found at '%s'"%self.ReleaseFolder
            return

        if os.path.isfile(self.ReleaseFolder):
            # Assume that we've been given the path to a release archive (e.g. /illumina/scratch/Isis/Builds/Development/Trunk/2.5.5/2.5.5 - x64.zip)
            ZipPath = self.ReleaseFolder
            FileName = os.path.split(ZipPath)[1]
            (Stub, Extension) = os.path.splitext(FileName)
            if Extension != ".zip":
                print "* Error: I don't know how to handle file %s"%ZipPath
                return
            VersionNumber = Stub.replace(r"_x64","").split()[0]
            self.ReleaseFolder = os.path.join(self.IsisRoot, VersionNumber)
            CommandLine = "unzip \"%s\" -d %s"%(ZipPath, self.IsisRoot) # Unzip will add the version# to the folder name automatically
            print CommandLine
            os.system(CommandLine)

        try:
            self.ChangeDirectoryName()
            self.CleanupUnusedFiles()
            if not self.CleanOnly:
                self.MakeRNATools()
                self.SafeDeleteAll(self.ReleaseFolder, skip=['Workflows'])
        finally:
            if os.path.exists(self.OriginalFolder):
                os.chdir(self.OriginalFolder)
        print "...Work complete!"

    def PrintUsageInfo(self):
        print "Usage info:"
        print "  -r ReleaseFolder"
        print "  -d DirectoryName (if you want to rename the release folder for a particular branch)"
        print "  -c only clean up files except for the IsisRNA_ToolsSrc folders"
        print
        print "Example usage:"
        print "python BuildIsisRNA_Tools.py -r /illumina/development/Isis/1.2.3"
        print "python BuildIsisRNA_Tools.py -r \"/illumina/scratch/Isis/Builds/Development/Trunk/2.5.18/2.5.18 - x64.zip\""
        print "python BuildIsisRNA_Tools.py -r \"/illumina/scratch/Isis/Builds/Development/Trunk/2.5.18/2.5.18 - x64.zip\" -d \"2.5.18.Enrichment\""

    def ParseCommandLine(self, Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:d:c")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                self.ReleaseFolder = os.path.abspath(Value) #trouble will arise for relative path due to os.chdir() in some functions.
            elif Option == "-d":
                self.DirectoryName = Value
            elif Option == "-c":
                self.CleanOnly = True
            else:
                print "* Error: Unhandled option!", Option

if __name__ == "__main__":
    FolderPreparer = Preparer()
    FolderPreparer.Main(sys.argv[1:])
