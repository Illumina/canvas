#!/usr/bin/env python
"""
UpdateRelease prepares a Linux release folder of Isas.
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
    "List the workflow names once and for all, to avoid magic strings"
    AmpliconDSWorker = "AmpliconDSWorker"
    ChIPSeqWorker = "ChIPSeqWorker"
    EnrichmentWorker = "EnrichmentWorker"
    ForensicsWorker = "ForenSeqWorker"
    GenerateFASTQWorker = "GenerateFASTQWorker"
    MetagenomicsWorker = "MetagenomicsWorker"
    MethylSeqWorker = "MethylSeqWorker"
    MitoWorker = "MitoWorker" 
    ResequencingWorker = "ResequencingWorker"
    RNAQuantWorker = "RNAQuantWorker"
    SmallPedigreeWorker = "SmallPedigreeWorker"
    SmallRNAWorker = "SmallRNAWorker"
    TargetedRnaSeqWorker = "TargetedRnaSeqWorker"
    TruSeqAmpliconWorker = "TruSeqAmpliconWorker"
    TumorNormalWorker = "TumorNormalWorker"
    WholeGenomeRnaSeqWorker = "WholeGenomeRnaSeqWorker"
    ZodiacRNAWorker = "ZodiacRNAWorker"
    VeriSeqPGSWorker = "VeriSeqPGSWorker"
    IsaacWorkflows = [EnrichmentWorker, ResequencingWorker]
    ReseqWorkflows = [EnrichmentWorker, VeriSeqPGSWorker,
                      ResequencingWorker, ChIPSeqWorker]
    MantaWorkflows = [ResequencingWorker, TumorNormalWorker, EnrichmentWorker, SmallPedigreeWorker]
    MantaRNAWorkflows = [WholeGenomeRnaSeqWorker, ZodiacRNAWorker]
    #NOTE: opt-in paradigm. If your workflow is not listed here it will be removed from the release
    Isas = [AmpliconDSWorker, ChIPSeqWorker, EnrichmentWorker,
        GenerateFASTQWorker, MetagenomicsWorker, MethylSeqWorker, MitoWorker,
        VeriSeqPGSWorker, ResequencingWorker,
        SmallPedigreeWorker, TargetedRnaSeqWorker, TruSeqAmpliconWorker, TumorNormalWorker]
    #Workflows in specialized Isas releases
    RNA = [GenerateFASTQWorker, TargetedRnaSeqWorker, WholeGenomeRnaSeqWorker, SmallRNAWorker, RNAQuantWorker, ZodiacRNAWorker]
    Forensics = [GenerateFASTQWorker, ForensicsWorker]
    All = None

class Preparer:
    def __init__(self):
        self.ReleaseFolder = None
        self.ReleaseType = "isas"
        self.DirectoryName = None
        self.IsasRoot = "/illumina/development/Isas"
        self.IsasRuntimePackages = "/illumina/sync/software/unofficial/Isis/packages"
        self.IsasInstallRoots = [self.IsasRoot, "/illumina/sync/software/unofficial/Isis", "/illumina/sync/software/unofficial/IsisRNA"]
        self.IsaacVersion = "03.15.11.19"
        self.bcl2fastqVersion = "v2.18.0.6"
        self.PumaVersion = "Metrics-1.0.10.0_PUMA-00.15.11.02"
        self.MonoVersion = "4.0.2" 
        self.PythonVersion = "2.7"
        self.PythonEnvironment =  "envs/20151030/"
        self.MarkDuplicatesVersion = "1.0.3"

        # ideally starling and strelka should have the same version
        # since they are taken from the same starka release
        self.StarlingVersion = "2.4.3"
        self.StrelkaVersion = "2.4.3"
        
        self.MantaFolder = '/illumina/development/manta/manta-0.29.1'
        self.RNAMantaFolder = '/illumina/sync/software/unofficial/IsisRNA/IsisRNATools/Manta.rnamaster.20151104'

        self.IncludeRPMBuild = 0
        self.IncludeRsync = 0
        self.ForceOverwrite = False
        self.Sge = False
        self.BuildIsasRNA_Tools = True
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
    def PurgeUnusedWorkflows(self):
        "If the workflow was not included in the list above, remove it from this release."
        WorkflowsFolder = os.path.join(self.ReleaseFolder, "Workflows")
        for WorkflowName in os.listdir(WorkflowsFolder):
            if WorkflowName not in Workflows.All:
                print "Removing unused workflow", WorkflowName
                shutil.rmtree(os.path.join(WorkflowsFolder, WorkflowName))
    def PurgeUnmanagedWindowsCode(self, TidyFolders):
        "Unmanaged code isn't portable.  Remove the unmanaged Windows executables; later we'll replace them with (unmanaged) Linux executables."
        Executables = "AlignmentResolver bgzip9 bowtie-build-debug bowtie-build bowtie-debug bowtie bowtie-inspect-debug bowtie-inspect bwa bwa-0.6.1 ExtractUnalignedReads gatk_to_gvcf MarkDuplicates samtools starling2 tabix".split()
        DLLs = "FileCompression drmaa Locfit PosixThreads".split()
        for WorkflowName in TidyFolders:
            WorkflowFolder = os.path.join(self.ReleaseFolder, "Workflows", WorkflowName)
            if not os.path.exists(WorkflowFolder):
                continue
            for Executable in Executables:
                for Extension in ("exe", "pdb"):
                    TempPath = os.path.join(WorkflowFolder, "%s.%s"%(Executable, Extension))
                    self.SafeDelete(TempPath)
            for DLL in DLLs:
                for Extension in ("dll", "exp", "lib", "pdb"):
                    TempPath = os.path.join(WorkflowFolder, "%s.%s"%(DLL, Extension))
                    self.SafeDelete(TempPath)
    def CopyLatest(self, SourcePath, targetWorkflows, OverrideTargetName = None):
        "Copy SourcePath to the each workflow in the TargetWorkflows list"
        # If we got a string, change it to a length-1 list of workflow names:
        if (type(targetWorkflows) == type("")):
            targetWorkflows = [targetWorkflows]
        if not os.path.exists(SourcePath):
            print "* Error: Not copying expected file %s"%SourcePath
            return
        if OverrideTargetName != None:
            FileName = OverrideTargetName
        else:
            FileName = os.path.split(SourcePath)[1]
        for Workflow in targetWorkflows:
            if not Workflow in Workflows.All: continue
            TargetFolder = os.path.join(self.ReleaseFolder, "Workflows", Workflow)
            if not os.path.exists(TargetFolder):
                print "Warning: Missing expected workflow folder at:", TargetFolder
                continue
            TargetPath = os.path.join(TargetFolder, FileName)
            if not os.path.isdir(SourcePath):
                if self.ForceOverwrite or not os.path.exists(TargetPath):
                    if not os.path.exists(os.path.dirname(TargetPath)):
                        os.makedirs(os.path.dirname(TargetPath))
                    print "Copy: %s"%TargetPath
                    shutil.copy2(SourcePath, TargetPath)#Note: copy2() copies file metadata as well
            else:
                self.MirrorFolder(SourcePath, TargetPath)
    def PurgeExisting(self, targetWorkflows, TargetName):
        "Purge TargetName in each workflow in the TargetWorkflows list"
        # If we got a string, change it to a length-1 list of workflow names:
        if (type(targetWorkflows) == type("")):
            targetWorkflows = [targetWorkflows]
        for Workflow in targetWorkflows:
            if not Workflow in Workflows.All: continue
            TargetFolder = os.path.join(self.ReleaseFolder, "Workflows", Workflow)
            if not os.path.exists(TargetFolder):
                print "Warning: Missing expected workflow folder at:", TargetFolder
                continue
            TargetPath = os.path.join(TargetFolder, TargetName)
            self.SafeDelete(TargetPath)
    def MirrorFolder(self, SourceFolder, TargetFolder, ExcludeFiles = None):
        "Copy a directory tree"
        if os.path.isfile(TargetFolder) or os.path.islink(TargetFolder):
            os.remove(TargetFolder)
        if not os.path.exists(TargetFolder):
            os.makedirs(TargetFolder)
        for FileName in os.listdir(SourceFolder):
            SourcePath = os.path.join(SourceFolder, FileName)
            if ExcludeFiles != None and ExcludeFiles.has_key(FileName):
                print "Skip: %s"%SourcePath
                continue
            TargetPath = os.path.join(TargetFolder, FileName)
            if os.path.isdir(SourcePath):
                self.MirrorFolder(SourcePath, TargetPath, ExcludeFiles)
                continue
            if self.ForceOverwrite or not os.path.exists(TargetPath):
                if os.path.islink(SourcePath):
                    if os.path.exists(TargetPath):
                        os.unlink(TargetPath)
                    linkto = os.readlink(SourcePath)
                    os.symlink(linkto, TargetPath)
                    print "Link: %s"%TargetPath
                else:
                    print "Copy: %s"%TargetPath
                    shutil.copy2(SourcePath, TargetPath)#Note: copy2() copies file metadata as well

    def PrepareLinuxRelease(self):
        print "Prepare Isas release folder at: %s"%self.ReleaseFolder
        self.PackagesFolder = os.path.join(self.IsasRoot, "packages")
        self.SourceFolder = os.path.join(self.IsasRoot, "source")
        # Make sure java is up-to-date:
        JavaFolder = os.path.join(self.ReleaseFolder, "Java")
        os.system("rm -rf %s"%JavaFolder)
        os.system("ln -s %s %s"%(os.path.join(self.PackagesFolder, "jdk1.6.0_32"), JavaFolder))
        self.PurgeUnusedWorkflows()
        # include Canvas directory (if it exists) when purging unmanaged code and updating pdb to mdb:
        TidyFolders = Workflows.All[:]
        for Workflow in Workflows.All:
            if os.path.exists(os.path.join(self.ReleaseFolder, "Workflows", Workflow, "Canvas")):
                TidyFolders.append(os.path.join(Workflow, "Canvas"))
        print "Tidy unmanaged code..."
        self.PurgeUnmanagedWindowsCode(TidyFolders)
        print "Copy / softlink files..."
        self.CopyPrerequisitesAlignment()
        self.CopyPrerequisitesVariantCallers()
        self.CopyPrerequisitesMisc()
        self.UpdateIsasConfig()
        self.CopyIsaac()
        self.CopyBcl2Fastq()
        self.CopyExpansionHunter()
        print "Convert pdb to mdb..."
        self.ConvertPdbToMdb(TidyFolders)
        self.CleanTopLevelFolder()
        self.MakeIsasScript()

    def MakeRNATools(self):
        for Workflow in Workflows.RNA:
            if not Workflow in Workflows.All: continue
            workDir = os.path.join(self.ReleaseFolder, "Workflows", Workflow)
            buildScript = os.path.join(workDir, "build_release.sh")
            basespaceBuildScript = os.path.join(workDir, "build_basespace_release.sh")
            if os.path.exists(buildScript):
                if self.BuildIsasRNA_Tools:
                    if self.ReleaseType == 'rna-basespace' and os.path.exists(basespaceBuildScript):
                        buildScript = basespaceBuildScript
                    print "Building RNA tools for ", Workflow
                    exitCode = os.system("/bin/sh " + buildScript)
                    if exitCode != 0:
                        print "Failed building RNA tools for", Workflow
                        sys.exit(-1)
                    print "Done building RNA tools for ", Workflow
                else:
                    rnaToolsSrc = os.path.join(workDir, 'IsasRNA_ToolsSrc')
                    print "Removing RNA tools source and build scripts for", Workflow
                    self.SafeDelete([buildScript, basespaceBuildScript, rnaToolsSrc])

    def ConvertPdbToMdb(self, TidyFolders):
        ExePath = os.path.join(self.IsasRoot, "packages", "mono-%s" %self.MonoVersion, "lib", "mono", "4.5", "pdb2mdb.exe")
        MonoPath = "/illumina/development/Isas/packages/mono-%s/bin/mono" % self.MonoVersion
        for Workflow in TidyFolders + [".."]:
            Folder = os.path.abspath(os.path.join(self.ReleaseFolder, "Workflows", Workflow))
            if not os.path.exists(Folder):
                continue
            os.chdir(Folder)    #make sure Folder is absolute path, otherwise trouble will arise here
            for FileName in os.listdir(Folder):
                (Stub, Extension) = os.path.splitext(FileName)
                if Extension != ".dll" and Extension != ".exe":
                    continue
                FilePath = os.path.join(Folder, FileName)
                PDBPath = os.path.join(Folder, "%s.pdb"%Stub)
                if not os.path.exists(PDBPath): continue
                CommandLine = "%s %s %s > /dev/null &> /dev/null"%(MonoPath, ExePath, FilePath)
                #print CommandLine 
                try:
                    os.system(CommandLine)
                except:
                    pass #traceback.print_exc()
                self.SafeDelete(PDBPath)
    def CleanTopLevelFolder(self):
        "Delete some files that are not needed by Isas.exe (.pdf files and Isas.vshost)"
        self.SafeDelete(glob.glob(os.path.join(self.ReleaseFolder, "*.pdb")))
        self.SafeDelete(glob.glob(os.path.join(self.ReleaseFolder, "Isas.vshost.*")))
    def CopyExpansionHunter(self):
        "Copy ExpansionHunter" 
        ExpansionHunter = "/illumina/development/readrecovery/readrecovery-Dumpster-1.6.2/RepeatExpansion/bin/expansionHunter"
        self.CopyLatest(ExpansionHunter, [Workflows.ResequencingWorker])
    def CopyIsaac(self):
        "Copy in the Isaac binary and related files"
        IsaacFolder = "/illumina/development/iSAAC/iSAAC-%s"%self.IsaacVersion
        for Workflow in Workflows.IsaacWorkflows:
            if not Workflow in Workflows.All: continue
            TargetRoot = os.path.join(self.ReleaseFolder, "Workflows", Workflow, "Isaac")
            # Copy bin/isaac-align and bin/isaac-sort-reference and bin/isaac-unpack-reference:
            TargetFolder = os.path.join(TargetRoot, "bin")
            if not os.path.exists(TargetFolder):
                os.makedirs(TargetFolder)
            for ExeName in ("isaac-align", "isaac-sort-reference", "isaac-unpack-reference"):
                SourcePath = os.path.join(IsaacFolder, "bin", ExeName)
                TargetPath = os.path.join(TargetFolder, ExeName)
                if self.ForceOverwrite or not os.path.exists(TargetPath):
                    print "Copy: %s"%TargetPath
                    shutil.copy2(SourcePath, TargetPath)#Note: copy2() copies file metadata as well
            # Copy other stuff:
            for FolderName in ("xsl", "css"):
                SourceFolder = os.path.join(IsaacFolder, "share", "iSAAC-%s"%self.IsaacVersion, FolderName)
                TargetFolder = os.path.join(TargetRoot, "share", "iSAAC-%s"%self.IsaacVersion, FolderName)
                self.MirrorFolder(SourceFolder, TargetFolder)
    def CopyBcl2Fastq(self):
        "Pull in the correct bcl2fastq version, including reporting-related files"
        Bcl2fastqFolder = "/illumina/development/bcl2fastq2/bcl2fastq2-%s/"%self.bcl2fastqVersion
        TargetFolder = os.path.join(self.ReleaseFolder, "bcl2fastq2")
        self.MirrorFolder(Bcl2fastqFolder, TargetFolder)
    def UpdateIsasConfig(self):
        "Repair a couple of settings in Isas.exe.config"
        FilePath = os.path.join(self.ReleaseFolder, "Isas.exe.config")
        File = open(FilePath, "rb")
        FileLines = File.readlines()
        File.close()
        File = open(FilePath, "wb")
        for FileLine in FileLines:
            if FileLine.find("<add key=\"TempFolder\"") != -1:
                File.write('    <add key="TempFolder" value="/tmp" />\n')
                continue
            if FileLine.find("<add key=\"GenomePath\"") != -1:
                File.write('    <add key="GenomePath" value="%s" />\n' % "/illumina/sync/software/unofficial/Isis/Genomes/")
                continue            
            if FileLine.find("<add key=\"Repository\"") != -1:
                continue
            if self.Sge:
                if FileLine.find("<add key=\"ClusterParallelEnvironmentName\"") != -1:
                    File.write("<add key=\"ClusterQueueName\" value=\"slice.q\"/>\n")
            File.write(FileLine)
        File.close()
    def MakeIsasScript(self):
        "Make script to start an Isas Analysis on linux"
        FilePath = os.path.join(self.ReleaseFolder, "Isas")
        File = open(FilePath, "wb")
        File.write('#!/usr/bin/env bash\n')
        File.write('ISAS_DIR="$( dirname "$( readlink -f "${BASH_SOURCE[0]}" )" )"\n')
        File.write('source ${ISAS_DIR}/UpdateEnvironmentMono ')
        if os.path.dirname(self.ReleaseFolder) in self.IsasInstallRoots:
            File.write('%s' % self.MonoVersion)
        else:
            File.write('%s' % os.path.join(self.IsasRuntimePackages, "mono-%s" % self.MonoVersion))
        python = "python-{0}".format(self.PythonVersion)
        if self.PythonEnvironment: python = os.path.join(python, self.PythonEnvironment)
        File.write(" %s\n" % os.path.join(self.IsasRuntimePackages, python))
        
        File.write('exec mono ${ISAS_DIR}/Isas.exe "$@"\n')
        File.close()
        Command = "chmod 755 \"%s\""%FilePath
        os.system(Command)

    def CopyPrerequisitesAlignment(self):
        self.CopyLatest(os.path.join(self.SourceFolder, "bwa-0.6.1-isis", "bwa"), Workflows.ReseqWorkflows + [Workflows.MitoWorker], "bwa-0.6.1")
        BwaWorkflows = Workflows.ReseqWorkflows + [Workflows.MitoWorker, Workflows.RNAQuantWorker, Workflows.WholeGenomeRnaSeqWorker, Workflows.ZodiacRNAWorker]
        self.CopyLatest(os.path.join(self.SourceFolder, "bwa-0.7.9a-isis", "bwa"), BwaWorkflows)
        BowtieBinDir = "/illumina/scripts/IsisRNA/20130207/IsisRNA_Tools/bin"
        BowtieWorkflows = [Workflows.SmallRNAWorker]
        self.CopyLatest(os.path.join(BowtieBinDir, "bowtie"), BowtieWorkflows)
        self.CopyLatest(os.path.join(BowtieBinDir, "bowtie-inspect"), BowtieWorkflows)
        self.CopyLatest(os.path.join(BowtieBinDir, "bowtie-build"), BowtieWorkflows)
        Bowtie2Workflows = [Workflows.MethylSeqWorker]
        self.CopyLatest(os.path.join(self.SourceFolder, "bowtie2-2.2.2", "bowtie2"), Bowtie2Workflows)	
        self.CopyLatest(os.path.join(self.SourceFolder, "bowtie2-2.2.2", "bowtie2-align-s"), Bowtie2Workflows)
        self.CopyLatest(os.path.join(self.SourceFolder, "scramble-1.13.10", "scramble"), [Workflows.ResequencingWorker])

    def CopyPrerequisitesVariantCallers(self):
        ReseqPlusTruseq = Workflows.ReseqWorkflows + [Workflows.TruSeqAmpliconWorker]
        self.CopyLatest(os.path.join(self.PackagesFolder, "bin", "gatk_to_gvcf"), ReseqPlusTruseq)
        self.FixMantaPermissions()
        self.CopyPrerequisitesPisces()
        self.CopyPrerequisitesStarling()
        self.CopyPrerequisitesStrelka()
        self.CopyPrerequisitesManta()

    def CopyPrerequisitesPisces(self):
        CopyWorkflows = [Workflows.ResequencingWorker, Workflows.TruSeqAmpliconWorker, Workflows.AmpliconDSWorker, Workflows.EnrichmentWorker,Workflows.TargetedRnaSeqWorker, Workflows.TumorNormalWorker]
        PiscesFolder = "/illumina/development/CallSomaticVariants/4.0.13.1"
        for Workflow in CopyWorkflows:
            if not Workflow in Workflows.All: continue
            TargetFolder = os.path.join(self.ReleaseFolder, "Workflows", Workflow)
            if not os.path.exists(TargetFolder):
                print "Warning: Missing expected workflow folder at:", TargetFolder
                continue
            TargetFolder = os.path.join(TargetFolder, "CallSomaticVariants")
            ExcludeFiles = { }
            self.MirrorFolder(PiscesFolder, TargetFolder, ExcludeFiles)        
           
    def CopyPrerequisitesManta(self):
        for Workflow in Workflows.MantaWorkflows + Workflows.MantaRNAWorkflows:
            if not Workflow in Workflows.All: continue
            WorkflowFolder = os.path.join(self.ReleaseFolder, "Workflows", Workflow)
            if not os.path.exists(WorkflowFolder):
                print "Warning: Missing expected workflow folder at:", WorkflowFolder
                continue
            TargetFolder = os.path.join(WorkflowFolder, "Manta")
            ExcludeFiles = { 
                "demo": 1, 
                "checkChromSet.pyc": 1, 
                "configBuildTimeInfo.pyc": 1,
                "configureOptions.pyc": 1,
                "configureUtil.pyc": 1,
                "estimateHardware.pyc": 1,
                "makeRunScript.pyc": 1,
                "mantaOptions.pyc": 1,
                "mantaWorkflow.pyc": 1,
                "workflowUtil.pyc": 1,
                "pyflowConfig.pyc": 1,
                "pyflow.pyc": 1,
                "pyflowTaskWrapper.pyc": 1 
            }
            if Workflow in Workflows.MantaWorkflows:
                self.MirrorFolder(self.MantaFolder, TargetFolder, ExcludeFiles)        
            else:
                self.MirrorFolder(self.RNAMantaFolder, TargetFolder, ExcludeFiles)        

    def CopyPrerequisitesStarling(self):
        ReseqPlusTruseq = Workflows.ReseqWorkflows + [Workflows.TruSeqAmpliconWorker, Workflows.WholeGenomeRnaSeqWorker, Workflows.SmallPedigreeWorker]
        StarkaFolder = "/illumina/development/STARKA/starka-%s" % self.StarlingVersion
        StarkaRnaFolder = "/illumina/sync/software/unofficial/IsisRNA/IsisRNATools/starka.247.20150826" #todo remove once starka-247 is merged into master

        # Todo: remove these two lines once all workflows are using starling pyflow
        self.CopyLatest(os.path.join(StarkaFolder, "libexec", "starling2"), ReseqPlusTruseq)
        self.CopyLatest(os.path.join(StarkaFolder, "share", "config", "model.json"), ReseqPlusTruseq)
        
        # Copy starling pyflow
        for Workflow in ReseqPlusTruseq:
            if not Workflow in Workflows.All: continue
            TargetFolder = os.path.join(self.ReleaseFolder, "Workflows", Workflow)
            if not os.path.exists(TargetFolder):
                print "Warning: Missing expected workflow folder at:", TargetFolder
                continue
            TargetFolder = os.path.join(TargetFolder, "Starling")
            ExcludeFiles = {
                "demo": 1,
                "strelka_config_isaac_default.ini": 1,
                "strelka_config_bwa_default.ini": 1,
                "runStarlingWorkflowDemo.bash": 1,
                "runStrelkaWorkflowDemo.bash": 1,
                "strelkaSiteSimulator": 1,"configureStrelkaNoiseWorkflow.py.ini": 1,
                "configureStrelkaNoiseWorkflow.py": 1,
                "configureStrelkaWorkflow.py": 1,
                "configureStrelkaWorkflow.py.ini": 1,
            }
            if Workflow == Workflows.WholeGenomeRnaSeqWorker:  #todo remove once starka-247 is merged into master        
                self.MirrorFolder(StarkaRnaFolder, TargetFolder, ExcludeFiles)        
            else:
                self.MirrorFolder(StarkaFolder, TargetFolder, ExcludeFiles)        

    def CopyPrerequisitesStrelka(self):
        StarkaFolder = "/illumina/development/STARKA/starka-%s" % self.StrelkaVersion
        TargetFolder = os.path.join(self.ReleaseFolder, "Workflows", Workflows.TumorNormalWorker)
        if not os.path.exists(TargetFolder):
            print "Warning: Missing expected workflow folder at:", TargetFolder
            return
        TargetFolder = os.path.join(TargetFolder, "Strelka")
        ExcludeFiles = {
            "demo": 1,
            "runStarlingWorkflowDemo.bash": 1,
            "runStrelkaWorkflowDemo.bash": 1,
            "configureStarlingWorkflow.py": 1,
            "configureStarlingWorkflow.py.ini": 1,
        }  
        self.MirrorFolder(StarkaFolder, TargetFolder, ExcludeFiles)

    def FixMantaPermissions(self):
        for worker in Workflows.MantaWorkflows + Workflows.MantaRNAWorkflows:
            Root = os.path.join(self.ReleaseFolder, "Workflows", worker, "Manta")
            if not os.path.isdir(Root):
                continue
            for (DirPath, DirNames, FileNames) in os.walk(Root):
                for FileName in FileNames:
                    if os.path.splitext(FileName)[1] == ".py":
                        FilePath = os.path.join(DirPath, FileName)
                        Command = "chmod 755 \"%s\""%FilePath
                        os.system(Command)
            Folder = os.path.join(Root, "libexec")
            for FileName in os.listdir(Folder):
                FilePath = os.path.join(Folder, FileName)
                Command = "chmod 755 \"%s\""%FilePath
                os.system(Command)

    def CopyPrerequisitesMisc(self):
        ReseqPlusTruseq = Workflows.ReseqWorkflows + [Workflows.TruSeqAmpliconWorker]
        TabixWorkflows = Workflows.All
        self.CopyLatest(os.path.join(self.PackagesFolder, "bin", "libFileCompression.so"), Workflows.All)
        # Also, get libFileCompression placed for Canvas:
        for WorkflowName in Workflows.All:
            TempCanvasFolder = os.path.join(self.ReleaseFolder, "Workflows", WorkflowName, "Canvas")
            if not os.path.exists(TempCanvasFolder):
                continue
            SoftLinkPath = os.path.join(TempCanvasFolder, "libFileCompression.so")
            if not os.path.exists(SoftLinkPath):
                Command = "ln -s ../libFileCompression.so %s"%SoftLinkPath
                print Command
                os.system(Command)
        # Also, get libFileCompression placed for Nirvana:
        for WorkflowName in Workflows.All:
            TempNirvanaFolder = os.path.join(self.ReleaseFolder, "Workflows", WorkflowName, "Nirvana")
            if not os.path.exists(TempNirvanaFolder):
                continue
            SoftLinkPath = os.path.join(TempNirvanaFolder, "libFileCompression.so")
            if not os.path.exists(SoftLinkPath):
                Command = "ln -s ../libFileCompression.so %s"%SoftLinkPath
                print Command
                os.system(Command)
        self.CopyLatest(os.path.join(self.PackagesFolder, "bin", "bgzip9"), ReseqPlusTruseq + [Workflows.MethylSeqWorker])
        self.CopyLatest(os.path.join(self.PackagesFolder, "bin", "bgzf_cat"), ReseqPlusTruseq)
        self.CopyLatest("/illumina/thirdparty/samtools/samtools-1.2/bin/samtools",
            Workflows.ReseqWorkflows + [Workflows.TruSeqAmpliconWorker, Workflows.SmallRNAWorker, Workflows.TargetedRnaSeqWorker,
            Workflows.ForensicsWorker, Workflows.AmpliconDSWorker, 
            Workflows.MethylSeqWorker, Workflows.TumorNormalWorker, Workflows.MitoWorker]),
        self.CopyLatest(os.path.join(self.SourceFolder, "tabix-0.2.6", "tabix"), TabixWorkflows)
        self.CopyLatest(os.path.join(self.SourceFolder, "MarkDuplicates", self.MarkDuplicatesVersion, "MarkDuplicates"),
                        ReseqPlusTruseq + [Workflows.MethylSeqWorker, Workflows.MitoWorker, Workflows.RNAQuantWorker, Workflows.WholeGenomeRnaSeqWorker, Workflows.ZodiacRNAWorker])
        self.CopyLatest(os.path.join(self.SourceFolder, "ExtractUnalignedReads-1.0.3", "ExtractUnalignedReads"),
                        Workflows.ReseqWorkflows + [Workflows.SmallRNAWorker])
        self.CopyLatest(os.path.join(self.SourceFolder, "Locfit", "liblocfit.so"), Workflows.TargetedRnaSeqWorker)
        self.PurgeExisting([Workflows.ResequencingWorker, Workflows.EnrichmentWorker, Workflows.TumorNormalWorker, Workflows.WholeGenomeRnaSeqWorker], "Puma")
        self.CopyLatest(os.path.join(self.PackagesFolder, "bin", "Puma", self.PumaVersion), [Workflows.ResequencingWorker, Workflows.EnrichmentWorker, Workflows.TumorNormalWorker, Workflows.WholeGenomeRnaSeqWorker], "Puma")
    def DoRPMBuild(self):
        print "Generating an RPM:"
        CommandLine = "/illumina/development/Isas/rpmbuild/build.sh --user=%s --version=%s"%(getpass.getuser(),
            os.path.split(self.ReleaseFolder)[1])
        print CommandLine
        os.system(CommandLine)
    def DoRsync(self):
        print "Syncing build to other clusters..."
        Clusters=["ussd-services.illumina.com","ukch-prd-lndt.chuk.illumina.com"]
        jobs=[]
        for cluster in Clusters:
            print "\tChecking for ssh key authentication to %s"%(cluster)
            try:
                CommandLine = "ssh -o BatchMode=yes -t -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no %s \"echo \"SSH encryption is setup. We can proceed\"\" &> /dev/null"%(cluster)
                #print CommandLine
                retcode = subprocess.call(CommandLine, shell=True)
                if retcode is not 0:
                    raise Exception()
            except Exception, e:
                print "\tError: key authentication for user %s is not set up to %s. No rsync to that cluster will be done!"%(getpass.getuser(), cluster)
                print "(See 'https://confluence.illumina.com/x/FpJX' for notes on how to set this up)"
                continue          
            CommandLine = "rsync -aqz --no-g %s %s:%s &> /dev/null"%(self.ReleaseFolder, cluster, self.IsasRoot)
            print "\t"+CommandLine
            jobs.append([CommandLine,subprocess.Popen(CommandLine, shell=True)])
        for job in jobs:
            returncode = job[1].wait()
            if returncode is not 0:
                print "failed command: %s"%(job[0])
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

        if self.ReleaseType == 'isas':
            Workflows.All = Workflows.Isas
        if self.ReleaseType == 'rna' or self.ReleaseType == 'rna-basespace':
            Workflows.All = Workflows.RNA
        if self.ReleaseType == 'forensics':
            Workflows.All = Workflows.Forensics
        if self.ReleaseType == 'all':
            Workflows.All = list(set(Workflows.Isas + Workflows.RNA + Workflows.Forensics))

        if os.path.isfile(self.ReleaseFolder):
            # Assume that we've been given the path to a release archive (e.g. /illumina/scratch/Isas/Builds/Development/Trunk/2.5.5/2.5.5 - x64.zip)
            ZipPath = self.ReleaseFolder
            FileName = os.path.split(ZipPath)[1]
            (Stub, Extension) = os.path.splitext(FileName)
            if Extension != ".zip":
                print "* Error: I don't know how to handle file %s"%ZipPath
                return
            VersionNumber = Stub.replace(r"_x64","").split()[0]
            self.ReleaseFolder = os.path.join(self.IsasRoot, VersionNumber)
            CommandLine = "unzip \"%s\" -d %s"%(ZipPath, self.IsasRoot) # Unzip will add the version# to the folder name automatically
            print CommandLine
            os.system(CommandLine)
        if os.path.dirname(self.ReleaseFolder) != self.IsasRoot:
            print "Note: Installing to non-standard location. No rpm builds or rsync will happen."
            self.IncludeRPMBuild = 0
            self.IncludeRsync = 0
        try:
            self.ChangeDirectoryName()
            self.PrepareLinuxRelease()
            self.MakeRNATools()
            if self.IncludeRPMBuild:
                self.DoRPMBuild()
            if self.IncludeRsync:
                self.DoRsync()
        finally:
            if os.path.exists(self.OriginalFolder):
                os.chdir(self.OriginalFolder)
        print "...Work complete!"

    def PrintUsageInfo(self):
        print "Usage info:"
        print "  -r ReleaseFolder"
        print "  -d DirectoryName (if you want to rename the release folder for a particular branch)"
        print "  -t ReleaseType (Default 'isas' includes all workflows in a normal Isas release. Other options: all, rna, rna-basespace, forensics"
        print "  -p: Also kick off an RPM build"
        print "  -s: Also rsync to the other clusters"
        print "  -q Add ClusterQueueName option to Isas.exe.config to enable SGE (multi-node) mode."   
        print "  -f force overwrite of existing files (not on by default!)"
        print "  -o turn off building IsasRNA_Tools"
        print
        print "Example usage:"
        print "python UpdateRelease.py -r /illumina/development/Isas/1.2.3"
        print "python UpdateRelease.py -r \"/illumina/scratch/Isas/Builds/Development/Trunk/2.5.18/2.5.18 - x64.zip\""
        print "python UpdateRelease.py -r \"/illumina/scratch/Isas/Builds/Development/Trunk/2.5.18/2.5.18 - x64.zip\" -d \"2.5.18.Enrichment\""

    def ParseCommandLine(self, Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:d:t:psfqo")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                self.ReleaseFolder = os.path.abspath(Value) #trouble will arise for relative path due to os.chdir() in some functions.
            elif Option == "-d":
                self.DirectoryName = Value
            elif Option == "-t":
                self.ReleaseType = Value.lower()
            elif Option == "-p":
                self.IncludeRPMBuild = 1
            elif Option == "-s":
                self.IncludeRsync = 1
            elif Option == "-f":
                self.ForceOverwrite = True
            elif Option == "-q":
                self.Sge = True
            elif Option == "-o":
                self.BuildIsasRNA_Tools = False
            else:
                print "* Error: Unhandled option!", Option

if __name__ == "__main__":
    FolderPreparer = Preparer()
    FolderPreparer.Main(sys.argv[1:])
