#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time

#import ROOT
from ROOT import *
############################################
#            Job steering                  #
############################################
from optparse import OptionParser
parser = OptionParser()
parser.add_option('--keeptar', action='store_true', dest='keeptar', default=False, help='keep old tarball for condor jobs (default = %default)')
parser.add_option("--model", dest="model", default = "TChiHH",help="SMS model", metavar="model")
parser.add_option("--inDir", dest="inDir", default = "/store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/",help="EOS output directory  (default = %default)", metavar="outDir")
parser.add_option("--outDir", dest="outDir", default = "/store/user/emacdona/BoostedH/",help="EOS output directory  (default = %default)", metavar="outDir")
(options, args) = parser.parse_args()

# -----------------------------------------------------------------
#Create CACHEDIR.TAG files on the fly to exclude output directories from condor tarball
# -----------------------------------------------------------------
def cachedir(DIR):
    if len(DIR)==0:
        return

    tagfile = DIR+"/CACHEDIR.TAG"
    if not os.path.isfile(tagfile):
        tag = open(tagfile,'w')
        tag.write("Signature: 8a477f597d28d172789f06886806bc55\n")
        tag.write("# This file is a cache directory tag.\n")
        tag.write("# For information about cache directory tags, see:\n")
        tag.write("#       http://www.brynosaurus.com/cachedir/")
        tag.close()

def condorize(command,tag,odir,CMSSWVER):

    print "--------------------"
    print "Launching phase space point:",tag


    # setup environment
    #change to a tmp dir
    os.chdir("tmp");
    startdir = os.getcwd();
    f1n = "tmp_%s.sh" %(tag);
    f1=open(f1n, 'w')
    f1.write("#!/bin/sh \n");
    f1.write("tar -xzf %s.tar.gz \n" % (CMSSWVER));
    f1.write("source /cvmfs/cms.cern.ch/cmsset_default.sh \n");
    f1.write("set SCRAM_ARCH=slc6_amd64_gcc530\n")
    f1.write("cd %s/src \n" %(CMSSWVER));
    f1.write("eval `scramv1 runtime -sh`\n")
    f1.write("scramv1 b ProjectRename\n")
    f1.write("scramv1 b \n")
    f1.write("eval `scramv1 runtime -sh`\n")
    ####Go to Executable directory
    f1.write("cd boostedHiggsPlusMET/src \n");
    f1.write(command+" \n")
    #### COPY OUTPUT   higgsCombineTChiHH1000_LSP1_2BoostedH.AsymptoticLimits.mH120.root
    #ALPHABET output
    f1.write("xrdcp -f ALPHABETMC2016_V17_%s.root  root://cmseos.fnal.gov/%s/ALPHABET/ALPHABETMC2016_V17_%s.root\n" %(tag,odir,tag))

    #datacard output
    f1.write("xrdcp -f %s.txt  root://cmseos.fnal.gov/%s/datacards/%s.txt\n" %(tag,odir,tag))

    #limits output
    f1.write("xrdcp -f higgsCombine%s.AsymptoticLimits.mH120.root  root://cmseos.fnal.gov/%s/limits/higgsCombine%s.AsymptoticLimits.mH120.root\n" %(tag,odir,tag))

    f1.write("rm -r *.py *.root *.txt *.tar.gz \n")
    f1.close();
    f2n = "tmp_%s.condor" % (tag);
    outtag = "out_%s_$(Cluster)" % (tag)

    f2=open(f2n, 'w')
    f2.write("universe = vanilla \n");
    f2.write("Executable = %s \n" % (f1n) );
    f2.write("Requirements = OpSys == \"LINUX\"&& (Arch != \"DUMMY\" ) \n");
    #f2.write("request_disk = 100MB \n");
    #f2.write("request_memory = 2100 \n");
    f2.write("Should_Transfer_Files = YES \n");
    f2.write("WhenToTransferOutput  = ON_EXIT \n");
    f2.write("Transfer_Input_Files = %s, %s.tar.gz \n" % (f1n,CMSSWVER));
    f2.write("Output = "+outtag+".stdout \n");
    f2.write("Error = "+outtag+".stderr \n");
    f2.write("Log = "+outtag+".log \n");
    f2.write("Notification    = Error \n");
    f2.write("x509userproxy = $ENV(X509_USER_PROXY) \n")
    f2.write("Queue 1 \n");
    f2.close();

    os.system("condor_submit %s" % (f2n));
    os.chdir("../.");
    #else:
    #	os.system("Qsub -l lnxfarm -o OutPut_LogFile%s -e %s" %(tag,command))
    #	print "Qsub -l lnxfarm -o OutPut_LogFile%s -e %s" %(tag,command)



if __name__ == '__main__':
    outDir = options.outDir

    # get some info from the OS
    CMSSWVER = os.getenv("CMSSW_VERSION")
    CMSSWBASE = os.getenv("CMSSW_BASE")
    print("CMSSWBASE: "+CMSSWBASE)
    print("CMSSWVER: "+CMSSWVER)
    # tar it up for usage

    if not os.path.exists('tmp'):
        os.makedirs('tmp')
        cachedir('tmp')

    if not options.keeptar:
        print "Make tar "
        os.system("tar --exclude-caches-all --exclude inputHistograms/fastsimSignalT*  -zcf tmp/"+CMSSWVER+".tar.gz -C "+CMSSWBASE+"/.. "+CMSSWVER)

    #f = TFile.Open("inputHistograms/fastsimSignalT1tttt/RA2bin_signal.root");
    localpath=options.inDir
    eosDir = "root://cmseos.fnal.gov/"+localpath
    filenames = next(os.walk("/eos/uscms/"+localpath))[2]

    models = []
    mHino=[]
    mLSPs=[]

    for f in filenames:
	parse=f.split("_")
	print parse
	#if not "proc" in parse[1]:continue
	if not "MC2016" in parse[6]: continue
	#if options.model==parse[2] :
	if options.model==parse[1]:
		models.append(parse[1])
		mHino.append(int(parse[4]))
		mLSPs.append(int(parse[5]))
		print("Mgo %d mLSP %d " %(int(parse[4]),int(parse[5])))
    for m in range(len(mHino)):
        #    for mLSP in mLSPs:
        command="../bin/ALPHABET 1 MC2016 %d %d 0" %(mHino[m],mLSPs[m]) #last line determines whether or not we run with veto
        command +="\n"
        command += "python QuickDataCardsABCDNorm_Higgsino.py ";
        command += "%i " % mHino[m];
        command += "%i " % mLSPs[m];
        #command += " --eos %s" % (eosDir);
        # os.system(command)
        print command
        tag="%s%d_LSP%d_2BoostedH" %(options.model,mHino[m],mLSPs[m])
        condorize( command, tag, outDir, CMSSWVER );
        time.sleep(0.05);


    # else:

    #     os.system('python analysisBuilderCondor.py -b --signal T1bbbb --mGo 1500 --mLSP 100 --realData --tag allBkgs');
    #     os.system('python analysisBuilderCondor.py -b --signal T1bbbb --mGo 1000 --mLSP 100 --realData --tag allBkgs');
    #     os.system('python analysisBuilderCondor.py -b --signal T1tttt --mGo 1500 --mLSP 800 --realData --tag allBkgs');
    #     os.system('python analysisBuilderCondor.py -b --signal T1tttt --mGo 1200 --mLSP 800 --realData --tag allBkgs');
    #     os.system('python analysisBuilderCondor.py -b --signal T1qqqq --mGo 1400 --mLSP 800 --realData --tag allBkgs');
    #     os.system('python analysisBuilderCondor.py -b --signal T1qqqq --mGo 1000 --mLSP 800 --realData --tag allBkgs');
