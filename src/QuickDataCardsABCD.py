import os
from ROOT import *
import sys
#f=TFile("ALPHABEThistosMC2016_AllMETBins.root", "READ");
f=TFile("ALPHABEThistosMC2016_ResolvedVetoBoosted.root", "READ");
#f=TFile("ALPHABEThistos_V17bkg.root", "READ");

ZregionA=f.Get("MET_doubletagSR_ZJets"); ZregionB=f.Get("MET_doubletagSB_ZJets");
ZregionC=f.Get("MET_antitagSR_ZJets"); ZregionD=f.Get("MET_antitagSB_ZJets");

WregionA=f.Get("MET_doubletagSR_WJets"); WregionB=f.Get("MET_doubletagSB_WJets");
WregionC=f.Get("MET_antitagSR_WJets"); WregionD=f.Get("MET_antitagSB_WJets");

TregionA=f.Get("MET_doubletagSR_TT"); TregionB=f.Get("MET_doubletagSB_TT");
TregionC=f.Get("MET_antitagSR_TT"); TregionD=f.Get("MET_antitagSB_TT");

QregionA=f.Get("MET_doubletagSR_QCD"); QregionB=f.Get("MET_doubletagSB_QCD");
QregionC=f.Get("MET_antitagSR_QCD"); QregionD=f.Get("MET_antitagSB_QCD");

scale=137000.0/35862.824;
ZregionA.Scale(scale); ZregionB.Scale(scale);
ZregionC.Scale(scale); ZregionD.Scale(scale);

WregionA.Scale(scale); WregionB.Scale(scale);
WregionC.Scale(scale); WregionD.Scale(scale);

TregionA.Scale(scale); TregionB.Scale(scale);
TregionC.Scale(scale); TregionD.Scale(scale);

QregionA.Scale(scale); QregionB.Scale(scale);
QregionC.Scale(scale); QregionD.Scale(scale);


hino=int(sys.argv[1])
TotalEvents=float(sys.argv[2])
#fsig=TFile("ALPHABEThistosMC2016_SignalAllMETBins.root", "READ");
fsig=TFile("ALPHABEThistosMC2016TChiHHResolvedVeto.root", "READ");
SignalA=fsig.Get("MET_doubletagSR_TChiHH%d" %(hino));
SignalB=fsig.Get("MET_doubletagSB_TChiHH%d" %(hino));
SignalC=fsig.Get("MET_antitagSR_TChiHH%d" %(hino));
SignalD=fsig.Get("MET_antitagSB_TChiHH%d" %(hino));
scale=137000.0/35862.824;
SignalA.Scale((0.5823329*0.5823329*scale)/TotalEvents);
SignalB.Scale((0.5823329*0.5823329*scale)/TotalEvents);
SignalC.Scale((0.5823329*0.5823329*scale)/TotalEvents);
SignalD.Scale((0.5823329*0.5823329*scale)/TotalEvents);

YieldsFile=TFile("OutputYields.root","READ");
TotalBackground=YieldsFile.Get("SignalRegionYields")
Predictions=YieldsFile.Get("BackgroundPred")
ControlB=YieldsFile.Get("ControlB")
ControlC=YieldsFile.Get("ControlC")
ControlD=YieldsFile.Get("ControlD")
#kappas=[0.82484,0.841152,0.932512,1.41155,1.0]
kappas=[ 0.689671,0.779551,0.914834,1.413,1.0] #sim/pred, these are hardcoded
fcard=open("SignalRegionTemplate.txt", 'r');
fcardout=open("SignalRegionTChiHH%d.txt" %hino , 'w')
for line in fcard:
	if "observation" in line:
		newline="observation %g %g %g %g %g"	 %(Predictions.GetBinContent(1)*kappas[0],Predictions.GetBinContent(2)*kappas[1],Predictions.GetBinContent(3)*kappas[2],Predictions.GetBinContent(4)*kappas[3], Predictions.GetBinContent(5)*kappas[4])
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline+" %g %g %g %g %g "  %(SignalA.GetBinContent(1),QregionA.GetBinContent(1),ZregionA.GetBinContent(1),WregionA.GetBinContent(1),TregionA.GetBinContent(1))
		newline=newline +"%g %g %g %g %g "  %(SignalA.GetBinContent(2),QregionA.GetBinContent(2),ZregionA.GetBinContent(2),WregionA.GetBinContent(2),TregionA.GetBinContent(2))
		newline=newline +"%g %g %g %g %g "  %(SignalA.GetBinContent(3),QregionA.GetBinContent(3),ZregionA.GetBinContent(3),WregionA.GetBinContent(3),TregionA.GetBinContent(3))
		newline=newline +"%g %g %g %g %g "  %(SignalA.GetBinContent(4),QregionA.GetBinContent(4),ZregionA.GetBinContent(4),WregionA.GetBinContent(4),TregionA.GetBinContent(4))
		newline=newline +"%g %g %g %g %g \n"%(SignalA.GetBinContent(5),QregionA.GetBinContent(5),ZregionA.GetBinContent(5),WregionA.GetBinContent(5),TregionA.GetBinContent(5))
		fcardout.write(newline);
		print newline
	else:
		fcardout.write(line);
fcardout.close()
fcard.close()
fcardout=open("ControlRegionBTChiHH%d.txt" %hino , 'w')
fcard=open("RegionBTemplate.txt", 'r');
for line in fcard:
	if "observation" in line:
		print line
		newline="observation %g %g %g %g %g" %(ControlB.GetBinContent(1), ControlB.GetBinContent(2),ControlB.GetBinContent(3), ControlB.GetBinContent(4), ControlB.GetBinContent(5))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g %g %g %g "   %(SignalB.GetBinContent(1),QregionB.GetBinContent(1),ZregionB.GetBinContent(1),WregionB.GetBinContent(1),TregionB.GetBinContent(1))
		newline=newline +"%g %g %g %g %g "   %(SignalB.GetBinContent(2),QregionB.GetBinContent(2),ZregionB.GetBinContent(2),WregionB.GetBinContent(2),TregionB.GetBinContent(2))
		newline=newline +"%g %g %g %g %g "   %(SignalB.GetBinContent(3),QregionB.GetBinContent(3),ZregionB.GetBinContent(3),WregionB.GetBinContent(3),TregionB.GetBinContent(3))
		newline=newline +"%g %g %g %g %g "   %(SignalB.GetBinContent(4),QregionB.GetBinContent(4),ZregionB.GetBinContent(4),WregionB.GetBinContent(4),TregionB.GetBinContent(4))
		newline=newline +"%g %g %g %g %g \n" %(SignalB.GetBinContent(5),QregionB.GetBinContent(5),ZregionB.GetBinContent(5),WregionB.GetBinContent(5),TregionB.GetBinContent(5))
		fcardout.write(newline);
		print newline
	else:
		fcardout.write(line);
fcardout.close()


fcardout=open("ControlRegionCTChiHH%d.txt" %hino , 'w')
fcard=open("RegionCTemplate.txt", 'r');
for line in fcard:
	if "observation" in line:
		print line
		newline="observation %g %g %g %g %g" %(ControlC.GetBinContent(1), ControlC.GetBinContent(2),ControlC.GetBinContent(3), ControlC.GetBinContent(4), ControlC.GetBinContent(5))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g %g %g %g "   %(SignalC.GetBinContent(1),QregionC.GetBinContent(1),ZregionC.GetBinContent(1),WregionC.GetBinContent(1),TregionC.GetBinContent(1))
		newline=newline +"%g %g %g %g %g "   %(SignalC.GetBinContent(2),QregionC.GetBinContent(2),ZregionC.GetBinContent(2),WregionC.GetBinContent(2),TregionC.GetBinContent(2))
		newline=newline +"%g %g %g %g %g "   %(SignalC.GetBinContent(3),QregionC.GetBinContent(3),ZregionC.GetBinContent(3),WregionC.GetBinContent(3),TregionC.GetBinContent(3))
		newline=newline +"%g %g %g %g %g "   %(SignalC.GetBinContent(4),QregionC.GetBinContent(4),ZregionC.GetBinContent(4),WregionC.GetBinContent(4),TregionC.GetBinContent(4))
		newline=newline +"%g %g %g %g %g \n" %(SignalC.GetBinContent(5),QregionC.GetBinContent(5),ZregionC.GetBinContent(5),WregionC.GetBinContent(5),TregionC.GetBinContent(5))
		fcardout.write(newline);
		print newline
	else:
		fcardout.write(line);
fcardout.close()

fcardout=open("ControlRegionDTChiHH%d.txt" %hino , 'w')
fcard=open("RegionDTemplate.txt", 'r');
for line in fcard:
	if "observation" in line:
		print line
		newline="observation %g %g %g %g %g" %(ControlD.GetBinContent(1), ControlD.GetBinContent(2),ControlD.GetBinContent(3), ControlD.GetBinContent(4), ControlD.GetBinContent(5))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g %g %g %g "   %(SignalD.GetBinContent(1),QregionD.GetBinContent(1),ZregionD.GetBinContent(1),WregionD.GetBinContent(1),TregionD.GetBinContent(1))
		newline=newline +"%g %g %g %g %g "   %(SignalD.GetBinContent(2),QregionD.GetBinContent(2),ZregionD.GetBinContent(2),WregionD.GetBinContent(2),TregionD.GetBinContent(2))
		newline=newline +"%g %g %g %g %g "   %(SignalD.GetBinContent(3),QregionD.GetBinContent(3),ZregionD.GetBinContent(3),WregionD.GetBinContent(3),TregionD.GetBinContent(3))
		newline=newline +"%g %g %g %g %g "   %(SignalD.GetBinContent(4),QregionD.GetBinContent(4),ZregionD.GetBinContent(4),WregionD.GetBinContent(4),TregionD.GetBinContent(4))
		newline=newline +"%g %g %g %g %g \n" %(SignalD.GetBinContent(5),QregionD.GetBinContent(5),ZregionD.GetBinContent(5),WregionD.GetBinContent(5),TregionD.GetBinContent(5))
		fcardout.write(newline);
		print newline
	else:
		fcardout.write(line);
fcardout.close()

##########################RESOLVED CARDS
fres=TFile("ALPHABEThistosMC2016_ResolvedSignalRegion.root", "READ");

ZresregionA  =fres.Get("MET_fourbSR_ZJets");  ZresregionB  =fres.Get("MET_fourbSB_ZJets");
Zres3bregionA=fres.Get("MET_threebSR_ZJets"); Zres3bregionB=fres.Get("MET_threebSB_ZJets");
ZresregionC  =fres.Get("MET_twobSR_ZJets");   ZresregionD  =fres.Get("MET_twobSB_ZJets");

WresregionA  =fres.Get("MET_fourbSR_WJets");  WresregionB  =fres.Get("MET_fourbSB_WJets");
Wres3bregionA=fres.Get("MET_threebSR_WJets"); Wres3bregionB=fres.Get("MET_threebSB_WJets");
WresregionC  =fres.Get("MET_twobSR_WJets");   WresregionD  =fres.Get("MET_twobSB_WJets");

TresregionA  =fres.Get("MET_fourbSR_TT");  TresregionB  =fres.Get("MET_fourbSB_TT");
Tres3bregionA=fres.Get("MET_threebSR_TT"); Tres3bregionB=fres.Get("MET_threebSB_TT");
TresregionC  =fres.Get("MET_twobSR_TT");   TresregionD  =fres.Get("MET_twobSB_TT");

QresregionA  =fres.Get("MET_fourbSR_QCD");  QresregionB  =fres.Get("MET_fourbSB_QCD");
Qres3bregionA=fres.Get("MET_threebSR_QCD"); Qres3bregionB=fres.Get("MET_threebSB_QCD");
QresregionC  =fres.Get("MET_twobSR_QCD");   QresregionD  =fres.Get("MET_twobSB_QCD");

ZresregionA.Scale(scale); ZresregionB.Scale(scale);
Zres3bregionA.Scale(scale); Zres3bregionB.Scale(scale);
ZresregionC.Scale(scale); ZresregionD.Scale(scale);
WresregionA.Scale(scale); WresregionB.Scale(scale);
Wres3bregionA.Scale(scale); Wres3bregionB.Scale(scale);
WresregionC.Scale(scale); WresregionD.Scale(scale);
TresregionA.Scale(scale); TresregionB.Scale(scale);
Tres3bregionA.Scale(scale); Tres3bregionB.Scale(scale);
TresregionC.Scale(scale); TresregionD.Scale(scale);

QresregionA.Scale(scale); QresregionB.Scale(scale);
Qres3bregionA.Scale(scale); Qres3bregionB.Scale(scale);
QresregionC.Scale(scale); QresregionD.Scale(scale);

YieldsFile4bResolved=TFile("OutputYieldsResolved4b.root","READ");
TotalBackgroundResolved4b=YieldsFile4bResolved.Get("SignalRegionYields")
PredictionsResolved4b=YieldsFile4bResolved.Get("BackgroundPred")
ControlBResolved4b=YieldsFile4bResolved.Get("ControlB")
ControlCResolved4b=YieldsFile4bResolved.Get("ControlC")
ControlDResolved4b=YieldsFile4bResolved.Get("ControlD")

fsigres=TFile("ALPHABEThistosMC2016TChiHHResolvedSearch.root", "READ");
SignalResA=fsigres.Get("MET_fourbSR_TChiHH%d" %(hino));
SignalResB=fsigres.Get("MET_fourbSB_TChiHH%d" %(hino));
SignalResA3b=fsigres.Get("MET_threebSR_TChiHH%d" %(hino));
SignalResB3b=fsigres.Get("MET_threebSB_TChiHH%d" %(hino));
SignalResC=fsigres.Get("MET_twobSR_TChiHH%d" %(hino));
SignalResD=fsigres.Get("MET_twobSB_TChiHH%d" %(hino));
SignalResA.Scale((0.5823329*0.5823329*scale)/TotalEvents);
SignalResB.Scale((0.5823329*0.5823329*scale)/TotalEvents);
SignalResA3b.Scale((0.5823329*0.5823329*scale)/TotalEvents);
SignalResB3b.Scale((0.5823329*0.5823329*scale)/TotalEvents);
SignalResC.Scale((0.5823329*0.5823329*scale)/TotalEvents);
SignalResD.Scale((0.5823329*0.5823329*scale)/TotalEvents);


fcard=open("SignalRegionResolvedTemplate.txt", 'r');
fcardout=open("SignalRegionResolved4bTChiHH%d.txt" %hino , 'w')
kappas4b=[0.817174,1.13345,0.709114,1.0] #these are hardcoded
print PredictionsResolved4b.GetBinContent(1)*kappas4b[0],QresregionA.GetBinContent(1)+ZresregionA.GetBinContent(1)+WresregionA.GetBinContent(1)+TresregionA.GetBinContent(1)
for line in fcard:
	if "observation" in line:
		print line
		newline="observation %g %g %g %g " %(PredictionsResolved4b.GetBinContent(1)*kappas4b[0],PredictionsResolved4b.GetBinContent(2)*kappas4b[1],PredictionsResolved4b.GetBinContent(3)*kappas4b[2],PredictionsResolved4b.GetBinContent(4)*kappas4b[3])
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g %g %g %g "  %(SignalResA.GetBinContent(1),QresregionA.GetBinContent(1),ZresregionA.GetBinContent(1),WresregionA.GetBinContent(1),TresregionA.GetBinContent(1))
		newline=newline +"%g %g %g %g %g "  %(SignalResA.GetBinContent(2),QresregionA.GetBinContent(2),ZresregionA.GetBinContent(2),WresregionA.GetBinContent(2),TresregionA.GetBinContent(2))
		newline=newline +"%g %g %g %g %g "  %(SignalResA.GetBinContent(3),QresregionA.GetBinContent(3),ZresregionA.GetBinContent(3),WresregionA.GetBinContent(3),TresregionA.GetBinContent(3))
		newline=newline +"%g %g %g %g %g \n"%(SignalResA.GetBinContent(4),QresregionA.GetBinContent(4),ZresregionA.GetBinContent(4),WresregionA.GetBinContent(4),TresregionA.GetBinContent(4))
		fcardout.write(newline);
	else:
		fcardout.write(line);
fcardout.close()
fcard.close()

YieldsFile3bResolved=TFile("OutputYieldsResolved3b.root","READ");
TotalBackgroundResolved3b=YieldsFile3bResolved.Get("SignalRegionYields")
PredictionsResolved3b=YieldsFile3bResolved.Get("BackgroundPred")
ControlBResolved3b=YieldsFile3bResolved.Get("ControlB")

fcard=open("SignalRegionResolvedTemplate3Btags.txt", 'r');
fcardout=open("SignalRegionResolved3bTChiHH%d.txt" %hino , 'w')
kappas3b=[0.932591,0.99398,0.824331,1.33238] #these are hardcoded

for line in fcard:
	if "observation" in line:
		print line
		newline="observation %g %g %g %g " %(PredictionsResolved3b.GetBinContent(1)*kappas3b[0],PredictionsResolved3b.GetBinContent(2)*kappas3b[1],PredictionsResolved3b.GetBinContent(3)*kappas3b[2],PredictionsResolved3b.GetBinContent(4)*kappas3b[3])
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g %g %g %g "  %(SignalResA3b.GetBinContent(1),Qres3bregionA.GetBinContent(1),Zres3bregionA.GetBinContent(1),Wres3bregionA.GetBinContent(1),Tres3bregionA.GetBinContent(1))
		newline=newline +"%g %g %g %g %g "  %(SignalResA3b.GetBinContent(2),Qres3bregionA.GetBinContent(2),Zres3bregionA.GetBinContent(2),Wres3bregionA.GetBinContent(2),Tres3bregionA.GetBinContent(2))
		newline=newline +"%g %g %g %g %g "  %(SignalResA3b.GetBinContent(3),Qres3bregionA.GetBinContent(3),Zres3bregionA.GetBinContent(3),Wres3bregionA.GetBinContent(3),Tres3bregionA.GetBinContent(3))
		newline=newline +"%g %g %g %g %g \n"%(SignalResA3b.GetBinContent(4),Qres3bregionA.GetBinContent(4),Zres3bregionA.GetBinContent(4),Wres3bregionA.GetBinContent(4),Tres3bregionA.GetBinContent(4))
		fcardout.write(newline);
		print newline
	else:
		fcardout.write(line);
fcardout.close()
fcard.close()

fcardout=open("ControlRegionResBTChiHH%d.txt" %hino , 'w')
fcard=open("RegionBTemplateRes.txt", 'r');
for line in fcard:
	if "observation" in line:
		print line
		newline="observation %g %g %g %g " %(ControlBResolved4b.GetBinContent(1), ControlBResolved4b.GetBinContent(2),ControlBResolved4b.GetBinContent(3),ControlBResolved4b.GetBinContent(4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g %g %g %g "  %(SignalResB.GetBinContent(1),QresregionB.GetBinContent(1),ZresregionB.GetBinContent(1),WresregionB.GetBinContent(1),TresregionB.GetBinContent(1))
		newline=newline +"%g %g %g %g %g "  %(SignalResB.GetBinContent(2),QresregionB.GetBinContent(2),ZresregionB.GetBinContent(2),WresregionB.GetBinContent(2),TresregionB.GetBinContent(2))
		newline=newline +"%g %g %g %g %g "  %(SignalResB.GetBinContent(3),QresregionB.GetBinContent(3),ZresregionB.GetBinContent(3),WresregionB.GetBinContent(3),TresregionB.GetBinContent(3))
		newline=newline +"%g %g %g %g %g \n"%(SignalResB.GetBinContent(4),QresregionB.GetBinContent(4),ZresregionB.GetBinContent(4),WresregionB.GetBinContent(4),TresregionB.GetBinContent(4))
		fcardout.write(newline);
		print newline
	else:
		fcardout.write(line);
fcardout.close()

fcardout=open("ControlRegionResCTChiHH%d.txt" %hino , 'w')
fcard=open("RegionCTemplateRes.txt", 'r');
for line in fcard:
	if "observation" in line:
		print line
		newline="observation %g %g %g %g " %(ControlCResolved4b.GetBinContent(1), ControlCResolved4b.GetBinContent(2),ControlCResolved4b.GetBinContent(3), ControlCResolved4b.GetBinContent(4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g %g %g %g "  %(SignalResC.GetBinContent(1),QresregionC.GetBinContent(1),ZresregionC.GetBinContent(1),WresregionC.GetBinContent(1),TresregionC.GetBinContent(1))
		newline=newline +"%g %g %g %g %g "  %(SignalResC.GetBinContent(2),QresregionC.GetBinContent(2),ZresregionC.GetBinContent(2),WresregionC.GetBinContent(2),TresregionC.GetBinContent(2))
		newline=newline +"%g %g %g %g %g "  %(SignalResC.GetBinContent(3),QresregionC.GetBinContent(3),ZresregionC.GetBinContent(3),WresregionC.GetBinContent(3),TresregionC.GetBinContent(3))
		newline=newline +"%g %g %g %g %g \n"%(SignalResC.GetBinContent(4),QresregionC.GetBinContent(4),ZresregionC.GetBinContent(4),WresregionC.GetBinContent(4),TresregionC.GetBinContent(4))
		fcardout.write(newline);
		print newline
	else:
		fcardout.write(line);
fcardout.close()
fcardout=open("ControlRegionResB3bTChiHH%d.txt" %hino , 'w')
fcard=open("RegionBTemplateRes3b.txt", 'r');
for line in fcard:
	if "observation" in line:
		print line
		newline="observation %g %g %g %g " %(ControlBResolved3b.GetBinContent(1), ControlBResolved3b.GetBinContent(2),ControlBResolved3b.GetBinContent(3), ControlBResolved3b.GetBinContent(4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g %g %g %g "  %(SignalResB3b.GetBinContent(1),Qres3bregionB.GetBinContent(1),Zres3bregionB.GetBinContent(1),Wres3bregionB.GetBinContent(1),Tres3bregionB.GetBinContent(1))
		newline=newline +"%g %g %g %g %g "  %(SignalResB3b.GetBinContent(2),Qres3bregionB.GetBinContent(2),Zres3bregionB.GetBinContent(2),Wres3bregionB.GetBinContent(2),Tres3bregionB.GetBinContent(2))
		newline=newline +"%g %g %g %g %g "  %(SignalResB3b.GetBinContent(3),Qres3bregionB.GetBinContent(3),Zres3bregionB.GetBinContent(3),Wres3bregionB.GetBinContent(3),Tres3bregionB.GetBinContent(3))
		newline=newline +"%g %g %g %g %g \n"%(SignalResB3b.GetBinContent(4),Qres3bregionB.GetBinContent(4),Zres3bregionB.GetBinContent(4),Wres3bregionB.GetBinContent(4),Tres3bregionB.GetBinContent(4))
		fcardout.write(newline);
		print newline
	else:
		fcardout.write(line);

fcardout.close()

fcardout=open("ControlRegionResDTChiHH%d.txt" %hino , 'w')
fcard=open("RegionDTemplateRes.txt", 'r');
for line in fcard:
	if "observation" in line:
		print line
		newline="observation %g %g %g %g " %(ControlDResolved4b.GetBinContent(1), ControlDResolved4b.GetBinContent(2),ControlDResolved4b.GetBinContent(3), ControlDResolved4b.GetBinContent(4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g %g %g %g "  %(SignalResD.GetBinContent(1),QresregionD.GetBinContent(1),ZresregionD.GetBinContent(1),WresregionD.GetBinContent(1),TresregionD.GetBinContent(1))
		newline=newline +"%g %g %g %g %g "  %(SignalResD.GetBinContent(2),QresregionD.GetBinContent(2),ZresregionD.GetBinContent(2),WresregionD.GetBinContent(2),TresregionD.GetBinContent(2))
		newline=newline +"%g %g %g %g %g "  %(SignalResD.GetBinContent(3),QresregionD.GetBinContent(3),ZresregionD.GetBinContent(3),WresregionD.GetBinContent(3),TresregionD.GetBinContent(3))
		newline=newline +"%g %g %g %g %g \n"%(SignalResD.GetBinContent(4),QresregionD.GetBinContent(4),ZresregionD.GetBinContent(4),WresregionD.GetBinContent(4),TresregionD.GetBinContent(4))
		fcardout.write(newline);
		print newline
	else:
		fcardout.write(line);
fcardout.close()
os.system("combineCards.py SignalRegionTChiHH%d.txt ControlRegionBTChiHH%d.txt ControlRegionCTChiHH%d.txt  ControlRegionDTChiHH%d.txt>TChiHH%d2BoostedH.txt "%(hino,hino,hino,hino,hino))
#os.system("combine -M Asymptotic -n TChiHH%d TChiHH%d2BoostedH.txt " %(hino,hino))
#os.system("combineCards.py SignalRegionTChiHH%d.txt ControlRegionBTChiHH%d.txt ControlRegionCTChiHH%d.txt  ControlRegionDTChiHH%d.txt>TChiHH%d.txt "%(hino,hino,hino,hino,hino))
os.system("combineCards.py SignalRegionResolved4bTChiHH%d.txt ControlRegionResBTChiHH%d.txt ControlRegionResCTChiHH%d.txt  ControlRegionResDTChiHH%d.txt>TChiHH%dResolved4b.txt "%(hino,hino,hino,hino,hino))
os.system("combineCards.py SignalRegionResolved3bTChiHH%d.txt ControlRegionResB3bTChiHH%d.txt >TChiHH%dResolved3b.txt "%(hino,hino,hino))
os.system("combineCards.py TChiHH%d2BoostedH.txt TChiHH%dResolved4b.txt  TChiHH%dResolved3b.txt > TChiHH%d.txt "%(hino,hino,hino,hino))

#os.system("combine -M Asymptotic -n TChiHH%d TChiHH%d.txt " %(hino,hino))
