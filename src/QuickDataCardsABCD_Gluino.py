import os
from ROOT import *
import sys

def binContent(histo,binNum):
	thisBinContent = histo.GetBinContent(binNum);
	if (thisBinContent>0):
	 	return thisBinContent;
  	else:
		 return 0.000001;

f=TFile("ALPHABET_V18_2BoostedH.root", "READ"); #BoostedOnly
fVeto=TFile("ALPHABETBoost_V17_resVeto_all.root", "READ"); #Boosted with resolved veto

sumRegionA2=f.Get("MET_doubletagSR_sum"); sumRegionB2=f.Get("MET_doubletagSB_sum");
sumRegionA1=f.Get("MET_tagSR_sum"); sumRegionB1=f.Get("MET_tagSB_sum");
sumRegionC=f.Get("MET_antitagSR_sum"); sumRegionD=f.Get("MET_antitagSB_sum");

scale=1.0;
sumRegionA2.Scale(scale); sumRegionB2.Scale(scale);
sumRegionA1.Scale(scale); sumRegionB1.Scale(scale);
sumRegionC.Scale(scale); sumRegionD.Scale(scale);

hino=int(sys.argv[1])
print("\n\n\n-----------------------------------------------------")
print("Looping over mass=%d now" %(hino))
print("-----------------------------------------------------\n")
f2=TFile("ALPHABET_V18_1Dsignal_2BoostedH.root", "READ")

SignalA2=f2.Get("MET_doubletagSR_T5qqqqZH_%d_1" %(hino));
SignalB2=f2.Get("MET_doubletagSB_T5qqqqZH_%d_1" %(hino));
SignalA1=f2.Get("MET_tagSR_T5qqqqZH_%d_1" %(hino));
SignalB1=f2.Get("MET_tagSB_T5qqqqZH_%d_1" %(hino));
SignalC=f2.Get("MET_antitagSR_T5qqqqZH_%d_1" %(hino));
SignalD=f2.Get("MET_antitagSB_T5qqqqZH_%d_1" %(hino));

scale=1.0;
SignalA2.Scale(scale)
SignalB2.Scale(scale)
SignalA1.Scale(scale)
SignalB1.Scale(scale)
SignalC.Scale(scale)
SignalD.Scale(scale)

Predictions2H=sumRegionB2.Clone("Predictions2H");
Predictions2H.Multiply(sumRegionC)
Predictions2H.Divide(sumRegionD)

Predictions1H=sumRegionB1.Clone("Predictions1H");
Predictions1H.Multiply(sumRegionC)
Predictions1H.Divide(sumRegionD)

kappas2H=[1.17,0.46,0.84] #sim/pred, will need to measure these
kappas1H=[1.04,0.98,0.86] #sim/pred, will need to measure these


fcard=open("SignalRegionTemplate.txt", 'r');
fcardout=open("SRPure_T5HH%d.txt" %hino , 'w')
for line in fcard:
	if "observation" in line:
		newline="observation %g %g %g" %(Predictions2H.GetBinContent(1)*kappas2H[0],Predictions2H.GetBinContent(2)*kappas2H[1],Predictions2H.GetBinContent(3)*kappas2H[2])
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline+" %g %g "  %(SignalA2.GetBinContent(1),sumRegionA2.GetBinContent(1))
		newline=newline +"%g %g "  %(SignalA2.GetBinContent(2),sumRegionA2.GetBinContent(2))
		newline=newline +"%g %g \n"%(SignalA2.GetBinContent(3)+SignalA2.GetBinContent(4),sumRegionA2.GetBinContent(3)+sumRegionA2.GetBinContent(4))
		fcardout.write(newline);
	else:
		fcardout.write(line);
fcardout.close()
fcard.close()

fcard=open("SignalRegion1Template.txt", 'r');
fcardout=open("SR1Pure_T5HH%d.txt" %(hino) , 'w')
for line in fcard:
	if "observation" in line:
		newline="observation %g %g %g" %(Predictions1H.GetBinContent(1)*kappas1H[0],Predictions1H.GetBinContent(2)*kappas1H[1],Predictions1H.GetBinContent(3)*kappas1H[2])
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline+" %g %g "  %(SignalA1.GetBinContent(1),sumRegionA1.GetBinContent(1))
		newline=newline +"%g %g "  %(SignalA1.GetBinContent(2),sumRegionA1.GetBinContent(2))
		newline=newline +"%g %g \n"%(SignalA1.GetBinContent(3)+SignalA1.GetBinContent(4),sumRegionA1.GetBinContent(3)+sumRegionA1.GetBinContent(4))
		fcardout.write(newline);
	else:
		fcardout.write(line);
fcardout.close()
fcard.close()

fcardout=open("CRBPure_T5HH%d.txt" %hino , 'w')
fcard=open("RegionBTemplate.txt", 'r');
for line in fcard:
	if "observation" in line:
		newline="observation %g %g %g" %(sumRegionB2.GetBinContent(1), sumRegionB2.GetBinContent(2),sumRegionB2.GetBinContent(3))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g "   %(SignalB2.GetBinContent(1),sumRegionB2.GetBinContent(1))
		newline=newline +"%g %g "   %(SignalB2.GetBinContent(2),sumRegionB2.GetBinContent(2))
		newline=newline +"%g %g \n" %(SignalB2.GetBinContent(3)+SignalB2.GetBinContent(4),sumRegionB2.GetBinContent(3)+sumRegionB2.GetBinContent(4))
		fcardout.write(newline);
	else:
		fcardout.write(line);
fcardout.close()

fcardout=open("CRB1Pure_T5HH%d.txt" %(hino) , 'w')
fcard=open("RegionB1Template.txt", 'r');
for line in fcard:
	if "observation" in line:
		newline="observation %g %g %g" %(sumRegionB1.GetBinContent(1), sumRegionB1.GetBinContent(2),sumRegionB1.GetBinContent(3))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g "   %(SignalB1.GetBinContent(1),sumRegionB1.GetBinContent(1))
		newline=newline +"%g %g "   %(SignalB1.GetBinContent(2),sumRegionB1.GetBinContent(2))
		newline=newline +"%g %g \n" %(SignalB1.GetBinContent(3)+SignalB1.GetBinContent(4),sumRegionB1.GetBinContent(3)+sumRegionB1.GetBinContent(4))
		fcardout.write(newline);
	else:
		fcardout.write(line);
fcardout.close()

fcardout=open("CRCPure_T5HH%d.txt" %hino , 'w')
fcard=open("RegionCTemplate.txt", 'r');
for line in fcard:
	if "observation" in line:
		newline="observation %g %g %g" %(sumRegionC.GetBinContent(1), sumRegionC.GetBinContent(2),sumRegionC.GetBinContent(3))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g "   %(SignalC.GetBinContent(1),sumRegionC.GetBinContent(1))
		newline=newline +"%g %g "   %(SignalC.GetBinContent(2),sumRegionC.GetBinContent(2))
		newline=newline +"%g %g \n" %(SignalC.GetBinContent(3)+SignalC.GetBinContent(4),sumRegionC.GetBinContent(3)+sumRegionC.GetBinContent(4))
		fcardout.write(newline);
	else:
		fcardout.write(line);
fcardout.close()


fcardout=open("CRDPure_T5HH%d.txt" %hino , 'w')
fcard=open("RegionDTemplate.txt", 'r');
for line in fcard:
	if "observation" in line:
		newline="observation %g %g %g" %(sumRegionD.GetBinContent(1), sumRegionD.GetBinContent(2),sumRegionD.GetBinContent(3))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g "   %(SignalD.GetBinContent(1),sumRegionD.GetBinContent(1))
		newline=newline +"%g %g "   %(SignalD.GetBinContent(2),sumRegionD.GetBinContent(2))
		newline=newline +"%g %g \n" %(SignalD.GetBinContent(3)+SignalD.GetBinContent(4),sumRegionD.GetBinContent(3)+sumRegionD.GetBinContent(4))
		fcardout.write(newline);
	else:
		fcardout.write(line);
fcardout.close()

'''
########################################################################
######################### Resolved Cards ###############################
########################################################################
fres=TFile("ALPHABETResMC2016_V18.root", "READ");

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

YieldsFileRes=TFile("OutputYields_res.root","READ");
TotalBackgroundRes4b=YieldsFileRes.Get("SRYields_4b")
PredictionsRes4b=YieldsFileRes.Get("BkgPred_4b")
ControlBRes4b=YieldsFileRes.Get("ControlB")
ControlCRes4b=YieldsFileRes.Get("ControlC")
ControlDRes4b=YieldsFileRes.Get("ControlD")

SignalResA=fres.Get("MET_fourbSR_TChiHH%d" %(hino));
SignalResB=fres.Get("MET_fourbSB_TChiHH%d" %(hino));
SignalResA3b=fres.Get("MET_threebSR_TChiHH%d" %(hino));
SignalResB3b=fres.Get("MET_threebSB_TChiHH%d" %(hino));
SignalResC=fres.Get("MET_twobSR_TChiHH%d" %(hino));
SignalResD=fres.Get("MET_twobSB_TChiHH%d" %(hino));

SignalResA.Scale(scale); SignalResA3b.Scale(scale);
SignalResB.Scale(scale); SignalResB3b.Scale(scale);
SignalResC.Scale(scale); SignalResD.Scale(scale);

fcard=open("SignalRegionResTemplate.txt", 'r');
fcardout=open("SignalRegionRes4bTChiHH%d.txt" %hino , 'w')
kappas4b=[0.817174,1.13345,0.709114,1.0] #these are hardcoded
print PredictionsRes4b.GetBinContent(1)*kappas4b[0],binContent(QresregionA,1)+binContent(ZresregionA,1)+binContent(WresregionA,1)+binContent(TresregionA,1)
for line in fcard:
	if "observation" in line:
		print line
		newline="observation %g %g %g %g " %(PredictionsRes4b.GetBinContent(1)*kappas4b[0],PredictionsRes4b.GetBinContent(2)*kappas4b[1],PredictionsRes4b.GetBinContent(3)*kappas4b[2],PredictionsRes4b.GetBinContent(4)*kappas4b[3])
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g %g %g %g "  %(binContent(SignalResA,1),binContent(QresregionA,1),binContent(ZresregionA,1),binContent(WresregionA,1),binContent(TresregionA,1))
		newline=newline +"%g %g %g %g %g "  %(binContent(SignalResA,2),binContent(QresregionA,2),binContent(ZresregionA,2),binContent(WresregionA,2),binContent(TresregionA,2))
		newline=newline +"%g %g %g %g %g "  %(binContent(SignalResA,3),binContent(QresregionA,3),binContent(ZresregionA,3),binContent(WresregionA,3),binContent(TresregionA,3))
		newline=newline +"%g %g %g %g %g \n"%(binContent(SignalResA,4)+binContent(SignalResA,5),binContent(QresregionA,4)+binContent(QresregionA,5),binContent(ZresregionA,4)+binContent(ZresregionA,5),binContent(WresregionA,4)+binContent(WresregionA,5),binContent(TresregionA,4)+binContent(TresregionA,5))
		fcardout.write(newline);
	else:
		fcardout.write(line);
fcardout.close()
fcard.close()

TotalBackgroundRes3b=YieldsFileRes.Get("SRYields_3b")
PredictionsRes3b=YieldsFileRes.Get("BkgPred_3b")
ControlBRes3b=YieldsFileRes.Get("ControlB1")

fcard=open("SignalRegionResTemplate.txt", 'r');
fcardout=open("SignalRegionRes3bTChiHH%d.txt" %hino , 'w')
kappas3b=[0.932591,0.99398,0.824331,1.33238] #these are hardcoded

for line in fcard:
	if "observation" in line:
		print line
		newline="observation %g %g %g %g " %(PredictionsRes3b.GetBinContent(1)*kappas3b[0],PredictionsRes3b.GetBinContent(2)*kappas3b[1],PredictionsRes3b.GetBinContent(3)*kappas3b[2],PredictionsRes3b.GetBinContent(4)*kappas3b[3])
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g %g %g %g "  %(binContent(SignalResA3b,1),binContent(Qres3bregionA,1),binContent(Zres3bregionA,1),binContent(Wres3bregionA,1),binContent(Tres3bregionA,1))
		newline=newline +"%g %g %g %g %g "  %(binContent(SignalResA3b,2),binContent(Qres3bregionA,2),binContent(Zres3bregionA,2),binContent(Wres3bregionA,2),binContent(Tres3bregionA,2))
		newline=newline +"%g %g %g %g %g "  %(binContent(SignalResA3b,3),binContent(Qres3bregionA,3),binContent(Zres3bregionA,3),binContent(Wres3bregionA,3),binContent(Tres3bregionA,3))
		newline=newline +"%g %g %g %g %g \n"%(binContent(SignalResA3b,4)+binContent(SignalResA3b,5),binContent(Qres3bregionA,4)+binContent(Qres3bregionA,5),binContent(Zres3bregionA,4)+binContent(Zres3bregionA,5),binContent(Wres3bregionA,4)+binContent(Wres3bregionA,5),binContent(Tres3bregionA,4)+binContent(Tres3bregionA,5))
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
		newline="observation %g %g %g %g " %(ControlBRes4b.GetBinContent(1), ControlBRes4b.GetBinContent(2),ControlBRes4b.GetBinContent(3),ControlBRes4b.GetBinContent(4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g %g %g %g "  %(binContent(SignalResB,1),binContent(QresregionB,1),binContent(ZresregionB,1),binContent(WresregionB,1),binContent(TresregionB,1))
		newline=newline +"%g %g %g %g %g "  %(binContent(SignalResB,2),binContent(QresregionB,2),binContent(ZresregionB,2),binContent(WresregionB,2),binContent(TresregionB,2))
		newline=newline +"%g %g %g %g %g "  %(binContent(SignalResB,3),binContent(QresregionB,3),binContent(ZresregionB,3),binContent(WresregionB,3),binContent(TresregionB,3))
		newline=newline +"%g %g %g %g %g \n"%(binContent(SignalResB,4)+binContent(SignalResB,5),binContent(QresregionB,4)+binContent(QresregionB,5),binContent(ZresregionB,4)+binContent(ZresregionB,5),binContent(WresregionB,4)+binContent(WresregionB,5),binContent(TresregionB,4)+binContent(TresregionB,5))
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
		newline="observation %g %g %g %g " %(ControlCRes4b.GetBinContent(1), ControlCRes4b.GetBinContent(2),ControlCRes4b.GetBinContent(3), ControlCRes4b.GetBinContent(4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g %g %g %g "  %(binContent(SignalResC,1),binContent(QresregionC,1),binContent(ZresregionC,1),binContent(WresregionC,1),binContent(TresregionC,1))
		newline=newline +"%g %g %g %g %g "  %(binContent(SignalResC,2),binContent(QresregionC,2),binContent(ZresregionC,2),binContent(WresregionC,2),binContent(TresregionC,2))
		newline=newline +"%g %g %g %g %g "  %(binContent(SignalResC,3),binContent(QresregionC,3),binContent(ZresregionC,3),binContent(WresregionC,3),binContent(TresregionC,3))
		newline=newline +"%g %g %g %g %g \n"%(binContent(SignalResC,4)+binContent(SignalResC,5),binContent(QresregionC,4)+binContent(QresregionC,5),binContent(ZresregionC,4)+binContent(ZresregionC,5),binContent(WresregionC,4)+binContent(WresregionC,5),binContent(TresregionC,4)+binContent(TresregionC,5))
		fcardout.write(newline);
		print newline
	else:
		fcardout.write(line);
fcardout.close()

fcardout=open("ControlRegionResB3bTChiHH%d.txt" %hino , 'w')
fcard=open("RegionBTemplateRes.txt", 'r');
for line in fcard:
	if "observation" in line:
		print line
		newline="observation %g %g %g %g " %(ControlBRes3b.GetBinContent(1), ControlBRes3b.GetBinContent(2),ControlBRes3b.GetBinContent(3), ControlBRes3b.GetBinContent(4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g %g %g %g "  %(binContent(SignalResB3b,1),binContent(Qres3bregionB,1),binContent(Zres3bregionB,1),binContent(Wres3bregionB,1),binContent(Tres3bregionB,1))
		newline=newline +"%g %g %g %g %g "  %(binContent(SignalResB3b,2),binContent(Qres3bregionB,2),binContent(Zres3bregionB,2),binContent(Wres3bregionB,2),binContent(Tres3bregionB,2))
		newline=newline +"%g %g %g %g %g "  %(binContent(SignalResB3b,3),binContent(Qres3bregionB,3),binContent(Zres3bregionB,3),binContent(Wres3bregionB,3),binContent(Tres3bregionB,3))
		newline=newline +"%g %g %g %g %g \n"%(binContent(SignalResB3b,4)+binContent(SignalResB3b,5),binContent(Qres3bregionB,4)+binContent(Qres3bregionB,5),binContent(Zres3bregionB,4)+binContent(Zres3bregionB,5),binContent(Wres3bregionB,4)+binContent(Wres3bregionB,5),binContent(Tres3bregionB,4)+binContent(Tres3bregionB,5))
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
		newline="observation %g %g %g %g " %(ControlDRes4b.GetBinContent(1), ControlDRes4b.GetBinContent(2),ControlDRes4b.GetBinContent(3), ControlDRes4b.GetBinContent(4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g %g %g %g "  %(binContent(SignalResD,1),binContent(QresregionD,1),binContent(ZresregionD,1),binContent(WresregionD,1),binContent(TresregionD,1))
		newline=newline +"%g %g %g %g %g "  %(binContent(SignalResD,2),binContent(QresregionD,2),binContent(ZresregionD,2),binContent(WresregionD,2),binContent(TresregionD,2))
		newline=newline +"%g %g %g %g %g "  %(binContent(SignalResD,3),binContent(QresregionD,3),binContent(ZresregionD,3),binContent(WresregionD,3),binContent(TresregionD,3))
		newline=newline +"%g %g %g %g %g \n"%(binContent(SignalResD,4)+binContent(SignalResD,5),binContent(QresregionD,4)+binContent(QresregionD,5),binContent(ZresregionD,4)+binContent(ZresregionD,5),binContent(WresregionD,4)+binContent(WresregionD,5),binContent(TresregionD,4)+binContent(TresregionD,5))
		fcardout.write(newline);
		print newline
	else:
		fcardout.write(line);
fcardout.close()
'''
#BoostedOnly
os.system("combineCards.py SRPure_T5HH%d.txt CRBPure_T5HH%d.txt CRCPure_T5HH%d.txt  CRDPure_T5HH%d.txt > T5HH%d_2BoostedHPure.txt "%(hino,hino,hino,hino,hino))
os.system("combineCards.py SR1Pure_T5HH%d.txt CRB1Pure_T5HH%d.txt CRCPure_T5HH%d.txt  CRDPure_T5HH%d.txt > T5HH%d_1BoostedHPure.txt "%(hino,hino,hino,hino,hino))
os.system("combineCards.py T5HH%d_1BoostedH.txt T5HH%d_2BoostedH.txt  > T5HH%d_BothBoostedHPure.txt "%(hino,hino,hino))
os.system("combine -M AsymptoticLimits -n T5HH%d_BothBoostedHPure T5HH%d_BothBoostedHPure.txt " %(hino,hino))

#ResolvedOnly
# os.system("combineCards.py SignalRegionRes4bTChiHH%d.txt ControlRegionResBTChiHH%d.txt ControlRegionResCTChiHH%d.txt  ControlRegionResDTChiHH%d.txt>TChiHH%dRes4b.txt "%(hino,hino,hino,hino,hino))
# os.system("combineCards.py SignalRegionRes3bTChiHH%d.txt ControlRegionResB3bTChiHH%d.txt >TChiHH%dRes3b.txt "%(hino,hino,hino))
# os.system("combineCards.py TChiHH%dRes4b.txt  TChiHH%dRes3b.txt > TChiHH%dRes.txt "%(hino,hino,hino))
# os.system("combine -M AsymptoticLimits -n TChiHH%dRes TChiHH%dRes.txt " %(hino,hino))

#Combo
# os.system("combineCards.py SignalRegionRes4bTChiHH%d.txt ControlRegionResBTChiHH%d.txt ControlRegionResCTChiHH%d.txt  ControlRegionResDTChiHH%d.txt>TChiHH%dRes4b.txt "%(hino,hino,hino,hino,hino))
# os.system("combineCards.py SignalRegionRes3bTChiHH%d.txt ControlRegionResB3bTChiHH%d.txt >TChiHH%dRes3b.txt "%(hino,hino,hino))
# os.system("combineCards.py SignalRegionTChiHH%d.txt ControlRegionBTChiHH%d.txt ControlRegionCTChiHH%d.txt  ControlRegionDTChiHH%d.txt > TChiHH%d_2BoostedH_veto.txt "%(hino,hino,hino,hino,hino))
# os.system("combineCards.py TChiHH%d_2BoostedH_veto.txt TChiHH%dRes4b.txt  TChiHH%dRes3b.txt > TChiHH%dCombo.txt "%(hino,hino,hino,hino))
# os.system("combine -M AsymptoticLimits -n TChiHH%dCombo TChiHH%dCombo.txt " %(hino,hino))
