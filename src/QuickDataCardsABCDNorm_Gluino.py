import os
from ROOT import *
import sys

def binContent(histo,binNum):
	thisBinContent = histo.GetBinContent(binNum);
	if (thisBinContent>0):
	 	return thisBinContent;
  	else:
		 return 0.000001;

f=TFile("ALPHABETBoost_V17_all.root", "READ"); #BoostedOnly
fVeto=TFile("ALPHABETBoost_V17_resVeto_all.root", "READ"); #Boosted with resolved veto


# f=TFile("ALPHABETInputTest.root", "READ"); #BoostedOnly


#f=TFile("ALPHABETBoostMC2016_V18.root", "READ"); #BoostedOnly
#f=TFile("ALPHABETBoostMC2016_V18_resVeto.root", "READ"); #Boosted with resolved veto

ZregionA=f.Get("MET_doubletagSR_ZJets"); ZregionB=f.Get("MET_doubletagSB_ZJets");
ZregionC=f.Get("MET_antitagSR_ZJets"); ZregionD=f.Get("MET_antitagSB_ZJets");
WregionA=f.Get("MET_doubletagSR_WJets"); WregionB=f.Get("MET_doubletagSB_WJets");
WregionC=f.Get("MET_antitagSR_WJets"); WregionD=f.Get("MET_antitagSB_WJets");
TregionA=f.Get("MET_doubletagSR_TT"); TregionB=f.Get("MET_doubletagSB_TT");
TregionC=f.Get("MET_antitagSR_TT"); TregionD=f.Get("MET_antitagSB_TT");
QregionA=f.Get("MET_doubletagSR_QCD"); QregionB=f.Get("MET_doubletagSB_QCD");
QregionC=f.Get("MET_antitagSR_QCD"); QregionD=f.Get("MET_antitagSB_QCD");

#for vetoed
sumregionA=fVeto.Get("MET_doubletagSR_sum"); sumregionB=fVeto.Get("MET_doubletagSB_sum");
sumregionC=fVeto.Get("MET_antitagSR_sum"); sumregionD=fVeto.Get("MET_antitagSB_sum");


scale=1.0;
ZregionA.Scale(scale); ZregionB.Scale(scale);
ZregionC.Scale(scale); ZregionD.Scale(scale);
WregionA.Scale(scale); WregionB.Scale(scale);
WregionC.Scale(scale); WregionD.Scale(scale);
TregionA.Scale(scale); TregionB.Scale(scale);
TregionC.Scale(scale); TregionD.Scale(scale);
QregionA.Scale(scale); QregionB.Scale(scale);
QregionC.Scale(scale); QregionD.Scale(scale);

hino=int(sys.argv[1])
# TotalEvents=float(sys.argv[2])
print("\n\n\n-----------------------------------------------------")
print("Looping over mass=%d now" %(hino))
print("-----------------------------------------------------\n")
SignalA=f.Get("MET_doubletagSR_T5HH%d" %(hino));
SignalB=f.Get("MET_doubletagSB_T5HH%d" %(hino));
SignalC=f.Get("MET_antitagSR_T5HH%d" %(hino));
SignalD=f.Get("MET_antitagSB_T5HH%d" %(hino));
scale=1.0;#137000.0/35862.824*4.0; #*4.0 for gluino is in ALPHABET now
# scale=4.0;#137000.0/35862.824*4.0; #*4.0 for gluino is in ALPHABET now
SignalA.Scale(scale)
SignalB.Scale(scale)
SignalC.Scale(scale)
SignalD.Scale(scale)

# YieldsFile=TFile("OutputYields.root", "READ"); #BoostedOnly
#YieldsFile=TFile("OutputYields_veto.root", "READ"); #Boosted with resolved veto

TotalBackground=ZregionA.Clone("TotalBackground");
TotalBackground.Add(WregionA);
TotalBackground.Add(TregionA);
TotalBackground.Add(QregionA);
ControlB=ZregionB.Clone("ControlB");
ControlB.Add(WregionB)
ControlB.Add(TregionB)
ControlB.Add(QregionB)
ControlC=ZregionC.Clone("ControlC")
ControlC.Add(WregionC)
ControlC.Add(TregionC)
ControlC.Add(QregionC)
ControlD=ZregionD.Clone("ControlD")
ControlD.Add(WregionD)
ControlD.Add(TregionD)
ControlD.Add(QregionD)

Predictions=ControlB.Clone("Predictions");
Predictions.Multiply(ControlC)
Predictions.Divide(ControlD)

#vetoed
SignalAVeto=fVeto.Get("MET_doubletagSR_T5HH%d" %(hino));
SignalBVeto=fVeto.Get("MET_doubletagSB_T5HH%d" %(hino));
SignalCVeto=fVeto.Get("MET_antitagSR_T5HH%d" %(hino));
SignalDVeto=fVeto.Get("MET_antitagSB_T5HH%d" %(hino));


TotalBackgroundVeto=sumregionA.Clone("TotalBackgroundVeto");
ControlBVeto=sumregionB.Clone("ControlBVeto");
ControlCVeto=sumregionC.Clone("ControlCVeto");
ControlDVeto=sumregionD.Clone("ControlDVeto");

PredictionsVeto=ControlBVeto.Clone("PredictionsVeto");
PredictionsVeto.Multiply(ControlCVeto)
PredictionsVeto.Divide(ControlDVeto)

bkgfrac=[0.85,0.1254,0.02367]
kappas=[1.0,1.0,1.0] #sim/pred, these are hardcoded
norm=ControlB.Integral(1,4)*ControlC.Integral(1,4)/ControlD.Integral(1,4)
normVeto=ControlBVeto.Integral(1,4)*ControlCVeto.Integral(1,4)/ControlDVeto.Integral(1,4)

for i in range(1,4):
	Predictions.SetBinContent(i, norm*bkgfrac[i-1]*kappas[0])
	#print ("zvv B C D %g %g %g " %(ZregionB.GetBinContent(i),ZregionC.GetBinContent(i), ZregionD.GetBinContent(i)))
	# print "B %g C %g D %g " %(ControlB.Integral(1,4), ControlC.Integral(1,4), ControlD.Integral(1,4))
	# print "Pred %g " %(Predictions.GetBinContent(i))
	PredictionsVeto.SetBinContent(i, normVeto*bkgfrac[i-1]*kappas[0])

#Predictions=f.Get("BackgroundPred")

#BoostedOnly
fcard=open("SignalRegionTemplateMergeBkg.txt", 'r');
fcardout=open("SignalRegionT5HH%d.txt" %hino , 'w')
for line in fcard:
	if "observation" in line:
		newline="observation %.3f %.3f %.3f" %(Predictions.GetBinContent(1),Predictions.GetBinContent(2),Predictions.GetBinContent(3)+Predictions.GetBinContent(4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline+" %g %g %g %g %g %g \n"  %(binContent(SignalA,1),norm, binContent(SignalA,2), norm,SignalA.GetBinContent(3)+SignalA.GetBinContent(4),norm)
		fcardout.write(newline);
		print newline
	else:
		fcardout.write(line);
fcardout.close()
fcard.close()

fcardout=open("ControlRegionBT5HH%d.txt" %hino , 'w')
fcard=open("RegionBTemplateScaleMergeBkg.txt", 'r');
for line in fcard:
	if "observation" in line:
		print line
		newline="observation %.3f " %(ControlB.Integral(1,4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g \n" %(SignalB.Integral(1,4),ControlB.Integral(1,4))
		fcardout.write(newline);
		print newline
	else:
		fcardout.write(line);
fcardout.close()

fcardout=open("ControlRegionCT5HH%d.txt" %hino , 'w')
fcard=open("RegionCTemplateScaleMergeBkg.txt", 'r');
for line in fcard:
	if "observation" in line:
		print line
		newline="observation %.3f " %(ControlC.Integral(1,4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g  \n" %(SignalC.Integral(1,4),ControlC.Integral(1,4))
		fcardout.write(newline);
		print newline
	else:
		fcardout.write(line);
fcardout.close()


fcardout=open("ControlRegionDT5HH%d.txt" %hino , 'w')
fcard=open("RegionDTemplateScaleMergeBkg.txt", 'r');
for line in fcard:
	if "observation" in line:
		print line
		newline="observation %.3f " %(ControlD.Integral(1,4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g \n" %(SignalD.Integral(1,4),ControlD.Integral(1,4))
		fcardout.write(newline);
		print newline
	else:
		fcardout.write(line);
fcardout.close()

#Boosted w/ veto
fcard=open("SignalRegionTemplateMergeBkg.txt", 'r');
fcardout=open("SignalRegionT5HH%d_veto.txt" %hino , 'w')
for line in fcard:
	if "observation" in line:
		newline="observation %.3f %.3f %.3f" %(PredictionsVeto.GetBinContent(1),PredictionsVeto.GetBinContent(2),PredictionsVeto.GetBinContent(3)+PredictionsVeto.GetBinContent(4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline+" %g %g %g %g %g %g \n"  %(binContent(SignalAVeto,1),normVeto, binContent(SignalAVeto,2), normVeto,SignalAVeto.GetBinContent(3)+SignalAVeto.GetBinContent(4),normVeto)
		fcardout.write(newline);
		print newline
	else:
		fcardout.write(line);
fcardout.close()
fcard.close()

fcardout=open("ControlRegionBT5HH%d_veto.txt" %hino , 'w')
fcard=open("RegionBTemplateScaleMergeBkg.txt", 'r');
for line in fcard:
	if "observation" in line:
		print line
		newline="observation %.3f " %(ControlBVeto.Integral(1,4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g \n" %(SignalBVeto.Integral(1,4),ControlBVeto.Integral(1,4))
		fcardout.write(newline);
		print newline
	else:
		fcardout.write(line);
fcardout.close()

fcardout=open("ControlRegionCT5HH%d_veto.txt" %hino , 'w')
fcard=open("RegionCTemplateScaleMergeBkg.txt", 'r');
for line in fcard:
	if "observation" in line:
		print line
		newline="observation %.3f " %(ControlCVeto.Integral(1,4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g  \n" %(SignalCVeto.Integral(1,4),ControlCVeto.Integral(1,4))
		fcardout.write(newline);
		print newline
	else:
		fcardout.write(line);
fcardout.close()


fcardout=open("ControlRegionDT5HH%d_veto.txt" %hino , 'w')
fcard=open("RegionDTemplateScaleMergeBkg.txt", 'r');
for line in fcard:
	if "observation" in line:
		print line
		newline="observation %.3f " %(ControlDVeto.Integral(1,4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g \n" %(SignalDVeto.Integral(1,4),ControlDVeto.Integral(1,4))
		fcardout.write(newline);
		print newline
	else:
		fcardout.write(line);
fcardout.close()

########################################################################
######################### Resolved Cards ###############################
########################################################################

fres=TFile("ALPHABETRes_V17_all.root", "READ");

RegARes4bLow=fres.Get("MET_fourbSRLow_sum")
RegARes4bHigh=fres.Get("MET_fourbSRHigh_sum")

RegBRes4bLow=fres.Get("MET_fourbSBLow_sum")
RegBRes4bHigh=fres.Get("MET_fourbSBHigh_sum")

RegARes3bLow=fres.Get("MET_threebSRLow_sum")
RegARes3bHigh=fres.Get("MET_threebSRHigh_sum")

RegBRes3bLow=fres.Get("MET_threebSBLow_sum")
RegBRes3bHigh=fres.Get("MET_threebSBHigh_sum")

RegCResLow=fres.Get("MET_twobSRLow_sum")
RegCResHigh=fres.Get("MET_twobSRHigh_sum")

RegDResLow=fres.Get("MET_twobSBLow_sum")
RegDResHigh=fres.Get("MET_twobSBHigh_sum")


PredictionsRes4bLow=RegBRes4bLow.Clone("PredictionsRes4bLow");
PredictionsRes4bLow.Multiply(RegCResLow)
PredictionsRes4bLow.Divide(RegDResLow)

PredictionsRes4bHigh=RegBRes4bHigh.Clone("PredictionsRes4bHigh");
PredictionsRes4bHigh.Multiply(RegCResHigh)
PredictionsRes4bHigh.Divide(RegDResHigh)

PredictionsRes3bLow=RegBRes3bLow.Clone("PredictionsRes3bLow");
PredictionsRes3bLow.Multiply(RegCResLow)
PredictionsRes3bLow.Divide(RegDResLow)

PredictionsRes3bHigh=RegBRes3bHigh.Clone("PredictionsRes3bHigh");
PredictionsRes3bHigh.Multiply(RegCResHigh)
PredictionsRes3bHigh.Divide(RegDResHigh)


SignalResA4bLow=fres.Get("MET_fourbSRLow_T5HH%d" %(hino));
SignalResA4bHigh=fres.Get("MET_fourbSRHigh_T5HH%d" %(hino));
SignalResB4bLow=fres.Get("MET_fourbSBLow_T5HH%d" %(hino));
SignalResB4bHigh=fres.Get("MET_fourbSBHigh_T5HH%d" %(hino));
SignalResA3bLow=fres.Get("MET_threebSRLow_T5HH%d" %(hino));
SignalResA3bHigh=fres.Get("MET_threebSRHigh_T5HH%d" %(hino));
SignalResB3bLow=fres.Get("MET_threebSBLow_T5HH%d" %(hino));
SignalResB3bHigh=fres.Get("MET_threebSBHigh_T5HH%d" %(hino));
SignalResCLow=fres.Get("MET_twobSRLow_T5HH%d" %(hino));
SignalResCHigh=fres.Get("MET_twobSRHigh_T5HH%d" %(hino));
SignalResDLow=fres.Get("MET_twobSBLow_T5HH%d" %(hino));
SignalResDHigh=fres.Get("MET_twobSBHigh_T5HH%d" %(hino));

# SignalResA.Scale(scale); SignalResA3b.Scale(scale);
# SignalResB.Scale(scale); SignalResB3b.Scale(scale);
# SignalResC.Scale(scale); SignalResD.Scale(scale);
kappas4b=[0.817174,1.13345,0.709114,1.0] #these are hardcoded
kappas3b=[0.932591,0.99398,0.824331,1.33238] #these are hardcoded


#4b low
fcard=open("SignalRegionResTemplate.txt", 'r');
fcardout=open("SignalRegionRes4bLowT5HH%d.txt" %hino , 'w')
#print PredictionsRes4bLow.GetBinContent(1)*kappas4b[0],binContent(RegARes4bLow,1)
for line in fcard:
	if "observation" in line:
		#print line
		newline="observation %g %g %g %g " %(PredictionsRes4bLow.GetBinContent(1)*kappas4b[0],PredictionsRes4bLow.GetBinContent(2)*kappas4b[1],PredictionsRes4bLow.GetBinContent(3)*kappas4b[2],PredictionsRes4bLow.GetBinContent(4)*kappas4b[3])
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g "  %(binContent(SignalResA4bLow,1),binContent(RegARes4bLow,1))
		newline=newline +"%g %g "  %(binContent(SignalResA4bLow,2),binContent(RegARes4bLow,2))
		newline=newline +"%g %g "  %(binContent(SignalResA4bLow,3),binContent(RegARes4bLow,3))
		newline=newline +"%g %g \n"%(binContent(SignalResA4bLow,4)+binContent(SignalResA4bLow,5),binContent(RegARes4bLow,4)+binContent(RegARes4bLow,5))
		fcardout.write(newline);
	else:
		fcardout.write(line);
fcardout.close()
fcard.close()

#4b high
fcard=open("SignalRegionResTemplate.txt", 'r');
fcardout=open("SignalRegionRes4bHighT5HH%d.txt" %hino , 'w')
#print PredictionsRes4bHigh.GetBinContent(1)*kappas4b[0],binContent(RegARes4bHigh,1)
for line in fcard:
	if "observation" in line:
		#print line
		newline="observation %g %g %g %g " %(PredictionsRes4bHigh.GetBinContent(1)*kappas4b[0],PredictionsRes4bHigh.GetBinContent(2)*kappas4b[1],PredictionsRes4bHigh.GetBinContent(3)*kappas4b[2],PredictionsRes4bHigh.GetBinContent(4)*kappas4b[3])
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g "  %(binContent(SignalResA4bHigh,1),binContent(RegARes4bHigh,1))
		newline=newline +"%g %g "  %(binContent(SignalResA4bHigh,2),binContent(RegARes4bHigh,2))
		newline=newline +"%g %g "  %(binContent(SignalResA4bHigh,3),binContent(RegARes4bHigh,3))
		newline=newline +"%g %g \n"%(binContent(SignalResA4bHigh,4)+binContent(SignalResA4bHigh,5),binContent(RegARes4bHigh,4)+binContent(RegARes4bHigh,5))
		fcardout.write(newline);
	else:
		fcardout.write(line);
fcardout.close()
fcard.close()

#3b low
fcard=open("SignalRegionResTemplate.txt", 'r');
fcardout=open("SignalRegionRes3bLowT5HH%d.txt" %hino , 'w')
#print PredictionsRes3bLow.GetBinContent(1)*kappas3b[0],binContent(RegARes3bLow,1)
for line in fcard:
	if "observation" in line:
		#print line
		newline="observation %g %g %g %g " %(PredictionsRes3bLow.GetBinContent(1)*kappas3b[0],PredictionsRes3bLow.GetBinContent(2)*kappas3b[1],PredictionsRes3bLow.GetBinContent(3)*kappas3b[2],PredictionsRes3bLow.GetBinContent(4)*kappas3b[3])
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g "  %(binContent(SignalResA3bLow,1),binContent(RegARes3bLow,1))
		newline=newline +"%g %g "  %(binContent(SignalResA3bLow,2),binContent(RegARes3bLow,2))
		newline=newline +"%g %g "  %(binContent(SignalResA3bLow,3),binContent(RegARes3bLow,3))
		newline=newline +"%g %g \n"%(binContent(SignalResA3bLow,4)+binContent(SignalResA3bLow,5),binContent(RegARes3bLow,4)+binContent(RegARes3bLow,5))
		fcardout.write(newline);
	else:
		fcardout.write(line);
fcardout.close()
fcard.close()

#3b high
fcard=open("SignalRegionResTemplate.txt", 'r');
fcardout=open("SignalRegionRes3bHighT5HH%d.txt" %hino , 'w')
#print PredictionsRes3bHigh.GetBinContent(1)*kappas3b[0],binContent(RegARes3bHigh,1)
for line in fcard:
	if "observation" in line:
		#print line
		newline="observation %g %g %g %g " %(PredictionsRes3bHigh.GetBinContent(1)*kappas3b[0],PredictionsRes3bHigh.GetBinContent(2)*kappas3b[1],PredictionsRes3bHigh.GetBinContent(3)*kappas3b[2],PredictionsRes3bHigh.GetBinContent(4)*kappas3b[3])
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g "  %(binContent(SignalResA3bHigh,1),binContent(RegARes3bHigh,1))
		newline=newline +"%g %g "  %(binContent(SignalResA3bHigh,2),binContent(RegARes3bHigh,2))
		newline=newline +"%g %g "  %(binContent(SignalResA3bHigh,3),binContent(RegARes3bHigh,3))
		newline=newline +"%g %g \n"%(binContent(SignalResA3bHigh,4)+binContent(SignalResA3bHigh,5),binContent(RegARes3bHigh,4)+binContent(RegARes3bHigh,5))
		fcardout.write(newline);
	else:
		fcardout.write(line);
fcardout.close()
fcard.close()

#Region B
#4b low
fcardout=open("ControlRegionResB4bLowT5HH%d.txt" %hino , 'w')
fcard=open("RegionBTemplateRes.txt", 'r');
for line in fcard:
	if "observation" in line:
		#print line
		newline="observation %g %g %g %g " %(RegBRes4bLow.GetBinContent(1), RegBRes4bLow.GetBinContent(2),RegBRes4bLow.GetBinContent(3),RegBRes4bLow.GetBinContent(4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g "  %(binContent(SignalResB4bLow,1),binContent(RegBRes4bLow,1))
		newline=newline +"%g %g "  %(binContent(SignalResB4bLow,2),binContent(RegBRes4bLow,2))
		newline=newline +"%g %g "  %(binContent(SignalResB4bLow,3),binContent(RegBRes4bLow,3))
		newline=newline +"%g %g \n"%(binContent(SignalResB4bLow,4)+binContent(SignalResB4bLow,5),binContent(RegBRes4bLow,4)+binContent(RegBRes4bLow,5))
		fcardout.write(newline);
		#print newline
	else:
		fcardout.write(line);
fcardout.close()

#4b high
fcardout=open("ControlRegionResB4bHighT5HH%d.txt" %hino , 'w')
fcard=open("RegionBTemplateRes.txt", 'r');
for line in fcard:
	if "observation" in line:
		#print line
		newline="observation %g %g %g %g " %(RegBRes4bHigh.GetBinContent(1), RegBRes4bHigh.GetBinContent(2),RegBRes4bHigh.GetBinContent(3),RegBRes4bHigh.GetBinContent(4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g "  %(binContent(SignalResB4bHigh,1),binContent(RegBRes4bHigh,1))
		newline=newline +"%g %g "  %(binContent(SignalResB4bHigh,2),binContent(RegBRes4bHigh,2))
		newline=newline +"%g %g "  %(binContent(SignalResB4bHigh,3),binContent(RegBRes4bHigh,3))
		newline=newline +"%g %g \n"%(binContent(SignalResB4bHigh,4)+binContent(SignalResB4bHigh,5),binContent(RegBRes4bHigh,4)+binContent(RegBRes4bHigh,5))
		fcardout.write(newline);
		#print newline
	else:
		fcardout.write(line);
fcardout.close()

#3b low
fcardout=open("ControlRegionResB3bLowT5HH%d.txt" %hino , 'w')
fcard=open("RegionBTemplateRes.txt", 'r');
for line in fcard:
	if "observation" in line:
		#print line
		newline="observation %g %g %g %g " %(RegBRes3bLow.GetBinContent(1), RegBRes3bLow.GetBinContent(2),RegBRes3bLow.GetBinContent(3),RegBRes3bLow.GetBinContent(4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g "  %(binContent(SignalResB3bLow,1),binContent(RegBRes3bLow,1))
		newline=newline +"%g %g "  %(binContent(SignalResB3bLow,2),binContent(RegBRes3bLow,2))
		newline=newline +"%g %g "  %(binContent(SignalResB3bLow,3),binContent(RegBRes3bLow,3))
		newline=newline +"%g %g \n"%(binContent(SignalResB3bLow,4)+binContent(SignalResB3bLow,5),binContent(RegBRes3bLow,4)+binContent(RegBRes3bLow,5))
		fcardout.write(newline);
		#print newline
	else:
		fcardout.write(line);
fcardout.close()

#3b high
fcardout=open("ControlRegionResB3bHighT5HH%d.txt" %hino , 'w')
fcard=open("RegionBTemplateRes.txt", 'r');
for line in fcard:
	if "observation" in line:
		#print line
		newline="observation %g %g %g %g " %(RegBRes3bHigh.GetBinContent(1), RegBRes3bHigh.GetBinContent(2),RegBRes3bHigh.GetBinContent(3),RegBRes3bHigh.GetBinContent(4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g "  %(binContent(SignalResB3bHigh,1),binContent(RegBRes3bHigh,1))
		newline=newline +"%g %g "  %(binContent(SignalResB3bHigh,2),binContent(RegBRes3bHigh,2))
		newline=newline +"%g %g "  %(binContent(SignalResB3bHigh,3),binContent(RegBRes3bHigh,3))
		newline=newline +"%g %g \n"%(binContent(SignalResB3bHigh,4)+binContent(SignalResB3bHigh,5),binContent(RegBRes3bHigh,4)+binContent(RegBRes3bHigh,5))
		fcardout.write(newline);
		#print newline
	else:
		fcardout.write(line);
fcardout.close()

#Region C
#Low
fcardout=open("ControlRegionResCLowT5HH%d.txt" %hino , 'w')
fcard=open("RegionCTemplateRes.txt", 'r');
for line in fcard:
	if "observation" in line:
		#print line
		newline="observation %g %g %g %g " %(RegCResLow.GetBinContent(1), RegCResLow.GetBinContent(2),RegCResLow.GetBinContent(3),RegCResLow.GetBinContent(4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g "  %(binContent(SignalResCLow,1),binContent(RegCResLow,1))
		newline=newline +"%g %g "  %(binContent(SignalResCLow,2),binContent(RegCResLow,2))
		newline=newline +"%g %g "  %(binContent(SignalResCLow,3),binContent(RegCResLow,3))
		newline=newline +"%g %g \n"%(binContent(SignalResCLow,4)+binContent(SignalResCLow,5),binContent(RegCResLow,4)+binContent(RegCResLow,5))
		fcardout.write(newline);
		#print newline
	else:
		fcardout.write(line);
fcardout.close()

#High
fcardout=open("ControlRegionResCHighT5HH%d.txt" %hino , 'w')
fcard=open("RegionCTemplateRes.txt", 'r');
for line in fcard:
	if "observation" in line:
		#print line
		newline="observation %g %g %g %g " %(RegCResHigh.GetBinContent(1), RegCResHigh.GetBinContent(2),RegCResHigh.GetBinContent(3),RegCResHigh.GetBinContent(4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g "  %(binContent(SignalResCHigh,1),binContent(RegCResHigh,1))
		newline=newline +"%g %g "  %(binContent(SignalResCHigh,2),binContent(RegCResHigh,2))
		newline=newline +"%g %g "  %(binContent(SignalResCHigh,3),binContent(RegCResHigh,3))
		newline=newline +"%g %g \n"%(binContent(SignalResCHigh,4)+binContent(SignalResCHigh,5),binContent(RegCResHigh,4)+binContent(RegCResHigh,5))
		fcardout.write(newline);
		#print newline
	else:
		fcardout.write(line);
fcardout.close()


#Region D
#Low
fcardout=open("ControlRegionResDLowT5HH%d.txt" %hino , 'w')
fcard=open("RegionDTemplateRes.txt", 'r');
for line in fcard:
	if "observation" in line:
		#print line
		newline="observation %g %g %g %g " %(RegDResLow.GetBinContent(1), RegDResLow.GetBinContent(2),RegDResLow.GetBinContent(3),RegDResLow.GetBinContent(4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g "  %(binContent(SignalResDLow,1),binContent(RegDResLow,1))
		newline=newline +"%g %g "  %(binContent(SignalResDLow,2),binContent(RegDResLow,2))
		newline=newline +"%g %g "  %(binContent(SignalResDLow,3),binContent(RegDResLow,3))
		newline=newline +"%g %g \n"%(binContent(SignalResDLow,4)+binContent(SignalResDLow,5),binContent(RegDResLow,4)+binContent(RegDResLow,5))
		fcardout.write(newline);
		#print newline
	else:
		fcardout.write(line);
fcardout.close()

#High
fcardout=open("ControlRegionResDHighT5HH%d.txt" %hino , 'w')
fcard=open("RegionDTemplateRes.txt", 'r');
for line in fcard:
	if "observation" in line:
		#print line
		newline="observation %g %g %g %g " %(RegDResHigh.GetBinContent(1), RegDResHigh.GetBinContent(2),RegDResHigh.GetBinContent(3),RegDResHigh.GetBinContent(4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g "  %(binContent(SignalResDHigh,1),binContent(RegDResHigh,1))
		newline=newline +"%g %g "  %(binContent(SignalResDHigh,2),binContent(RegDResHigh,2))
		newline=newline +"%g %g "  %(binContent(SignalResDHigh,3),binContent(RegDResHigh,3))
		newline=newline +"%g %g \n"%(binContent(SignalResDHigh,4)+binContent(SignalResDHigh,5),binContent(RegDResHigh,4)+binContent(RegDResHigh,5))
		fcardout.write(newline);
		#print newline
	else:
		fcardout.write(line);
fcardout.close()



#BoostedOnly
os.system("combineCards.py SignalRegionT5HH%d.txt ControlRegionBT5HH%d.txt ControlRegionCT5HH%d.txt  ControlRegionDT5HH%d.txt > T5HH%d_2BoostedH.txt "%(hino,hino,hino,hino,hino))
os.system("combine -M AsymptoticLimits -n T5HH%d_2BoostedH T5HH%d_2BoostedH.txt " %(hino,hino))

#ResolvedOnly
os.system("combineCards.py SignalRegionRes4bHighT5HH%d.txt ControlRegionResB4bHighT5HH%d.txt ControlRegionResCHighT5HH%d.txt  ControlRegionResDHighT5HH%d.txt>T5HH%dRes4bHigh.txt "%(hino,hino,hino,hino,hino))
os.system("combineCards.py SignalRegionRes4bLowT5HH%d.txt ControlRegionResB4bLowT5HH%d.txt ControlRegionResCLowT5HH%d.txt  ControlRegionResDLowT5HH%d.txt>T5HH%dRes4bLow.txt "%(hino,hino,hino,hino,hino))
os.system("combineCards.py SignalRegionRes3bHighT5HH%d.txt ControlRegionResB3bHighT5HH%d.txt ControlRegionResCHighT5HH%d.txt  ControlRegionResDHighT5HH%d.txt>T5HH%dRes3bHigh.txt "%(hino,hino,hino,hino,hino))
os.system("combineCards.py SignalRegionRes3bLowT5HH%d.txt ControlRegionResB3bLowT5HH%d.txt ControlRegionResCLowT5HH%d.txt  ControlRegionResDLowT5HH%d.txt>T5HH%dRes3bLow.txt "%(hino,hino,hino,hino,hino))
os.system("combineCards.py T5HH%dRes4bHigh.txt T5HH%dRes4bLow.txt  T5HH%dRes3bHigh.txt T5HH%dRes3bLow.txt > T5HH%dRes.txt "%(hino,hino,hino,hino,hino))
os.system("combine -M AsymptoticLimits -n T5HH%dRes T5HH%dRes.txt " %(hino,hino))

#Combo
#os.system("combineCards.py SignalRegionRes4bT5HH%d.txt ControlRegionResBT5HH%d.txt ControlRegionResCT5HH%d.txt  ControlRegionResDT5HH%d.txt>T5HH%dRes4b.txt "%(hino,hino,hino,hino,hino))
#os.system("combineCards.py SignalRegionRes3bT5HH%d.txt ControlRegionResB3bT5HH%d.txt >T5HH%dRes3b.txt "%(hino,hino,hino))
os.system("combineCards.py SignalRegionT5HH%d_veto.txt ControlRegionBT5HH%d_veto.txt ControlRegionCT5HH%d_veto.txt  ControlRegionDT5HH%d_veto.txt > T5HH%d_2BoostedH_veto.txt "%(hino,hino,hino,hino,hino))
os.system("combineCards.py T5HH%d_2BoostedH_veto.txt T5HH%dRes.txt > T5HH%dCombo.txt "%(hino,hino,hino))
os.system("combine -M AsymptoticLimits -n T5HH%dCombo T5HH%dCombo.txt " %(hino,hino))
