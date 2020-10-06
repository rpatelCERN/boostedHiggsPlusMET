import os
from ROOT import *
import sys

def binContent(histo,binNum):
	thisBinContent = histo.GetBinContent(binNum);
	if (thisBinContent>0):
	 	return thisBinContent;
  	else:
		 return 0.000001;

# f=TFile("ALPHABET_V18_2BoostedH.root", "READ"); #BoostedOnly
f=TFile("ALPHABET_0l_MET3.root", "READ"); #BoostedOnly
fVeto=TFile("ALPHABETBoost_V17_resVeto_all.root", "READ"); #Boosted with resolved veto

sumRegionA2=f.Get("MET_doubletagSR_sum"); sumRegionA1=f.Get("MET_tagSR_sum");
sumRegionB2=f.Get("MET_doubletagSB_sum"); sumRegionB1=f.Get("MET_tagSB_sum");
#sumRegionB2=f.Get("MET_doubletagSB2_sum"); sumRegionB1=f.Get("MET_tagSB2_sum");
sumRegionC=f.Get("MET_antitagSR_sum"); sumRegionD=f.Get("MET_antitagSB_sum");
#sumRegionD=f.Get("MET_antitagSB2_sum");

#for vetoed
vetoSumRegionA2=fVeto.Get("MET_doubletagSR_sum"); vetoSumRegionB2=fVeto.Get("MET_doubletagSB_sum");
vetoSumRegionA1=fVeto.Get("MET_tagSR_sum"); vetoSumRegionB1=fVeto.Get("MET_tagSB_sum");
vetoSumRegionC=fVeto.Get("MET_antitagSR_sum"); vetoSumRegionD=fVeto.Get("MET_antitagSB_sum");

scale=1.0;
sumRegionA2.Scale(scale); sumRegionB2.Scale(scale);
sumRegionA1.Scale(scale); sumRegionB1.Scale(scale);
sumRegionC.Scale(scale); sumRegionD.Scale(scale);

hino=int(sys.argv[1])
LSP=int(sys.argv[2])
# TotalEvents=float(sys.argv[2])
print("\n\n\n-----------------------------------------------------")
print("Looping over mass=%d, mLSP=%d now" %(hino,LSP))
print("-----------------------------------------------------\n")
f2=TFile("ALPHABET_sig_MET3.root", "READ")
# f2=TFile("ALPHABET_V18_1Dsignal_2BoostedH.root", "READ")
# f2=TFile("ALPHABET_V18_2Dsignal_2BoostedH.root", "READ")

SignalA2=f2.Get("MET_doubletagSR_TChiHH_%d_%d" %(hino,LSP));
SignalB2=f2.Get("MET_doubletagSB_TChiHH_%d_%d" %(hino,LSP));
#SignalB2=f2.Get("MET_doubletagSB2_TChiHH_%d_%d" %(hino,LSP));
SignalA1=f2.Get("MET_tagSR_TChiHH_%d_%d" %(hino,LSP));
SignalB1=f2.Get("MET_tagSB_TChiHH_%d_%d" %(hino,LSP));
#SignalB1=f2.Get("MET_tagSB2_TChiHH_%d_%d" %(hino,LSP));
SignalC=f2.Get("MET_antitagSR_TChiHH_%d_%d" %(hino,LSP));
SignalD=f2.Get("MET_antitagSB_TChiHH_%d_%d" %(hino,LSP));
#SignalD=f2.Get("MET_antitagSB2_TChiHH_%d_%d" %(hino,LSP));

# SignalA2=f2.Get("MET_doubletagSR_TChiHH%d_LSP%d" %(hino,LSP));
# SignalB2=f2.Get("MET_doubletagSB_TChiHH%d_LSP%d" %(hino,LSP));
# SignalA1=f2.Get("MET_tagSR_TChiHH%d_LSP%d" %(hino,LSP));
# SignalB1=f2.Get("MET_tagSB_TChiHH%d_LSP%d" %(hino,LSP));
# SignalC=f2.Get("MET_antitagSR_TChiHH%d_LSP%d" %(hino,LSP));
# SignalD=f2.Get("MET_antitagSB_TChiHH%d_LSP%d" %(hino,LSP));

scale=1.0;#137000.0/35862.824;
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

#vetoed
SignalA2Veto=fVeto.Get("MET_doubletagSR_TChiHH%d" %(hino));
SignalB2Veto=fVeto.Get("MET_doubletagSB_TChiHH%d" %(hino));
SignalA1Veto=fVeto.Get("MET_tagSR_TChiHH%d" %(hino));
SignalB1Veto=fVeto.Get("MET_tagSB_TChiHH%d" %(hino));
SignalCVeto=fVeto.Get("MET_antitagSR_TChiHH%d" %(hino));
SignalDVeto=fVeto.Get("MET_antitagSB_TChiHH%d" %(hino));

sumRegionA2Veto=vetoSumRegionA2.Clone("sumRegionA2Veto");
sumRegionB2Veto=vetoSumRegionB2.Clone("sumRegionB2Veto");
sumRegionA1Veto=vetoSumRegionA1.Clone("sumRegionA1Veto");
sumRegionB1Veto=vetoSumRegionB1.Clone("sumRegionB1Veto");
sumRegionCVeto=vetoSumRegionC.Clone("sumRegionCVeto");
sumRegionDVeto=vetoSumRegionD.Clone("sumRegionDVeto");

Predictions2HVeto=sumRegionB2Veto.Clone("Predictions2HVeto");
Predictions2HVeto.Multiply(sumRegionCVeto)
Predictions2HVeto.Divide(sumRegionDVeto)

Predictions1HVeto=sumRegionB1Veto.Clone("Predictions1HVeto");
Predictions1HVeto.Multiply(sumRegionCVeto)
Predictions1HVeto.Divide(sumRegionDVeto)

# bkgfrac2H=[0.851,0.119,0.031] #V18
# bkgfrac1H=[0.851,0.119,0.031] #V18
bkgfrac2H=[0.763,0.190,0.048] #TEST
bkgfrac1H=[0.763,0.190,0.048] #TEST
# kappas2H = 1.08 #V18
# kappas1H = 1.03 #V18
kappas2H = 1.00 #V18
kappas1H = 1.00 #V18
# kappas2H = 1.51 #TEST
# kappas1H = 1.29 #TEST
# kappas2H = 1.0 #V18, original SB, removing kappas
# kappas1H = 1.0 #V18, original SB, removing kappas
#kappas2H = 1.66 #V18, SB2
#kappas1H = 1.19 #V18, SB2

# bkgfrac2H=[0.857,0.118,0.023] #V17
# bkgfrac1H=[0.843,0.130,0.027] #V17
# kappas2H=[1.1,1.1,1.1] #V17
# kappas1H=[1.1,1.1,1.1] #V17

norm2H=sumRegionB2.Integral(1,4)*sumRegionC.Integral(1,4)/sumRegionD.Integral(1,4)
norm2HVeto=sumRegionB2Veto.Integral(1,4)*sumRegionCVeto.Integral(1,4)/sumRegionDVeto.Integral(1,4)

norm1H=sumRegionB1.Integral(1,4)*sumRegionC.Integral(1,4)/sumRegionD.Integral(1,4)
norm1HVeto=sumRegionB1Veto.Integral(1,4)*sumRegionCVeto.Integral(1,4)/sumRegionDVeto.Integral(1,4)

for i in range(1,4):
	Predictions2H.SetBinContent(i, norm2H*bkgfrac2H[i-1]*kappas2H)
	Predictions1H.SetBinContent(i, norm1H*bkgfrac1H[i-1]*kappas1H)

	Predictions2HVeto.SetBinContent(i, norm2HVeto*bkgfrac2H[i-1]*kappas2H)
	Predictions1HVeto.SetBinContent(i, norm1HVeto*bkgfrac1H[i-1]*kappas1H)


#BoostedOnly
fcard=open("SignalRegionTemplateMergeBkg.txt", 'r');
fcardout=open("SRMerge_TChiHH%d_LSP%d.txt" %(hino,LSP) , 'w')
#fcardout=open("SRMerge_TChiHH%d_LSP%d_SB2.txt" %(hino,LSP) , 'w')
for line in fcard:
	if "observation" in line:
		newline="observation %.3f %.3f %.3f" %(Predictions2H.GetBinContent(1),Predictions2H.GetBinContent(2),Predictions2H.GetBinContent(3)+Predictions2H.GetBinContent(4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline+" %g %g %g %g %g %g \n"  %(binContent(SignalA2,1),norm2H, binContent(SignalA2,2), norm2H,SignalA2.GetBinContent(3)+SignalA2.GetBinContent(4),norm2H)
		fcardout.write(newline);
		# print newline
	else:
		fcardout.write(line);
fcardout.close()
fcard.close()


fcard=open("SignalRegion1TemplateMergeBkg.txt", 'r');
fcardout=open("SR1Merge_TChiHH%d_LSP%d.txt" %(hino,LSP) , 'w')
#fcardout=open("SR1Merge_TChiHH%d_LSP%d_SB2.txt" %(hino,LSP) , 'w')
for line in fcard:
	if "observation" in line:
		newline="observation %.3f %.3f %.3f" %(Predictions1H.GetBinContent(1),Predictions1H.GetBinContent(2),Predictions1H.GetBinContent(3)+Predictions1H.GetBinContent(4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline+" %g %g %g %g %g %g \n"  %(binContent(SignalA1,1),norm1H, binContent(SignalA1,2), norm1H,SignalA1.GetBinContent(3)+SignalA1.GetBinContent(4),norm1H)
		fcardout.write(newline);
		# print newline
	else:
		fcardout.write(line);
fcardout.close()
fcard.close()

fcardout=open("CRBMerge_TChiHH%d_LSP%d.txt" %(hino,LSP) , 'w')
#fcardout=open("CRBMerge_TChiHH%d_LSP%d_SB2.txt" %(hino,LSP) , 'w')
fcard=open("RegionBTemplateScaleMergeBkg.txt", 'r');
for line in fcard:
	if "observation" in line:
		# print line
		newline="observation %.3f " %(sumRegionB2.Integral(1,4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g \n" %(SignalB2.Integral(1,4),sumRegionB2.Integral(1,4))
		fcardout.write(newline);
		# print newline
	else:
		fcardout.write(line);
fcardout.close()

fcardout=open("CRB1Merge_TChiHH%d_LSP%d.txt" %(hino,LSP) , 'w')
#fcardout=open("CRB1Merge_TChiHH%d_LSP%d_SB2.txt" %(hino,LSP) , 'w')
fcard=open("RegionB1TemplateScaleMergeBkg.txt", 'r');
for line in fcard:
	if "observation" in line:
		# print line
		newline="observation %.3f " %(sumRegionB1.Integral(1,4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g \n" %(SignalB1.Integral(1,4),sumRegionB1.Integral(1,4))
		fcardout.write(newline);
		# print newline
	else:
		fcardout.write(line);
fcardout.close()

fcardout=open("CRCMerge_TChiHH%d_LSP%d.txt" %(hino,LSP) , 'w')
#fcardout=open("CRCMerge_TChiHH%d_LSP%d_SB2.txt" %(hino,LSP) , 'w')
fcard=open("RegionCTemplateScaleMergeBkg.txt", 'r');
for line in fcard:
	if "observation" in line:
		# print line
		newline="observation %.3f " %(sumRegionC.Integral(1,4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g  \n" %(SignalC.Integral(1,4),sumRegionC.Integral(1,4))
		fcardout.write(newline);
		# print newline
	else:
		fcardout.write(line);
fcardout.close()


fcardout=open("CRDMerge_TChiHH%d_LSP%d.txt" %(hino,LSP) , 'w')
#fcardout=open("CRDMerge_TChiHH%d_LSP%d_SB2.txt" %(hino,LSP) , 'w')
fcard=open("RegionDTemplateScaleMergeBkg.txt", 'r');
for line in fcard:
	if "observation" in line:
		# print line
		newline="observation %.3f " %(sumRegionD.Integral(1,4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g \n" %(SignalD.Integral(1,4),sumRegionD.Integral(1,4))
		fcardout.write(newline);
		# print newline
	else:
		fcardout.write(line);
fcardout.close()

#Boosted w/ veto
fcard=open("SignalRegionVetoTemplateMergeBkg.txt", 'r');
fcardout=open("SRMerge_TChiHH%d_LSP%d_veto.txt" %(hino,LSP) , 'w')
for line in fcard:
	if "observation" in line:
		newline="observation %.3f %.3f %.3f" %(Predictions2HVeto.GetBinContent(1),Predictions2HVeto.GetBinContent(2),Predictions2HVeto.GetBinContent(3)+Predictions2HVeto.GetBinContent(4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline+" %g %g %g %g %g %g \n"  %(binContent(SignalA2Veto,1),norm2HVeto, binContent(SignalA2Veto,2), norm2HVeto,SignalA2Veto.GetBinContent(3)+SignalA2Veto.GetBinContent(4),norm2HVeto)
		fcardout.write(newline);
		# print newline
	else:
		fcardout.write(line);
fcardout.close()
fcard.close()

fcard=open("SignalRegionVeto1TemplateMergeBkg.txt", 'r');
fcardout=open("SR1Merge_TChiHH%d_LSP%d_veto.txt" %(hino,LSP) , 'w')
for line in fcard:
	if "observation" in line:
		newline="observation %.3f %.3f %.3f" %(Predictions1HVeto.GetBinContent(1),Predictions1HVeto.GetBinContent(2),Predictions1HVeto.GetBinContent(3)+Predictions1HVeto.GetBinContent(4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline+" %g %g %g %g %g %g \n"  %(binContent(SignalA1Veto,1),norm1HVeto, binContent(SignalA1Veto,2), norm1HVeto,SignalA1Veto.GetBinContent(3)+SignalA1Veto.GetBinContent(4),norm1HVeto)
		fcardout.write(newline);
		# print newline
	else:
		fcardout.write(line);
fcardout.close()
fcard.close()

fcardout=open("CRBMerge_TChiHH%d_LSP%d_veto.txt" %(hino,LSP) , 'w')
fcard=open("RegionBTemplateScaleMergeBkg.txt", 'r');
for line in fcard:
	if "observation" in line:
		# print line
		newline="observation %.3f " %(sumRegionB2Veto.Integral(1,4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g \n" %(SignalB2Veto.Integral(1,4),sumRegionB2Veto.Integral(1,4))
		fcardout.write(newline);
		# print newline
	else:
		fcardout.write(line);
fcardout.close()


fcardout=open("CRB1Merge_TChiHH%d_LSP%d_veto.txt" %(hino,LSP) , 'w')
fcard=open("RegionB1TemplateScaleMergeBkg.txt", 'r');
for line in fcard:
	if "observation" in line:
		# print line
		newline="observation %.3f " %(sumRegionB1Veto.Integral(1,4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g \n" %(SignalB1Veto.Integral(1,4),sumRegionB1Veto.Integral(1,4))
		fcardout.write(newline);
		# print newline
	else:
		fcardout.write(line);
fcardout.close()


fcardout=open("CRCMerge_TChiHH%d_LSP%d_veto.txt" %(hino,LSP) , 'w')
fcard=open("RegionCTemplateScaleMergeBkg.txt", 'r');
for line in fcard:
	if "observation" in line:
		# print line
		newline="observation %.3f " %(sumRegionCVeto.Integral(1,4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g  \n" %(SignalCVeto.Integral(1,4),sumRegionCVeto.Integral(1,4))
		fcardout.write(newline);
		# print newline
	else:
		fcardout.write(line);
fcardout.close()


fcardout=open("CRDMerge_TChiHH%d_LSP%d_veto.txt" %(hino,LSP) , 'w')
fcard=open("RegionDTemplateScaleMergeBkg.txt", 'r');
for line in fcard:
	if "observation" in line:
		# print line
		newline="observation %.3f " %(sumRegionDVeto.Integral(1,4))
		fcardout.write(newline);
	elif "rate" in line and not "rateParam" in line:
		newline="rate "
		newline=newline +"%g %g \n" %(SignalDVeto.Integral(1,4),sumRegionDVeto.Integral(1,4))
		fcardout.write(newline);
		# print newline
	else:
		fcardout.write(line);
fcardout.close()

########################################################################
######################### Resolved Cards ###############################
########################################################################
'''
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

SignalResA=fres.Get("MET_fourbSR_T5HH%d" %(hino));
SignalResB=fres.Get("MET_fourbSB_T5HH%d" %(hino));
SignalResA3b=fres.Get("MET_threebSR_T5HH%d" %(hino));
SignalResB3b=fres.Get("MET_threebSB_T5HH%d" %(hino));
SignalResC=fres.Get("MET_twobSR_T5HH%d" %(hino));
SignalResD=fres.Get("MET_twobSB_T5HH%d" %(hino));

SignalResA.Scale(scale); SignalResA3b.Scale(scale);
SignalResB.Scale(scale); SignalResB3b.Scale(scale);
SignalResC.Scale(scale); SignalResD.Scale(scale);

fcard=open("SignalRegionResTemplateScale.txt", 'r');
fcardout=open("SignalRegionRes4bT5HH%d.txt" %hino , 'w')
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

fcard=open("SignalRegionResTemplateScale.txt", 'r');
fcardout=open("SignalRegionRes3bT5HH%d.txt" %hino , 'w')
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

fcardout=open("ControlRegionResBT5HH%d.txt" %hino , 'w')
fcard=open("RegionBTemplateScaleRes.txt", 'r');
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

fcardout=open("ControlRegionResCT5HH%d.txt" %hino , 'w')
fcard=open("RegionCTemplateScaleRes.txt", 'r');
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

fcardout=open("ControlRegionResB3bT5HH%d.txt" %hino , 'w')
fcard=open("RegionBTemplateScaleRes.txt", 'r');
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

fcardout=open("ControlRegionResDT5HH%d.txt" %hino , 'w')
fcard=open("RegionDTemplateScaleRes.txt", 'r');
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
#os.system("combineCards.py SRMerge_TChiHH%d_LSP%d_SB2.txt CRBMerge_TChiHH%d_LSP%d_SB2.txt CRCMerge_TChiHH%d_LSP%d_SB2.txt  CRDMerge_TChiHH%d_LSP%d_SB2.txt > TChiHH%d_LSP%d_2BoostedHMerge_SB2.txt "%(hino,LSP,hino,LSP,hino,LSP,hino,LSP,hino,LSP))
#os.system("combineCards.py SR1Merge_TChiHH%d_LSP%d_SB2.txt CRB1Merge_TChiHH%d_LSP%d_SB2.txt CRCMerge_TChiHH%d_LSP%d_SB2.txt  CRDMerge_TChiHH%d_LSP%d_SB2.txt > TChiHH%d_LSP%d_1BoostedHMerge_SB2.txt "%(hino,LSP,hino,LSP,hino,LSP,hino,LSP,hino,LSP))
#os.system("combineCards.py TChiHH%d_LSP%d_1BoostedHMerge_SB2.txt TChiHH%d_LSP%d_2BoostedHMerge_SB2.txt  > TChiHH%d_LSP%d_BothBoostedHMerge_SB2.txt "%(hino,LSP,hino,LSP,hino,LSP))
#os.system("combine -M AsymptoticLimits -n TChiHH%d_LSP%d_BothBoostedHMerge_SB2 TChiHH%d_LSP%d_BothBoostedHMerge_SB2.txt " %(hino,LSP,hino,LSP))


os.system("combineCards.py SRMerge_TChiHH%d_LSP%d.txt CRBMerge_TChiHH%d_LSP%d.txt CRCMerge_TChiHH%d_LSP%d.txt  CRDMerge_TChiHH%d_LSP%d.txt > TChiHH%d_LSP%d_2H_MET3.txt "%(hino,LSP,hino,LSP,hino,LSP,hino,LSP,hino,LSP))
os.system("combineCards.py SR1Merge_TChiHH%d_LSP%d.txt CRB1Merge_TChiHH%d_LSP%d.txt CRCMerge_TChiHH%d_LSP%d.txt  CRDMerge_TChiHH%d_LSP%d.txt > TChiHH%d_LSP%d_1H_MET3.txt "%(hino,LSP,hino,LSP,hino,LSP,hino,LSP,hino,LSP))
os.system("combineCards.py TChiHH%d_LSP%d_1H_MET3.txt TChiHH%d_LSP%d_2H_MET3.txt  > TChiHH%d_LSP%d_MET3.txt "%(hino,LSP,hino,LSP,hino,LSP))
os.system("combine -M AsymptoticLimits -n TChiHH%d_LSP%d_MET3 TChiHH%d_LSP%d_MET3.txt " %(hino,LSP,hino,LSP))

# os.system("combine -M AsymptoticLimits -n TChiHH%d_LSP%d_2BoostedH TChiHH%d_LSP%d_2BoostedH.txt " %(hino,LSP,hino,LSP))
# os.system("combine -M AsymptoticLimits -n TChiHH%d_LSP%d_1BoostedH TChiHH%d_LSP%d_1BoostedH.txt " %(hino,LSP,hino,LSP))

#BoostedOnlyVeto
#os.system("combineCards.py SignalRegionTChiHH%d_LSP%d_veto.txt ControlRegionBTChiHH%d_LSP%d_veto.txt ControlRegionCTChiHH%d_LSP%d_veto.txt  ControlRegionDTChiHH%d_LSP%d_veto.txt > TChiHH%d_LSP%d_2BoostedH_veto.txt "%(hino,LSP,hino,LSP,hino,LSP,hino,LSP,hino,LSP))
#os.system("combineCards.py SignalRegion1TChiHH%d_LSP%d_veto.txt ControlRegionB1TChiHH%d_LSP%d_veto.txt ControlRegionCTChiHH%d_LSP%d_veto.txt  ControlRegionDTChiHH%d_LSP%d_veto.txt > TChiHH%d_LSP%d_1BoostedH_veto.txt "%(hino,LSP,hino,LSP,hino,LSP,hino,LSP,hino,LSP))
#os.system("combineCards.py TChiHH%d_LSP%d_1BoostedH_veto.txt TChiHH%d_LSP%d_2BoostedH_veto.txt  > TChiHH%d_LSP%d_BothBoostedH_veto.txt "%(hino,LSP,hino,LSP,hino,LSP))
# os.system("combine -M AsymptoticLimits -n TChiHH%d_LSP%d_BothBoostedH_veto TChiHH%d_LSP%d_BothBoostedH_veto.txt " %(hino,LSP,hino,LSP))
# os.system("combine -M AsymptoticLimits -n TChiHH%d_LSP%d_2BoostedH_veto TChiHH%d_LSP%d_2BoostedH_veto.txt " %(hino,LSP,hino,LSP))
# os.system("combine -M AsymptoticLimits -n TChiHH%d_LSP%d_1BoostedH_veto TChiHH%d_LSP%d_1BoostedH_veto.txt " %(hino,LSP,hino,LSP))


#ResolvedOnly
# os.system("combineCards.py SignalRegionRes4bT5HH%d.txt ControlRegionResBT5HH%d.txt ControlRegionResCT5HH%d.txt  ControlRegionResDT5HH%d.txt>T5HH%dRes4b.txt "%(hino,hino,hino,hino,hino))
# os.system("combineCards.py SignalRegionRes3bT5HH%d.txt ControlRegionResB3bT5HH%d.txt >T5HH%dRes3b.txt "%(hino,hino,hino))
# os.system("combineCards.py T5HH%dRes4b.txt  T5HH%dRes3b.txt > T5HH%dRes.txt "%(hino,hino,hino))
# os.system("combine -M AsymptoticLimits -n TChiHH%dResOnly datacard-TChiHH_mChi-%d_mLSP-0_Tune_2016_resolved.txt " %(hino,hino))

#Combo
#os.system("combineCards.py SignalRegionRes4bT5HH%d.txt ControlRegionResBT5HH%d.txt ControlRegionResCT5HH%d.txt  ControlRegionResDT5HH%d.txt>T5HH%dRes4b.txt "%(hino,hino,hino,hino,hino))
#os.system("combineCards.py SignalRegionRes3bT5HH%d.txt ControlRegionResB3bT5HH%d.txt >T5HH%dRes3b.txt "%(hino,hino,hino))
#os.system("combineCards.py SignalRegionT5HH%d.txt ControlRegionBT5HH%d.txt ControlRegionCT5HH%d.txt  ControlRegionDT5HH%d.txt > T5HH%d_2BoostedH_veto.txt "%(hino,hino,hino,hino,hino))

#Using Jaebak's datacards + BoostedOnlyVeto
#os.system("combineCards.py TChiHH%d_LSP%d_BothBoostedH_veto.txt datacard-TChiHH_mChi-%d_mLSP-0_Tune_2016_resolved.txt > TChiHH%d_LSP%d_Combo.txt "%(hino,LSP,hino,hino,LSP))
#os.system("combine -M AsymptoticLimits -n TChiHH%d_LSP%d_Combo TChiHH%d_LSP%d_Combo.txt " %(hino,LSP,hino,LSP))
