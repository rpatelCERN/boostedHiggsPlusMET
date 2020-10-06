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
# TotalEvents=float(sys.argv[2])
print("\n\n\n-----------------------------------------------------")
print("Looping over mass=%d now" %(hino))
print("-----------------------------------------------------\n")
f2=TFile("ALPHABET_sig_MET3.root", "READ")
# f2=TFile("ALPHABET_V18_1Dsignal_2BoostedH.root", "READ")

SignalA2=f2.Get("MET_doubletagSR_T5qqqqZH_%d_1" %(hino));
SignalB2=f2.Get("MET_doubletagSB_T5qqqqZH_%d_1" %(hino));
#SignalB2=f2.Get("MET_doubletagSB2_T5qqqqZH_%d_1" %(hino));
SignalA1=f2.Get("MET_tagSR_T5qqqqZH_%d_1" %(hino));
SignalB1=f2.Get("MET_tagSB_T5qqqqZH_%d_1" %(hino));
#SignalB1=f2.Get("MET_tagSB2_T5qqqqZH_%d_1" %(hino));
SignalC=f2.Get("MET_antitagSR_T5qqqqZH_%d_1" %(hino));
#SignalD=f2.Get("MET_antitagSB2_T5qqqqZH_%d_1" %(hino));
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

#vetoed
SignalA2Veto=fVeto.Get("MET_doubletagSR_T5HH%d" %(hino));
SignalB2Veto=fVeto.Get("MET_doubletagSB_T5HH%d" %(hino));
SignalA1Veto=fVeto.Get("MET_tagSR_T5HH%d" %(hino));
SignalB1Veto=fVeto.Get("MET_tagSB_T5HH%d" %(hino));
SignalCVeto=fVeto.Get("MET_antitagSR_T5HH%d" %(hino));
SignalDVeto=fVeto.Get("MET_antitagSB_T5HH%d" %(hino));

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
# kappas2H = 1.0 #V18, original SB, removing kappas
# kappas1H = 1.0 #V18, original SB, removing kappas
# kappas2H = 1.51 #TEST
# kappas1H = 1.29 #TEST
kappas2H = 1.00 #V18
kappas1H = 1.00 #V18
#kappas2H = 1.66 #V18, SB
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

#Predictions=f.Get("BackgroundPred")

#BoostedOnly
fcard=open("SignalRegionTemplateMergeBkg.txt", 'r');
fcardout=open("SRMerge_T5HH%d.txt" %hino , 'w')
#fcardout=open("SRMerge_T5HH%d_SB2.txt" %hino , 'w')
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
fcardout=open("SR1Merge_T5HH%d.txt" %(hino) , 'w')
#fcardout=open("SR1Merge_T5HH%d_SB2.txt" %(hino) , 'w')
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

fcardout=open("CRBMerge_T5HH%d.txt" %hino , 'w')
#fcardout=open("CRBMerge_T5HH%d_SB2.txt" %hino , 'w')
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

fcardout=open("CRB1Merge_T5HH%d.txt" %(hino) , 'w')
#fcardout=open("CRB1Merge_T5HH%d_SB2.txt" %(hino) , 'w')
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

fcardout=open("CRCMerge_T5HH%d.txt" %hino , 'w')
#fcardout=open("CRCMerge_T5HH%d_SB2.txt" %hino , 'w')
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


fcardout=open("CRDMerge_T5HH%d.txt" %hino , 'w')
#fcardout=open("CRDMerge_T5HH%d_SB2.txt" %hino , 'w')
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
fcardout=open("SRMerge_T5HH%d_veto.txt" %hino , 'w')
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
fcardout=open("SR1Merge_T5HH%d_veto.txt" %(hino) , 'w')
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

fcardout=open("CRBMerge_T5HH%d_veto.txt" %hino , 'w')
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


fcardout=open("CRB1Merge_T5HH%d_veto.txt" %(hino) , 'w')
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

fcardout=open("CRCMerge_T5HH%d_veto.txt" %hino , 'w')
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


fcardout=open("CRDMerge_T5HH%d_veto.txt" %hino , 'w')
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

'''

#BoostedOnly
#os.system("combineCards.py SRMerge_T5HH%d_SB2.txt CRBMerge_T5HH%d_SB2.txt CRCMerge_T5HH%d_SB2.txt  CRDMerge_T5HH%d_SB2.txt > T5HH%d_2BoostedHMerge_SB2.txt "%(hino,hino,hino,hino,hino))
#os.system("combineCards.py SR1Merge_T5HH%d_SB2.txt CRB1Merge_T5HH%d_SB2.txt CRCMerge_T5HH%d_SB2.txt  CRDMerge_T5HH%d_SB2.txt > T5HH%d_1BoostedHMerge_SB2.txt "%(hino,hino,hino,hino,hino))
#os.system("combineCards.py T5HH%d_1BoostedHMerge_SB2.txt T5HH%d_2BoostedHMerge_SB2.txt  > T5HH%d_BothBoostedHMerge_SB2.txt "%(hino,hino,hino))
#os.system("combine -M AsymptoticLimits -n T5HH%d_BothBoostedHMerge_SB2 T5HH%d_BothBoostedHMerge_SB2.txt " %(hino,hino))


os.system("combineCards.py SRMerge_T5HH%d.txt CRBMerge_T5HH%d.txt CRCMerge_T5HH%d.txt  CRDMerge_T5HH%d.txt > T5HH%d_2H_MET3.txt "%(hino,hino,hino,hino,hino))
os.system("combineCards.py SR1Merge_T5HH%d.txt CRB1Merge_T5HH%d.txt CRCMerge_T5HH%d.txt  CRDMerge_T5HH%d.txt > T5HH%d_1H_MET3.txt "%(hino,hino,hino,hino,hino))
os.system("combineCards.py T5HH%d_1H_MET3.txt T5HH%d_2H_MET3.txt  > T5HH%d_MET3.txt "%(hino,hino,hino))
os.system("combine -M AsymptoticLimits -n T5HH%d_MET3 T5HH%d_MET3.txt " %(hino,hino))

# os.system("combine -M AsymptoticLimits -n T5HH%d_2BoostedH T5HH%d_2BoostedH.txt " %(hino,hino))
# os.system("combine -M AsymptoticLimits -n T5HH%d_1BoostedH T5HH%d_1BoostedH.txt " %(hino,hino))

#ResolvedOnly
# os.system("combineCards.py SignalRegionRes4bHighT5HH%d.txt ControlRegionResB4bHighT5HH%d.txt ControlRegionResCHighT5HH%d.txt  ControlRegionResDHighT5HH%d.txt>T5HH%dRes4bHigh.txt "%(hino,hino,hino,hino,hino))
# os.system("combineCards.py SignalRegionRes4bLowT5HH%d.txt ControlRegionResB4bLowT5HH%d.txt ControlRegionResCLowT5HH%d.txt  ControlRegionResDLowT5HH%d.txt>T5HH%dRes4bLow.txt "%(hino,hino,hino,hino,hino))
# os.system("combineCards.py SignalRegionRes3bHighT5HH%d.txt ControlRegionResB3bHighT5HH%d.txt ControlRegionResCHighT5HH%d.txt  ControlRegionResDHighT5HH%d.txt>T5HH%dRes3bHigh.txt "%(hino,hino,hino,hino,hino))
# os.system("combineCards.py SignalRegionRes3bLowT5HH%d.txt ControlRegionResB3bLowT5HH%d.txt ControlRegionResCLowT5HH%d.txt  ControlRegionResDLowT5HH%d.txt>T5HH%dRes3bLow.txt "%(hino,hino,hino,hino,hino))
# os.system("combineCards.py T5HH%dRes4bHigh.txt T5HH%dRes4bLow.txt  T5HH%dRes3bHigh.txt T5HH%dRes3bLow.txt > T5HH%dRes.txt "%(hino,hino,hino,hino,hino))
# os.system("combine -M AsymptoticLimits -n T5HH%dRes T5HH%dRes.txt " %(hino,hino))
#
# #Combo
# #os.system("combineCards.py SignalRegionRes4bT5HH%d.txt ControlRegionResBT5HH%d.txt ControlRegionResCT5HH%d.txt  ControlRegionResDT5HH%d.txt>T5HH%dRes4b.txt "%(hino,hino,hino,hino,hino))
# #os.system("combineCards.py SignalRegionRes3bT5HH%d.txt ControlRegionResB3bT5HH%d.txt >T5HH%dRes3b.txt "%(hino,hino,hino))
# os.system("combineCards.py SignalRegionT5HH%d_veto.txt ControlRegionBT5HH%d_veto.txt ControlRegionCT5HH%d_veto.txt  ControlRegionDT5HH%d_veto.txt > T5HH%d_2BoostedH_veto.txt "%(hino,hino,hino,hino,hino))
# os.system("combineCards.py T5HH%d_2BoostedH_veto.txt T5HH%dRes.txt > T5HH%dCombo.txt "%(hino,hino,hino))
# os.system("combine -M AsymptoticLimits -n T5HH%dCombo T5HH%dCombo.txt " %(hino,hino))
