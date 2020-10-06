import os
import glob
import math
import array
import sys
import time
import ROOT
from array import array
from ROOT import gROOT

# define the palette for z axis
NRGBs = 5
NCont = 255
stops = array("d",[0.00, 0.34, 0.61, 0.84, 1.00])
red= array("d",[0.50, 0.50, 1.00, 1.00, 1.00])
green = array("d",[ 0.50, 1.00, 1.00, 0.60, 0.50])
blue = array("d",[1.00, 1.00, 0.50, 0.40, 0.50])
#FI = ROOT.TColor.CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont)
#ROOT.gStyle.SetNumberContours(NCont)
#MyPalette = []
#for i in range(0,100): MyPalette.append(FI+i)
ROOT.gStyle.SetPalette(112);


def MakeExpectedSignificancePlot(vmx, vmy, vobs):
	xparticle = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{2}}}";
  	yparticle = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}";
  	title = ";m_{"+xparticle+"} [GeV];m_{"+yparticle+"} [GeV]; Expected Significance";
	g=ROOT.TGraph2D()
	g.SetName(title)

	for m in range(len(vmx)):
		g.SetPoint(g.GetN(), vmx[m], vmy[m], vobs[m])
	g.SetNpx(128)
	g.SetNpy(160)

  	the_max = 0
  	for i in range(0,g.GetN()):
		z = g.GetZ()[i]
		if (z>the_max):
			the_max = z

  	if (the_max > 6):
	  	the_max = 6

  	g.SetMinimum(0.)
  	g.SetMaximum(the_max)
  	g.GetHistogram().SetTitle(title)

  	c = ROOT.TCanvas("can","can",1800,1200)
  	ltitle = ROOT.TLatex(c.GetLeftMargin(), 1.-0.5*c.GetTopMargin(),"#font[62]{CMS}#scale[0.76]{#font[52]{ Supplementary}}");
  	rtitle = ROOT.TLatex(1.-c.GetRightMargin(), 1.-0.5*c.GetTopMargin(),"#scale[0.8]{137 fb^{-1} (13 TeV)}");
  	ltitle.SetNDC();
  	rtitle.SetNDC();
  	ltitle.SetTextAlign(12);
  	rtitle.SetTextAlign(32);
  	model = ROOT.TLatex(0.64,0.88,"pp#rightarrow #tilde{#chi}^{0}_{1} #tilde{#chi}^{0}_{1} #rightarrow HH #tilde{G} #tilde{G}");

  	g.Draw("colz");

  	for z in range(0,the_max):
		if (z==5):
			style=1
			width=2
		else:
			style=2
			width=1
		# DrawContours(g, 1, style, width, 0, z)
		DrawContours(g)

  	model.Draw("same");
  	ltitle.Draw("same");
  	rtitle.Draw("same");

  	c.Print("TChiHH2D_sigexp.pdf");
  	c.Print("TChiHH2D_sigexp.root");
  	return g;

def DrawContours(g2):
	#OK I am oversimplifying this, then I'll add shit back in
	graph = ROOT.TGraph()
	style=2
	h_2D = g2.GetHistogram();
	l = g2.GetContourList(1.0);

  	max_points = -1;
  	g=[];
  	for i in range(0,len(l)):
		g.append(l[i]);
		n_points = g[i].GetN();
		if (n_points > max_points):
			old_graph = graph.Clone();
			if (i>0):
				ReverseGraph(g[i]);
			graph = joinGraphs(old_graph, g[i]);
			graph.RemovePoint(graph.GetN()-1);
			max_points = n_points;
		g[i].Draw("L same");
	return graph;


def FixGraph(graph):
	glu_lsp = 225;
  	np = graph.GetN();
	iniglu=0;inilsp=0; endglu=0; endlsp=0;
  	graph.GetPoint(0, ROOT.Double(iniglu), ROOT.Double(inilsp));
  	graph.GetPoint(np-1, ROOT.Double(endglu), ROOT.Double(endlsp));

  	#Reversing graph if printed towards decreasing mgluino
  	if (inilsp < endlsp):
                ReverseGraph(graph);
                endglu = iniglu;
                endlsp = inilsp;

  	# Adding a point so that it goes down to mLSP = 0, but not for WZ,SOS
  	if (endlsp<30):
                graph.SetPoint(graph.GetN(), endglu, 0);
                np+=1
 	ReverseGraph(graph);
  	# Adding a point at mLSP = 0, and removing points beyond the diagonal
  	for point in range(0,np):
  		# for(int point(0); point < np; point++){
		mglu=0; mlsp=0;
                graph.GetPoint(point, ROOT.Double(mglu), ROOT.Double(mlsp));
                if (mlsp > mglu-glu_lsp-5):
                        while (point<=graph.GetN() and point!=0):
                                graph.RemovePoint(graph.GetN()-1);
                                np-=1
			break


def ReverseGraph(graph):
	np = graph.GetN();
  	mglus=[]; mlsps=[];
  	for point in range(0,np):
  		# for(int point(np-1); point >= 0; point--){
		mglu=0;mlsp=0;
                graph.GetPoint(point, ROOT.Double(mglu), ROOT.Double(mlsp));
                mglus.append(mglu);
                mlsps.append(mlsp);
        # for(int point(0); point < np; point++)
  	for point in range(0,np):
                graph.SetPoint(point, mglus[point], mlsps[point]);


def joinGraphs(graph1, graph2):
	graph = ROOT.TGraph();
  	for point in range(0,graph1.GetN()):
  		# for(int point(0); point < graph1->GetN(); point++) {
                mglu = 0; mlsp=0;
                graph1.GetPoint(point, ROOT.Double(mglu), ROOT.Double(mlsp));
                graph.SetPoint(graph.GetN(), mglu, mlsp);
  	# Points in graph1
  	for point in range(0,graph2.GetN()):
  		# for(int point(0); point < graph2->GetN(); point++) {
                mglu = 0; mlsp=0;
                graph2.GetPoint(point, ROOT.Double(mglu), ROOT.Double(mlsp));
                graph.SetPoint(graph.GetN(), mglu, mlsp);
  	# Points in graph1
  	graph1.GetPoint(0, ROOT.Double(mglu), ROOT.Double(mlsp));
  	graph.SetPoint(graph.GetN(), mglu, mlsp);
  	gname = graph1.GetName(); gname += graph2.GetName();
  	graph.SetName(gname);
  	return graph;

def higgsino2DCrossSection(hig_mass): #this is in pb
	if (hig_mass==127):
		xsec = 1.44725; xsec_unc = 0.0395277;
	elif (hig_mass==150):
		xsec = 0.71514; xsec_unc = 0.0421496;
	elif (hig_mass==175):
		xsec = 0.419059; xsec_unc = 0.0453279;
	elif (hig_mass==200):
		xsec = 0.244213; xsec_unc = 0.047925;
	elif (hig_mass==225):
		xsec = 0.156286; xsec_unc = 0.0502876;
	elif (hig_mass==250):
		xsec = 0.104252; xsec_unc = 0.0526169;
	elif (hig_mass==275):
		xsec = 0.0719125; xsec_unc = 0.0549666;
	elif (hig_mass==300):
		xsec = 0.0509994; xsec_unc = 0.0572762;
	elif (hig_mass==325):
		xsec = 0.0369715; xsec_unc = 0.0590317;
	elif (hig_mass==350):
		xsec = 0.0273286; xsec_unc = 0.0607766;
	elif (hig_mass==375):
		xsec = 0.0205429; xsec_unc = 0.0625031;
	elif (hig_mass==400):
		xsec = 0.0156691; xsec_unc = 0.0642085;
	elif (hig_mass==425):
		xsec = 0.0120965; xsec_unc = 0.0657801;
	elif (hig_mass==450):
		xsec = 0.00944017; xsec_unc = 0.0674544;
	elif (hig_mass==475):
		xsec = 0.00743587; xsec_unc = 0.0686033;
	elif (hig_mass==500):
		xsec = 0.00590757; xsec_unc = 0.0699909;
	elif (hig_mass==525):
		xsec = 0.00469101; xsec_unc = 0.0713704;
	elif (hig_mass==550):
		xsec = 0.0038167; xsec_unc = 0.0722834;
	elif (hig_mass==575):
		xsec = 0.003073; xsec_unc = 0.0739957;
	elif (hig_mass==600):
		xsec = 0.00253015; xsec_unc = 0.0754291;
	elif (hig_mass==625):
		xsec = 0.00206136; xsec_unc = 0.0763466;
	elif (hig_mass==650):
		xsec = 0.00171418; xsec_unc = 0.0775695;
	elif (hig_mass==675):
		xsec = 0.00140934; xsec_unc = 0.0783375;
	elif (hig_mass==700):
		xsec = 0.00118113; xsec_unc = 0.0796388;
	elif (hig_mass==725):
		xsec = 0.000979349; xsec_unc = 0.0809883;
	elif (hig_mass==750):
		xsec = 0.000826366; xsec_unc = 0.081879;
	elif (hig_mass==775):
		xsec = 0.000690208; xsec_unc = 0.0842049;
	elif (hig_mass==800):
		xsec = 0.000586211; xsec_unc = 0.0862527;
	elif (hig_mass==825):
		xsec = 0.00049277; xsec_unc = 0.0864444;
	elif (hig_mass==850):
		xsec = 0.000420556; xsec_unc = 0.085742;
	elif (hig_mass==875):
		xsec = 0.000358734; xsec_unc = 0.0889174;
	elif (hig_mass==900):
		xsec = 0.000305935; xsec_unc = 0.0912439;
	elif (hig_mass==925):
		xsec = 0.000260948; xsec_unc = 0.091372;
	elif (hig_mass==950):
		xsec = 0.00022285; xsec_unc = 0.0919538;
	elif (hig_mass==975):
		xsec = 0.000189681; xsec_unc = 0.0938108;
	elif (hig_mass==1000):
		xsec = 0.00016428; xsec_unc = 0.0954285;
	elif (hig_mass==1025):
		xsec = 0.000142206; xsec_unc = 0.0957231;
	elif (hig_mass==1050):
		xsec = 0.000120971; xsec_unc = 0.0968997;
	elif (hig_mass==1075):
		xsec = 0.000105301; xsec_unc = 0.0979041;
	elif (hig_mass==1100):
		xsec = 9.12469e-05; xsec_unc = 0.0964142;
	elif (hig_mass==1125):
		xsec = 7.9765e-05; xsec_unc = 0.099902;
	elif (hig_mass==1150):
		xsec = 6.78234e-05; xsec_unc = 0.101061;
	elif (hig_mass==1175):
		xsec = 5.9016e-05; xsec_unc = 0.102051;
	elif (hig_mass==1200):
		xsec = 5.16263e-05; xsec_unc = 0.102499;
	elif (hig_mass==1225):
		xsec = 4.5147e-05; xsec_unc = 0.10403;
	elif (hig_mass==1250):
		xsec = 3.88343e-05; xsec_unc = 0.105206;
	elif (hig_mass==1275):
		xsec = 3.41304e-05; xsec_unc = 0.10619;
	elif (hig_mass==1300):
		xsec = 2.99353e-05; xsec_unc = 0.10783;
	elif (hig_mass==1325):
		xsec = 2.63637e-05; xsec_unc = 0.108024;
	elif (hig_mass==1350):
		xsec = 2.26779e-05; xsec_unc = 0.109016;
	elif (hig_mass==1375):
		xsec = 1.99318e-05; xsec_unc = 0.109822;
	elif (hig_mass==1400):
		xsec = 1.75031e-05; xsec_unc = 0.111631;
	elif (hig_mass==1425):
		xsec = 1.53974e-05; xsec_unc = 0.111417;
	elif (hig_mass==1450):
		xsec = 1.3245e-05; xsec_unc = 0.112313;
	elif (hig_mass==1475):
		xsec = 1.16416e-05; xsec_unc = 0.113058;
	else:
		xsec = 0; xsec_unc = 0;
	return xsec, xsec_unc
	#return xsec


def MakeObservedSignificancePlot(vmx,vmy,vobs):
	#SetupSignedColors();
	xparticle = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{2}}}";
	yparticle = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}";
	title = ";m_{"+xparticle+"} [GeV];m_{"+yparticle+"} [GeV];Observed Significance";
	g=ROOT.TGraph2D()
	g.SetName(title)
	for m in range(len(vobs)):
		g.SetPoint(g.GetN(), vmx[m], vmy[m], vobs[m])
	g.SetNpx(128)
	g.SetNpy(160)
	the_max = 3
	g.SetMinimum(-the_max);
	g.SetMaximum(the_max);
	g.SetNpx(128)
	g.SetNpy(160)

	g.GetHistogram().SetTitle(title);
	g.GetHistogram().SetTickLength(0, "Z");
	c = ROOT.TCanvas("can","can", 800, 900);
	c.cd();

	ltitle = ROOT.TLatex(c.GetLeftMargin(), 1.-0.5*c.GetTopMargin(),"#font[62]{CMS}#scale[0.76]{#font[52]{ Supplementary}}");
	rtitle = ROOT.TLatex(1.-c.GetRightMargin(), 1.-0.5*c.GetTopMargin(),"#scale[0.8]{137 fb^{-1} (13 TeV)}");
	ltitle.SetNDC();
	rtitle.SetNDC();
	ltitle.SetTextAlign(12);
	rtitle.SetTextAlign(32);
	model = ROOT.TLatex(0.64,0.88,"pp#rightarrow #tilde{#chi}^{0}_{1} #tilde{#chi}^{0}_{1} #rightarrow HH #tilde{G} #tilde{G}");

	g.Draw("colz");

	lines=[];
	max_list = [0.0,0.5,1.0,1.5,2.0,2.5,3.0]
	for z in max_list:
		style = int(round(min(5, 2*abs(z)+1)));
		width = max(1, int(round(3-abs(z))));
		# DrawContours(g, 1, style, width, 0, z);
		DrawContours(g);
		x1 = 1-c.GetRightMargin()+0.0047;
		x2 = 1-c.GetRightMargin()+0.05;
		ybot = c.GetBottomMargin();
		ytop = 1-c.GetTopMargin();
		zpos = ybot+(ytop-ybot)*(the_max+z)/(2*the_max);
		this_line = ROOT.TLine(x1, zpos, x2, zpos)
		lines.append(this_line);
		lines[-1].SetLineColor(1);
		lines[-1].SetLineStyle(style);
		lines[-1].SetLineWidth(width);
		lines[-1].SetNDC(True);
		if (z != 0):
			DrawContours(g);
			# DrawContours(g, 1, style, width, 0, -z);
			zneg = ybot+(ytop-ybot)*(the_max-z)/(2*the_max);
			this_line = ROOT.TLine(x1, zneg, x2, zneg);
			lines.append(this_line);
			lines[-1].SetLineColor(1);
			lines[-1].SetLineStyle(style);
			lines[-1].SetLineWidth(width);
			lines[-1].SetNDC(True);


	for l in lines:
		l.Draw("same");

	model.Draw("same");
	ltitle.Draw("same");
	rtitle.Draw("same");

	c.Print("TChiHH2D_sigobs.pdf");
	c.Print("TChiHH2D_sigobs.root");
	return g;


def MakeLimitPlot(vmx,vmy,vlim,vobs,vobs1sig,vexp,vmx2,vmy2):
	#SetupColors();
	xparticle = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{2}}}";
	yparticle = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}";
	title = ";m_{"+xparticle+"} [GeV];m_{"+yparticle+"} [GeV];95% CL upper limit on cross section [pb]";

	glim=ROOT.TGraph2D()
	glim.SetName(title)
	for m in range(len(vobs)):
		glim.SetPoint(glim.GetN(), vmx[m], vmy[m], vobs[m])
	# for m in range(len(vlim)):
	# 	glim.SetPoint(glim.GetN(), vmx[m], vmy[m], vlim[m])

	gobs=ROOT.TGraph2D()
	gobs.SetName("Observed Limit")
	for m in range(len(vobs)):
		gobs.SetPoint(gobs.GetN(), vmx[m], vmy[m], vobs[m])

	gobs1sig=ROOT.TGraph2D()
	gobs1sig.SetName("Observed 1#sigma Limit")
	for m in range(len(vobs1sig)):
		gobs1sig.SetPoint(gobs1sig.GetN(), vmx2[m], vmy2[m], vobs1sig[m])

	gexp=ROOT.TGraph2D()
	gexp.SetName("Expected Limit")
	for m in range(len(vexp)):
		gexp.SetPoint(gexp.GetN(), vmx[m], vmy[m], vexp[m])

	# gup=ROOT.TGraph2D()
	# gup.SetName("Expected +1#sigma Limit")
	# for m in range(len(vup)):
	# 	gup.SetPoint(gup.GetN(), vmx[m], vmy[m], vup[m])
	#
	# gdown=ROOT.TGraph2D()
	# gdown.SetName("Expected -1#sigma Limit")
	# for m in range(len(vdown)):
	# 	gdown.SetPoint(gdown.GetN(), vmx[m], vmy[m], vdown[m])

	glim.SetMinimum(0.00001);
	glim.SetMaximum(2);
	glim.SetNpx(128)
	glim.SetNpy(160)
	glim.SetTitle(title);
	num_smooth_=0
	# cup = DrawContours(gup, 2, 2, 5, num_smooth_);
	# cdown = DrawContours(gdown, 2, 2, 5, num_smooth_);
	# cexp = DrawContours(gexp, 2, 1, 5, num_smooth_, 1.);
	# cobs1sig = DrawContours(gobs1sig, 1, 2, 5, num_smooth_, 1.);
	# cobs = DrawContours(gobs, 1, 1, 5, num_smooth_, 1.);

	# cexp = gexp.GetContourList(1.0);
	# cobs1sig = gobs1sig.GetContourList(1.0);
	# cobs = gobs.GetContourList(1.0);
	# cexp.Draw("l same")

	l = ROOT.TLegend(ROOT.gStyle.GetPadLeftMargin(), 1.-2.*ROOT.gStyle.GetPadTopMargin(),1.-ROOT.gStyle.GetPadRightMargin(), 1.-ROOT.gStyle.GetPadTopMargin());
	l.SetNColumns(2);
	l.SetTextSize(0.05);
	l.SetBorderSize(0);
	l.AddEntry(cexp, "Expected", "l");
	l.AddEntry(cobs, "Observed", "l");

	c = ROOT.TCanvas("","",1200, 950);
	c.cd();
	c.SetRightMargin(0.15);
	ltitle = ROOT.TLatex(c.GetLeftMargin(), 1.-0.5*c.GetTopMargin(),"#font[62]{CMS}#scale[0.76]{#font[52]{ Supplementary}}");
	rtitle = ROOT.TLatex(1.-c.GetRightMargin(), 1.-0.5*c.GetTopMargin(),"#scale[0.8]{137 fb^{-1} (13 TeV)}");
	ltitle.SetNDC();
	rtitle.SetNDC();
	ltitle.SetTextAlign(12);
	rtitle.SetTextAlign(32);

	c.SetTopMargin(2.*c.GetTopMargin());
	c.SetLeftMargin(0.12);
	c.SetLogz();
	glim.Draw("colz");
	#ObsLim=glim.GetContourList(1.0);

	l.Draw("same");
	ltitle.Draw("same");
	rtitle.Draw("same");

	filebase = "TChiHH2D_limit_scan";
	if (num_smooth_>0):
		filebase += "_smooth";
		filebase += num_smooth_;

	#gPad.Update();
	glim.GetZaxis().SetTitleOffset(1.3);
	glim.GetXaxis().SetRangeUser(0,1200);
	c.Print(filebase+".pdf");
	file = ROOT.TFile(filebase+".root", "recreate");
	glim.GetHistogram().Write("TChiHH2D_ObservedExcludedXsec");
	cobs.Write("TChiHH2D_ObservedLimit");
	cobs1sig.Write("TChiHH2D_ObservedLimit1Sig");
	cexp.Write("TChiHH2D_ExpectedLimit");
	#ObsLim.Write("ObsLim")
	# cup.Write(("TChiHH2D_ExpectedLimitUp");
	# cdown.Write("TChiHH2D_ExpectedLimitDown");
	# if (False): #Significances not saved together with limits as per recommendations
	# 	hsigobs.Write("ObservedSignificance");
	# 	hsigexp.Write("ExpectedSignificance");
	file.Close();
	print("\nSaved limit curves in " + filebase + ".root\n")

def setStyleCOLZ(histo,canv):
	# set z axis
	histo.GetZaxis().SetLabelFont(42)
	histo.GetZaxis().SetTitleFont(42)
	histo.GetZaxis().SetLabelSize(0.035)
	histo.GetZaxis().SetTitleSize(0.035)
	histo.SetMinimum(9.0E-05)
	histo.SetMaximum(10.0)

	canv.cd()
	#ROOT.gStyle.SetPalette(100, MyPalette);
	ROOT.gStyle.SetPalette(112);
	histo.Draw("colz")

	ROOT.gPad.Update()
	palette = histo.GetListOfFunctions().FindObject("palette")
	palette.SetX1NDC(1.-0.18)
	palette.SetY1NDC(0.14)
	palette.SetX2NDC(1.-0.13)
	palette.SetY2NDC(1.-0.08)
	palette.SetLabelFont(42)
	palette.SetLabelSize(0.035)
	print "I did things"
	return histo, canv


def setStyle(histo,canv):
	# canvas style
	ROOT.gStyle.SetOptStat(0)
	ROOT.gStyle.SetOptTitle(0)

	canv.SetLogz()
	canv.SetTickx(1)
	canv.SetTicky(1)

	canv.SetRightMargin(0.19)
	canv.SetTopMargin(0.08)
	canv.SetLeftMargin(0.14)
	canv.SetBottomMargin(0.14)

	# set x axis
	xparticle = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{2}}}";
	yparticle = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}";
	title = ";m_{"+xparticle+"} [GeV];m_{"+yparticle+"} [GeV];95% CL upper limit on cross section [pb]";
	histo.GetXaxis().SetLabelFont(42)
	histo.GetXaxis().SetLabelSize(0.035)
	histo.GetXaxis().SetTitleFont(42)
	histo.GetXaxis().SetTitleSize(0.05)
	histo.GetXaxis().SetTitleOffset(1.2)
	histo.GetXaxis().SetTitle(xparticle)
	#self.emptyHisto.GetXaxis().CenterTitle(True)

	# set y axis
	histo.GetYaxis().SetLabelFont(42)
	histo.GetYaxis().SetLabelSize(0.035)
	histo.GetYaxis().SetTitleFont(42)
	histo.GetYaxis().SetTitleSize(0.05)
	histo.GetYaxis().SetTitleOffset(1.3)
	histo.GetYaxis().SetTitle(yparticle)
	#self.emptyHisto.GetYaxis().CenterTitle(True)
	return histo, canv


def DrawCont(g_obs,histo_obs,histo_obs1sig,histo_exp,histo_exp1sig,histo_exp2sig,canv):
	# setStyle(histo_obs,canv)
	# setStyleCOLZ(histo_obs,canv)
	histo_obs.Draw()
	#DrawObsArea(g_obs,canv)
	DrawLines(histo_obs,histo_obs1sig,histo_exp,histo_exp1sig,histo_exp2sig,canv)
	# DrawText()
	# DrawLegend()


def DrawObsArea(g_obs,canv): #pass this observed histo
	# add points to observed to close area
	# this will disappea
	g_obs.SetPoint(g_obs.GetN(), 1300,-1300)
	g_obs.SetPoint(g_obs.GetN(), -1300,-1300)
	# observed

	trasparentColor = ROOT.gROOT.GetColor(15)
	trasparentColor.SetAlpha(0.5)
	g_obs.SetFillColor(15)
	g_obs.SetLineStyle(1)
	# DRAW AREAS
	g_obs.Draw("FSAME")


def DrawLines(g_obs,g_obs1sig,g_exp,g_exp1sig,g_exp2sig,canv):
	# observed
	histo_obs = DrawContours(g_obs)
	histo_obs.SetLineColor(1)
	histo_obs.SetLineStyle(1)
	histo_obs.SetLineWidth(4)
	# observed  1sigma
	#histo_obs1sig = DrawContours(g_obs1sig)
	#histo_obs1sig.SetLineColor(1)
	#histo_obs1sig.SetLineStyle(1)
	#histo_obs1sig.SetLineWidth(2)

	# expected 1sigma 2 sigma
	# histo_exp1sig = DrawContours(g_exp1sig)
	# histo_exp1sig.SetLineColor(2)
	# histo_exp1sig.SetLineStyle(7)
	# histo_exp1sig.SetLineWidth(2)
	#
	# histo_exp2sig = DrawContours(g_exp2sig)
	# histo_exp2sig.SetLineColor(2)
	# histo_exp2sig.SetLineStyle(2)
	# histo_exp2sig.SetLineWidth(2)

	# expected
	# histo_exp = DrawContours(g_exp)
	# histo_exp.SetLineColor(2)
	# histo_exp.SetLineStyle(7)
	# histo_exp.SetLineWidth(4)

	# DRAW LINES
	# histo_exp.Draw("l SAME")
	# histo_exp1sig.Draw("l SAME")
	# histo_exp2sig.Draw("l SAME")
	histo_obs.Draw("l SAME")
	#histo_obs1sig.Draw("l SAME")
