#include "tdrstyle.C"
#include "CMS_lumi.cc"
void MakeABCDRes(){
    setTDRStyle();
    TFile*f=new TFile("ALPHABETResMC2016_V18.root", "READ");

    TH1F*h_A_ZJets = (TH1F*)f->Get("MET_fourbSR_ZJets"); TH1F*h_B_ZJets = (TH1F*)f->Get("MET_fourbSB_ZJets");
    TH1F*h_A1_ZJets = (TH1F*)f->Get("MET_threebSR_ZJets"); TH1F*h_B1_ZJets = (TH1F*)f->Get("MET_threebSB_ZJets");
    TH1F*h_C_ZJets= (TH1F*)f->Get("MET_twobSR_ZJets"); TH1F*h_D_ZJets = (TH1F*)f->Get("MET_twobSB_ZJets");
    TH1F*h_A_WJets = (TH1F*)f->Get("MET_fourbSR_WJets"); TH1F*h_B_WJets = (TH1F*)f->Get("MET_fourbSB_WJets");
    TH1F*h_A1_WJets = (TH1F*)f->Get("MET_threebSR_WJets"); TH1F*h_B1_WJets = (TH1F*)f->Get("MET_threebSB_WJets");
    TH1F*h_C_WJets= (TH1F*)f->Get("MET_twobSR_WJets"); TH1F*h_D_WJets = (TH1F*)f->Get("MET_twobSB_WJets");
    TH1F*h_A_QCD = (TH1F*)f->Get("MET_fourbSR_QCD"); TH1F*h_B_QCD = (TH1F*)f->Get("MET_fourbSB_QCD");
    TH1F*h_A1_QCD = (TH1F*)f->Get("MET_threebSR_QCD"); TH1F*h_B1_QCD = (TH1F*)f->Get("MET_threebSB_QCD");
    TH1F*h_C_QCD= (TH1F*)f->Get("MET_twobSR_QCD"); TH1F*h_D_QCD = (TH1F*)f->Get("MET_twobSB_QCD");
    TH1F*h_A_TT = (TH1F*)f->Get("MET_fourbSR_TT"); TH1F*h_B_TT = (TH1F*)f->Get("MET_fourbSB_TT");
    TH1F*h_A1_TT = (TH1F*)f->Get("MET_threebSR_TT"); TH1F*h_B1_TT = (TH1F*)f->Get("MET_threebSB_TT");
    TH1F*h_C_TT= (TH1F*)f->Get("MET_twobSR_TT"); TH1F*h_D_TT = (TH1F*)f->Get("MET_twobSB_TT");

    TH1F*TotalBkg=(TH1F*)h_A_ZJets->Clone("TotalBkg");
    TH1F*TotalBkg3b=(TH1F*)h_A1_ZJets->Clone("TotalBkg3b");
    double scale=137000.0/35862.824;
    // double scale=1.;
    TH1D*ControlB=(TH1D*)h_B_ZJets->Clone("ControlB");
    TH1D*ControlB1=(TH1D*)h_B1_ZJets->Clone("ControlB1");
    TH1D*ControlC=(TH1D*)h_C_ZJets->Clone("ControlC");
    TH1D*ControlD=(TH1D*)h_D_ZJets->Clone("ControlD");

    TotalBkg->Add(h_A_WJets); TotalBkg3b->Add(h_A1_WJets);
    ControlB->Add(h_B_WJets); ControlB1->Add(h_B1_WJets);
    ControlC->Add(h_C_WJets); ControlD->Add(h_D_WJets);

    TotalBkg->Add(h_A_QCD); TotalBkg3b->Add(h_A1_QCD);
    ControlB->Add(h_B_QCD); ControlB1->Add(h_B1_QCD);
    ControlC->Add(h_C_QCD); ControlD->Add(h_D_QCD);

    TotalBkg->Add(h_A_TT); TotalBkg3b->Add(h_A1_TT);
    ControlB->Add(h_B_TT); ControlB1->Add(h_B1_TT);
    ControlC->Add(h_C_TT); ControlD->Add(h_D_TT);

    TotalBkg->Scale(scale); TotalBkg3b->Scale(scale);
    ControlB->Scale(scale); ControlB1->Scale(scale);
    ControlC->Scale(scale); ControlD->Scale(scale);

    TH1D*Pred=(TH1D*)ControlB->Clone("Pred");
    Pred->Multiply(ControlC);
    Pred->Divide(ControlD);

    TH1D*Pred3b=(TH1D*)ControlB1->Clone("Pred3b");
    Pred3b->Multiply(ControlC);
    Pred3b->Divide(ControlD);

    THStack * BackgroundStack = new THStack("hs","");
    h_A_ZJets->Scale(scale); h_A_WJets->Scale(scale);
    h_A_QCD->Scale(scale); h_A_TT->Scale(scale);
    h_A_QCD->SetFillColor(kGray); h_A_TT->SetFillColor(kCyan);
    h_A_ZJets->SetFillColor(kGreen+2); h_A_WJets->SetFillColor(kBlue);
    BackgroundStack->Add(h_A_ZJets); BackgroundStack->Add(h_A_WJets);
    BackgroundStack->Add(h_A_TT); BackgroundStack->Add(h_A_QCD);

    THStack * BackgroundStack3b = new THStack("hs3b","");
    h_A1_ZJets->Scale(scale); h_A1_WJets->Scale(scale);
    h_A1_QCD->Scale(scale); h_A1_TT->Scale(scale);
    h_A1_QCD->SetFillColor(kGray); h_A1_TT->SetFillColor(kCyan);
    h_A1_ZJets->SetFillColor(kGreen+2); h_A1_WJets->SetFillColor(kBlue);
    BackgroundStack3b->Add(h_A1_ZJets); BackgroundStack3b->Add(h_A1_WJets);
    BackgroundStack3b->Add(h_A1_TT); BackgroundStack3b->Add(h_A1_QCD);

    double W = 800; double H = 600;
    double T = 0.08*H; double B = 0.12*H;
    double L = 0.12*W; double R = 0.06*W;
    double up_height = 0.8;
    double dw_correction = 1.18;
    double font_size_dw  = 0.1;
    double dw_height = (1. - up_height) * dw_correction;
    double dw_height_offset = 0.02;

    TCanvas * can_h = new TCanvas("can_h","Closure", 50, 50, W, H);
    can_h->SetFillColor(0); can_h->SetBorderMode(0);
    can_h->SetFrameFillStyle(0); can_h->SetFrameBorderMode(0);
    can_h->SetLeftMargin(L/W); can_h->SetRightMargin(R/W);
    can_h->SetTopMargin(T/H); can_h->SetBottomMargin(B/H);
    can_h->SetTickx(0); can_h->SetTicky(0);
    TPad *pad1 = new TPad("pad1", "top pad" , 0.0, 0.25, 1.0, 1.0);
    TPad *pad2 = new TPad("pad2", "bottom pad", 0.0, 0.0, 1.0, 0.25);
    pad1->SetTickx(0); pad1->SetTicky(0);
    pad1->SetPad(0., 1 - up_height, 1., 1.00);
    pad1->SetFrameFillColor(0); pad1->SetFillColor(0);
    pad1->SetTopMargin(0.08);
    pad1->SetLeftMargin(0.12);
    pad1->SetRightMargin(0.06);
    // pad1->SetLogy(logy)
    pad1->Draw();

    can_h->Update();
    TotalBkg->GetYaxis()->SetRangeUser(0,90);
    Pred->GetYaxis()->SetRangeUser(0,45);
    Pred->SetMarkerStyle(kFullCircle);
    h_B_ZJets->SetFillColor(kBlack);
    Pred->SetTitle(";MET [GeV];Events 137/fb");
    Pred->GetYaxis()->SetTitleOffset(0.7);
    Pred->Draw("pe");
    BackgroundStack->Draw("histsame");
    Pred->Draw("pesame");

    TLegend* legend = new TLegend(0.56,0.68,0.90,0.90);
    legend->SetNColumns(2);
    legend->AddEntry(h_A_QCD, "QCD", "f");
    legend->AddEntry(h_A_TT, "TT", "f");
    legend->AddEntry(h_A_WJets, "WJets", "f");
    legend->AddEntry(h_A_ZJets, "ZJets", "f");
    legend->AddEntry(Pred, "Pred: B*C/D", "pe");

    legend->SetBorderSize(0);
    legend->Draw();
    writeExtraText = true;
    extraText  = "  Simulation";
    lumi_sqrtS = "137 fb^{-1}(13 TeV)";
    CMS_lumi(can_h,0,0);

    TFile*output = new TFile("OutputYields_res.root","RECREATE");
    TotalBkg->Write("SRYields_4b");
    TotalBkg3b->Write("SRYields_3b");
    Pred->Write("BkgPred_4b");
    Pred3b->Write("BkgPred_3b");
    ControlB->Write("ControlB");
    ControlB1->Write("ControlB1"); //This is B for 3b case
    ControlC->Write("ControlC");
    ControlD->Write("ControlD");
    can_h->Write("Fourb");

    // TotalBkg3b->GetYaxis()->SetRangeUser(0,20);
    // Pred3b->SetMarkerStyle(kFullCircle);
    // h_B1_ZJets->SetFillColor(kBlack);
    // Pred3b->SetTitle(";MET [GeV];Events 137/fb");
    // Pred3b->GetYaxis()->SetTitleOffset(0.7);
    // Pred3b->Draw("pe");
    // BackgroundStack3b->Draw("histsame");
    // Pred3b->Draw("pesame");
    // legend->Draw();
    // writeExtraText = true;
    // extraText  = "  Simulation";
    // lumi_sqrtS = "137 fb^{-1}(13 TeV)";
    // CMS_lumi(can_h,0,0);
    // can_h->Write("Threeb");
}
