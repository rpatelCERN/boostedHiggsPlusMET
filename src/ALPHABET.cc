#include "TString.h"
#include "TChain.h"
#include "TH1F.h"
#include "TROOT.h"
#include "THStack.h"
#include "TPad.h"

#include <vector>
#include <map>
#include <iostream>
#include <assert.h>

#include "plotterUtils.cc"
#include "skimSamples.cc"
#include "definitions.cc"
#include "RA2bTree.cc"
#include "ALPHABET.h"
#include "TriggerEfficiencySextet.cc"

using namespace std;
using namespace alphabet;

int main(int argc, char** argv){

    int region(0);
    bool looseCuts(false);
    int MAX_EVENTS(99999999);

    if( argc >= 2 ){
        region = atoi(argv[1]);
        if( argc >= 3 )
            looseCuts = atoi(argv[2]);
        if( argc >= 4 )
            MAX_EVENTS = atoi(argv[3]);
    }else
        cout << "Running over the default (signal) region ... " << endl;

    gROOT->ProcessLine(".L tdrstyle.C");
    gROOT->ProcessLine("setTDRStyle()");

    skimSamples* skims_;
    if( region == 0 )
        skims_ = new skimSamples(skimSamples::kSignal);
    else if( region == 1 )
        skims_ = new skimSamples(skimSamples::kSLm);
    else if( region == 2 )
        skims_ = new skimSamples(skimSamples::kSLe);
    else if( region == 3 )
        skims_ = new skimSamples(skimSamples::kLowDphi);
    else
        assert(1);

    typedef bool(*cuts)(RA2bTree*);
    vector<cuts> baselineCuts;

    if( looseCuts ){
        baselineCuts.push_back(*FiltersCut<RA2bTree>);
        if( region == 3 ){
            baselineCuts.push_back(*lowDPhiCuts<RA2bTree>);
        }else{
            baselineCuts.push_back(*DeltaPhiCuts<RA2bTree>);
        }
        if( region == 1 ){
            baselineCuts.push_back(*singleMuCut<RA2bTree>);
        }
        if( region == 2 ){
            baselineCuts.push_back(*singleEleCut<RA2bTree>);
        }
        baselineCuts.push_back(*METHTlooseCut<RA2bTree>);
        baselineCuts.push_back(*AK8MultCut<RA2bTree>);
    }
    else{
        if( region == 0 ){
            baselineCuts.push_back(*baselineCut<RA2bTree>);
        }
        else if( region == 1){
            baselineCuts.push_back(*singleMuBaselineCut<RA2bTree>);
        }
        else if( region == 2){
            baselineCuts.push_back(*singleEleBaselineCut<RA2bTree>);
        }
        else if( region == 3){
            baselineCuts.push_back(*lowDphiBaselineCut<RA2bTree>);
        }
        else
            assert(1);
    }

    skimSamples skims = *skims_;

    typedef plot<RA2bTree> plot;

    double mJbins[4]={50.,85.,135.,250.};
    double METbins[6]={150.,200.,300.,500.,700.,1000.};

    vector<vector<plot> > plots;

    /*
    for( int i = 0 ; i < numMETbins ; i++ ) {
        TString tag="_";
        tag+=lowestMET+i*binWidth;
        vector<plot> plotsTemp;
        plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_tagSR"+tag,"m_{J} [GeV]",3,mJbins));
        plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_tagSB"+tag,"m_{J} [GeV]",3,mJbins));

        plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_antitagSR"+tag,"m_{J} [GeV]",3,mJbins));
        plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_antitagSB"+tag,"m_{J} [GeV]",3,mJbins));

        plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_doubletagSR"+tag,"m_{J} [GeV]",3,mJbins));
        plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_doubletagSB"+tag,"m_{J} [GeV]",3,mJbins));

        plots.push_back(plotsTemp);
    }
    */
    //vector<plot> tempPlots;
    // plot MET_Plot(*fillMET<RA2bTree>,"MET","m_{J} [GeV]",3,100,700);
    plot MET_Plot(*fillMET<RA2bTree>,"MET","m_{J} [GeV]",5,METbins);
    plot J1_pt(*fillLeadingJetPt<RA2bTree>,"J1_pt","p_{T,J} [GeV]",50,300.,1300.);
    plot J2_pt(*fillSubLeadingJetPt<RA2bTree>,"J2_pt","p_{T,J} [GeV]",50,300.,1300.);
    plot J1_m(*fillLeadingJetMass<RA2bTree>,"J1_m","m_{J} [GeV]",50,50.,250.);
    plot J2_m(*fillSubLeadingJetMass<RA2bTree>,"J2_m","m_{J} [GeV]",50,50.,250.);
    plot ClosestMass(*fillClosestJetMass<RA2bTree>,"ClosestMass","m_{J} [GeV]",3,mJbins);
    plot FarthestMass(*fillFarthestJetMass<RA2bTree>,"FarthestMass","m_{J} [GeV]",3,mJbins);

    plot DeltaRMaxPlot(*fillDeltaRMax<RA2bTree>,"DeltaRMax","deltaR_max [cm]",50,0,5.0);

    plot J1_doubleB_lowMET(*fillLeadingBBtagLow<RA2bTree>,"J1_DoubleB_lowMET","double-B discriminator, MET < 300 GeV",50,-1.0,1.0);
    plot J2_doubleB_lowMET(*fillSubLeadingBBtagLow<RA2bTree>,"J2_DoubleB_lowMET","double-B discriminator, MET < 300 GeV",50,-1.0,1.0);

    plot J1_doubleB_highMET(*fillLeadingBBtagHigh<RA2bTree>,"J1_DoubleB_highMET","double-B discriminator, MET > 300 GeV",50,-1.0,1.0);
    plot J2_doubleB_highMET(*fillSubLeadingBBtagHigh<RA2bTree>,"J2_DoubleB_highMET","double-B discriminator, MET > 300 GeV",50,-1.0,1.0);

    vector<TH1F*> EmilysDeltaRMaxPlotsDouble;
    vector<TH1F*> EmilysDeltaRMaxPlotsSingle;

    vector<plot> doubletagSRPlots;
    doubletagSRPlots.push_back(plot(DeltaRMaxPlot));
    doubletagSRPlots.push_back(plot(J1_pt));
    doubletagSRPlots.push_back(plot(J2_pt));
    doubletagSRPlots.push_back(plot(J1_m));
    doubletagSRPlots.push_back(plot(J2_m));
    doubletagSRPlots.push_back(plot(ClosestMass));
    doubletagSRPlots.push_back(plot(FarthestMass));

    doubletagSRPlots.push_back(plot(J1_doubleB_lowMET));
    doubletagSRPlots.push_back(plot(J2_doubleB_lowMET));
    doubletagSRPlots.push_back(plot(J1_doubleB_highMET));
    doubletagSRPlots.push_back(plot(J2_doubleB_highMET));
    doubletagSRPlots.push_back(plot(MET_Plot));



    vector<plot> doubletagSBPlots;
    doubletagSBPlots.push_back(plot(DeltaRMaxPlot));
    doubletagSBPlots.push_back(plot(J1_pt));
    doubletagSBPlots.push_back(plot(J2_pt));
    doubletagSBPlots.push_back(plot(J1_m));
    doubletagSBPlots.push_back(plot(J2_m));
    doubletagSBPlots.push_back(plot(ClosestMass));
    doubletagSBPlots.push_back(plot(FarthestMass));

    doubletagSBPlots.push_back(plot(J1_doubleB_lowMET));
    doubletagSBPlots.push_back(plot(J2_doubleB_lowMET));
    doubletagSBPlots.push_back(plot(J1_doubleB_highMET));
    doubletagSBPlots.push_back(plot(J2_doubleB_highMET));
    doubletagSBPlots.push_back(plot(MET_Plot));


    vector<plot> tagSRPlots;
    tagSRPlots.push_back(plot(DeltaRMaxPlot));
    tagSRPlots.push_back(plot(J1_pt));
    tagSRPlots.push_back(plot(J2_pt));
    tagSRPlots.push_back(plot(J1_m));
    tagSRPlots.push_back(plot(J2_m));
    tagSRPlots.push_back(plot(ClosestMass));
    tagSRPlots.push_back(plot(FarthestMass));

    tagSRPlots.push_back(plot(J1_doubleB_lowMET));
    tagSRPlots.push_back(plot(J2_doubleB_lowMET));
    tagSRPlots.push_back(plot(J1_doubleB_highMET));
    tagSRPlots.push_back(plot(J2_doubleB_highMET));
    tagSRPlots.push_back(plot(MET_Plot));


    vector<plot> tagSBPlots;
    tagSBPlots.push_back(plot(DeltaRMaxPlot));
    tagSBPlots.push_back(plot(J1_pt));
    tagSBPlots.push_back(plot(J2_pt));
    tagSBPlots.push_back(plot(J1_m));
    tagSBPlots.push_back(plot(J2_m));
    tagSBPlots.push_back(plot(ClosestMass));
    tagSBPlots.push_back(plot(FarthestMass));

    tagSBPlots.push_back(plot(J1_doubleB_lowMET));
    tagSBPlots.push_back(plot(J2_doubleB_lowMET));
    tagSBPlots.push_back(plot(J1_doubleB_highMET));
    tagSBPlots.push_back(plot(J2_doubleB_highMET));
    tagSBPlots.push_back(plot(MET_Plot));


    vector<plot> antitagSRPlots;
    antitagSRPlots.push_back(plot(MET_Plot));
    antitagSRPlots.push_back(plot(DeltaRMaxPlot));
    antitagSRPlots.push_back(plot(J1_pt));
    antitagSRPlots.push_back(plot(J2_pt));
    antitagSRPlots.push_back(plot(J1_m));
    antitagSRPlots.push_back(plot(J2_m));
    antitagSRPlots.push_back(plot(ClosestMass));
    antitagSRPlots.push_back(plot(FarthestMass));

    antitagSRPlots.push_back(plot(J1_doubleB_lowMET));
    antitagSRPlots.push_back(plot(J2_doubleB_lowMET));
    antitagSRPlots.push_back(plot(J1_doubleB_highMET));
    antitagSRPlots.push_back(plot(J2_doubleB_highMET));

    vector<plot> antitagSBPlots;
    antitagSBPlots.push_back(plot(MET_Plot));
    antitagSBPlots.push_back(plot(DeltaRMaxPlot));
    antitagSBPlots.push_back(plot(J1_pt));
    antitagSBPlots.push_back(plot(J2_pt));
    antitagSBPlots.push_back(plot(J1_m));
    antitagSBPlots.push_back(plot(J2_m));
    antitagSBPlots.push_back(plot(ClosestMass));
    antitagSBPlots.push_back(plot(FarthestMass));

    antitagSBPlots.push_back(plot(J1_doubleB_lowMET));
    antitagSBPlots.push_back(plot(J2_doubleB_lowMET));
    antitagSBPlots.push_back(plot(J1_doubleB_highMET));
    antitagSBPlots.push_back(plot(J2_doubleB_highMET));



    // background MC samples - 0 lepton regions
    std::cout<<"Ntuples size: "<< skims.ntuples.size()<<std::endl;
    for( int iSample = 0 ; iSample < skims.ntuples.size() ; iSample++){

        RA2bTree* ntuple = skims.ntuples[iSample];

        // for( int iBin = 0 ; iBin < numMETbins ; iBin++){
        //     for( int iPlot = 0 ; iPlot < plots[iBin].size() ; iPlot++){
        //         plots[iBin][iPlot].addNtuple(ntuple,skims.sampleName[iSample]);
        //         plots[iBin][iPlot].setFillColor(ntuple,skims.fillColor[iSample]);
        //     }
        // }
        TString histoNameDoubleTag = "DeltaRMax_doubletagSR_"+skims.sampleName[iSample];
        TString histoNameSingleTag = "DeltaRMax_tagSR_"+skims.sampleName[iSample];
        TH1F * DeltaRMaxPlotDouble = new TH1F(histoNameDoubleTag,"deltaR_max [cm]",50, 0, 5.0);
        TH1F * DeltaRMaxPlotSingle = new TH1F(histoNameSingleTag,"deltaR_max [cm]",50, 0, 5.0);

        EmilysDeltaRMaxPlotsDouble.push_back(DeltaRMaxPlotDouble);
        EmilysDeltaRMaxPlotsSingle.push_back(DeltaRMaxPlotSingle);

        for( int i = 0 ; i < doubletagSRPlots.size() ; i++ ){
            doubletagSRPlots[i].addNtuple(ntuple,"doubletagSR_"+skims.sampleName[iSample]);
            doubletagSRPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
        }
        for( int i = 0 ; i < doubletagSBPlots.size() ; i++ ){
            doubletagSBPlots[i].addNtuple(ntuple,"doubletagSB_"+skims.sampleName[iSample]);
            doubletagSBPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
        }
        for( int i = 0 ; i < tagSRPlots.size() ; i++ ){
            tagSRPlots[i].addNtuple(ntuple,"tagSR_"+skims.sampleName[iSample]);
            tagSRPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
        }
        for( int i = 0 ; i < tagSBPlots.size() ; i++ ){
            tagSBPlots[i].addNtuple(ntuple,"tagSB_"+skims.sampleName[iSample]);
            tagSBPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
        }
        for( int i = 0 ; i < antitagSRPlots.size() ; i++ ){
            antitagSRPlots[i].addNtuple(ntuple,"antitagSR_"+skims.sampleName[iSample]);
            antitagSRPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
        }
        for( int i = 0 ; i < antitagSBPlots.size() ; i++ ){
            antitagSBPlots[i].addNtuple(ntuple,"antitagSB_"+skims.sampleName[iSample]);
            antitagSBPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
        }

        int numEvents = ntuple->fChain->GetEntries();
        TH1F* readThis = 0;
        ntuple->fChain->GetFile()->GetObject("nEventProc", readThis);
        int nEvents_signal = readThis->GetBinContent(1);
        // int nEvents_signal = 1;

        // file->GetObject("hpx", readThis);
        // std::cout<<"Supposedly num of preprocessed events: "<<nEvents_signal<<std::endl;
        ntupleBranchStatus<RA2bTree>(ntuple);
        int bin = -1;
        double weight=0.;
        float trigWeight=1.0;
        bool passBaseline;
        double jetMass1,jetMass2;
        TString filename;
        for( int iEvt = 0 ; iEvt < min(MAX_EVENTS,numEvents) ; iEvt++ ){
            ntuple->GetEntry(iEvt);
            if( iEvt % 10000 == 0 ) cout << skims.sampleName[iSample] << ": " << iEvt << "/" << min(MAX_EVENTS,numEvents) << endl;

            if(region==0){
                std::vector<double> EfficiencyCenterUpDown = Eff_MetMhtSextetReal_CenterUpDown(ntuple->HT, ntuple->MHT, ntuple->NJets);
                trigWeight=EfficiencyCenterUpDown[0];
            }
            else if( region == 1 ){
                trigWeight=singleMuonTrigWeights(ntuple);
            }else if( region == 2 ){
                trigWeight=singleElectronTrigWeights(ntuple);
            }else if( region == 3 ){
                trigWeight=lowDphiTrigWeights(ntuple);
            }

            // passBaseline=true;
            // for( auto baselineCut : baselineCuts ){
            //     passBaseline&=baselineCut(ntuple);
            // }

            // if( ! passBaseline ) continue;
            if( ! boostedBaselineCut(ntuple) ) continue;


            filename = ntuple->fChain->GetFile()->GetName();

            if( ( filename.Contains("SingleLept") || filename.Contains("DiLept") ) && ntuple->madHT>600. )continue;
            bin = -1;
            // weight = ntuple->Weight*lumi*trigWeight*customPUweights(ntuple); //this might be a problem
            weight = ntuple->Weight*lumi;

            if ( filename.Contains("TChiHH") ) {
              // std::cout<<"weight: "<<ntuple->Weight<<", lumi: "<<lumi<<", TMath: "<< TMath::Power(0.57,2)<< ", nEvents: "<<nEvents_signal;

              weight= ntuple->Weight*lumi*TMath::Power(0.57,2)/nEvents_signal;
              // cout<<",    New weight: "<< weight<<std::endl;
              // std::cout<<"Supposedly num of preprocessed events: "<<nEvents_signal<<std::endl;
            }
            // weight = ntuple->Weight*lumi; //this might be a problem
            // std::cout<<"weight: "<<weight<<std::endl;

            //if( skims.sampleName[iSample] == "TT" ){
            //    weight *= ISRweights(ntuple);
            //}

            // //Commenting this out for now
            // for( int iBin = 0 ; iBin < numMETbins ; iBin++ ){
            //     if( ntuple->MET > lowestMET ){
            //         if( ntuple->MET > numMETbins*(binWidth-1)+lowestMET )
            //             bin = numMETbins-1;
            //         else
            //             bin = int((ntuple->MET-lowestMET)/binWidth);
            //     }
            // }
            // if( bin < 0 ) continue;


            if( doubletagSRCut(ntuple) ){
              // std::cout<<"Does it make it here?\n";
              EmilysDeltaRMaxPlotsDouble.back()->Fill( deltaRMax(ntuple) );
                // plots[bin][4].fill(ntuple,weight);
                for( int i = 0 ; i < doubletagSRPlots.size() ; i++ ) doubletagSRPlots[i].fill (ntuple,weight);
            }
            else if( doubletagSBCut( ntuple ) ){
                // plots[bin][5].fill(ntuple,weight);
                for( int i = 0 ; i < doubletagSBPlots.size() ; i++ ) doubletagSBPlots[i].fill (ntuple,weight);
            }
            else if( tagSRCut( ntuple ) ){
                // plots[bin][0].fill(ntuple,weight);
                EmilysDeltaRMaxPlotsSingle.back()->Fill( deltaRMax(ntuple) );
                for( int i = 0 ; i < tagSRPlots.size() ; i++ ) tagSRPlots[i].fill (ntuple,weight);
            }

            else if( tagSBCut( ntuple ) ){
                // plots[bin][1].fill(ntuple,weight);
                for( int i = 0 ; i < tagSBPlots.size() ; i++ ) tagSBPlots[i].fill (ntuple,weight);
            }
            else if( antitagSRCut( ntuple ) ){
                // plots[bin][2].fill(ntuple,weight);
                for( int i = 0 ; i < antitagSRPlots.size() ; i++ ) antitagSRPlots[i].fill (ntuple,weight);
            }
            else if( antitagSBCut( ntuple ) ){
                // plots[bin][3].fill(ntuple,weight);
                for( int i = 0 ; i < antitagSBPlots.size() ; i++ ) antitagSBPlots[i].fill (ntuple,weight);
            } // end if-else-if block for tagging regions


        } // end event loop

    } // end sample loop


    /*
    // data
    RA2bTree* ntuple = skims.dataNtuple;

    for( int iBin = 0 ; iBin < numMETbins ; iBin++){
        for( int iPlot = 0 ; iPlot < plots[iBin].size() ; iPlot++){
            plots[iBin][iPlot].addDataNtuple(ntuple,"data");
        }
    }
    for( int i = 0 ; i < doubletagSRPlots.size() ; i++ ){
        doubletagSRPlots[i].addDataNtuple(ntuple,"doubletagSR_data");
    }
    for( int i = 0 ; i < doubletagSBPlots.size() ; i++ ){
        doubletagSBPlots[i].addDataNtuple(ntuple,"doubletagSB_data");
    }
    for( int i = 0 ; i < tagSRPlots.size() ; i++ ){
        tagSRPlots[i].addDataNtuple(ntuple,"tagSR_data");
    }
    for( int i = 0 ; i < tagSBPlots.size() ; i++ ){
        tagSBPlots[i].addDataNtuple(ntuple,"tagSB_data");
    }
    for( int i = 0 ; i < antitagSRPlots.size() ; i++ ){
        antitagSRPlots[i].addDataNtuple(ntuple,"antitagSR_data");
    }
    for( int i = 0 ; i < antitagSBPlots.size() ; i++ ){
        antitagSBPlots[i].addDataNtuple(ntuple,"antitagSB_data");
    }

    int numEvents = ntuple->fChain->GetEntries();
    ntupleBranchStatus<RA2bTree>(ntuple);
    bool passBaseline;
    double jetMass1,jetMass2;
    for( int iEvt = 0 ; iEvt < min(MAX_EVENTS,numEvents) ; iEvt++ ){
        ntuple->GetEntry(iEvt);
        if( iEvt % 100000 == 0 ) cout << "data: " << iEvt << "/" << min(MAX_EVENTS,numEvents) << endl;

        passBaseline=true;
        for( auto baselineCut : baselineCuts ){
            passBaseline&=baselineCut(ntuple);
        }
        if( ! passBaseline ) continue;

        if( region == 0 ){
            if( !signalTriggerCut(ntuple) ) continue;
        }else if( region == 1){
            if( !singleMuTriggerCut(ntuple) ) continue;
        }else if( region == 2){
            if( !singleEleTriggerCut(ntuple) ) continue;
        }else if( region == 3 ){
            if( !lowDphiTriggerCut(ntuple) ) continue;
        }

        int bin = -1;
        for( int iBin = 0 ; iBin < numMETbins ; iBin++ ){
            if( ntuple->MET > lowestMET ){
                if( ntuple->MET > lowestMET+binWidth*(numMETbins-1) )
                    bin = numMETbins-1;
                else
                    bin = int((ntuple->MET-lowestMET)/binWidth);
            }
        }
        if( bin < 0 ) continue;

        UInt_t runNum = ntuple->RunNum;
        UInt_t lumiBlk = ntuple->LumiBlockNum;
        ULong64_t evtNum  = ntuple->EvtNum;

        if( doubletagSRCut( ntuple ) && ntuple->MET>300 && ntuple->HT>600){
          std::cout<<"Double tag: Run Num: "<<runNum<<", lumi block: "<<lumiBlk<<", event num: "<<evtNum;
          std::cout<<", Lead jet BB: "<< fillLeadingBBtag(ntuple)<<", sublead BB: "<<fillSubLeadingBBtag(ntuple)<<std::endl;
            if( region != 0 ){
                plots[bin][4].fillData(ntuple);
                for( int i = 0 ; i < doubletagSRPlots.size() ; i++ )
                    doubletagSRPlots[i].fillData(ntuple);
            }
        }else if( doubletagSBCut( ntuple )  && ntuple->MET>300 && ntuple->HT>600){
            plots[bin][5].fillData(ntuple);
            for( int i = 0 ; i < doubletagSBPlots.size() ; i++ )
                doubletagSBPlots[i].fillData(ntuple);
        }else if( tagSRCut( ntuple )  && ntuple->MET>300 && ntuple->HT>600){
          std::cout<<"Single tag: Run Num: "<<runNum<<", lumi block: "<<lumiBlk<<", event num: "<<evtNum;
          std::cout<<", Lead jet BB: "<< fillLeadingBBtag(ntuple)<<", sublead BB: "<<fillSubLeadingBBtag(ntuple)<<std::endl;

            if( region != 0 ){
                plots[bin][0].fillData(ntuple);
                for( int i = 0 ; i < tagSRPlots.size() ; i++ )
                    tagSRPlots[i].fillData(ntuple);
            }
        }else if( tagSBCut( ntuple ) ){
            plots[bin][1].fillData(ntuple);
            for( int i = 0 ; i < tagSBPlots.size() ; i++ )
                tagSBPlots[i].fillData(ntuple);
        }else if( antitagSRCut( ntuple ) ){
                plots[bin][2].fillData(ntuple);
                for( int i = 0 ; i < antitagSRPlots.size() ; i++ )
                    antitagSRPlots[i].fillData(ntuple);
        }else if( antitagSBCut( ntuple ) ){
            plots[bin][3].fillData(ntuple);
            for( int i = 0 ; i < antitagSBPlots.size() ; i++ )
                antitagSBPlots[i].fillData(ntuple);
        }// end if-else-if block for tagging regions
    }// end event loop
    */

    TFile* outputFile;
    // TString regionName;
    // TString cutName="";
    // if( looseCuts )
    //     cutName="_looseCuts";
    // if( region == 0 )
    //     regionName="";
    // if( region == 1 )
    //     regionName="_singleMu";
    // if( region == 2 )
    //     regionName="_singleEle";
    // if( region == 3 )
    //     regionName="_lowDphi";
        // outputFile = new TFile("ALPHABEThistos"+cutName+regionName+".root","RECREATE");
        // outputFile = new TFile("ALPHABEThistos_TChiHH_V12.root","RECREATE");
        // outputFile = new TFile("ALPHABEThistos_V12_METonly_signalOnly.root","RECREATE");
        // outputFile = new TFile("ALPHABEThistos_V12_signal.root","RECREATE");
        outputFile = new TFile("ALPHABEThistos_V12_lowB.root","RECREATE");
        // outputFile = new TFile("ALPHABEThistos_V12_medB.root","RECREATE");

    // for( int iBin = 0 ; iBin < numMETbins; iBin++){
    //     for( int iPlot = 0 ; iPlot < plots[iBin].size() ; iPlot++){
    //         outputFile->cd();
    //         plots[iBin][iPlot].buildSum();
    //         plots[iBin][iPlot].Write();
    //         plots[iBin][iPlot].sum->Write();
    //
    //     }
    // }

    std::cout<<"size: "<<doubletagSRPlots.size()<<std::endl;
    for( int i = 0 ; i < doubletagSRPlots.size() ; i++ ){
        outputFile->cd();
        // doubletagSRPlots[i].buildSum("doubletagSR");
        doubletagSRPlots[i].Write();
        // doubletagSRPlots[i].sum->Write();
    }
    for( int i = 0 ; i < EmilysDeltaRMaxPlotsDouble.size() ; i++ ){
        outputFile->cd();
        EmilysDeltaRMaxPlotsDouble[i]->Write();
    }

    for( int i = 0 ; i < doubletagSBPlots.size() ; i++ ){
        outputFile->cd();
        // doubletagSBPlots[i].buildSum("doubletagSB");
        doubletagSBPlots[i].Write();
        // doubletagSBPlots[i].sum->Write();
    }
    for( int i = 0 ; i < tagSRPlots.size() ; i++ ){
        outputFile->cd();
        // tagSRPlots[i].buildSum("tagSR");
        tagSRPlots[i].Write();
        // tagSRPlots[i].sum->Write();
    }

    for( int i = 0 ; i < EmilysDeltaRMaxPlotsSingle.size() ; i++ ){
        outputFile->cd();
        EmilysDeltaRMaxPlotsSingle[i]->Write();
    }

    for( int i = 0 ; i < tagSBPlots.size() ; i++ ){
        outputFile->cd();
        // tagSBPlots[i].buildSum("tagSB");
        tagSBPlots[i].Write();
        // tagSBPlots[i].sum->Write();
    }
    for( int i = 0 ; i < antitagSRPlots.size() ; i++ ){
        outputFile->cd();
        // antitagSRPlots[i].buildSum("antitagSR");
        antitagSRPlots[i].Write();
        // antitagSRPlots[i].sum->Write();
    }
    for( int i = 0 ; i < antitagSBPlots.size() ; i++ ){
        outputFile->cd();
        // antitagSBPlots[i].buildSum("antitagSB");
        antitagSBPlots[i].Write();
        // antitagSBPlots[i].sum->Write();
    }

    outputFile->Close();

}
