#!/bin/bash

DIRECTORY=ResToBoost/
# DIRECTORY=BoostToRes/

cd $DIRECTORY

cd Boost/
rm -f Bkg.root TChiHH.root T5HH.root
hadd -f Bkg.root QCD_MET.root TT_MET.root WJets_MET.root ZJets_MET.root
hadd -f TChiHH.root TChiHH*_MET.root
hadd -f T5HH.root T5HH*_MET.root

cd ../DiAK8/
rm -f Bkg.root TChiHH.root T5HH.root
hadd -f Bkg.root QCD_MET.root TT_MET.root WJets_MET.root ZJets_MET.root
hadd -f TChiHH.root TChiHH*_MET.root
hadd -f T5HH.root T5HH*_MET.root

cd ../Res4b_lowR/
rm -f Bkg.root TChiHH.root T5HH.root
hadd -f Bkg.root QCD_MET.root TT_MET.root WJets_MET.root ZJets_MET.root
hadd -f TChiHH.root TChiHH*_MET.root
hadd -f T5HH.root T5HH*_MET.root

cd ../Res4b_highR/
rm -f Bkg.root TChiHH.root T5HH.root
hadd -f Bkg.root QCD_MET.root TT_MET.root WJets_MET.root ZJets_MET.root
hadd -f TChiHH.root TChiHH*_MET.root
hadd -f T5HH.root T5HH*_MET.root

cd ../Res3b_lowR/
rm -f Bkg.root TChiHH.root T5HH.root
hadd -f Bkg.root QCD_MET.root TT_MET.root WJets_MET.root ZJets_MET.root
hadd -f TChiHH.root TChiHH*_MET.root
hadd -f T5HH.root T5HH*_MET.root

cd ../Res3b_highR/
rm -f Bkg.root TChiHH.root T5HH.root
hadd -f Bkg.root QCD_MET.root TT_MET.root WJets_MET.root ZJets_MET.root
hadd -f TChiHH.root TChiHH*_MET.root
hadd -f T5HH.root T5HH*_MET.root

cd ../..
