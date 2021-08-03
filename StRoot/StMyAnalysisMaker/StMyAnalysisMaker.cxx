
#include "StMyAnalysisMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StThreeVectorF.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include <TH1F.h>
#include "Stiostream.h"
//#include "fstream.h"
#include <math.h>
#include <TMath.h>

//StEmc
#include "StEmcClusterCollection.h"
#include "StEmcPoint.h"
#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StEmcUtil/others/emcDetectorName.h"
#include "StEmcADCtoEMaker/StBemcData.h"
#include "StEmcADCtoEMaker/StEmcADCtoEMaker.h"

#include "StEmcRawMaker/StBemcRaw.h"
#include "StEmcRawMaker/StBemcTables.h"
#include "StEmcRawMaker/StEmcRawMaker.h"
#include "StEmcRawMaker/defines.h"

#include "tables/St_emcStatus_Table.h"
#include "tables/St_smdStatus_Table.h"

#include "StMuDSTMaker/COMMON/StMuEmcCollection.h"
#include "StEmcCollection.h"
#include "StEmcCluster.h"
#include "StMuDSTMaker/COMMON/StMuEmcPoint.h"
#include "StEmcUtil/projection/StEmcPosition.h"
#include "StEmcUtil/filters/StEmcFilter.h"

#include "StEmcRawHit.h"
#include "StEmcModule.h"
#include "StEmcDetector.h"
#include "StEmcClusterCollection.h"
#include "StDaqLib/EMC/StEmcDecoder.h"

#include "StBTofHeader.h"
#include "TF1.h"
#include "TFile.h"
#include "StarClassLibrary/StThreeVectorD.hh"
#include "StarClassLibrary/StElectron.hh"
#include "StarClassLibrary/StPhysicalHelix.hh"
#include "tables/St_emcStatus_Table.h"
#include "PhysicalConstants.h"
#include "StMcEventTypes.hh"
#include "StMcEvent.hh"
#include "StMcEventMaker/StMcEventMaker.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "phys_constants.h"
#include "TMath.h"
#include "StGlobals.hh"
#include "StGlobalTrack.h"
#include "StPhysicalHelixD.hh"
#include "StDaqLib/TRG/trgStructures2005.h"
#include "StBTofUtil/tofPathLength.hh"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TTreeHelper.h"
#include "TBranch.h"


ClassImp(StMyAnalysisMaker)

//                              0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 10 , 11 , 12 , 13 , 14
Float_t StMyAnalysisMaker::p_low[14] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8};
Float_t StMyAnalysisMaker::p_up[14]  = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0};
//-----------------------------------------------------------------------------
StMyAnalysisMaker::StMyAnalysisMaker(const char* name, StPicoDstMaker *picoMaker, /*const*/ char* outName)
  : StMaker(name)
{
  mPicoDstMaker = picoMaker;
  mPicoDst = 0;
  mOutName = outName;
  mOutName = mOutName + ".root";
}

/*----------------------------------------------------------------------------------------------------------------------*/
StMyAnalysisMaker::~StMyAnalysisMaker()
{/* */}

Int_t StMyAnalysisMaker::Init()
{

  href_vz = new TH1F("h_ref_vz","refmult_vz",1000,0.,1000.);
  hvz_b = new TH1F("h_vz_b","vz_dis_b",1000,-150,150);
  hvzvpdvz_b =new TH2F("h_vz_vpd_b","vz_vs_vpd_b",1000,-150,150,1000,-150,150);
  hvr_b = new TH2F("h_vr_b","vy_vs_vx_b",1000,-10,10,1000,-10,10);
  htofvsref_b = new TH2F("htofvsref_b","ref_vs_tof",2000,0.0,2000.0,2000,0.0,2000.0);
  htofmatchvsref_b=new TH2F("htofmatchvsref_b","",1500,0.0,1500.0,1500,0.0,1500.0);

  href = new TH1F("h_ref","refmult_dis",1000,0.,1000.);
  hvz = new TH1F("h_vz","vz_dis",1000,-100,100);
  hvzvpdvz =new TH2F("h_vz_vpd","vz_vs_vpd",1000,-150,150,1000,-150,150);
  hvr = new TH2F("h_vr","vy_vs_vx",1000,-10,10,1000,-10,10);
  htofvsref = new TH2F("htofvsref","ref_vs_tof",2000,0.0,2000.0,2000,0.0,2000.0);
  htofmatchvsref=new TH2F("htofmatchvsref","",1500,0.0,1500.0,1500,0.0,1500.0);

  hbetavsp =new TH2F("hbetavsp","beta_vs_p",1000,0,6,1000,0,10);
  hmassvsp =new TH2F("hmassvsp","m2_vs_p*q",3000,-6.,6.,2000,-0.2,15.8);
  hdedxvsp =new TH2F("hdedxvsp","dedx_vs_p*q",4000,-6.,6.,2000,0.,50.);
  h_eta_phi =new TH2F("h_etaphi","eta_vs_phi",1000,-6.3,6.3,300,-1.5,1.5);
  h_eta_phi_before =new TH2F("h_etaphi_before","eta_vs_phi w/o cut",1000,-6.3,6.3,300,-1.8,1.8);

  h_counter = new TH1F("h_counter","event_counter",58,-8,50.);

  h_pt = new TH1F("h_pt","",1000,0.,10.);
  h_eta_b= new TH1F("h_eta_b","",1000,-2.,2.);
  h_eta= new TH1F("h_eta","",1000,-2.,2.);
  h_nhitfit=new TH1F("h_nhitfit","",80,-0.5,79.5);
  h_nhitmax=new TH1F("h_nhitmax","",80,-0.5,79.5);
  h_nhitratio=new TH1F("h_nhitratio","",1000,-0.5,1.5);
  h_dca = new TH1F("h_dca","",1000,0.,5.);
  h_phi = new TH1F("h_phi","",1000,-6.28,6.28);

  h_nsigmapip=new TH2F("h_nsigmapip","",200,-5.,5.,200,-5.,5.);
  h_nsigmapik=new TH2F("h_nsigmapik","",200,-5.,5.,200,-5.,5.);
  h_nsigmakp=new TH2F("h_nsigmakp","",200,-5.,5.,200,-5.,5.);
  TString name_pip, name_pik, name_kp;
  TString name_pip_p, name_pik_p, name_kp_p;
  for(int i=0;i<14;i++)
  {
    name_pip = Form("m_h_nsigmapip_%d",i);
    name_pik = Form("m_h_nsigmapik_%d",i);
    name_kp = Form("m_h_nsigmakp_%d",i);
    name_pip_p = Form("m_h_nsigmapip, %1.1f < p < %1.1f",StMyAnalysisMaker::p_low[i],StMyAnalysisMaker::p_up[i]);
    name_pik_p = Form("m_h_nsigmapik, %1.1f < p < %1.1f",StMyAnalysisMaker::p_low[i],StMyAnalysisMaker::p_up[i]);
    name_kp_p = Form("m_h_nsigmakp, %1.1f < p < %1.1f",StMyAnalysisMaker::p_low[i],StMyAnalysisMaker::p_up[i]);
    m_h_nsigmapip[i]=new TH2F(name_pip.Data(),name_pip_p.Data(),200,-5.,5.,200,-5.,5.);
    m_h_nsigmapik[i]=new TH2F(name_pik.Data(),name_pik_p.Data(),200,-5.,5.,200,-5.,5.);
    m_h_nsigmakp[i]=new TH2F(name_kp.Data(),name_kp_p.Data(),200,-5.,5.,200,-5.,5.);
  }
  h_nsigmae=new TH1F("hnsigmae","",500,-10.,10.);

  h_test_nsigmapi = new TH1D("h_test_nsigmapi", "", 1000, -50.0, 50.0);
  h_test_nsigmapr = new TH1D("h_test_nsigmapr", "", 1000, -50.0, 50.0);

  h_runidvstofmult_b = new TProfile("runidvstofmult_b", "", 90000, 22031041, 22121041,"");
  h_runidvsrefmult_b = new TProfile("runidvsrefmult_b", "", 90000, 22031041, 22121041,"");

  h_runidvstofmult = new TProfile("runidvstofmult", "", 90000, 22031041, 22121041,"");
  h_runidvsrefmult = new TProfile("runidvsrefmult", "", 90000, 22031041, 22121041,"");

  return kStOK;
}

//-----------------------------------------------------------------------------
Int_t StMyAnalysisMaker::Finish() {

  TFile *f = new TFile(mOutName.Data(),"RECREATE");
  f->cd();
  href_vz->Write(), hvz_b->Write(), hvzvpdvz_b->Write(), hvr_b->Write(), htofvsref_b->Write(),htofmatchvsref_b->Write();
  href->Write(), hvz->Write(), hvzvpdvz->Write(), hvr->Write(), htofvsref->Write(), htofmatchvsref->Write();
  h_runidvstofmult_b->Write();
  h_runidvsrefmult_b->Write();
  h_runidvstofmult->Write();
  h_runidvsrefmult->Write();
  hbetavsp->Write();
  hmassvsp->Write();
  hdedxvsp->Write();
  h_eta_phi_before->Write();
  h_eta_phi->Write();

  h_pt->Write(), h_eta_b->Write(), h_eta->Write(), h_nhitfit->Write(),h_nhitmax->Write(), h_nhitratio->Write(), h_dca->Write(),h_phi->Write();
  h_nsigmapip->Write(), h_nsigmapik->Write(), h_nsigmakp->Write(), h_nsigmae->Write();

  h_test_nsigmapi->Write();
  h_test_nsigmapr->Write();
  for(int i=0;i<14;i++)
  {
    m_h_nsigmapip[i]->Write();
    m_h_nsigmapik[i]->Write();
    m_h_nsigmakp[i]->Write();
  }
  return kStOK;
}

//-----------------------------------------------------------------------------
void StMyAnalysisMaker::Clear(Option_t *opt) {
}

//-----------------------------------------------------------------------------
Int_t StMyAnalysisMaker::Make() {
  if(!mPicoDstMaker) {
    return kStWarn;
  }

  mPicoDst = mPicoDstMaker->picoDst();
  if(!mPicoDst) {
    return kStWarn;
  }
  if(!mPicoDst)return 0;
  //Load event
  mPicoEvent = (StPicoEvent*)mPicoDst->event();
  if(!mPicoEvent){
    return 0;
  }
  //select minbias trigger events
  vertexPos = mPicoEvent->primaryVertex();
  //get run index
  int runId = mPicoEvent->runId();
  // Event cut
  float tofMult = mPicoEvent->btofTrayMultiplicity();
  float tofmatch= mPicoEvent->nBTOFMatch();
  float refMult = mPicoEvent->refMult();
  hvz_b->Fill(vertexPos.Z());
  hvr_b->Fill(vertexPos.X(),vertexPos.Y());
  hvzvpdvz_b->Fill(vertexPos.Z(),mPicoEvent->vzVpd());
  htofvsref_b->Fill(refMult,tofMult);
  htofmatchvsref_b->Fill(refMult,tofmatch);
  h_runidvstofmult_b->Fill(runId,tofMult);
  h_runidvsrefmult_b->Fill(runId,refMult);
  if( TMath::Abs(vertexPos.Z()) > 70 ) return 0;
  if( ((vertexPos.X()-0)*(vertexPos.X()-0)+(vertexPos.Y()-0)*(vertexPos.Y()-0)) >= 4 ) return 0;
  if(TMath::Abs(vertexPos.Z()) < 10 ){
      href_vz->Fill(mPicoEvent->refMult());
  }

  href->Fill(mPicoEvent->refMult());
  hvz->Fill(vertexPos.Z());
  hvr->Fill(vertexPos.X(),vertexPos.Y());
  hvzvpdvz->Fill(vertexPos.Z(),mPicoEvent->vzVpd());
  htofvsref->Fill(refMult,tofMult);
  htofmatchvsref->Fill(refMult,tofmatch);

  h_runidvstofmult->Fill(runId,tofMult);
  h_runidvsrefmult->Fill(runId,refMult);
  const int ntracks=mPicoDst->numberOfTracks();
  float mField = mPicoEvent->bField();
  for(int i=0; i<ntracks; i++)
  {
      mPicoTrack = (StPicoTrack*)mPicoDst->track(i);
      if(!mPicoTrack) continue;
      if(!mPicoTrack->isPrimary()) continue;
      momentum = mPicoTrack->pMom();
      StPicoPhysicalHelix helix = mPicoTrack->helix(mField);
      float dca=mPicoTrack->gDCA(vertexPos).Mag();

      h_pt->Fill(momentum.Perp());
      h_eta_b->Fill(momentum.PseudoRapidity());
      h_nhitfit->Fill(mPicoTrack->nHitsFit());
      h_nhitmax->Fill(mPicoTrack->nHitsMax());
      h_nhitratio->Fill((float)mPicoTrack->nHitsFit()/(float)mPicoTrack->nHitsMax());
      h_dca->Fill(dca);
      h_phi->Fill(momentum.Phi());
      h_eta_phi_before->Fill(momentum.Phi(),momentum.PseudoRapidity());
      // if(fabs(momentum.PseudoRapidity()) > 1.0) continue;
      if(momentum.Perp() < 0.15) continue;
      if(momentum.Mag() > 10.0) continue;
      if((Float_t)mPicoTrack->nHitsFit()/(Float_t)mPicoTrack->nHitsMax() < 0.52) continue;
      if(mPicoTrack->nHitsFit()<15) continue;
      if(fabs(dca) > 1.0) continue;
      h_eta_phi->Fill(momentum.Phi(),momentum.PseudoRapidity());
      h_eta->Fill(momentum.PseudoRapidity());
      int tofIndex = mPicoTrack->bTofPidTraitsIndex();
      Int_t   btofMatchFlag =  0;
      Float_t btofYLocal    =  -999;
      float tof = 0, L=0, beta=0.0, mass2= 0.0;
      if(tofIndex>=0) {
          StPicoBTofPidTraits *tofPid = mPicoDst->btofPidTraits(tofIndex);
          btofMatchFlag = tofPid->btofMatchFlag();
          btofYLocal    = tofPid->btofYLocal();
          if(tofPid) {
              beta = tofPid->btofBeta();
              tof = tofPid->btof();
              if(beta<1e-4) {
                  TVector3 btofHitPos_ = tofPid->btofHitPos();
                  const StThreeVectorF *btofHitPos = new StThreeVectorF(btofHitPos_.X(),btofHitPos_.Y(),btofHitPos_.Z());
                  const StThreeVectorF *vertexPos_ = new StThreeVectorF(vertexPos.X(), vertexPos.Y(), vertexPos.Z());
                  L = tofPathLength(vertexPos_, btofHitPos, helix.curvature());
                  if(tof>0) beta = L/(tof*(C_C_LIGHT/1.e9));
                  else beta = -1;
              }
          }
      }
      bool isGoodTof = btofMatchFlag >0 && beta > 0 && fabs(btofYLocal) < 1.8;
      if(isGoodTof) mass2 = momentum.Mag()*momentum.Mag()*(1./pow(beta,2)-1); else mass2 = -999;

      if(TMath::Abs(beta)>1e-5)  hbetavsp->Fill(momentum.Mag(), 1/beta);
      else hbetavsp->Fill(momentum.Mag(), 0);
      hmassvsp->Fill(momentum.Mag()/mPicoTrack->charge(), mass2);
      hdedxvsp->Fill(momentum.Mag()/mPicoTrack->charge(),mPicoTrack->dEdx());

      h_nsigmapip->Fill(mPicoTrack->nSigmaPion(),mPicoTrack->nSigmaProton());
      h_nsigmapik->Fill(mPicoTrack->nSigmaPion(),mPicoTrack->nSigmaKaon());
      h_nsigmakp->Fill(mPicoTrack->nSigmaKaon(),mPicoTrack->nSigmaProton());
      for(int i=0;i<14;i++)
      {
        if(momentum.Mag()>=StMyAnalysisMaker::p_low[i] && momentum.Mag()<=StMyAnalysisMaker::p_up[i])
        {
          m_h_nsigmapip[i]->Fill(mPicoTrack->nSigmaPion(),mPicoTrack->nSigmaProton());
          m_h_nsigmapik[i]->Fill(mPicoTrack->nSigmaPion(),mPicoTrack->nSigmaKaon());
          m_h_nsigmakp[i]->Fill(mPicoTrack->nSigmaKaon(),mPicoTrack->nSigmaProton());
        }
      }
      h_nsigmae->Fill(mPicoTrack->nSigmaElectron());

      if(momentum.Mag() > 0.2 && momentum.Mag() < 0.4){
          h_test_nsigmapi->Fill(mPicoTrack->nSigmaPion());
          h_test_nsigmapr->Fill(mPicoTrack->nSigmaProton());
      }

  }
  return kStOK;
}

//--------------------------------------------------------------
