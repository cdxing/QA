#ifndef STAR_picoMaker
#define STAR_picoMaker
#include "StMaker.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TNtuple.h"
#include "TF1.h"
#include "StarClassLibrary/StThreeVectorF.hh"
#include "StarClassLibrary/StThreeVectorD.hh"
#include "StarClassLibrary/StLorentzVectorF.hh"
#include "StarClassLibrary/StLorentzVectorD.hh"
#include "StarClassLibrary/StHelixD.hh"
#include "TString.h"
#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StEmcUtil/filters/StEmcFilter.h"
#include "StEmcUtil/projection/StEmcPosition.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StEmcRawHit.h"
#include "TVector3.h"
#include <map>


class StPicoDst;
class StPicoTrack;
class StPicoEvent;
class StPicoDstMaker;
class TH1F;
class TH2F;
class TProfile;

class StMyAnalysisMaker : public StMaker {
  public:
    StMyAnalysisMaker(const char *name, StPicoDstMaker *picoMaker, /*const*/ char *outName);
    virtual ~StMyAnalysisMaker();

    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();
    int GetRunIndex(int);
    static Float_t p_low[14];
    static Float_t p_up[14];

  private:
    StPicoDstMaker *mPicoDstMaker;
    StPicoDst      *mPicoDst;
    StPicoEvent    *mPicoEvent;
    StPicoTrack    *mPicoTrack;

    TString    mOutName;
    TString mOutputName;

    TH1F *href_vz, *hvz_b;
    TH1F *href, *hvz;
    TH2F *hvzvpdvz_b, *hvr_b, *hvzvpdvz, *hvr, *hmassvsp, *hdedxvsp, *htofvsref_b, *htofvsref;
    TH2F *htofmatchvsref, *htofmatchvsref_b;
    TH2F *hbetavsp;
    TH2F *h_eta_phi;
    TH2F *h_eta_phi_before;
    TH1F *h_counter;

    TH1F *h_pt, *h_eta_b, *h_eta, *h_nhitfit, *h_nhitmax, *h_nhitratio, *h_dca, *h_phi;
    TH2F *h_nsigmapip, *h_nsigmapik, *h_nsigmakp; TH1F *h_nsigmae;
    TH2F *m_h_nsigmapip[20], *m_h_nsigmapik[20], *m_h_nsigmakp[20];

    TH1D *h_test_nsigmapi;
    TH1D *h_test_nsigmapr;

    TProfile *h_runidvstofmult_b, *h_runidvsrefmult_b;
    TProfile *h_runidvstofmult, *h_runidvsrefmult;

    TVector3 vertexPos;
    TVector3 momentum;

    ClassDef(StMyAnalysisMaker,1)
};

#endif
