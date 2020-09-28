// -*- C++ -*-
//
// Package:    myAnalyzer/MyVtx
// Class:      MyVtx
//
/**\class MyVtx MyVtx.cc myAnalyzer/MyVtx/plugins/MyVtx.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Kevin Stenson
//         Created:  Wed, 31 Jul 2019 20:32:37 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// math
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

// reco vertex
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

// simulated track
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"

// pile-up
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// vertexing
#include "RecoVertex/PrimaryVertexProducer/interface/TrackFilterForPVFinding.h"

// simulated vertex
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
// associator
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/VertexAssociation/interface/calculateVertexSharedTracks.h"
#include "SimDataFormats/Associations/interface/VertexToTrackingVertexAssociator.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include <numeric>

#include "TMath.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class MyVtx : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  typedef math::XYZTLorentzVector LorentzVector;

  // auxiliary class holding simulated vertices
  struct simPrimaryVertex {
    simPrimaryVertex(double x1, double y1, double z1)
        :x(x1), y(y1), z(z1),
         pt(0), ptsq(0), pt_rec(0), ptsq_rec(0), trksph(0), trksph_rec(0), closest_vertex_distance_z(-1.),
         nGenTrk(0),
         num_matched_reco_tracks(0),
         average_match_quality(0.0) {
      ptot.setPx(0);
      ptot.setPy(0);
      ptot.setPz(0);
      ptot.setE(0);
      ptot_rec.setPx(0);
      ptot_rec.setPy(0);
      ptot_rec.setPz(0);
      ptot_rec.setE(0);
      p4 = LorentzVector(0, 0, 0, 0);
      r = sqrt(x*x + y*y);
    };
    double x, y, z, r;
    HepMC::FourVector ptot;
    HepMC::FourVector ptot_rec;
    LorentzVector p4;
    double pt;
    double ptsq;
    double pt_rec;
    double ptsq_rec;
    double trksph;
    double trksph_rec;
    double closest_vertex_distance_z;
    int nGenTrk;
    int num_matched_reco_tracks;
    float average_match_quality;
    EncodedEventId eventId;
    TrackingVertexRef sim_vertex;
    std::vector<const reco::Vertex *> rec_vertices;
  };

  // auxiliary class holding reconstructed vertices
  struct recoPrimaryVertex {
    enum VertexProperties {
      NONE = 0,
      MATCHED = 1,
      DUPLICATE = 2,
      MERGED = 4
    };
    recoPrimaryVertex(double x1, double y1, double z1)
      :x(x1), y(y1), z(z1), chi2prob(0),
       pt(0), ptsq(0), trkmet(0), trkmet2(0),
       trkd0(0), trkdz(0), trksph(0), trkmass(0), closest_vertex_distance_z(-1.), purity(-1.),
       nRecoTrk(0),
       num_matched_sim_tracks(0),matchnum(-1),
       kind_of_vertex(0),
       recVtx(nullptr) {
      r = sqrt(x*x + y*y);
    };
    double x, y, z, r;
    double chi2prob;
    double pt, ptsq, trkmet, trkmet2;
    double trkd0, trkdz, trksph, trkmass;
    double closest_vertex_distance_z;
    double purity; // calculated and assigned in calculatePurityAndFillHistograms
    int nRecoTrk;
    int num_matched_sim_tracks;
    int matchnum;
    int kind_of_vertex;
    std::vector<const TrackingVertex *> sim_vertices;
    std::vector<const simPrimaryVertex *> sim_vertices_internal;
    std::vector<unsigned int> sim_vertices_num_shared_tracks;
    const reco::Vertex *recVtx;
    reco::VertexBaseRef recVtxRef;
  };

   public:
      explicit MyVtx(const edm::ParameterSet&);
      ~MyVtx();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  void resetSimPVAssociation(std::vector<simPrimaryVertex>&);
  void matchSim2RecoVertices(std::vector<simPrimaryVertex>&,
                             const reco::VertexSimToRecoCollection&);
  void matchReco2SimVertices(std::vector<recoPrimaryVertex>&,
                             const reco::VertexRecoToSimCollection&,
                             const std::vector<simPrimaryVertex>&);
  bool matchRecoTrack2SimSignal(const reco::TrackBaseRef&);

  std::vector<MyVtx::simPrimaryVertex> getSimPVs(
      const edm::Handle<TrackingVertexCollection>&);

  std::vector<MyVtx::recoPrimaryVertex> getRecoPVs(
      const edm::Handle<edm::View<reco::Vertex>>&);

  // ----------member data ---------------------------
  bool verbose_;
  bool use_only_charged_tracks_;

  const reco::RecoToSimCollection *r2s_;
  const reco::SimToRecoCollection *s2r_;

  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > vecPileupSummaryInfoToken_;
  edm::EDGetTokenT<edm::View<reco::Vertex>> reco_vertex_collection_token_;
  edm::EDGetTokenT<TrackingParticleCollection> trackingParticleCollectionToken_;
  edm::EDGetTokenT<TrackingVertexCollection> trackingVertexCollectionToken_;
  edm::EDGetTokenT<reco::SimToRecoCollection> simToRecoAssociationToken_;
  edm::EDGetTokenT<reco::RecoToSimCollection> recoToSimAssociationToken_;
  edm::EDGetTokenT<reco::VertexToTrackingVertexAssociator> vertexAssociatorToken_;
  edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file

  TTree* rvtx;
  TTree* svtx;

  int n_rvtx, n_trk_rvtx[300], match_rvtx[300], matchsig_rvtx[300];
  float purity_rvtx[300], x_rvtx[300], y_rvtx[300], z_rvtx[300], prob_rvtx[300], trkpt2_rvtx[300], trkpt_rvtx[300];
  float trkmet_rvtx[300], trkmet2_rvtx[300], trkd0_rvtx[300], trkdz_rvtx[300], trksph_rvtx[300], trkmass_rvtx[300];
  int n_svtx, n_trk_svtx[300], n_match_trk_svtx[300], match_svtx[300], matchsig_svtx[300];
  float x_svtx[300], y_svtx[300], z_svtx[300], trkpt2_svtx[300], trkpt_svtx[300];
  float trkmet_svtx[300], trkmet2_svtx[300], trksph_svtx[300], trkmass_svtx[300];
  float x_rec_svtx[300], y_rec_svtx[300], z_rec_svtx[300], trkpt2_rec_svtx[300], trkpt_rec_svtx[300];
  float trkmet_rec_svtx[300], trkmet2_rec_svtx[300], trksph_rec_svtx[300], trkmass_rec_svtx[300];

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MyVtx::MyVtx(const edm::ParameterSet& iConfig)
  : verbose_(iConfig.getUntrackedParameter<bool>("verbose", false)),
    use_only_charged_tracks_(iConfig.getUntrackedParameter<bool>(
          "use_only_charged_tracks", true)),
    vecPileupSummaryInfoToken_(consumes<std::vector<PileupSummaryInfo> >(
          edm::InputTag(std::string("addPileupInfo")))),
    trackingParticleCollectionToken_(consumes<TrackingParticleCollection>(
          iConfig.getUntrackedParameter<edm::InputTag>("trackingParticleCollection"))),
    trackingVertexCollectionToken_(consumes<TrackingVertexCollection>(
          iConfig.getUntrackedParameter<edm::InputTag>("trackingVertexCollection"))),
    simToRecoAssociationToken_(consumes<reco::SimToRecoCollection>(
          iConfig.getUntrackedParameter<edm::InputTag>("trackAssociatorMap"))),
    recoToSimAssociationToken_(consumes<reco::RecoToSimCollection>(
          iConfig.getUntrackedParameter<edm::InputTag>("trackAssociatorMap"))),
    vertexAssociatorToken_(consumes<reco::VertexToTrackingVertexAssociator>(
	   iConfig.getUntrackedParameter<edm::InputTag>("vertexAssociator"))),
    tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks")))
{
  
  reco_vertex_collection_token_ = edm::EDGetTokenT<edm::View<reco::Vertex>>(
	     consumes<edm::View<reco::Vertex>>(
             iConfig.getUntrackedParameter<edm::InputTag>("vertexRecoCollection")));
  
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> fs;
   rvtx = fs->make<TTree>("rvtx","tree of reco vertices");
   rvtx->Branch("n_rvtx", &n_rvtx, "n_rvtx/I");
   rvtx->Branch("match_rvtx", match_rvtx, "match_rvtx[n_rvtx]/I");
   rvtx->Branch("matchsig_rvtx", matchsig_rvtx, "matchsig_rvtx[n_rvtx]/I");
   rvtx->Branch("purity_rvtx", purity_rvtx, "purity_rvtx[n_rvtx]/F");
   rvtx->Branch("n_trk_rvtx", n_trk_rvtx, "n_trk_rvtx[n_rvtx]/I");
   rvtx->Branch("x_rvtx", x_rvtx, "x_rvtx[n_rvtx]/F");
   rvtx->Branch("y_rvtx", y_rvtx, "y_rvtx[n_rvtx]/F");
   rvtx->Branch("z_rvtx", z_rvtx, "z_rvtx[n_rvtx]/F");
   rvtx->Branch("prob_rvtx", prob_rvtx, "prob_rvtx[n_rvtx]/F");
   rvtx->Branch("trkpt2_rvtx", trkpt2_rvtx, "trkpt2_rvtx[n_rvtx]/F");
   rvtx->Branch("trkpt_rvtx", trkpt_rvtx, "trkpt_rvtx[n_rvtx]/F");
   rvtx->Branch("trkd0_rvtx", trkd0_rvtx, "trkd0_rvtx[n_rvtx]/F");
   rvtx->Branch("trkdz_rvtx", trkdz_rvtx, "trkdz_rvtx[n_rvtx]/F");
   rvtx->Branch("trkmet_rvtx", trkmet_rvtx, "trkmet_rvtx[n_rvtx]/F");
   rvtx->Branch("trkmet2_rvtx", trkmet2_rvtx, "trkmet2_rvtx[n_rvtx]/F");
   rvtx->Branch("trksph_rvtx", trksph_rvtx, "trksph_rvtx[n_rvtx]/F");
   rvtx->Branch("trkmass_rvtx", trkmass_rvtx, "trkmass_rvtx[n_rvtx]/F");
   
   svtx = fs->make<TTree>("svtx","tree of sim vertices");
   svtx->Branch("n_svtx", &n_svtx, "n_svtx/I");
   svtx->Branch("match_svtx", match_svtx, "match_svtx[n_svtx]/I");
   svtx->Branch("matchsig_svtx", matchsig_svtx, "matchsig_svtx[n_svtx]/I");
   svtx->Branch("n_trk_svtx", n_trk_svtx, "n_trk_svtx[n_svtx]/I");
   svtx->Branch("n_match_trk_svtx", n_match_trk_svtx, "n_match_trk_svtx[n_svtx]/I");
   svtx->Branch("x_svtx", x_svtx, "x_svtx[n_svtx]/F");
   svtx->Branch("y_svtx", y_svtx, "y_svtx[n_svtx]/F");
   svtx->Branch("z_svtx", z_svtx, "z_svtx[n_svtx]/F");
   svtx->Branch("trkpt2_svtx", trkpt2_svtx, "trkpt2_svtx[n_svtx]/F");
   svtx->Branch("trkpt_svtx", trkpt_svtx, "trkpt_svtx[n_svtx]/F");
   svtx->Branch("trkmet_svtx", trkmet_svtx, "trkmet_svtx[n_svtx]/F");
   svtx->Branch("trkmet2_svtx", trkmet2_svtx, "trkmet2_svtx[n_svtx]/F");
   svtx->Branch("trksph_svtx", trksph_svtx, "trksph_svtx[n_svtx]/F");
   svtx->Branch("trkmass_svtx", trkmass_svtx, "trkmass_svtx[n_svtx]/F");
   svtx->Branch("trkpt2_rec_svtx", trkpt2_rec_svtx, "trkpt2_rec_svtx[n_svtx]/F");
   svtx->Branch("trkpt_rec_svtx", trkpt_rec_svtx, "trkpt_rec_svtx[n_svtx]/F");
   svtx->Branch("trkmet_rec_svtx", trkmet_rec_svtx, "trkmet_rec_svtx[n_svtx]/F");
   svtx->Branch("trkmet2_rec_svtx", trkmet2_rec_svtx, "trkmet2_rec_svtx[n_svtx]/F");
   svtx->Branch("trksph_rec_svtx", trksph_rec_svtx, "trksph_rec_svtx[n_svtx]/F");
   svtx->Branch("trkmass_rec_svtx", trkmass_rec_svtx, "trkmass_rec_svtx[n_svtx]/F");

   
  //   histo = fs->make<TH1I>("charge" , "Charges" , 2 , -1 , 1 );

}


MyVtx::~MyVtx()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  rvtx->Write();
  svtx->Write();

}


//
// member functions
//
bool MyVtx::matchRecoTrack2SimSignal(const reco::TrackBaseRef& recoTrack) {
  auto found = r2s_->find(recoTrack);

  // reco track not matched to any TP
  if(found == r2s_->end())
    return false;

  // reco track matched to some TP from signal vertex
  for(const auto& tp: found->val) {
    if(tp.first->eventId().bunchCrossing() == 0 && tp.first->eventId().event() == 0)
      return true;
  }

  // reco track not matched to any TP from signal vertex
  return false;
}

/* Extract information form TrackingParticles/TrackingVertex and fill
 * the helper class simPrimaryVertex with proper generation-level
 * information */
std::vector<MyVtx::simPrimaryVertex>
MyVtx::getSimPVs(
    const edm::Handle<TrackingVertexCollection>& tVC) {
  std::vector<MyVtx::simPrimaryVertex> simpv;
  int current_event = -1;

  if (verbose_) {
    std::cout << "getSimPVs TrackingVertexCollection " << std::endl;
  }

  for (TrackingVertexCollection::const_iterator v = tVC->begin();
       v != tVC->end(); ++v) {
//    if (verbose_) {
//      std::cout << "BunchX.EventId: " << v->eventId().bunchCrossing() << "."
//                << (v->eventId()).event() << " Position: " << v->position()
//                << " G4/HepMC Vertices: " << v->g4Vertices().size() << "/"
//                << v->genVertices().size()
//                << "   t = " << v->position().t() * 1.e12
//                << "    == 0:" << (v->position().t() > 0) << std::endl;
//      for (TrackingVertex::g4v_iterator gv = v->g4Vertices_begin();
//           gv != v->g4Vertices_end(); gv++) {
//        std::cout << *gv << std::endl;
//      }
//      std::cout << "----------" << std::endl;
//
//    }  // end of verbose_ session

    // I'd rather change this and select only vertices that come from
    // BX=0.  We should keep only the first vertex from all the events
    // at BX=0.
    if (v->eventId().bunchCrossing() != 0) continue;
    if (v->eventId().event() != current_event) {
      current_event = v->eventId().event();
    } else {
      continue;
    }
    // TODO(rovere) is this really necessary?
    if (fabs(v->position().z()) > 1000) continue;  // skip funny junk vertices

    // could be a new vertex, check  all primaries found so far to avoid
    // multiple entries
    simPrimaryVertex sv(v->position().x(), v->position().y(),
                        v->position().z());
    sv.eventId = v->eventId();
    sv.sim_vertex = TrackingVertexRef(tVC, std::distance(tVC->begin(), v));

    // TODO(rovere) maybe get rid of this old logic completely ... ?
    simPrimaryVertex* vp = nullptr;  // will become non-NULL if a vertex
                                  // is found and then point to it
    for (std::vector<simPrimaryVertex>::iterator v0 = simpv.begin();
         v0 != simpv.end(); v0++) {
      if ((sv.eventId == v0->eventId) && (fabs(sv.x - v0->x) < 1e-5) &&
          (fabs(sv.y - v0->y) < 1e-5) && (fabs(sv.z - v0->z) < 1e-5)) {
        vp = &(*v0);
        break;
      }
    }
    if (!vp) {
      // this is a new vertex, add it to the list of sim-vertices
      simpv.push_back(sv);
      vp = &simpv.back();
      if (verbose_) {
        std::cout << "this is a new vertex " << sv.eventId.event() << "   "
                  << sv.x << " " << sv.y << " " << sv.z << std::endl;
      }
    } else {
      if (verbose_) {
        std::cout << "this is not a new vertex" << sv.x << " " << sv.y << " "
                  << sv.z << std::endl;
      }
    }

    // Loop over daughter track(s) as Tracking Particles
    for (TrackingVertex::tp_iterator iTP = v->daughterTracks_begin();
         iTP != v->daughterTracks_end(); ++iTP) {
      auto momentum = (*(*iTP)).momentum();
      const reco::Track* matched_best_reco_track = nullptr;
      double match_quality = -1;
      if (use_only_charged_tracks_ && (**iTP).charge() == 0)
          continue;
      if (s2r_->find(*iTP) != s2r_->end()) {
        matched_best_reco_track = (*s2r_)[*iTP][0].first.get();
        match_quality = (*s2r_)[*iTP][0].second;
      }
//      if (verbose_) {
//        std::cout << "  Daughter momentum:      " << momentum;
//        std::cout << "  Daughter type     " << (*(*iTP)).pdgId();
//        std::cout << "  matched: " << (matched_best_reco_track != nullptr);
//        std::cout << "  match-quality: " << match_quality;
//        std::cout << std::endl;
//      }
      vp->ptot.setPx(vp->ptot.x() + momentum.x());
      vp->ptot.setPy(vp->ptot.y() + momentum.y());
      vp->ptot.setPz(vp->ptot.z() + momentum.z());
      vp->ptot.setE(vp->ptot.e() + (**iTP).energy());
      vp->pt += (**iTP).pt();
      vp->ptsq += ((**iTP).pt() * (**iTP).pt());
      // TODO(rovere) only select charged sim-particles? If so, maybe
      // put it as a configuration parameter?
      if (matched_best_reco_track) {
        vp->num_matched_reco_tracks++;
        vp->average_match_quality += match_quality;
	vp->ptot_rec.setPx(vp->ptot_rec.x() + momentum.x());
	vp->ptot_rec.setPy(vp->ptot_rec.y() + momentum.y());
	vp->ptot_rec.setPz(vp->ptot_rec.z() + momentum.z());
	vp->ptot_rec.setE(vp->ptot_rec.e() + (**iTP).energy());
	vp->pt_rec += (**iTP).pt();
	vp->ptsq_rec += ((**iTP).pt() * (**iTP).pt());	
      }
      // TODO(rovere) get rid of cuts on sim-tracks
      // TODO(rovere) be consistent between simulated tracks and
      // reconstructed tracks selection
      // count relevant particles
      if (((**iTP).pt() > 0.2) && (fabs((**iTP).eta()) < 3.8) &&
          (**iTP).charge() != 0) {
        vp->nGenTrk++;
      }
    }  // End of for loop on daughters sim-particles
    if (vp->num_matched_reco_tracks)
      vp->average_match_quality /=
          static_cast<float>(vp->num_matched_reco_tracks);
    if (verbose_) {
      std::cout << "matched reco tracks: " << vp->num_matched_reco_tracks
		<< ", gen trks in acceptance: " << vp->nGenTrk
		<< ", ratio: " << vp->num_matched_reco_tracks / static_cast<float>(vp->nGenTrk)
                << " with average quality: " << vp->average_match_quality
                << std::endl;
    }
  }  // End of for loop on tracking vertices

  if (verbose_) {
    std::cout << "------- MyVtx simPVs from TrackingVertices -------" << std::endl;
    for (std::vector<simPrimaryVertex>::iterator v0 = simpv.begin();
         v0 != simpv.end(); v0++) {
      std::cout << "z=" << v0->z << "  event=" << v0->eventId.event()
                << std::endl;
    }
    std::cout << "-----------------------------------------------" << std::endl;
  }  // End of for summary on discovered simulated primary vertices.

  // In case of no simulated vertices, break here
  if(simpv.empty())
    return simpv;

  // Now compute the closest distance in z between all simulated vertex
  // first initialize
  auto prev_z = simpv.back().z;
  for(simPrimaryVertex& vsim: simpv) {
    vsim.closest_vertex_distance_z = std::abs(vsim.z - prev_z);
    prev_z = vsim.z;
  }
  // then calculate
  for (std::vector<simPrimaryVertex>::iterator vsim = simpv.begin();
       vsim != simpv.end(); vsim++) {
    std::vector<simPrimaryVertex>::iterator vsim2 = vsim;
    vsim2++;
    for (; vsim2 != simpv.end(); vsim2++) {
      double distance = std::abs(vsim->z - vsim2->z);
      // need both to be complete
      vsim->closest_vertex_distance_z = std::min(vsim->closest_vertex_distance_z, distance);
      vsim2->closest_vertex_distance_z = std::min(vsim2->closest_vertex_distance_z, distance);
    }
  }
  return simpv;
}

/* Extract information form recoVertex and fill the helper class
 * recoPrimaryVertex with proper reco-level information */
std::vector<MyVtx::recoPrimaryVertex>
MyVtx::getRecoPVs(
    const edm::Handle<edm::View<reco::Vertex>>& tVC) {
  std::vector<MyVtx::recoPrimaryVertex> recopv;

  if (verbose_) {
    std::cout << "getRecoPVs VertexCollection " << std::endl;
  }

  for (auto v = tVC->begin(); v != tVC->end(); ++v) {
    if (verbose_) {
      std::cout << " Position: " << v->position() << std::endl;
    }

    // Skip junk vertices
    if (fabs(v->z()) > 1000) continue;
    if (v->isFake() || !v->isValid()) continue;

    recoPrimaryVertex sv(v->position().x(), v->position().y(),
                         v->position().z());
    sv.recVtx = &(*v);
    sv.recVtxRef = reco::VertexBaseRef(tVC, std::distance(tVC->begin(), v));
    double pvchi2 = v->chi2();
    double pvndof = v->ndof();
    sv.chi2prob = ChiSquaredProbability(pvchi2,pvndof);
    // this is a new vertex, add it to the list of reco-vertices
    recopv.push_back(sv);
    MyVtx::recoPrimaryVertex* vp = &recopv.back();

    // Loop over daughter track(s)
    double pvpx = 0;
    double pvpy = 0;
    double pvd0 = 0;
    double pvdz = 0;
    double pve = 0;
    double pvpz = 0;
    for (auto iTrack = v->tracks_begin(); iTrack != v->tracks_end(); ++iTrack) {
      auto momentum = (*(*iTrack)).innerMomentum();
      // TODO(rovere) better handle the pixelVertices, whose tracks
      // do not have the innerMomentum defined. This is a temporary
      // hack to overcome this problem.
      if (momentum.mag2() == 0)
        momentum = (*(*iTrack)).momentum();
//      if (verbose_) {
//        std::cout << "  Daughter momentum:      " << momentum << std::endl;
//      }
      vp->pt += std::sqrt(momentum.perp2());
      vp->ptsq += (momentum.perp2());
      vp->nRecoTrk++;
      pvpx += momentum.x();
      pvpy += momentum.y();
      pvpz += momentum.z();
      pve += std::sqrt(pvpx*pvpx + pvpy*pvpy + pvpz*pvpz);
      auto matched = r2s_->find(*iTrack);
      if(matched != r2s_->end()) {
        vp->num_matched_sim_tracks++;
      }

    }  // End of for loop on daughters reconstructed tracks
    vp->trkmet2 = pvpx*pvpx + pvpy*pvpy;
    vp->trkmet = std::sqrt(vp->trkmet2);
    vp->trkmass = std::sqrt(pve*pve - (pvpx*pvpx + pvpy*pvpy + pvpz*pvpz));
    vp->trkdz = pvdz;
    vp->trkd0 = pvd0;
  }    // End of for loop on reco vertices

  if (verbose_) {
    std::cout << "------- MyVtx recoPVs from VertexCollection -------" << std::endl;
    for (std::vector<recoPrimaryVertex>::iterator v0 = recopv.begin();
         v0 != recopv.end(); v0++) {
      std::cout << "z=" << v0->z << std::endl;
    }
    std::cout << "-----------------------------------------------" << std::endl;
  }  // End of for summary on reconstructed primary vertices.

  // In case of no reco vertices, break here
  if(recopv.empty())
    return recopv;

  // Now compute the closest distance in z between all reconstructed vertex
  // first initialize
  auto prev_z = recopv.back().z;
  for(recoPrimaryVertex& vreco: recopv) {
    vreco.closest_vertex_distance_z = std::abs(vreco.z - prev_z);
    prev_z = vreco.z;
  }
  for (std::vector<recoPrimaryVertex>::iterator vreco = recopv.begin();
       vreco != recopv.end(); vreco++) {
    std::vector<recoPrimaryVertex>::iterator vreco2 = vreco;
    vreco2++;
    for (; vreco2 != recopv.end(); vreco2++) {
      double distance = std::abs(vreco->z - vreco2->z);
      // need both to be complete
      vreco->closest_vertex_distance_z = std::min(vreco->closest_vertex_distance_z, distance);
      vreco2->closest_vertex_distance_z = std::min(vreco2->closest_vertex_distance_z, distance);
    }
  }
  return recopv;
}

void MyVtx::resetSimPVAssociation(
    std::vector<simPrimaryVertex> & simpv) {
  for (auto & v : simpv) {
    v.rec_vertices.clear();
  }
}

// ------------ method called to produce the data  ------------
void MyVtx::matchSim2RecoVertices(
    std::vector<simPrimaryVertex>& simpv,
    const reco::VertexSimToRecoCollection& vertex_s2r) {
  if (verbose_) {
    std::cout << "MyVtx::matchSim2RecoVertices " << std::endl;
  }
  for (std::vector<simPrimaryVertex>::iterator vsim = simpv.begin();
       vsim != simpv.end(); vsim++) {

    auto matched = vertex_s2r.find(vsim->sim_vertex);
    if(matched != vertex_s2r.end()) {
      for(const auto vertexRefQuality: matched->val) {
        vsim->rec_vertices.push_back(&(*(vertexRefQuality.first)));
      }
    }

    if (verbose_) {
      if (!vsim->rec_vertices.empty()) {
	if (vsim->rec_vertices.size()>1) {
	  std::cout << "Found two matching vertices for genVtx" << std::endl;
	}
        for (auto const& v : vsim->rec_vertices) {
          std::cout << "Found a matching vertex for genVtx "
                    << vsim->z << " at " << v->z()
                    << " with sign: " << fabs(v->z() - vsim->z) / v->zError()
                    << std::endl;
        }
      } else {
        std::cout << "No matching vertex for " << vsim->z << std::endl;
      }
    }
  }  // end for loop on simulated vertices
  if (verbose_) {
    std::cout << "Done with matching sim vertices" << std::endl;
  }
}

void MyVtx::matchReco2SimVertices(
    std::vector<recoPrimaryVertex>& recopv,
    const reco::VertexRecoToSimCollection& vertex_r2s,
    const std::vector<simPrimaryVertex>& simpv) {
  for (std::vector<recoPrimaryVertex>::iterator vrec = recopv.begin();
       vrec != recopv.end(); vrec++) {

    auto matched = vertex_r2s.find(vrec->recVtxRef);
    if(matched != vertex_r2s.end()) {
      for(const auto vertexRefQuality: matched->val) {
        const auto tvPtr = &(*(vertexRefQuality.first));
        vrec->sim_vertices.push_back(tvPtr);
      }

      for(const TrackingVertex *tv: vrec->sim_vertices) {
        // Set pointers to internal simVertex objects
	int ivtx = -1;
        for(const auto& vv: simpv) {
	  ivtx++;
          if (&(*(vv.sim_vertex)) == tv) {
	    vrec->matchnum = ivtx;
            vrec->sim_vertices_internal.push_back(&vv);
            continue;
          }
        }

        // Calculate number of shared tracks
        vrec->sim_vertices_num_shared_tracks.push_back(calculateVertexSharedTracks(*(vrec->recVtx), *tv, *r2s_));
      }
    }

    if (verbose_) {
      for (auto v : vrec->sim_vertices) {
        std::cout << "Found a matching vertex for reco: " << vrec->z
                  << " at gen:" << v->position().z() << " with sign: "
                  << fabs(vrec->z - v->position().z()) / vrec->recVtx->zError()
                  << std::endl;
      }
    }
  }  // end for loop on reconstructed vertices
}

// ------------ method called for each event  ------------
void
MyVtx::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using std::vector;
  using std::cout;
  using std::endl;
  using edm::Handle;
  using edm::View;
  using namespace reco;
  using namespace edm;

  std::vector<float> pileUpInfo_z;

  edm::Handle<std::vector<PileupSummaryInfo> > puinfoH;
  if (iEvent.getByToken(vecPileupSummaryInfoToken_, puinfoH)) {
    for (auto const& pu_info : *puinfoH.product()) {
//      if(do_generic_sim_plots_) {
//        mes_["root_folder"]["GenVtx_vs_BX"]
//          ->Fill(pu_info.getBunchCrossing(), pu_info.getPU_NumInteractions());
//      }
      if (pu_info.getBunchCrossing() == 0) {
        pileUpInfo_z = pu_info.getPU_zpositions();
        if (verbose_) {
          for (auto const& p : pileUpInfo_z) {
            std::cout << "PileUpInfo on Z vertex: " << p << std::endl;
          }
        }
        break;
      }
    }
  }

  edm::Handle<TrackingParticleCollection> TPCollectionH;
  iEvent.getByToken(trackingParticleCollectionToken_, TPCollectionH);
  if (!TPCollectionH.isValid()) 
    edm::LogWarning("MyVtx")
      << "TPCollectionH is not valid";

  edm::Handle<TrackingVertexCollection> TVCollectionH;
  iEvent.getByToken(trackingVertexCollectionToken_, TVCollectionH);
  if (!TVCollectionH.isValid())
    edm::LogWarning("MyVtx")
      << "TVCollectionH is not valid";

  edm::Handle<reco::SimToRecoCollection> simToRecoH;
  iEvent.getByToken(simToRecoAssociationToken_, simToRecoH);
  if ( simToRecoH.isValid() )
    s2r_ = simToRecoH.product();
  else
    edm::LogWarning("MyVtx") << "simToRecoH is not valid";

  edm::Handle<reco::RecoToSimCollection> recoToSimH;
  iEvent.getByToken(recoToSimAssociationToken_, recoToSimH);
  if ( recoToSimH.isValid() )
    r2s_ = recoToSimH.product();
  else
    edm::LogWarning("MyVtx") << "recoToSimH is not valid";

  // Vertex associator
  edm::Handle<reco::VertexToTrackingVertexAssociator> vertexAssociatorH;
  iEvent.getByToken(vertexAssociatorToken_, vertexAssociatorH);
  if (!vertexAssociatorH.isValid()) {
    edm::LogWarning("MyVtx") << "vertexAssociatorH is not valid";
    return;
  }
  const reco::VertexToTrackingVertexAssociator& vertexAssociator = *(vertexAssociatorH.product());

  std::vector<simPrimaryVertex> simpv;  // a list of simulated primary MC vertices
  simpv = getSimPVs(TVCollectionH);
  int num_pileup_vertices = simpv.size();

  std::vector<recoPrimaryVertex> recopv;  // a list of reconstructed primary MC vertices
  edm::Handle<edm::View<reco::Vertex>> RecoVCollectionH;
  iEvent.getByToken(reco_vertex_collection_token_, RecoVCollectionH);
  if (!RecoVCollectionH.isValid())
    edm::LogWarning("MyVtx") << "RecoVCollectionH is not valid";
  
  reco::VertexRecoToSimCollection vertex_r2s = vertexAssociator.associateRecoToSim(RecoVCollectionH, TVCollectionH);
  reco::VertexSimToRecoCollection vertex_s2r = vertexAssociator.associateSimToReco(RecoVCollectionH, TVCollectionH);

  resetSimPVAssociation(simpv);
  matchSim2RecoVertices(simpv, vertex_s2r);
  recopv = getRecoPVs(RecoVCollectionH);
  matchReco2SimVertices(recopv, vertex_r2s, simpv);

  int num_total_gen_vertices_assoc2reco = 0;
  int num_total_reco_vertices_assoc2gen = 0;
  int num_total_gen_vertices_multiassoc2reco = 0;
  int num_total_reco_vertices_multiassoc2gen = 0;
  int num_total_reco_vertices_duplicate = 0;
  int genpv_position_in_reco_collection = -1;
  for (auto const& v : simpv) {
    float mistag = 1.;
    if (v.eventId.event() == 0) {
      if (!RecoVCollectionH->empty() &&
	  std::find(v.rec_vertices.begin(), v.rec_vertices.end(),
		    &((*RecoVCollectionH.product())[0])) != v.rec_vertices.end()) {
	mistag = 0.;
      }
      if (verbose_) std::cout << "mistag = " << mistag << std::endl;
//      mes_[label]["MisTagRate"]->Fill(mistag);
//      mes_[label]["MisTagRate_vs_PU"]->Fill(simpv.size(), mistag);
//      mes_[label]["MisTagRate_vs_sum-pt2"]->Fill(v.ptsq, mistag);
//      mes_[label]["MisTagRate_vs_Z"]->Fill(v.z, mistag);
//      mes_[label]["MisTagRate_vs_R"]->Fill(v.r, mistag);
//      mes_[label]["MisTagRate_vs_NumTracks"]->Fill(v.nGenTrk, mistag);
        // Now check at which location the Simulated PV has been
        // reconstructed in the primary vertex collection
        // at-hand. Mark it with fake index -1 if it was not
        // reconstructed at all.

      auto iv = (*RecoVCollectionH.product()).begin();
      for (int pv_position_in_reco_collection = 0;
	   iv != (*RecoVCollectionH.product()).end();
	   ++pv_position_in_reco_collection, ++iv) {
	if (std::find(v.rec_vertices.begin(), v.rec_vertices.end(),
		      &(*iv)) != v.rec_vertices.end()) {
	  //	  mes_[label]["TruePVLocationIndex"]->Fill(pv_position_in_reco_collection);
	  const bool genPVMatchedToRecoPV = (pv_position_in_reco_collection == 0);
	  //	  mes_[label]["TruePVLocationIndexCumulative"]->Fill(genPVMatchedToRecoPV ? 0 : 1);

	  //	  fillRecoAssociatedGenPVHistograms(label, v, genPVMatchedToRecoPV);
	  if(genPVMatchedToRecoPV) {
	    auto pv = recopv[0];
	    assert(pv.recVtx == &(*iv));
	    //	    fillResolutionAndPullHistograms(label, num_pileup_vertices, pv, true);
	  }
	  genpv_position_in_reco_collection = pv_position_in_reco_collection;
	  if (verbose_) std::cout << "genpv_position_in_reco_collection = " << genpv_position_in_reco_collection << std::endl;
	  break;
	}
      }

      // If we reached the end, it means that the Simulated PV has not
      // been associated to any reconstructed vertex: mark it as
      // missing in the reconstructed vertex collection using the fake
      // index -1.
      if (iv == (*RecoVCollectionH.product()).end()) {
	//	mes_[label]["TruePVLocationIndex"]->Fill(-1.);
	//	mes_[label]["TruePVLocationIndexCumulative"]->Fill(-1.);
      }
    }
    
    if (!v.rec_vertices.empty()) num_total_gen_vertices_assoc2reco++;
    if (v.rec_vertices.size() > 1) num_total_gen_vertices_multiassoc2reco++;
    // No need to N-tplicate the Gen-related cumulative histograms:
    // fill them only at the first iteration
    //    fillRecoAssociatedGenVertexHistograms(label, v);
  }
    //  calculatePurityAndFillHistograms(label, recopv, genpv_position_in_reco_collection, signal_is_highest_pt);

  std::vector<double> vtx_sumpt_sigmatched;
  std::vector<double> vtx_sumpt2_sigmatched;
  vtx_sumpt_sigmatched.reserve(recopv.size());
  vtx_sumpt2_sigmatched.reserve(recopv.size());

  // Calculate purity
  for(auto& v: recopv) {
    double sumpt_all = 0;
    double sumpt_sigmatched = 0;
    double sumpt2_sigmatched = 0;
    const reco::Vertex *vertex = v.recVtx;
    for(auto iTrack = vertex->tracks_begin(); iTrack != vertex->tracks_end(); ++iTrack) {
      double pt = (*iTrack)->pt();
      sumpt_all += pt;
      if(matchRecoTrack2SimSignal(*iTrack)) {
        sumpt_sigmatched += pt;
        sumpt2_sigmatched += pt*pt;
      }
    }
    v.purity = sumpt_sigmatched / sumpt_all;

    vtx_sumpt_sigmatched.push_back(sumpt_sigmatched);
    vtx_sumpt2_sigmatched.push_back(sumpt2_sigmatched);
  }

  const double vtxAll_sumpt_sigmatched = std::accumulate(vtx_sumpt_sigmatched.begin(), vtx_sumpt_sigmatched.end(), 0.0);
  const double vtxNot0_sumpt_sigmatched = vtxAll_sumpt_sigmatched - vtx_sumpt_sigmatched[0];
  const double missing = vtxNot0_sumpt_sigmatched / vtxAll_sumpt_sigmatched;

    //  mes_[label]["GenAllAssoc2Reco_NumVertices"]->Fill(simpv.size(), simpv.size());
    //  mes_[label]["GenAllAssoc2RecoMatched_NumVertices"]->Fill(simpv.size(), num_total_gen_vertices_assoc2reco);
    //  mes_[label]["GenAllAssoc2RecoMultiMatched_NumVertices"]->Fill(simpv.size(), num_total_gen_vertices_multiassoc2reco);
  for (auto & v : recopv) {
    //    fillGenAssociatedRecoVertexHistograms(label, num_pileup_vertices, v);
    if (!v.sim_vertices.empty()) {
      num_total_reco_vertices_assoc2gen++;
      if (v.sim_vertices_internal[0]->rec_vertices.size() > 1) {
	num_total_reco_vertices_duplicate++;
      }
    }
    if (v.sim_vertices.size() > 1) num_total_reco_vertices_multiassoc2gen++;
  }
//  mes_[label]["RecoAllAssoc2Gen_NumVertices"]->Fill(recopv.size(), recopv.size());
//  mes_[label]["RecoAllAssoc2GenMatched_NumVertices"]->Fill(recopv.size(), num_total_reco_vertices_assoc2gen);
//  mes_[label]["RecoAllAssoc2GenMultiMatched_NumVertices"]->Fill(recopv.size(), num_total_reco_vertices_multiassoc2gen);
//  mes_[label]["RecoAllAssoc2MultiMatchedGen_NumVertices"]->Fill(recopv.size(), num_total_reco_vertices_duplicate);
//  mes_[label]["RecoVtx_vs_GenVtx"]->Fill(simpv.size(), recopv.size());
//  mes_[label]["MatchedRecoVtx_vs_GenVtx"]->Fill(simpv.size(), num_total_reco_vertices_assoc2gen);
  
  n_rvtx = recopv.size();
  n_svtx = num_pileup_vertices;
  for (int i = 0; i<n_svtx; ++i) {
    match_svtx[i] = -1;
  }
  for (int i = 0; i<n_rvtx; ++i) {
    match_rvtx[i] = recopv[i].matchnum;
    if (recopv[i].matchnum >= 0) {
      match_svtx[match_rvtx[i]] = i;
    }
    matchsig_rvtx[i] = 0;
    purity_rvtx[i] = recopv[i].purity;
    n_trk_rvtx[i] = recopv[i].nRecoTrk;
    x_rvtx[i] = recopv[i].x;
    y_rvtx[i] = recopv[i].y;
    z_rvtx[i] = recopv[i].z;
    prob_rvtx[i] = recopv[i].chi2prob;
    trkpt2_rvtx[i] = recopv[i].ptsq;
    trkpt_rvtx[i] = recopv[i].pt;
    trkd0_rvtx[i] = recopv[i].trkd0;
    trkdz_rvtx[i] = recopv[i].trkdz;
    trkmet_rvtx[i] = recopv[i].trkmet;
    trkmet2_rvtx[i] = recopv[i].trkmet2;
    trksph_rvtx[i] = recopv[i].trksph;
    trkmass_rvtx[i] = recopv[i].trkmass;
  }
  rvtx->Fill();
  std::cout << "recopv size = " << n_rvtx << ", simpv size = " << n_svtx << std::endl;
  for (int i = 0; i<n_svtx; ++i) {
    matchsig_svtx[i] = 0;
    n_trk_svtx[i] = simpv[i].nGenTrk;
    n_match_trk_svtx[i] = simpv[i].num_matched_reco_tracks;
    x_svtx[i] = simpv[i].x;
    y_svtx[i] = simpv[i].y;
    z_svtx[i] = simpv[i].z;
    trkpt2_svtx[i] = simpv[i].ptsq;
    trkpt_svtx[i] = simpv[i].pt;
    trkmet_svtx[i] = simpv[i].ptot.perp();
    trkmet2_svtx[i] = simpv[i].ptot.perp2();
    trksph_svtx[i] = simpv[i].trksph;
    trkmass_svtx[i] = simpv[i].ptot.m();
    trkpt2_rec_svtx[i] = simpv[i].ptsq_rec;
    trkpt_rec_svtx[i] = simpv[i].pt_rec;
    trkmet_rec_svtx[i] = simpv[i].ptot_rec.perp();
    trkmet2_rec_svtx[i] = simpv[i].ptot_rec.perp2();
    trksph_rec_svtx[i] = simpv[i].trksph_rec;
    trkmass_rec_svtx[i] = simpv[i].ptot_rec.m();
  }
  
  svtx->Fill();

//  for(const auto& track : iEvent.get(tracksToken_) ) {
//      // do something with track parameters, e.g, plot the charge.
//    //    int charge = track.charge();
//    //       histo->Fill( track.charge() );
//  }

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void
MyVtx::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
MyVtx::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MyVtx::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyVtx);
