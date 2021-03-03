#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h" 

#include <TTree.h>
#include <TLorentzVector.h>
 
#include <vector>
#include <map>
#include <any>
#include <set>
#include <memory>
#include <cmath>
#include <iostream>
#include <string>
#include <cstdlib>
using std::vector;
using std::string;
using std::map;
using std::unordered_map;

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/View.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

template <class T> string typeStr(){return typeid(T).name();}
template <> string typeStr<bool>(){return "bool";}
template <> string typeStr<int>(){return "int";}
template <> string typeStr<float>(){return "float";}

class hgcalfinecalontupler: public edm::one::EDAnalyzer<edm::one::SharedResources>  {
    public:
        explicit hgcalfinecalontupler(const edm::ParameterSet&);
        ~hgcalfinecalontupler() {}
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    private:
        void beginJob() override;
        void doBeginRun_(const edm::Run&, const edm::EventSetup&) override {}
        void analyze(const edm::Event&, const edm::EventSetup&) override;
        void doEndRun_(const edm::Run&, const edm::EventSetup&) override {}
        void endJob() override {}

        edm::Service<TFileService> fs;
        TTree* tree_;
        hgcal::RecHitTools hgcalRecHitToolInstance_ ;
        edm::EDGetTokenT<edm::View<PCaloHit>> hgcalEEHitsToken_;
        edm::EDGetTokenT<edm::View<PCaloHit>> hgcalHEfrontHitsToken_;
        edm::EDGetTokenT<edm::View<PCaloHit>> hgcalHEbackHitsToken_;
        edm::EDGetTokenT<edm::SimTrackContainer> tokenSimTracks;
        edm::EDGetTokenT<edm::SimVertexContainer> tokenSimVertices;

        unordered_map<string, std::any> vectors_;
        unordered_map<string, std::any> scalars_;

        template <class T> void makeBranch(string name){
            vectors_.emplace(name, vector<T>());
            tree_->Branch(
                name.c_str(),
                ("vector<" + typeStr<T>() +">").c_str(),
                std::any_cast< vector<T> >(&vectors_[name]),
                32000, 0
                );
            }

        /* Return branch as vector by reference */
        template <class T> vector<T>& getBranch(string name){
            return std::any_cast< vector<T>& >(vectors_[name]);
            }

        /* Add a value to a branch */
        template <class T> void fill(string name, T value){
            auto it = vectors_.find(name);
            if (it != vectors_.end()) {
                std::any_cast< vector<T> >(&(it->second))->push_back(value);
                }
            else {
                throw cms::Exception("Unknown")
                    << "Cannot fill branch " << name
                    << "; no such branch";
                }
            }

        template <class T> void clearBranch(string name){
            std::any_cast< vector<T>& >(vectors_[name]).clear();
            }

        template <class T> void makeScalarBranch(string name){
            T value = 0;
            scalars_.emplace(name, value);
            tree_->Branch(name.c_str(), std::any_cast<T>(&scalars_[name]));
            }

        template <class T> void fillScalar(string name, T value){
            *(std::any_cast<T>(&scalars_[name])) = value;
            }

        /*
        Clears the vectors stored in the branches.
        Unfortunately can't figure out how to do this without the manual type check
        */
        void clear(){
            for( auto& element : vectors_ ) {
                if (auto* ptr = std::any_cast<vector<bool>>(&(element.second))) {
                    ptr->clear();
                    }
                else if (auto* ptr = std::any_cast<vector<int>>(&(element.second))) {
                    ptr->clear();
                    }
                else if (auto* ptr = std::any_cast<vector<float>>(&(element.second))) {
                    ptr->clear();
                    }
                else {
                    throw cms::Exception("Unknown")
                        << "Tried to clear "
                        << element.first
                        << ", but type is unknown";
                    }
            for( auto& element : scalars_ ) {
                if (auto* ptr = std::any_cast<bool>(&(element.second))) {
                    *ptr = false;
                    }
                else if (auto* ptr = std::any_cast<int>(&(element.second))) {
                    *ptr = 0;
                    }
                else if (auto* ptr = std::any_cast<float>(&(element.second))) {
                    *ptr = 0.f;
                    }
                else {
                    throw cms::Exception("Unknown")
                        << "Tried to clear "
                        << element.first
                        << ", but type is unknown";
                    }
                }
                }
            }
    };

hgcalfinecalontupler::hgcalfinecalontupler(const edm::ParameterSet& iConfig) : 
    tree_(NULL),
    hgcalEEHitsToken_(consumes<edm::View<PCaloHit>>(edm::InputTag("g4SimHits", "HGCHitsEE"))),
    hgcalHEfrontHitsToken_(consumes<edm::View<PCaloHit>>(edm::InputTag("g4SimHits", "HGCHitsHEfront"))),
    hgcalHEbackHitsToken_(consumes<edm::View<PCaloHit>>(edm::InputTag("g4SimHits", "HGCHitsHEback"))),
    tokenSimTracks(consumes<edm::SimTrackContainer>(edm::InputTag("g4SimHits"))),
    tokenSimVertices(consumes<edm::SimVertexContainer>(edm::InputTag("g4SimHits")))
    {
        usesResource("TFileService");
        }

void hgcalfinecalontupler::beginJob() {
    tree_ = fs->make<TTree>("tree", "tree");
    makeScalarBranch<int>("event_number");

    makeBranch<float>("simtrack_x");
    makeBranch<float>("simtrack_y");
    makeBranch<float>("simtrack_z");
    makeBranch<float>("simtrack_pt");
    makeBranch<float>("simtrack_eta");
    makeBranch<float>("simtrack_phi");
    makeBranch<float>("simtrack_energy");
    makeBranch<float>("simtrack_mass");
    makeBranch<int>("simtrack_trackid");
    makeBranch<int>("simtrack_vertexindex");
    makeBranch<int>("simtrack_pdgid");
    makeBranch<bool>("simtrack_crossedboundary");
    makeBranch<float>("simtrack_boundary_x");
    makeBranch<float>("simtrack_boundary_y");
    makeBranch<float>("simtrack_boundary_z");
    makeBranch<float>("simtrack_boundary_t");
    makeBranch<float>("simtrack_boundary_pt");
    makeBranch<float>("simtrack_boundary_eta");
    makeBranch<float>("simtrack_boundary_phi");
    makeBranch<float>("simtrack_boundary_energy");
    makeBranch<float>("simtrack_boundary_mass");
    makeBranch<int>("simtrack_vertex_id");
    makeBranch<float>("simtrack_vertex_x");
    makeBranch<float>("simtrack_vertex_y");
    makeBranch<float>("simtrack_vertex_z");
    makeBranch<float>("simtrack_vertex_t");
    makeBranch<int>("simtrack_vertex_processtype");
    makeBranch<int>("simtrack_parenttrackid");
    makeBranch<bool>("simtrack_noparent");
    makeBranch<bool>("simtrack_parentexists");
    makeBranch<bool>("simtrack_hashits");

    makeBranch<int>("simhit_detid");
    makeBranch<float>("simhit_x");
    makeBranch<float>("simhit_y");
    makeBranch<float>("simhit_z");
    makeBranch<int>("simhit_layer");
    makeBranch<float>("simhit_energy");
    makeBranch<int>("simhit_pdgid");
    makeBranch<float>("simhit_emenergy");
    makeBranch<float>("simhit_time");
    makeBranch<int>("simhit_trackid");
    makeBranch<bool>("simhit_inEE");
    makeBranch<bool>("simhit_inHSi");
    makeBranch<bool>("simhit_inHsc");
    }

void hgcalfinecalontupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    clear();
    edm::ESHandle<CaloGeometry> geom;
    iSetup.get<CaloGeometryRecord>().get(geom);
    hgcalRecHitToolInstance_.setGeometry(*geom);

    fillScalar<int>("event_number", iEvent.id().event());

    edm::Handle<edm::SimVertexContainer> handleSimVertices;
    iEvent.getByLabel("g4SimHits", handleSimVertices);
    edm::Handle<edm::SimTrackContainer> handleSimTracks;
    iEvent.getByLabel("g4SimHits", handleSimTracks);

    // Fill a set with all available track ids to check for broken parentage
    std::set<int> all_track_ids;
    std::set<int> all_track_ids_that_crossed_boundary;
    unordered_map<int, int> track_id_to_pdgid;
    for(const auto& track : *(handleSimTracks.product())){
        all_track_ids.insert(track.trackId());
        track_id_to_pdgid.emplace(track.trackId(), track.type());
        if (track.crossedBoundary()) all_track_ids_that_crossed_boundary.insert(track.trackId());
        }
    std::set<int> all_track_ids_with_hits;

    // Store the hits
    std::vector<edm::EDGetTokenT<edm::View<PCaloHit>>> tokens = {
        hgcalEEHitsToken_,
        hgcalHEfrontHitsToken_,
        hgcalHEbackHitsToken_
        };
    for (edm::EDGetTokenT<edm::View<PCaloHit>> token : tokens ) {
        edm::Handle< edm::View<PCaloHit> > handle;
        iEvent.getByToken(token, handle);
        for (auto const & hit : handle->ptrs()) {
            DetId id = hit->id();
            fill<int>("simhit_detid", id.rawId());
            GlobalPoint position = hgcalRecHitToolInstance_.getPosition(id);
            fill<float>("simhit_x", position.x());
            fill<float>("simhit_y", position.y());
            fill<float>("simhit_z", position.z());
            fill<int>("simhit_layer", hgcalRecHitToolInstance_.getLayer(id));
            fill<float>("simhit_energy", hit->energy());
            fill<float>("simhit_emenergy", hit->energyEM());
            fill<float>("simhit_time", hit->time());
            fill<int>("simhit_trackid", hit->geantTrackId());
            fill<bool>("simhit_inEE", id.det() == DetId::HGCalEE);
            fill<bool>("simhit_inHSi", id.det() == DetId::HGCalHSi);
            fill<bool>("simhit_inHsc", id.det() == DetId::HGCalHSc);
            // Check whether the parent track exists, and if so get its pdgid
            if (!all_track_ids.count(hit->geantTrackId())){
                edm::LogError("DoFineCalo")
                    << "Event " << iEvent.id().event()
                    << ": Hit " << id.rawId()
                    << " has parent " << hit->geantTrackId()
                    << ", which has NOT been saved!";
                fill<int>("simhit_pdgid", 0);
                }
            else{
                fill<int>("simhit_pdgid", track_id_to_pdgid[hit->geantTrackId()]);
                }
            // Check whether the parent track crossed the boundary
            // (it must by definition for finecalo volumes)
            if (!all_track_ids_that_crossed_boundary.count(hit->geantTrackId())){
                edm::LogError("DoFineCalo")
                    << "Event " << iEvent.id().event()
                    << ": Hit " << id.rawId()
                    << " has parent " << hit->geantTrackId()
                    << ", which did NOT cross the boundary!";
                }
            all_track_ids_with_hits.insert(hit->geantTrackId());
            }
        }

    // Store the tracks
    for(const auto& track : *(handleSimTracks.product())){
        fill<int>("simtrack_trackid", track.trackId());
        fill<float>("simtrack_x", track.trackerSurfacePosition().X());
        fill<float>("simtrack_y", track.trackerSurfacePosition().Y());
        fill<float>("simtrack_z", track.trackerSurfacePosition().Z());
        const math::XYZTLorentzVectorD momentum = track.momentum();
        fill<float>("simtrack_pt", momentum.Pt());
        fill<float>("simtrack_eta", momentum.Eta());
        fill<float>("simtrack_phi", momentum.Phi());
        fill<float>("simtrack_energy", momentum.E());
        fill<float>("simtrack_mass", momentum.M());
        fill<int>("simtrack_vertexindex", track.vertIndex());
        fill<int>("simtrack_pdgid", track.type());
        fill<bool>("simtrack_crossedboundary", track.crossedBoundary());
        if (track.crossedBoundary()){
            const math::XYZTLorentzVectorF boundaryPosition = track.getPositionAtBoundary();
            fill<float>("simtrack_boundary_x", boundaryPosition.X());
            fill<float>("simtrack_boundary_y", boundaryPosition.Y());
            fill<float>("simtrack_boundary_z", boundaryPosition.Z());
            fill<float>("simtrack_boundary_t", boundaryPosition.T());
            const math::XYZTLorentzVectorF boundaryMomentum = track.getMomentumAtBoundary();
            fill<float>("simtrack_boundary_pt", boundaryMomentum.Pt());
            fill<float>("simtrack_boundary_eta", boundaryMomentum.Eta());
            fill<float>("simtrack_boundary_phi", boundaryMomentum.Phi());
            fill<float>("simtrack_boundary_energy", boundaryMomentum.E());
            fill<float>("simtrack_boundary_mass", boundaryMomentum.M());
            }
        else{
            // Better to have defaults for selections with uproot
            fill<float>("simtrack_boundary_x", 0.);
            fill<float>("simtrack_boundary_y", 0.);
            fill<float>("simtrack_boundary_z", 0.);
            fill<float>("simtrack_boundary_t", 0.);
            fill<float>("simtrack_boundary_pt", 0.);
            fill<float>("simtrack_boundary_eta", 0.);
            fill<float>("simtrack_boundary_phi", 0.);
            fill<float>("simtrack_boundary_energy", 0.);
            fill<float>("simtrack_boundary_mass", 0.);
            }
        fill<bool>("simtrack_hashits", all_track_ids_with_hits.count(track.trackId()));
        SimVertex vertex = handleSimVertices.product()->at(track.vertIndex());
        fill<int>("simtrack_vertex_id", vertex.vertexId());
        const math::XYZTLorentzVectorD position = vertex.position();
        fill<float>("simtrack_vertex_x", position.X());
        fill<float>("simtrack_vertex_y", position.Y());
        fill<float>("simtrack_vertex_z", position.Z());
        fill<float>("simtrack_vertex_t", position.T());
        fill<int>("simtrack_vertex_processtype", vertex.processType());
        fill<int>("simtrack_parenttrackid", vertex.parentIndex());
        fill<bool>("simtrack_noparent", vertex.noParent());
        // Check whether the parent track exists
        bool parentExists = (vertex.noParent() or all_track_ids.count(vertex.parentIndex()));
        fill<bool>("simtrack_parentexists", parentExists);
        if (!parentExists){
            edm::LogError("DoFineCalo")
                << "Event " << iEvent.id().event()
                << ": Track " << track.trackId()
                << " has parent " << vertex.parentIndex()
                << ", which has NOT been saved!";
            }
        }

    tree_->Fill();
    }

void hgcalfinecalontupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("SimTrackTag", edm::InputTag("g4SimHits"));
    desc.add<edm::InputTag>("SimVertexTag", edm::InputTag("g4SimHits"));
    descriptions.add("hgcalfinecalontupler", desc);
    }

DEFINE_FWK_MODULE(hgcalfinecalontupler);