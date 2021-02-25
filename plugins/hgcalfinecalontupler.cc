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

struct Hits{
    vector<unsigned int> detid_;
    vector<float> x_;
    vector<float> y_;
    vector<float> z_;
    vector<float> eta_;
    vector<float> phi_;
    vector<float> pt_;
    vector<int> zside_;
    vector<float> siThickness_;
    vector<float> radiusToSide_;
    vector<int> siThickIndex_;
    vector<int> layer_;
    vector<int> layerWithOffset_;
    vector<std::pair<int, int>> wafer_;
    vector<std::pair<int, int>> cell_;
    vector<bool> isHalfCell_;
    vector<bool> isSilicon_;
    vector<bool> inEE_;
    vector<bool> inHsi_;
    vector<bool> inHsc_;
    vector<float> energy_;
    vector<float> emFraction_;
    vector<float> time_;
    vector<int> trackId_;
    vector<int> depth_;
    
    void linkToTree(TTree* tree, std::string name){
        tree->Branch((name + "_detid").c_str(), "vector<unsigned int>", &detid_, 32000, 0);
        tree->Branch((name + "_x").c_str(), "vector<float>", &x_, 32000, 0);
        tree->Branch((name + "_y").c_str(), "vector<float>", &y_, 32000, 0);
        tree->Branch((name + "_z").c_str(), "vector<float>", &z_, 32000, 0);
        tree->Branch((name + "_eta").c_str(), "vector<float>", &eta_, 32000, 0);
        tree->Branch((name + "_phi").c_str(), "vector<float>", &phi_, 32000, 0);
        tree->Branch((name + "_pt").c_str(), "vector<float>", &pt_, 32000, 0);
        tree->Branch((name + "_zside").c_str(), "vector<int>", &zside_, 32000, 0);
        tree->Branch((name + "_siThickness").c_str(), "vector<float>", &siThickness_, 32000, 0);
        tree->Branch((name + "_radiusToSide").c_str(), "vector<float>", &radiusToSide_, 32000, 0);
        tree->Branch((name + "_siThickIndex").c_str(), "vector<int>", &siThickIndex_, 32000, 0);
        tree->Branch((name + "_layer").c_str(), "vector<int>", &layer_, 32000, 0);
        tree->Branch((name + "_layerWithOffset").c_str(), "vector<int>", &layerWithOffset_, 32000, 0);
        tree->Branch((name + "_wafer").c_str(), "vector<pair<int, int>>", &wafer_, 32000, 0);
        tree->Branch((name + "_cell").c_str(), "vector<pair<int, int>>", &cell_, 32000, 0);
        tree->Branch((name + "_isHalfCell").c_str(), "vector<bool>", &isHalfCell_, 32000, 0);
        tree->Branch((name + "_isSilicon").c_str(), "vector<bool>", &isSilicon_, 32000, 0);
        tree->Branch((name + "_inEE").c_str(), "vector<bool>", &inEE_, 32000, 0);
        tree->Branch((name + "_inHsi").c_str(), "vector<bool>", &inHsi_, 32000, 0);
        tree->Branch((name + "_inHsc").c_str(), "vector<bool>", &inHsc_, 32000, 0);
        tree->Branch((name + "_energy").c_str(), "vector<float>", &energy_, 32000, 0);
        tree->Branch((name + "_emFraction").c_str(), "vector<float>", &emFraction_, 32000, 0);
        tree->Branch((name + "_time").c_str(), "vector<float>", &time_, 32000, 0);
        tree->Branch((name + "_trackId").c_str(), "vector<int>", &trackId_, 32000, 0);
        tree->Branch((name + "_depth").c_str(), "vector<int>", &depth_, 32000, 0);
        }

    void fill(const PCaloHit & hit, const hgcal::RecHitTools * hgcalRecHitToolInstance){
        DetId id = hit.id();
        detid_.push_back(id.rawId());

        GlobalPoint position = hgcalRecHitToolInstance->getPosition(id);
        x_.push_back(position.x());
        y_.push_back(position.y());
        z_.push_back(position.z());

        zside_.push_back(hgcalRecHitToolInstance->zside(id));
        siThickness_.push_back(hgcalRecHitToolInstance->getSiThickness(id));
        siThickIndex_.push_back(hgcalRecHitToolInstance->getSiThickIndex(id));
        layer_.push_back(hgcalRecHitToolInstance->getLayer(id));
        layerWithOffset_.push_back(hgcalRecHitToolInstance->getLayerWithOffset(id));
        isHalfCell_.push_back(hgcalRecHitToolInstance->isHalfCell(id));
        isSilicon_.push_back(hgcalRecHitToolInstance->isSilicon(id));
        eta_.push_back(hgcalRecHitToolInstance->getEta(id));
        phi_.push_back(hgcalRecHitToolInstance->getPhi(id));
        pt_.push_back(hgcalRecHitToolInstance->getPt(id, hit.energy()));

        bool inEE = false;
        bool inHsi = false;
        bool inHsc = false;
        if (id.det() == DetId::HGCalEE){
            inEE = true;
            }
        else if (id.det() == DetId::HGCalHSi){
            inHsi = true;
            }
        else if (id.det() == DetId::HGCalHSc){
            inHsc = true;
            }
        inEE_.push_back(inEE);
        inHsi_.push_back(inHsi);
        inHsc_.push_back(inHsc);

        if (inEE || inHsi){
            wafer_.push_back(hgcalRecHitToolInstance->getWafer(id));
            cell_.push_back(hgcalRecHitToolInstance->getCell(id));
            radiusToSide_.push_back(hgcalRecHitToolInstance->getRadiusToSide(id));
            }
        else{
            wafer_.push_back(std::pair<int, int>(-1, -1));
            cell_.push_back(std::pair<int, int>(-1, -1));
            radiusToSide_.push_back(-1.);
            }

        energy_.push_back(hit.energy());
        emFraction_.push_back(hit.energyEM() / hit.energy());
        time_.push_back(hit.time());
        trackId_.push_back(hit.geantTrackId());
        depth_.push_back(hit.depth());
        }
    };

struct Tracks {
    vector<float> x_;
    vector<float> y_;
    vector<float> z_;
    vector<TLorentzVector> momentum_;
    vector<int> trackId_;
    vector<int> vertexIndex_;
    vector<int> pdgid_;
    vector<bool> crossedBoundary_;
    vector<int> idAtBoundary_;
    vector<float> xAtBoundary_;
    vector<float> yAtBoundary_;
    vector<float> zAtBoundary_;
    vector<TLorentzVector> momentumAtBoundary_;
    // vector<TLorentzVector> correctedMomentumAtBoundary_;

    void fill(const SimTrack & track){
        x_.push_back(track.trackerSurfacePosition().X());
        y_.push_back(track.trackerSurfacePosition().Y());
        z_.push_back(track.trackerSurfacePosition().Z());
        momentum_.push_back(TLorentzVector(
            track.momentum().Px(), track.momentum().Py(), track.momentum().Pz(), track.momentum().E()
            ));
        trackId_.push_back(track.trackId());
        vertexIndex_.push_back(track.vertIndex());
        pdgid_.push_back(track.type());

        crossedBoundary_.push_back(track.crossedBoundary());
        if (track.crossedBoundary()){
            idAtBoundary_.push_back(track.getIDAtBoundary());
            xAtBoundary_.push_back(track.getPositionAtBoundary().X());
            yAtBoundary_.push_back(track.getPositionAtBoundary().Y());
            zAtBoundary_.push_back(track.getPositionAtBoundary().Z());
            momentumAtBoundary_.push_back(TLorentzVector(
                track.getMomentumAtBoundary().Px(), track.getMomentumAtBoundary().Py(), track.getMomentumAtBoundary().Pz(), track.getMomentumAtBoundary().E()
                ));
            // correctedMomentumAtBoundary_.push_back(TLorentzVector(
            //     track.getCorrectedMomentumAtBoundary().Px(), track.getCorrectedMomentumAtBoundary().Py(),
            //     track.getCorrectedMomentumAtBoundary().Pz(), track.getCorrectedMomentumAtBoundary().E()
            //     ));
            }
        else{
            idAtBoundary_.push_back(-1);
            xAtBoundary_.push_back(-1);
            yAtBoundary_.push_back(-1);
            zAtBoundary_.push_back(-1);
            momentumAtBoundary_.push_back(TLorentzVector(0,0,0,0));
            // correctedMomentumAtBoundary_.push_back(TLorentzVector(0,0,0,0));
            }
        }

    void linkToTree(TTree* tree, std::string name){
        tree->Branch((name + "_x").c_str(), "vector<float>", &x_, 32000, 0);
        tree->Branch((name + "_y").c_str(), "vector<float>", &y_, 32000, 0);
        tree->Branch((name + "_z").c_str(), "vector<float>", &z_, 32000, 0);
        tree->Branch((name + "_momentum").c_str(), "vector<TLorentzVector>", &momentum_, 32000, 0);
        tree->Branch((name + "_trackId").c_str(), "vector<int>", &trackId_, 32000, 0);
        tree->Branch((name + "_vertexIndex").c_str(), "vector<int>", &vertexIndex_, 32000, 0);
        tree->Branch((name + "_pdgid").c_str(), "vector<int>", &pdgid_, 32000, 0);
        tree->Branch((name + "_crossedBoundary").c_str(), "vector<bool>", &crossedBoundary_, 32000, 0);
        tree->Branch((name + "_idAtBoundary").c_str(), "vector<int>", &idAtBoundary_, 32000, 0);
        tree->Branch((name + "_xAtBoundary").c_str(), "vector<float>", &xAtBoundary_, 32000, 0);
        tree->Branch((name + "_yAtBoundary").c_str(), "vector<float>", &yAtBoundary_, 32000, 0);
        tree->Branch((name + "_zAtBoundary").c_str(), "vector<float>", &zAtBoundary_, 32000, 0);
        tree->Branch((name + "_momentumAtBoundary").c_str(), "vector<TLorentzVector>", &momentumAtBoundary_, 32000, 0);
        // tree->Branch((name + "_correctedMomentumAtBoundary").c_str(), "vector<TLorentzVector>", &correctedMomentumAtBoundary_, 32000, 0);
        }
    };

struct Vertices {
    vector<float> x_;
    vector<float> y_;
    vector<float> z_;
    vector<float> t_;
    vector<int> id_;
    vector<int> parentTrackId_;
    vector<int> processType_;
    vector<bool> noParent_;

    void fill(const SimVertex & vertex){
        x_.push_back(vertex.position().X());
        y_.push_back(vertex.position().Y());
        z_.push_back(vertex.position().Z());
        t_.push_back(vertex.position().T());
        id_.push_back(vertex.vertexId());
        parentTrackId_.push_back(vertex.parentIndex());
        processType_.push_back(vertex.processType());
        noParent_.push_back(vertex.noParent());
        }

    void linkToTree(TTree* tree, std::string name){
        tree->Branch((name + "_x").c_str(), "vector<float>", &x_, 32000, 0);
        tree->Branch((name + "_y").c_str(), "vector<float>", &y_, 32000, 0);
        tree->Branch((name + "_z").c_str(), "vector<float>", &z_, 32000, 0);
        tree->Branch((name + "_t").c_str(), "vector<float>", &t_, 32000, 0);
        tree->Branch((name + "_id").c_str(), "vector<int>", &id_, 32000, 0);
        tree->Branch((name + "_parentTrackId").c_str(), "vector<int>", &parentTrackId_, 32000, 0);
        tree->Branch((name + "_processType").c_str(), "vector<int>", &processType_, 32000, 0);
        tree->Branch((name + "_noParent").c_str(), "vector<bool>", &noParent_, 32000, 0);
        }
    };

struct Entry {
    int event_number;
    Hits hits;
    Tracks tracks;
    Vertices vertices;
    void linkToTree(TTree* tree){
        tree->Branch("event_number", &event_number);
        hits.linkToTree(tree, "simhit");
        tracks.linkToTree(tree, "simtrack");
        vertices.linkToTree(tree, "simvertex");
        }
    };


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
        Entry entry;
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
    makeBranch<int>("simtrack_idatboundary");
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
    for(const auto& track : *(handleSimTracks.product())){
        all_track_ids.insert(track.trackId());
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
            // Check whether the parent track exists
            if (!all_track_ids.count(hit->geantTrackId())){
                edm::LogError("DoFineCalo")
                    << "Event " << iEvent.id().event()
                    << ": Hit " << id.rawId()
                    << " has parent " << hit->geantTrackId()
                    << ", which has NOT been saved!";
                }
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