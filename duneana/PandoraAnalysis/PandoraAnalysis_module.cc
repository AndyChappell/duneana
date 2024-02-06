////////////////////////////////////////////////////////////////////////
// Class:       pandoraAnalysis
// Plugin Type: analyzer (art v3_05_01)
// File:        pandoraAnalysis_module.cc
//
// Generated at Fri Aug  7 15:01:35 2020 by Maria Brigida Brunetti using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/ArtDataHelper/TrackUtils.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "larsim/MCCheater/BackTrackerService.h"

#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dunereco/AnaUtils/DUNEAnaTrackUtils.h"
#include "dunereco/AnaUtils/DUNEAnaShowerUtils.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include <TTree.h>
#include <vector>
#include <string>
#include <TH1D.h>
#include <TLorentzVector.h>


namespace test {
    class pandoraAnalysis;
}


class test::pandoraAnalysis : public art::EDAnalyzer {
    public:
        explicit pandoraAnalysis(fhicl::ParameterSet const& p);
        // The compiler-generated destructor is fine for non-base
        // classes without bare pointers or other resource use.

        // Plugins should not be copied or assigned.
        pandoraAnalysis(pandoraAnalysis const&) = delete;
        pandoraAnalysis(pandoraAnalysis&&) = delete;
        pandoraAnalysis& operator=(pandoraAnalysis const&) = delete;
        pandoraAnalysis& operator=(pandoraAnalysis&&) = delete;

        // Required functions.
        void analyze(art::Event const& e) override;

        // Selected optional functions.
        void beginJob() override;
        void endJob() override;

        void reset(bool deepClean=false);

    private:

        // Declare member data here.

        TTree *fTree;

        unsigned int fEventID;
        unsigned int fRunID;
        unsigned int fSubRunID;

        unsigned int fNMCParticles;
        unsigned int fNMCReconstructible;
        unsigned int fNPFParticles;

        static const int kNMaxMCParticles = 20000;
        static const int kNMaxPFParticles = 2000;
        static const int kNMaxPFPClusters = 100;
        static const int kNViews = 3;


        bool fMCIsPrimary[kNMaxMCParticles];
        int fMCParticlePdgCode[kNMaxMCParticles];
        double fMCParticleTrueEnergy [kNMaxMCParticles];
        double fMCParticleTrueMom [kNMaxMCParticles]; 
        int fMCParticleTrackID[kNMaxMCParticles]; 
        int fMCParticleParentTrackID[kNMaxMCParticles]; 
        double fMCParticleStartPositionX[kNMaxMCParticles];
        double fMCParticleStartPositionY[kNMaxMCParticles];
        double fMCParticleStartPositionZ[kNMaxMCParticles];
        double fMCParticleStartPositionT[kNMaxMCParticles];
        double fMCParticleStartMomentumX[kNMaxMCParticles];
        double fMCParticleStartMomentumY[kNMaxMCParticles];
        double fMCParticleStartMomentumZ[kNMaxMCParticles];
        double fMCParticleStartMomentumE[kNMaxMCParticles];
        double fMCParticleEndPositionX[kNMaxMCParticles];
        double fMCParticleEndPositionY[kNMaxMCParticles];
        double fMCParticleEndPositionZ[kNMaxMCParticles];
        double fMCParticleEndPositionT[kNMaxMCParticles];
        double fMCParticleEndMomentumX[kNMaxMCParticles];
        double fMCParticleEndMomentumY[kNMaxMCParticles];
        double fMCParticleEndMomentumZ[kNMaxMCParticles];
        double fMCParticleEndMomentumE[kNMaxMCParticles];
        int fMCParticleNHits[kNMaxMCParticles];
        int fMCParticleNHitsView[kNMaxMCParticles][kNViews];

        int fPFPID[kNMaxPFParticles];
        bool fPFPIsPrimary[kNMaxPFParticles];
        int fPFPTrueParticleMatchedID[kNMaxPFParticles];
        int fPFPTrueParticleMatchedPosition[kNMaxPFParticles];
        int fPFPParentID[kNMaxPFParticles];
        int fPFPPdgCode[kNMaxPFParticles];
        int fPFPNChildren[kNMaxPFParticles];
        int fPFPNClusters[kNMaxPFParticles];
        int fPFPNHits[kNMaxPFParticles];
        int fPFPNHitsView[kNMaxPFParticles][kNViews];
        int fPFPNSharedTrueParticleHits[kNMaxPFParticles];
        int fPFPNSharedTrueParticleHitsView[kNMaxPFParticles][kNViews];
        bool fPFPIsTrack[kNMaxPFParticles];
        bool fPFPIsShower[kNMaxPFParticles];


        int fPFPTrackID[kNMaxPFParticles];
        double fPFPTrackLength[kNMaxPFParticles];
        double fPFPTrackStartX[kNMaxPFParticles];
        double fPFPTrackStartY[kNMaxPFParticles];
        double fPFPTrackStartZ[kNMaxPFParticles];
        double fPFPTrackVertexX[kNMaxPFParticles];
        double fPFPTrackVertexY[kNMaxPFParticles];
        double fPFPTrackVertexZ[kNMaxPFParticles];
        double fPFPTrackEndX[kNMaxPFParticles];
        double fPFPTrackEndY[kNMaxPFParticles];
        double fPFPTrackEndZ[kNMaxPFParticles];
        double fPFPTrackTheta[kNMaxPFParticles];
        double fPFPTrackPhi[kNMaxPFParticles];
        double fPFPTrackZenithAngle[kNMaxPFParticles];
        double fPFPTrackAzimuthAngle[kNMaxPFParticles];
        double fPFPTrackStartDirectionX[kNMaxPFParticles];
        double fPFPTrackStartDirectionY[kNMaxPFParticles];
        double fPFPTrackStartDirectionZ[kNMaxPFParticles];
        double fPFPTrackEndDirectionX[kNMaxPFParticles];
        double fPFPTrackEndDirectionY[kNMaxPFParticles];
        double fPFPTrackEndDirectionZ[kNMaxPFParticles];

        int    fPFPShowerID[kNMaxPFParticles];
        double fPFPShowerEnergy[kNMaxPFParticles];
        double fPFPShowerEnergyToTrueEnergyRatio[kNMaxPFParticles];
        double fPFPShowerEnergyToTrueMomentumRatio[kNMaxPFParticles];
        double fPFPShowerDirectionX[kNMaxPFParticles];
        double fPFPShowerDirectionY[kNMaxPFParticles];
        double fPFPShowerDirectionZ[kNMaxPFParticles];
        double fPFPShowerStartX[kNMaxPFParticles];
        double fPFPShowerStartY[kNMaxPFParticles];
        double fPFPShowerStartZ[kNMaxPFParticles];
        double fPFPShowerLength[kNMaxPFParticles];
        double fPFPShowerOpenAngle[kNMaxPFParticles];

        double fPFPCompleteness[kNMaxMCParticles];
        double fPFPCompletenessView[kNMaxMCParticles][kNViews];
        double fPFPPurity[kNMaxMCParticles];
        double fPFPPurityView[kNMaxMCParticles][kNViews];

        std::string fTruthLabel;
        std::string fHitLabel;
        std::string fTrackLabel;
        std::string fShowerLabel;
        std::string fPFParticleLabel;
        bool fRollUpUnsavedIDs;

        const geo::Geometry* fGeom;
};


test::pandoraAnalysis::pandoraAnalysis(fhicl::ParameterSet const& p)
    : EDAnalyzer{p}
{
    fTruthLabel = p.get<std::string>("TruthLabel");
    fHitLabel = p.get<std::string>("HitLabel");
    fPFParticleLabel = p.get<std::string>("PFParticleLabel");
    fTrackLabel = p.get<std::string>("TrackLabel");
    fShowerLabel = p.get<std::string>("ShowerLabel");
    fGeom    = &*art::ServiceHandle<geo::Geometry>();
    fRollUpUnsavedIDs = p.get<bool>("RollUpUnsavedIDs"); 
} 

void test::pandoraAnalysis::analyze(art::Event const& e)
{
    reset();
    fEventID = e.id().event();
    fRunID = e.id().run();
    fSubRunID = e.id().subRun();
    auto const clockData{art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e)};
    std::cout << "=============== EVENT ID " << fEventID << " == RUN ID " << fRunID << " == SUBRUN ID " << fSubRunID << " ================" << std::endl;

    //Get all hits
    std::vector<art::Ptr<recob::Hit> > allHits;
    auto hitHandle{e.getHandle<std::vector<recob::Hit>>(fHitLabel)};
    if (hitHandle)
        art::fill_ptr_vector(allHits, hitHandle);

    //Fill MC particle to hits map
    std::map<int,int> trueMCHits, trueMCHitsU, trueMCHitsV, trueMCHitsW;
    std::map<simb::MCParticle, bool> reconstructibleMC;
    std::map<int, simb::MCParticle> g4ToMCMap;
    std::map<simb::MCParticle, int> mcToPosMap;
    for (const auto& hit: allHits)
    {
        TruthMatchUtils::G4ID g4ID(TruthMatchUtils::TrueParticleID(clockData, hit, fRollUpUnsavedIDs));
        if (TruthMatchUtils::Valid(g4ID))
        {
            ++trueMCHits[g4ID];
            if (hit->View() == 0)
                ++trueMCHitsU[g4ID];
            if (hit->View() == 1)
                ++trueMCHitsV[g4ID];
            if (hit->View() == 2)
                ++trueMCHitsW[g4ID];
        }
    }

    if (!e.isRealData())
    {
        art::ValidHandle<std::vector<simb::MCParticle>> mcParticles{e.getValidHandle<std::vector<simb::MCParticle>>(fTruthLabel)};
        if (mcParticles.isValid())
        {
            for(unsigned int i = 0; i < mcParticles->size(); ++i)
            {
                const simb::MCParticle &mc{mcParticles->at(i)};
                unsigned int nGoodViews{0};
                nGoodViews += trueMCHitsU[mc.TrackId()] >= 5 ? 1 : 0;
                nGoodViews += trueMCHitsV[mc.TrackId()] >= 5 ? 1 : 0;
                nGoodViews += trueMCHitsW[mc.TrackId()] >= 5 ? 1 : 0;
                if (trueMCHits[mc.TrackId()] > 0 && nGoodViews >= 2)
                    reconstructibleMC[mc] = true;
            }
            fNMCParticles = reconstructibleMC.size();

            unsigned int i{0};
            for(const auto &[mc, dummy] : reconstructibleMC)
            {
                g4ToMCMap[mc.TrackId()] = mc;
                mcToPosMap[mc] = i;
                fMCParticleTrueEnergy[i] = mc.E();
                fMCParticleTrueMom[i] = mc.P();
                fMCParticlePdgCode[i] = mc.PdgCode();
                fMCParticleTrackID[i] = mc.TrackId();
                fMCParticleParentTrackID[i] = mc.Mother();
                fMCIsPrimary[i] = mc.Process() == "primary";;
                fMCParticleNHits[i] = trueMCHits[mc.TrackId()];
                fMCParticleStartPositionX[i] = mc.Position().X();
                fMCParticleStartPositionY[i] = mc.Position().Y();
                fMCParticleStartPositionZ[i] = mc.Position().Z();
                fMCParticleStartPositionT[i] = mc.Position().T();
                fMCParticleEndPositionX[i] = mc.EndPosition().X();
                fMCParticleEndPositionY[i] = mc.EndPosition().Y();
                fMCParticleEndPositionZ[i] = mc.EndPosition().Z();
                fMCParticleEndPositionT[i] = mc.EndPosition().T();
                fMCParticleStartMomentumX[i] = mc.Momentum().X();
                fMCParticleStartMomentumY[i] = mc.Momentum().Y();
                fMCParticleStartMomentumZ[i] = mc.Momentum().Z();
                fMCParticleStartMomentumE[i] = mc.Momentum().E();
                fMCParticleEndMomentumX[i] = mc.EndMomentum().X();
                fMCParticleEndMomentumY[i] = mc.EndMomentum().Y();
                fMCParticleEndMomentumZ[i] = mc.EndMomentum().Z();
                fMCParticleEndMomentumE[i] = mc.EndMomentum().E();
                ++i;
            }
        }
    }

    const std::vector<art::Ptr<recob::PFParticle>> pfParticles{dune_ana::DUNEAnaEventUtils::GetPFParticles(e, fPFParticleLabel)};
    fNPFParticles = pfParticles.size();
    if (!fNPFParticles)
    {
        std::cout << "No PFParticles found!" << std::endl;
        return;
    }

    std::vector<art::Ptr<recob::Cluster>> clusters;
    auto clusterHandle{e.getHandle<std::vector<recob::Cluster> >(fPFParticleLabel)};
    if (clusterHandle)
        art::fill_ptr_vector(clusters, clusterHandle);
    art::FindManyP<recob::Cluster> clusterParticleAssoc(pfParticles, e, fPFParticleLabel);

    auto trackHandle{e.getHandle<std::vector<recob::Track> >(fTrackLabel)};
    if (!trackHandle)
    {
        std::cout<<"Unable to find std::vector<recob::Track> with module label: " << fTrackLabel << std::endl;
        return;
    }
    std::vector<art::Ptr<recob::Track>> tracks;
    art::fill_ptr_vector(tracks, trackHandle);

    unsigned int p{0};
    for(const art::Ptr<recob::PFParticle> &pfp: pfParticles)
    {
        std::vector<art::Ptr<recob::Hit>> pfpHits;
        if(dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfp, e, fPFParticleLabel, fTrackLabel))
        {
            art::Ptr<recob::Track> track{dune_ana::DUNEAnaPFParticleUtils::GetTrack(pfp, e, fPFParticleLabel, fTrackLabel)};
            pfpHits = dune_ana::DUNEAnaTrackUtils::GetHits(track, e, fTrackLabel);
        }
        else if (dune_ana::DUNEAnaPFParticleUtils::IsShower(pfp, e, fPFParticleLabel, fShowerLabel))
        {
            art::Ptr<recob::Shower> shower = dune_ana::DUNEAnaPFParticleUtils::GetShower(pfp, e, fPFParticleLabel, fShowerLabel);
            pfpHits = dune_ana::DUNEAnaShowerUtils::GetHits(shower, e, fShowerLabel);
        }

        if(!e.isRealData())
        {
            TruthMatchUtils::G4ID g4ID(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, pfpHits,fRollUpUnsavedIDs));
            if (TruthMatchUtils::Valid(g4ID))
            {
                // Only consider PFPs where the matched MC was deemed reconstructible
                if (g4ToMCMap.find(g4ID) == g4ToMCMap.end())
                    continue;
                fPFPTrueParticleMatchedID[p] = g4ID;
                fPFPTrueParticleMatchedPosition[p] = mcToPosMap[g4ToMCMap[g4ID]];
                fPFPNHits[p] = pfpHits.size();
            }
        }

        fPFPID[p] = p;
        fPFPIsPrimary[p] = pfp->IsPrimary();
        fPFPPdgCode[p] = pfp->PdgCode();
        fPFPNChildren[p] = pfp->NumDaughters();
        fPFPParentID[p] = pfp->IsPrimary() ? -1 : pfp->Parent();

        std::vector<art::Ptr<recob::Cluster>> pfpClusters{clusterParticleAssoc.at(pfp.key())};
        fPFPNClusters[p] = pfpClusters.size();

        if(dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfp, e, fPFParticleLabel, fTrackLabel))
        {
            art::Ptr<recob::Track> track{dune_ana::DUNEAnaPFParticleUtils::GetTrack(pfp, e, fPFParticleLabel, fTrackLabel)};
            fPFPIsTrack[p] = true;
            fPFPTrackID[p] = track->ID();
            fPFPTrackLength[p] = track->Length();
            fPFPTrackStartX[p] = track->Start().X();
            fPFPTrackStartY[p] = track->Start().Y();
            fPFPTrackStartZ[p] = track->Start().Z();
            fPFPTrackVertexX[p] = track->Vertex().X();
            fPFPTrackVertexY[p] = track->Vertex().Y();
            fPFPTrackVertexZ[p] = track->Vertex().Z();
            fPFPTrackEndX[p] = track->End().X();
            fPFPTrackEndY[p] = track->End().Y();
            fPFPTrackEndZ[p] = track->End().Z();
            fPFPTrackTheta[p] = track->Theta();
            fPFPTrackPhi[p] = track->Phi();
            fPFPTrackZenithAngle[p] = track->ZenithAngle();
            fPFPTrackAzimuthAngle[p] = track->AzimuthAngle();
            fPFPTrackStartDirectionX[p] = track->StartDirection().X();
            fPFPTrackStartDirectionY[p] = track->StartDirection().Y();
            fPFPTrackStartDirectionZ[p] = track->StartDirection().Z();
            fPFPTrackEndDirectionX[p] = track->EndDirection().X();
            fPFPTrackEndDirectionY[p] = track->EndDirection().Y();
            fPFPTrackEndDirectionZ[p] = track->EndDirection().Z();
        }
        else if (dune_ana::DUNEAnaPFParticleUtils::IsShower(pfp, e, fPFParticleLabel, fShowerLabel))
        {
            art::Ptr<recob::Shower> shower = dune_ana::DUNEAnaPFParticleUtils::GetShower(pfp, e, fPFParticleLabel, fShowerLabel);
            fPFPIsShower[p] = true;
            fPFPShowerID[p] = shower->ID();
            //const std::vector<double> &energy{shower->Energy()};
            //fPFPShowerEnergy[p] = energy[2] / 1000.;
            fPFPShowerEnergy[p] = 0; // placeholder
            fPFPShowerEnergyToTrueEnergyRatio[p] = -1;
            fPFPShowerEnergyToTrueMomentumRatio[p] = -1;
            fPFPShowerDirectionX[p] = shower->Direction().X();
            fPFPShowerDirectionY[p] = shower->Direction().Y();
            fPFPShowerDirectionZ[p] = shower->Direction().Z();
            fPFPShowerStartX[p] = shower->ShowerStart().X();
            fPFPShowerStartY[p] = shower->ShowerStart().Y();
            fPFPShowerStartZ[p] = shower->ShowerStart().Z();
            fPFPShowerLength[p] = shower->Length();
            fPFPShowerOpenAngle[p] = shower->OpenAngle();
            if(!e.isRealData())
            {
                const int pos{fPFPTrueParticleMatchedPosition[p]};
                fPFPShowerEnergyToTrueEnergyRatio[p] = fPFPShowerEnergy[p] / fMCParticleTrueEnergy[pos];
                fPFPShowerEnergyToTrueMomentumRatio[p] = fPFPShowerEnergy[p] / fMCParticleTrueMom[pos];
            }
        }

        if(!e.isRealData())
        {
            for (const auto& hit: pfpHits)
            {
                TruthMatchUtils::G4ID g4ID(TruthMatchUtils::TrueParticleID(clockData, hit, fRollUpUnsavedIDs));
                if (TruthMatchUtils::Valid(g4ID) && g4ID == fPFPTrueParticleMatchedID[p])
                    ++fPFPNSharedTrueParticleHits[p];
            }

            fPFPPurity[p] = fPFPNHits[p] ? fPFPNSharedTrueParticleHits[p] / fPFPNHits[p] : 0.f;
            if (fPFPTrueParticleMatchedPosition[p]< 999999)
            {
                 fPFPCompleteness[p] = fMCParticleNHits[fPFPTrueParticleMatchedPosition[p]] ?
                     fPFPNSharedTrueParticleHits[p] / fMCParticleNHits[fPFPTrueParticleMatchedPosition[p]] : 0.f;
            }
            else
            {
                fPFPCompleteness[p] = 0.f;
            }

        }
        ++p;
    }
    fTree->Fill();
}

void test::pandoraAnalysis::beginJob()
{

    //deep clean the variables
    reset(true);
    // Implementation of optional member function here.
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("pandoraOutput","Pandora Output Tree");

    //Event branches
    fTree->Branch("eventID",&fEventID,"eventID/i");
    fTree->Branch("runID",&fRunID,"runID/i");
    fTree->Branch("subrunID",&fSubRunID,"subrunID/i");
    fTree->Branch("nMCParticles",&fNMCReconstructible,"nMCParticles/i");
    fTree->Branch("nPFParticles",&fNPFParticles,"nPFParticles/i");

    //MC truth branches
    //fTree->Branch("mcIsMCPrimary",&fMCIsPrimary);
    fTree->Branch("mcIsMCPrimary",&fMCIsPrimary,"MCIsPrimary[nMCParticles]/O");
    fTree->Branch("mcParticlePdgCode",&fMCParticlePdgCode,"MCParticlePdgCode[nMCParticles]/I");
    fTree->Branch("mcParticleTrueEnergy",&fMCParticleTrueEnergy,"MCParticleTrueEnergy[nMCParticles]/D");
    fTree->Branch("mcParticleTrackID",&fMCParticleTrackID,"MCParticleTrackID[nMCParticles]/I");
    fTree->Branch("mcParticleParentTrackID",&fMCParticleParentTrackID,"MCParticleParentTrackID[nMCParticles]/I");
    fTree->Branch("mcParticleStartPositionX",&fMCParticleStartPositionX,"MCParticleStartPositionX[nMCParticles]/D");
    fTree->Branch("mcParticleStartPositionY",&fMCParticleStartPositionY,"MCParticleStartPositionY[nMCParticles]/D");
    fTree->Branch("mcParticleStartPositionZ",&fMCParticleStartPositionZ,"MCParticleStartPositionZ[nMCParticles]/D");
    fTree->Branch("mcParticleStartPositionT",&fMCParticleStartPositionT,"MCParticleStartPositionT[nMCParticles]/D");
    fTree->Branch("mcParticleStartMomentumX",&fMCParticleStartMomentumX, "MCParticleStartMomentumX[nMCParticles]/D");
    fTree->Branch("mcParticleStartMomentumY",&fMCParticleStartMomentumY, "MCParticleStartMomentumY[nMCParticles]/D");
    fTree->Branch("mcParticleStartMomentumZ",&fMCParticleStartMomentumZ, "MCParticleStartMomentumZ[nMCParticles]/D");
    fTree->Branch("mcParticleStartMomentumE",&fMCParticleStartMomentumE, "MCParticleStartMomentumE[nMCParticles]/D");
    fTree->Branch("mcParticleEndPositionX",&fMCParticleEndPositionX, "MCParticleEndPositionX[nMCParticles]/D");
    fTree->Branch("mcParticleEndPositionY",&fMCParticleEndPositionY, "MCParticleEndPositionY[nMCParticles]/D");
    fTree->Branch("mcParticleEndPositionZ",&fMCParticleEndPositionZ, "MCParticleEndPositionZ[nMCParticles]/D");
    fTree->Branch("mcParticleEndPositionT",&fMCParticleEndPositionT, "MCParticleEndPositionT[nMCParticles]/D");
    fTree->Branch("mcParticleEndMomentumX",&fMCParticleEndMomentumX, "MCParticleEndMomentumX[nMCParticles]/D");
    fTree->Branch("mcParticleEndMomentumY",&fMCParticleEndMomentumY, "MCParticleEndMomentumY[nMCParticles]/D");
    fTree->Branch("mcParticleEndMomentumZ",&fMCParticleEndMomentumZ, "MCParticleEndMomentumZ[nMCParticles]/D");
    fTree->Branch("mcParticleEndMomentumE",&fMCParticleEndMomentumE, "MCParticleEndMomentumE[nMCParticles]/D");
    fTree->Branch("mcParticleNHits", &fMCParticleNHits, "MCParticleNHits[nMCParticles]/I");
    fTree->Branch("mcParticleNHitsView", &fMCParticleNHitsView, "MCParticleNHitsView[nMCParticles][3]/I");

    //PFP branches
    fTree->Branch("pfpTrueParticleMatchedID",&fPFPTrueParticleMatchedID,"PFPTrueParticleMatchedID[nPFParticles]/I");
    fTree->Branch("pfpTrueParticleMatchedPosition",&fPFPTrueParticleMatchedPosition,"PFPTrueParticleMatchedPosition[nPFParticles]/I");
    fTree->Branch("pfpIsPrimary",&fPFPIsPrimary,"PFPIsPrimary[nPFParticles]/O");
    fTree->Branch("pfpID",&fPFPID, "PFPID[nPFParticles]/I");
    fTree->Branch("pfpParentID",&fPFPParentID, "PFPParentID[nPFParticles]/I");
    fTree->Branch("pfpPdgCode",&fPFPPdgCode, "PFPPdgCode[nPFParticles]/I");
    fTree->Branch("pfpNChildren",&fPFPNChildren,"PFPNChildren[nPFParticles]/I");
    fTree->Branch("pfpNClusters",&fPFPNClusters,"PFPNClusters[nPFParticles]/I");
    fTree->Branch("pfpNHits",&fPFPNHits,"PFPNHits[nPFParticles]/I");
    fTree->Branch("pfpNHitsView",&fPFPNHitsView,"PFPNHitsView[nPFParticles][3]/I");
    fTree->Branch("pfpNSharedTrueParticleHits",&fPFPNSharedTrueParticleHits,"PFPNSharedTrueParticleHits[nPFParticles]/I");
    fTree->Branch("pfpNSharedTrueParticleHitsView",&fPFPNSharedTrueParticleHitsView,"PFPNSharedTrueParticleHitsView[nPFParticles][3]/I");
    fTree->Branch("pfpIsTrack",&fPFPIsTrack,"PFPIsTrack[nPFParticles]/O");
    fTree->Branch("pfpIsShower",&fPFPIsShower,"PFPIsShower[nPFParticles]/O");
    fTree->Branch("pfpTrackID", &fPFPTrackID,"PFPNClusters[nPFParticles]/I");
    fTree->Branch("pfpTrackLength",&fPFPTrackLength,"PFPTrackLength[nPFParticles]/D");
    fTree->Branch("pfpTrackStartX",&fPFPTrackStartX,"PFPTrackStartX[nPFParticles]/D");
    fTree->Branch("pfpTrackStartY",&fPFPTrackStartY,"PFPTrackStartY[nPFParticles]/D");
    fTree->Branch("pfpTrackStartZ",&fPFPTrackStartZ,"PFPTrackStartZ[nPFParticles]/D");
    fTree->Branch("pfpTrackVertexX",&fPFPTrackVertexX,"PFPTrackVertexX[nPFParticles]/D");
    fTree->Branch("pfpTrackVertexY",&fPFPTrackVertexY,"PFPTrackVertexY[nPFParticles]/D");
    fTree->Branch("pfpTrackVertexZ",&fPFPTrackVertexZ,"PFPTrackVertexZ[nPFParticles]/D");
    fTree->Branch("pfpTrackEndX",&fPFPTrackEndX,"PFPTrackEndX[nPFParticles]/D");
    fTree->Branch("pfpTrackEndY",&fPFPTrackEndY,"PFPTrackEndY[nPFParticles]/D");
    fTree->Branch("pfpTrackEndZ",&fPFPTrackEndZ,"PFPTrackEndZ[nPFParticles]/D");
    fTree->Branch("pfpTrackTheta",&fPFPTrackTheta,"PFPTrackTheta[nPFParticles]/D");
    fTree->Branch("pfpTrackPhi",&fPFPTrackPhi,"PFPTrackPhi[nPFParticles]/D");
    fTree->Branch("pfpTrackZenithAngle",&fPFPTrackZenithAngle,"PFPTrackZenithAngle[nPFParticles]/D");
    fTree->Branch("pfpTrackAzimuthAngle",&fPFPTrackAzimuthAngle,"PFPTrackAzimuthAngle[nPFParticles]/D");
    fTree->Branch("pfpTrackStartDirectionX",&fPFPTrackStartDirectionX,"PFPTrackStartDirectionX[nPFParticles]/D");
    fTree->Branch("pfpTrackStartDirectionY",&fPFPTrackStartDirectionY,"PFPTrackStartDirectionY[nPFParticles]/D");
    fTree->Branch("pfpTrackStartDirectionZ",&fPFPTrackStartDirectionZ,"PFPTrackStartDirectionZ[nPFParticles]/D");
    fTree->Branch("pfpTrackEndDirectionX",&fPFPTrackEndDirectionX,"PFPTrackEndDirectionX[nPFParticles]/D");
    fTree->Branch("pfpTrackEndDirectionY",&fPFPTrackEndDirectionY,"PFPTrackEndDirectionY[nPFParticles]/D");
    fTree->Branch("pfpTrackEndDirectionZ",&fPFPTrackEndDirectionZ,"PFPTrackEndDirectionZ[nPFParticles]/D");

    fTree->Branch("pfpShowerID",&fPFPShowerID,"PFPShowerID[nPFParticles]/I");
    fTree->Branch("pfpShowerEnergy",&fPFPShowerEnergy,"PFPShowerEnergy[nPFParticles]/D");
    fTree->Branch("pfpShowerEnergyToTrueEnergyRatio",&fPFPShowerEnergyToTrueEnergyRatio,"PFPShowerEnergyToTrueEnergyRatio[nPFParticles]/D");
    fTree->Branch("pfpShowerEnergyToTrueMomentumRatio",&fPFPShowerEnergyToTrueMomentumRatio,"PFPShowerEnergyToTrueMomentumRatio[nPFParticles]/D");
    fTree->Branch("pfpShowerDirectionX",&fPFPShowerDirectionX,"PFPShowerDirectionX[nPFParticles]/D");
    fTree->Branch("pfpShowerDirectionY",&fPFPShowerDirectionY,"PFPShowerDirectionY[nPFParticles]/D");
    fTree->Branch("pfpShowerDirectionZ",&fPFPShowerDirectionZ,"PFPShowerDirectionZ[nPFParticles]/D");
    fTree->Branch("pfpShowerStartX",&fPFPShowerStartX,"PFPShowerStartX[nPFParticles]/D");
    fTree->Branch("pfpShowerStartY",&fPFPShowerStartY,"PFPShowerStartY[nPFParticles]/D");
    fTree->Branch("pfpShowerStartZ",&fPFPShowerStartZ,"PFPShowerStartZ[nPFParticles]/D");
    fTree->Branch("pfpShowerLength",&fPFPShowerLength,"PFPShowerLength[nPFParticles]/D");
    fTree->Branch("pfpShowerOpeningAngle",&fPFPShowerOpenAngle,"PFPShowerOpenAngle[nPFParticles]/D");

    fTree->Branch("pfpCompleteness", &fPFPCompleteness, "PFPCompleteness[nPFParticles]/D");
    fTree->Branch("pfpCompletenessView", &fPFPCompletenessView, "PFPCompletenessView[nPFParticles][3]/D");
    fTree->Branch("pfpPurity", &fPFPPurity, "PFPPurity[nPFParticles]/D");
    fTree->Branch("pfpPurityView", &fPFPPurityView, "PFPPurityView[nPFParticles][3]/D");
}

void test::pandoraAnalysis::endJob()
{
    // Implementation of optional member function here.
}

void test::pandoraAnalysis::reset(bool deepClean)
{
    for(unsigned int iMc=0; iMc<(deepClean ? kNMaxMCParticles : fNMCParticles); iMc++){
        fMCIsPrimary[iMc]=0;
        fMCParticlePdgCode[iMc]=0;
        fMCParticleTrueEnergy[iMc]=-999999;
        fMCParticleTrackID[iMc]=999999;
        fMCParticleParentTrackID[iMc]=999999;
        fMCParticleStartPositionX[iMc]=999999;
        fMCParticleStartPositionY[iMc]=999999;
        fMCParticleStartPositionZ[iMc]=999999;
        fMCParticleStartPositionT[iMc]=999999;
        fMCParticleStartMomentumX[iMc]=999999;
        fMCParticleStartMomentumY[iMc]=999999;
        fMCParticleStartMomentumZ[iMc]=999999;
        fMCParticleStartMomentumE[iMc]=999999;
        fMCParticleEndPositionX[iMc]=999999;
        fMCParticleEndPositionY[iMc]=999999;
        fMCParticleEndPositionZ[iMc]=999999;
        fMCParticleEndPositionT[iMc]=999999;
        fMCParticleEndMomentumX[iMc]=999999;
        fMCParticleEndMomentumY[iMc]=999999;
        fMCParticleEndMomentumZ[iMc]=999999;
        fMCParticleEndMomentumE[iMc]=999999;
        fMCParticleNHits[iMc]=999999;
        fMCParticleNHitsView[iMc][0]=999999;
        fMCParticleNHitsView[iMc][1]=999999;
        fMCParticleNHitsView[iMc][2]=999999;
    }
    fNMCParticles = 0;
    fNMCReconstructible = 0;

    for(unsigned int iPfp=0; iPfp<(deepClean ? kNMaxPFParticles : fNPFParticles); iPfp++){

        fPFPID[iPfp]=999999;
        fPFPTrueParticleMatchedID[iPfp]=999999;
        fPFPTrueParticleMatchedPosition[iPfp]=999999;

        fPFPNHits[iPfp]=999999;
        fPFPNSharedTrueParticleHits[iPfp]=0;

        fPFPNClusters[iPfp]=999999;
        fPFPIsTrack[iPfp]=0;
        fPFPIsShower[iPfp]=0;

        fPFPTrackID[iPfp]=999999;
        fPFPTrackLength[iPfp]=999999;
        fPFPTrackStartX[iPfp]=999999;
        fPFPTrackStartY[iPfp]=999999;
        fPFPTrackStartZ[iPfp]=999999;
        fPFPTrackVertexX[iPfp]=999999;
        fPFPTrackVertexY[iPfp]=999999;
        fPFPTrackVertexZ[iPfp]=999999;
        fPFPTrackEndX[iPfp]=999999;
        fPFPTrackEndY[iPfp]=999999;
        fPFPTrackEndZ[iPfp]=999999;
        fPFPTrackTheta[iPfp]=999999;
        fPFPTrackPhi[iPfp]=999999;
        fPFPTrackZenithAngle[iPfp]=999999;
        fPFPTrackAzimuthAngle[iPfp]=999999;
        fPFPTrackStartDirectionX[iPfp]=999999;
        fPFPTrackStartDirectionY[iPfp]=999999;
        fPFPTrackStartDirectionZ[iPfp]=999999;
        fPFPTrackEndDirectionX[iPfp]=999999;
        fPFPTrackEndDirectionY[iPfp]=999999;
        fPFPTrackEndDirectionZ[iPfp]=999999;

        fPFPShowerID[iPfp]=999999;
        fPFPShowerEnergy[iPfp]=999999;
        fPFPShowerEnergyToTrueEnergyRatio[iPfp]=999999;
        fPFPShowerEnergyToTrueMomentumRatio[iPfp]=999999;
        fPFPShowerDirectionX[iPfp]=999999;
        fPFPShowerDirectionY[iPfp]=999999;
        fPFPShowerDirectionZ[iPfp]=999999;
        fPFPShowerStartX[iPfp]=999999;
        fPFPShowerStartY[iPfp]=999999;
        fPFPShowerStartZ[iPfp]=999999;
        fPFPShowerLength[iPfp]=999999;
        fPFPShowerOpenAngle[iPfp]=999999;

        fPFPCompleteness[iPfp]=999999;
        fPFPPurity[iPfp]=999999;

        for (unsigned int iView = 0; iView < kNViews; iView++)
        {
            fPFPNHitsView[iPfp][iView] = 999999;
            fPFPNSharedTrueParticleHitsView[iPfp][iView]=0;
            fPFPPurityView[iPfp][iView]=999999;
            fPFPCompletenessView[iPfp][iView]=999999;
        }
    }
    fNPFParticles = 0;

    return;
}

DEFINE_ART_MODULE(test::pandoraAnalysis)
