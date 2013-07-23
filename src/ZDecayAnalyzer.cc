// -*- C++ -*-
//
// Package:    ZDecayAnalyzer
// Class:      ZDecayAnalyzer
//
/**\class ZDecayAnalyzer ZDecayAnalyzer.cc MyAnalysis/ZDecayAnalyzer/src/ZDecayAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Dinko Ferencek,8 R-004,+41227676479,
//         Created:  Wed Jul 17 21:47:50 CEST 2013
// $Id: ZDecayAnalyzer.cc,v 1.1 2013/07/18 00:24:27 ferencek Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TVector3.h"

//
// class declaration
//

class ZDecayAnalyzer : public edm::EDAnalyzer {
   public:
      explicit ZDecayAnalyzer(const edm::ParameterSet&);
      ~ZDecayAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
      const edm::InputTag    genParticleTag;
      const int              bosonPdgId;
      const std::vector<int> bosonDecayProdPdgIds;

      edm::Service<TFileService> fs;

      TH1D *h1_BosonPt;
      TH1D *h1_BosonEta;
      TH1D *h1_BosonMass;
      TH1D *h1_BosonMassFromDecayProd;

      TH1D *h1_BosonRestFramePt;
      TH1D *h1_BosonRestFrameMass;
      TH1D *h1_BosonRestFrameMassFromDecayProd;

      TH1D *h1_CosThetaStar;
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
ZDecayAnalyzer::ZDecayAnalyzer(const edm::ParameterSet& iConfig) :

  genParticleTag(iConfig.getParameter<edm::InputTag>("GenParticleTag")),
  bosonPdgId(iConfig.getParameter<int>("BosonPdgId")),
  bosonDecayProdPdgIds(iConfig.getParameter<std::vector<int> >("BosonDecayProdPdgIds"))

{
   //now do what ever initialization is needed
   int ptBins=1000;
   double ptMin=0., ptMax=1000.;
   int etaBins=200;
   double etaMin=-5., etaMax=5.;
   int massBins=150;
   int massMin=0., massMax=150.;
   int thetaBins=50;
   int thetaMin=0., thetaMax=1.;

   h1_BosonPt = fs->make<TH1D>("h1_BosonPt",";p_{T} [GeV];",ptBins,ptMin,ptMax);
   h1_BosonPt->Sumw2();
   h1_BosonPt->SetDefaultSumw2(kTRUE);
   h1_BosonEta               = fs->make<TH1D>("h1_BosonEta",";#eta;",etaBins,etaMin,etaMax);
   h1_BosonMass              = fs->make<TH1D>("h1_BosonMass",";m [GeV];",massBins,massMin,massMax);
   h1_BosonMassFromDecayProd = fs->make<TH1D>("h1_BosonMassFromDecayProd",";m [GeV];",massBins,massMin,massMax);

   h1_BosonRestFramePt                = fs->make<TH1D>("h1_BosonRestFramePt",";p_{T} [GeV];",ptBins,ptMin,ptMax);
   h1_BosonRestFrameMass              = fs->make<TH1D>("h1_BosonRestFrameMass",";m [GeV];",massBins,massMin,massMax);
   h1_BosonRestFrameMassFromDecayProd = fs->make<TH1D>("h1_BosonRestFrameMassFromDecayProd",";m [GeV];",massBins,massMin,massMax);

   h1_CosThetaStar = fs->make<TH1D>("h1_CosThetaStar",";|cos#theta*|;",thetaBins,thetaMin,thetaMax);
}


ZDecayAnalyzer::~ZDecayAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ZDecayAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   edm::Handle<reco::GenParticleCollection> genParticles;
   iEvent.getByLabel(genParticleTag,genParticles);

   // loop over GenParticles and select bosons
   for(reco::GenParticleCollection::const_iterator it = genParticles->begin(); it != genParticles->end(); ++it)
   {
     if( it->status() == 2 ) break; // to speed things up (only works with Pythia6)
     if( abs(it->pdgId()) == abs(bosonPdgId) && it->status() == 3 )
     {
       bool decayProductsFound = false;
       // vector of pointers to bosob decay products
       std::vector<const reco::Candidate*> bosonDecayProducts;

       for(unsigned i=0; i<it->numberOfDaughters(); ++i)
       {
         if( it->daughter(i)->status()==2 ) continue; // only care about status=3 daughters
         for(std::vector<int>::const_iterator pdgIdIt = bosonDecayProdPdgIds.begin(); pdgIdIt != bosonDecayProdPdgIds.end(); ++pdgIdIt)
         {
           if( abs(it->daughter(i)->pdgId()) == abs(*pdgIdIt) )
           {
             decayProductsFound = true;
             bosonDecayProducts.push_back(it->daughter(i));
           }
         }
       }
       if( decayProductsFound && bosonDecayProducts.size()<2 )
       {
         edm::LogError("TooFewDecayProducts") << "Fewer than two boson decay products found.";
         return;
       }
       if( decayProductsFound && bosonDecayProducts.size()>2 )
       {
         edm::LogError("TooManyDecayProducts") << "More than two boson decay products found.";
         return;
       }
       if( !decayProductsFound ) return;

       h1_BosonPt->Fill( it->pt() );
       h1_BosonEta->Fill( it->eta() );
       h1_BosonMass->Fill( it->mass() );

       TLorentzVector v_boson, v_prod1, v_prod2, v_prod_sum;
       v_boson.SetPtEtaPhiM( it->pt(), it->eta(), it->phi(), it->mass() );
       v_prod1.SetPtEtaPhiM( bosonDecayProducts[0]->pt(), bosonDecayProducts[0]->eta(), bosonDecayProducts[0]->phi(), bosonDecayProducts[0]->mass() );
       v_prod2.SetPtEtaPhiM( bosonDecayProducts[1]->pt(), bosonDecayProducts[1]->eta(), bosonDecayProducts[1]->phi(), bosonDecayProducts[1]->mass() );

       v_prod_sum = v_prod1 + v_prod2;

       h1_BosonMassFromDecayProd->Fill( v_prod_sum.M() );

       TVector3 p_boson = v_boson.Vect();

       TVector3 boost = v_boson.BoostVector();

       v_boson.Boost(-boost.x(),-boost.y(),-boost.z());
       v_prod1.Boost(-boost.x(),-boost.y(),-boost.z());
       v_prod2.Boost(-boost.x(),-boost.y(),-boost.z());

       v_prod_sum = v_prod1 + v_prod2;

       h1_BosonRestFramePt->Fill( v_boson.Pt() );
       h1_BosonRestFrameMass->Fill( v_boson.M() );
       h1_BosonRestFrameMassFromDecayProd->Fill( v_prod_sum.M() );


       TVector3 p_prod1 = v_prod1.Vect();

       double cosThetaStar = (p_prod1.Dot(p_boson))/(p_prod1.Mag()*p_boson.Mag());

       h1_CosThetaStar->Fill(fabs(cosThetaStar));
     }
   }
}


// ------------ method called once each job just before starting event loop  ------------
void
ZDecayAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
ZDecayAnalyzer::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void
ZDecayAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
ZDecayAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
ZDecayAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
ZDecayAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZDecayAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZDecayAnalyzer);
