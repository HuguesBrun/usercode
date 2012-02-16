// -*- C++ -*-
//
// Package:    ReadTrigger
// Class:      ReadTrigger
// 
/**\class ReadTrigger ReadTrigger.cc hugues/ReadTrigger/src/ReadTrigger.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Hugues Louis Brun,40 4-B16,+41227671548,
//         Created:  Mon Feb 13 18:31:25 CET 2012
// $Id$
//
//


// system include files
#include <memory>
#include "TTree.h"
#include "TFile.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/Provenance/interface/ParameterSetID.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/ESHandle.h"

// recup the trigger for the RECO event
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
//
// class declaration
//

class ReadTrigger : public edm::EDAnalyzer {
   public:
      explicit ReadTrigger(const edm::ParameterSet&);
      ~ReadTrigger();

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
      edm::InputTag triggerResultsTag_;
      HLTConfigProvider hltConfig_;
      //edm::TriggerNames triggerNames_;
      std::string outputFile_;
      int eventNumber;
      int bit0, bit1, bit2, bit3, bit4, bit5, bit6, bit7, bit8, bit9, bit10, bit11, bit12, bit13;

	TTree *myTree_;
	TFile *myFile;

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
ReadTrigger::ReadTrigger(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   triggerResultsTag_ = iConfig.getParameter<edm::InputTag> ("hltProducer");
   outputFile_			       = iConfig.getParameter<std::string>("outputFile");
   myFile = TFile::Open(outputFile_.c_str(),"RECREATE");
}


ReadTrigger::~ReadTrigger()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ReadTrigger::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
//	std::cout << iEvent.id().event() << std::endl;
	edm::Handle<edm::TriggerResults> trigResults;

      iEvent.getByLabel(triggerResultsTag_,trigResults);
      const edm::TriggerNames & triggerNames = iEvent.triggerNames(*trigResults);
      std::vector<std::string>  hltNames_;
	hltNames_=triggerNames.triggerNames();
	 unsigned int n(hltNames_.size());
         bit0 = 0;
         bit1 = 0;
         bit2 = 0;
         bit3 =0;
         bit4 = 0;
         bit5 = 0;
         bit6 = 0;
         bit7 = 0;
         bit8 = 0;
         bit9 = 0;
         bit10 = 0;
         bit11 = 0;
         bit12 = 0;
         bit13 = 0;

//	for (unsigned int i =0 ; i < n; i++){
//		std::cout << hltNames_[i] << std::endl;
//	}
      //if (trigResults.isValid()) std::cout << "coucou " << std::endl;
  //     const unsigned int n(hltNames_.size());
   //    for (unsigned int i=0; i!=n; ++i) {
		
     //  }
       eventNumber = iEvent.id().event();
       bit0 = trigResults->accept(0);    
       bit1 = trigResults->accept(1);    
       bit2 = trigResults->accept(2);    
       bit3 = trigResults->accept(3);    
       bit4 = trigResults->accept(4);    
       bit5 = trigResults->accept(5);    
       bit6 = trigResults->accept(6);    
       bit7 = trigResults->accept(7);    
       bit8 = trigResults->accept(8);    
       bit9 = trigResults->accept(9);    
       bit10 = trigResults->accept(10);    
       bit11 = trigResults->accept(11);    
       bit12 = trigResults->accept(12);    
       bit13 = trigResults->accept(13);   

	myTree_->Fill(); 

/*	std::cout << " " << hltNames_[0] << "=" << bit0;
	std::cout << " " << hltNames_[1] << "=" << bit1;
	std::cout << " " << hltNames_[2] << "=" << bit2;
	std::cout << " " << hltNames_[3] << "=" << bit3;
	std::cout << " " << hltNames_[4] << "=" << bit4;
	std::cout << " " << hltNames_[5] << "=" << bit5;
	std::cout << " " << hltNames_[6] << "=" << bit6;
	std::cout << " " << hltNames_[7] << "=" << bit7;
	std::cout << " " << hltNames_[8] << "=" << bit8;
	std::cout << " " << hltNames_[9] << "=" << bit9;
	std::cout << " " << hltNames_[10] << "=" << bit10;
	std::cout << " " << hltNames_[11] << "=" << bit11;
	std::cout << " " << hltNames_[12] << "=" << bit12;
	std::cout << " " << hltNames_[13] << "=" << bit13;*/
}


// ------------ method called once each job just before starting event loop  ------------
void 
ReadTrigger::beginJob()
{
	myTree_ = new TTree("triggerResults","");
        myTree_->Branch("eventNumber",&eventNumber,"eventNumber/I");
	myTree_->Branch("bit0",&bit0,"bit0/I");
	myTree_->Branch("bit1",&bit1,"bit1/I");
	myTree_->Branch("bit2",&bit2,"bit2/I");
	myTree_->Branch("bit3",&bit3,"bit3/I");
	myTree_->Branch("bit4",&bit4,"bit4/I");
	myTree_->Branch("bit5",&bit5,"bit5/I");
	myTree_->Branch("bit6",&bit6,"bit6/I");
	myTree_->Branch("bit7",&bit7,"bit7/I");
	myTree_->Branch("bit8",&bit8,"bit8/I");
	myTree_->Branch("bit9",&bit9,"bit9/I");
	myTree_->Branch("bit10",&bit10,"bit10/I");
	myTree_->Branch("bit11",&bit11,"bit11/I");
	myTree_->Branch("bit12",&bit12,"bit12/I");
	myTree_->Branch("bit13",&bit13,"bit13/I");
}


// ------------ method called once each job just after ending the event loop  ------------
void 
ReadTrigger::endJob() 
{
myFile->cd();
myTree_->Write();
myFile->Close();
}

// ------------ method called when starting to processes a run  ------------
void 
ReadTrigger::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
    //hltConfig_.init(iRun,iSetup,triggerResultsTag_.process(),true);
}

// ------------ method called when ending the processing of a run  ------------
void 
ReadTrigger::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ReadTrigger::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ReadTrigger::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ReadTrigger::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ReadTrigger);
