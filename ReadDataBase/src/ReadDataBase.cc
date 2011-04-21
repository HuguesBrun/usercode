// -*- C++ -*-
//
// Package:    ReadDataBase
// Class:      ReadDataBase
// 
/**\class ReadDataBase ReadDataBase.cc hugues/ReadDataBase/src/ReadDataBase.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Hugues Brun
//         Created:  Wed Apr 20 15:43:04 CEST 2011
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CondFormats/EcalObjects/interface/EcalClusterEnergyCorrectionParameters.h"
#include "CondFormats/DataRecord/interface/EcalClusterEnergyCorrectionParametersRcd.h"



//
// class declaration
//

class ReadDataBase : public edm::EDAnalyzer {
   public:
      explicit ReadDataBase(const edm::ParameterSet&);
      ~ReadDataBase();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
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
ReadDataBase::ReadDataBase(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


ReadDataBase::~ReadDataBase()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
ReadDataBase::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

                 ESHandle<EcalClusterEnergyCorrectionParameters> theParams;
		 iSetup.get<EcalClusterEnergyCorrectionParametersRcd>().get(theParams);
                 const EcalClusterEnergyCorrectionParameters *params = theParams.product();

		
		 int offset = 20;
		 cout << "2+ offset" << (params->params())[2 + offset] << endl;
		 cout << "3+ offset" << (params->params())[3 + offset] << endl;
		 cout << "4+ offset" << (params->params())[4 + offset] << endl;
		 cout << "5+ offset" << (params->params())[5 + offset] << endl;
		 cout << "6+ offset" << (params->params())[6 + offset] << endl;
		 cout << "7+ offset" << (params->params())[7 + offset] << endl;
		 cout << "8+ offset" << (params->params())[8 + offset] << endl;
}


// ------------ method called once each job just before starting event loop  ------------
void 
ReadDataBase::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ReadDataBase::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(ReadDataBase);
