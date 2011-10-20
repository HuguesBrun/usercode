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
// $Id: ReadDataBase.cc,v 1.1 2011/04/21 06:48:41 hbrun Exp $
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

#include "CondFormats/EcalObjects/interface/EcalClusterCrackCorrParameters.h"
#include "CondFormats/DataRecord/interface/EcalClusterCrackCorrParametersRcd.h"

//#include "CondFormats/EcalObjects/interface/EcalClusterLocalContCorrParameters.h"
//#include "CondFormats/DataRecord/interface/EcalClusterLocalContCorrParametersRcd.h"

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

	edm::ESHandle<EcalClusterCrackCorrParameters> theParamsLocal;
	iSetup.get<EcalClusterCrackCorrParametersRcd>().get(theParamsLocal);
	const EcalClusterCrackCorrParameters *myParamsLocal = theParamsLocal.product();


   ESHandle<EcalClusterEnergyCorrectionParameters> theParams;
   iSetup.get<EcalClusterEnergyCorrectionParametersRcd>().get(theParams);
   const EcalClusterEnergyCorrectionParameters *params = theParams.product();

		
		 int offset = 0;
/*		cout << "offset= " << offset << endl;
		for(int j = 2 ; j <= 20 ; j++ )
		{
		 cout << j << "+ offset= " << setprecision(10) << fixed << (params->params())[j + offset] << endl;
		}

		cout << endl << endl;
		 offset = 20;
		cout << "offset= " << offset << endl;
		for(int j = 2 ; j <= 20 ; j++ )
		{
		 cout << j << "+ offset= " << (params->params())[j + offset] << endl;
		}
*/
	cout << endl;
	cout << "#################################################" << endl;
	cout << "C(eta)" << endl;
	cout << "#################################################" << endl;
	for( int j = 0 ; j <= 1 ; j++ )
	{
		cout << setprecision(10) << fixed << (params->params())[j + offset] << ", ";
	}

	cout << endl;
	cout << "#################################################" << endl;
  cout << "f(brem), EB" << endl;
  cout << "#################################################" << endl;
	for( int j = 2 ; j <= 8 ; j++ )
  {
    cout << setprecision(10) << fixed << (params->params())[j + offset] << ", ";
  }

	cout << endl;
	cout << "#################################################" << endl;
  cout << "f(et, eta), EB" << endl;
  cout << "#################################################" << endl;
  for( int j = 9 ; j <= 19 ; j++ )
  {
    cout << setprecision(10) << fixed << (params->params())[j + offset] << ", ";
  }

	offset = 20;
	cout << endl;
	cout << "#################################################" << endl;
  cout << "f(brem), EE" << endl;
  cout << "#################################################" << endl;
	for( int j = 2 ; j <= 8 ; j++ )
  {
    cout << setprecision(10) << fixed << (params->params())[j + offset] << ", ";
  }

	cout << endl;
	cout << "#################################################" << endl;
  cout << "f(et, eta), EE" << endl;
  cout << "#################################################" << endl;
  for( int j = 9 ; j <= 19 ; j++ )
  {
    cout << setprecision(10) << fixed << (params->params())[j + offset] << ", ";
  }
	cout << endl << endl;


	cout << endl;
	cout << endl;
	cout << "#################################################" << endl;
	cout << "#################################################" << endl;
	cout << "#################################################" << endl;
	cout << endl;
	cout << "#################################################" << endl;
  cout << "coefficients for eta side of an intermodule gap closer to the interaction point" << endl;
	cout << "#################################################" << endl;
	offset = 0;
	for (int k=0; k!=5; ++k) cout << "(myParams->params())[" << k+offset << "]= " << (myParamsLocal->params())[k+offset] << endl;;
	cout << endl;
	cout << "#################################################" << endl;
  cout << "coefficients for the other eta side" << endl;
	cout << "#################################################" << endl;
	offset = 5;
	for (int k=0; k!=5; ++k) cout << "(myParams->params())[" << k+offset << "]= " << (myParamsLocal->params())[k+offset] << endl;;
	cout << endl;
	cout << "#################################################" << endl;
	cout << "coefficients for one phi side of a SM" << endl;
	cout << "#################################################" << endl;
	offset = 10;
	for (int k=0; k!=5; ++k) cout << "(myParams->params())[" << k+offset << "]= " << (myParamsLocal->params())[k+offset] << endl;;
	cout << endl;
	cout << "#################################################" << endl;
	cout << "coefficients for the other side" << endl;
	cout << "#################################################" << endl;
	offset = 15;
	for (int k=0; k!=5; ++k) cout << "(myParams->params())[" << k+offset << "]= " << (myParamsLocal->params())[k+offset] << endl;;



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
