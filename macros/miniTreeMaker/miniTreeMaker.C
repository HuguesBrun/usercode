#include "miniTreeMaker.h"

void beginMacro(){

	doHLT                    = true;
  	doMC                     = false;
  	doJetMC                  = false;
  	doMETMC                  = false;
  	doPDFInfo                = true;
  	doSignalMuMuGamma        = false;
  	doSignalTopTop           = false;
  	doPhotonConversionMC     = false;
  	doBeamSpot               = true;
  	doPrimaryVertex          = true;
  	doZeePrimaryVertex       = false;
  	doTrack                  = true;
  	doJet                    = true;
  	doMuon                   = true;
  	doElectron               = true;
  	doPhoton                 = true;
  	doCluster                = true;
  	doPhotonConversion       = true;
  	doMET                    = true;
  	doBardak                 = false;
  	doPhotonVertexCorrection = false;
  	doPhotonIsolation        = true;


	mcParticles = new TClonesArray("TRootMCParticle", 0);
	genJets = new TClonesArray("TRootParticle", 0);
	genMETs = new TClonesArray("TRootParticle", 0);
	mcPhotons = new TClonesArray("TRootMCPhoton", 0);
	vertices = new TClonesArray("TRootVertex", 0);
	zeeVertices = new TClonesArray("TRootVertex", 0);
	tracks = new TClonesArray("TRootTrack", 0);
	jets = new TClonesArray("TRootJet", 0);
	pflowjets = new TClonesArray("TRootJet", 0);
	sisconejets = new TClonesArray("TRootJet", 0);
	muons = new TClonesArray("TRootMuon", 0);
	electrons = new TClonesArray("TRootElectron", 0);
	photons = new TClonesArray("TRootPhoton", 0);
	clusters = new TClonesArray("TRootCluster", 0);
	superClusters = new TClonesArray("TRootSuperCluster", 0);
	conversionTracks = new TClonesArray("TRootTrack", 0);
	met = new TClonesArray("TRootMET", 0);


  inputEventTree->SetBranchAddress("Event", &event, &event_br);
  inputEventTree->SetBranchStatus("Event", 1);
 
  inputRunTree->SetBranchAddress("runInfos", &runInfos, &run_br);
  inputRunTree->SetBranchStatus("runInfos", 1);
  if(doMC)
    {
      inputEventTree->SetBranchAddress("MCParticles", &mcParticles, &mcParticles_br);
      inputEventTree->SetBranchStatus("MCParticles", 1);
    }
  
  if(doJetMC)
    {
      inputEventTree->SetBranchAddress("genJets", &genJets, &genJets_br);
      inputEventTree->SetBranchStatus("genJets", 1);
    }
  
  if(doMETMC)
    {
      inputEventTree->SetBranchAddress("genMETs", &genMETs, &genMETs_br);
      inputEventTree->SetBranchStatus("genMETs", 1);
    }
	
  if(doSignalMuMuGamma)
    {
      inputEventTree->SetBranchAddress("MuMuGamma", &mcMuMuGammaEvent, &mcSignalMuMuGamma_br);
      inputEventTree->SetBranchStatus("MuMuGamma", 1);
    }
  
  if(doSignalTopTop)
    {
      inputEventTree->SetBranchAddress("rootMCTopTop", &mcTopTopEvent, &mcTopTopEvent_br);
      inputEventTree->SetBranchStatus("rootMCTopTop", 1);
    }
	
  if(doPhotonConversionMC)
    {
      inputEventTree->SetBranchAddress("MCPhotons", &mcPhotons, &mcPhotons_br);
      inputEventTree->SetBranchStatus("MCPhotons", 1);
    }

  if(doBeamSpot)
    {
      inputEventTree->SetBranchAddress("BeamSpot", &beamSpot, &beamSpot_br);
      inputEventTree->SetBranchStatus("BeamSpot", 1);
    }

  if(doPrimaryVertex)
    {
      inputEventTree->SetBranchAddress("Vertices", &vertices, &vertices_br);
      inputEventTree->SetBranchStatus("Vertices", 1);
    }

  if(doZeePrimaryVertex)
    {
      inputEventTree->SetBranchAddress("ZeeVertices", &zeeVertices, &zeeVertices_br);
      inputEventTree->SetBranchStatus("ZeeVertices", 1);
    }
  
  if(doTrack)
    {
      inputEventTree->SetBranchAddress("Tracks", &tracks, &tracks_br);
      inputEventTree->SetBranchStatus("Tracks", 1);
    }


  if(doJet)	{
    inputEventTree->SetBranchAddress("ak5CaloJets", &jets, &jets_br);
    inputEventTree->SetBranchStatus("ak5CaloJets", 1);
  }


  if(doMuon)
    {
      inputEventTree->SetBranchAddress("Muons", &muons, &muons_br);
      inputEventTree->SetBranchStatus("Muons", 1);
    }
  
  if(doElectron)
    {
      inputEventTree->SetBranchAddress("Electrons", &electrons, &electrons_br);
      inputEventTree->SetBranchStatus("Electrons", 1);
    }
  
  if(doPhoton)
    {
      inputEventTree->SetBranchAddress("Photons", &photons, &photons_br);
      inputEventTree->SetBranchStatus("Photons", 1);
    }
  
  if(doCluster)
    {
      inputEventTree->SetBranchAddress("BasicClusters", &clusters, &clusters_br);
      inputEventTree->SetBranchStatus("BasicClusters", 1);
      
      inputEventTree->SetBranchAddress("SuperClusters", &superClusters, &superClusters_br);
      inputEventTree->SetBranchStatus("SuperClusters", 1);
    }
  

  if(doPhotonConversion)
    {
      inputEventTree->SetBranchAddress("ConversionTracks", &conversionTracks, &conversions_br);
      inputEventTree->SetBranchStatus("ConversionTracks", 1);
    }

  if(doMET)
    {
      inputEventTree->SetBranchAddress("MET", &met, &met_br);
      inputEventTree->SetBranchStatus("MET", 1);
    }
  
  if(doBardak)
    {
      inputEventTree->SetBranchAddress("bardak", &bardak, &bardak_br);
      inputEventTree->SetBranchStatus("bardak", 1);
    }
  





	myTree_ = new TTree("photons","");
               myTree_->Branch("pho_iEvent",&pho_iEvent,"pho_iEvent/I");
                myTree_->Branch("pho_EventNumber",&pho_EventNumber,"pho_EventNumber/I");
                myTree_->Branch("pho_RunNumber",&pho_RunNumber,"pho_RunNumber/I");
		myTree_->Branch("pho_lumiSection",&pho_lumiSection,"pho_lumiSection/I");
                myTree_->Branch("pho_isEB",&pho_isEB,"pho_isEB/I");
                myTree_->Branch("pho_isEE",&pho_isEE,"pho_isEE/I");
                myTree_->Branch("pho_isEEP",&pho_isEEP,"pho_isEEP/I");
                myTree_->Branch("pho_isEEM",&pho_isEEM,"pho_isEEM/I");
                myTree_->Branch("pho_isAlsoElectron",&pho_isAlsoElectron,"pho_isAlsoElectron/I");
                myTree_->Branch("pho_ElectronClassification",&pho_ElectronClassification,"pho_ElectronClassification/I");
                myTree_->Branch("pho_hasPixelSeed",&pho_hasPixelSeed,"pho_hasPixelSeed/I");
                myTree_->Branch("pho_et",&pho_et,"pho_et/F");
                myTree_->Branch("pho_energy",&pho_energy,"pho_energy/F");
                myTree_->Branch("pho_px",&pho_px,"pho_px/F");
                myTree_->Branch("pho_py",&pho_py,"pho_py/F");
                myTree_->Branch("pho_pz",&pho_pz,"pho_pz/F");
                myTree_->Branch("pho_eta",&pho_eta,"pho_eta/F");
                myTree_->Branch("pho_phi",&pho_phi,"pho_phi/F");
                myTree_->Branch("pho_s4",&pho_s4,"pho_s4/F");
                myTree_->Branch("pho_s9",&pho_s9,"pho_s9/F");
                myTree_->Branch("pho_eMax",&pho_eMax,"pho_eMax/F");
                myTree_->Branch("pho_e2nd",&pho_e2nd,"pho_e2nd/F");
                myTree_->Branch("pho_r19",&pho_r19,"pho_r19/F");
                myTree_->Branch("pho_r9",&pho_r9,"pho_r9/F");
                myTree_->Branch("pho_S9overS9minusS1S2",&pho_S9overS9minusS1S2,"pho_S9overS9minusS1S2/F");
                myTree_->Branch("pho_cEE",&pho_cEE,"pho_cEE/F");
                myTree_->Branch("pho_cEP",&pho_cEP,"pho_cEP/F");
                myTree_->Branch("pho_cPP",&pho_cPP,"pho_cPP/F");
                myTree_->Branch("pho_phiwidth",&pho_phiwidth,"pho_phiwidth/F");
                myTree_->Branch("pho_etawidth",&pho_etawidth,"pho_etawidth/F");
                myTree_->Branch("pho_sigmaphi",&pho_sigmaphi,"pho_sigmaphi/F");
                myTree_->Branch("pho_sigmaeta",&pho_sigmaeta,"pho_sigmaeta/F");
                myTree_->Branch("pho_sigmaIetaIeta",&pho_sigmaIetaIeta,"pho_sigmaIetaIeta/F");
                myTree_->Branch("pho_sigmaEtaEta",&pho_sigmaEtaEta,"pho_sigmaEtaEta/F");
                myTree_->Branch("pho_hoe",&pho_hoe,"pho_hoe/F");
                myTree_->Branch("pho_IsoEcalRechit",&pho_IsoEcalRechit,"pho_IsoEcalRechit/F");
                myTree_->Branch("pho_IsoHcalRechit",&pho_IsoHcalRechit,"pho_IsoHcalRechit/F");
                myTree_->Branch("pho_IsoSolidTrkCone",&pho_IsoSolidTrkCone,"pho_IsoSolidTrkCone/F");
                myTree_->Branch("pho_IsoHollowTrkCone",&pho_IsoHollowTrkCone,"pho_IsoHollowTrkCone/F");
                myTree_->Branch("pho_IsoEcalRechit03",&pho_IsoEcalRechit03,"pho_IsoEcalRechit03/F");
                myTree_->Branch("pho_IsoHcalRechit03",&pho_IsoHcalRechit03,"pho_IsoHcalRechit03/F");
                myTree_->Branch("pho_IsoSolidTrkCone03",&pho_IsoSolidTrkCone03,"pho_IsoSolidTrkCone03/F");
                myTree_->Branch("pho_IsoHollowTrkCone03",&pho_IsoHollowTrkCone03,"pho_IsoHollowTrkCone03/F");
                myTree_->Branch("pho_esRatio",&pho_esRatio,"pho_esRatio/F");
                myTree_->Branch("pho_convNTracks",&pho_convNTracks,"pho_convNTracks/I");
                myTree_->Branch("pho_ptoverjetpt",&pho_ptoverjetpt,"pho_ptoverjetpt/F");
                myTree_->Branch("pho_DrSCclosest",&pho_DrSCclosest,"pho_DrSCclosest/F");
                myTree_->Branch("pho_DrTrkclosest",&pho_DrTrkclosest,"pho_DrTrkclosest/F");
                myTree_->Branch("pho_ptoverjetpt_pt2",&pho_ptoverjetpt_pt2,"pho_ptoverjetpt_pt2/F");
                myTree_->Branch("pho_DrTrkclosest_pt2",&pho_DrTrkclosest_pt2,"pho_DrTrkclosest_pt2/F");
                myTree_->Branch("pho_ptoverjetpt_pt5",&pho_ptoverjetpt_pt5,"pho_ptoverjetpt_pt5/F");
                myTree_->Branch("pho_DrTrkclosest_pt5",&pho_DrTrkclosest_pt5,"pho_DrTrkclosest_pt5/F");
                myTree_->Branch("pho_ptoverjetpt_pt10",&pho_ptoverjetpt_pt10,"pho_ptoverjetpt_pt10/F");
                myTree_->Branch("pho_DrTrkclosest_pt10",&pho_DrTrkclosest_pt10,"pho_DrTrkclosest_pt10/F");
                myTree_->Branch("pho_ptoverjetpt_pt15",&pho_ptoverjetpt_pt15,"pho_ptoverjetpt_pt15/F");
                myTree_->Branch("pho_ptoverjetpt_PFlow",&pho_ptoverjetpt_PFlow,"pho_ptoverjetpt_PFlow/F");
                myTree_->Branch("pho_ptoverjetpt_PFlow_pt2",&pho_ptoverjetpt_PFlow_pt2,"pho_ptoverjetpt_PFlow_pt2/F");
                myTree_->Branch("pho_ptoverjetpt_PFlow_pt5",&pho_ptoverjetpt_PFlow_pt5,"pho_ptoverjetpt_PFlow_pt5/F");
                myTree_->Branch("pho_ptoverjetpt_PFlow_pt10",&pho_ptoverjetpt_PFlow_pt10,"pho_ptoverjetpt_PFlow_pt10/F");
                myTree_->Branch("pho_ptoverjetpt_sisCone",&pho_ptoverjetpt_sisCone,"pho_ptoverjetpt_sisCone/F");
                myTree_->Branch("pho_ptoverjetpt_sisCone_pt2",&pho_ptoverjetpt_sisCone_pt2,"pho_ptoverjetpt_sisCone_pt2/F");
                myTree_->Branch("pho_ptoverjetpt_sisCone_pt5",&pho_ptoverjetpt_sisCone_pt5,"pho_ptoverjetpt_sisCone_pt5/F");
                myTree_->Branch("pho_ptoverjetpt_sisCone_pt10",&pho_ptoverjetpt_sisCone_pt10,"pho_ptoverjetpt_sisCone_pt10/F");
		myTree_->Branch("pho_transverseMomentumToJetDirection",&pho_transverseMomentumToJetDirection,"pho_transverseMomentumToJetDirection/F");
		myTree_->Branch("pho_transverseMomentumToJetDirection_pt2",&pho_transverseMomentumToJetDirection_pt2,"pho_transverseMomentumToJetDirection_pt2/F");
		myTree_->Branch("pho_transverseMomentumToJetDirection_pt5",&pho_transverseMomentumToJetDirection_pt5,"pho_transverseMomentumToJetDirection_pt5/F");
		myTree_->Branch("pho_transverseMomentumToJetDirection_pt10",&pho_transverseMomentumToJetDirection_pt10,"pho_transverseMomentumToJetDirection_pt10/F");
		myTree_->Branch("pho_transverseMomentumToJetDirection_pt20",&pho_transverseMomentumToJetDirection_pt20,"pho_transverseMomentumToJetDirection_pt20/F");
		myTree_->Branch("pho_transverseToJetRatio",&pho_transverseToJetRatio,"pho_transverseToJetRatio/F");
		myTree_->Branch("pho_transverseToJetRatio_pt2",&pho_transverseToJetRatio_pt2,"pho_transverseToJetRatio_pt2/F");
		myTree_->Branch("pho_transverseToJetRatio_pt5",&pho_transverseToJetRatio_pt5,"pho_transverseToJetRatio_pt5/F");
		myTree_->Branch("pho_transverseToJetRatio_pt10",&pho_transverseToJetRatio_pt10,"pho_transverseToJetRatio_pt10/F");
		myTree_->Branch("pho_transverseToJetRatio_pt20",&pho_transverseToJetRatio_pt20,"pho_transverseToJetRatio_pt20/F");
                myTree_->Branch("pho_isAfterCut1",&pho_isAfterCut1,"pho_isAfterCut1/I");
                myTree_->Branch("pho_isAfterCut2",&pho_isAfterCut2,"pho_isAfterCut2/I");
                myTree_->Branch("pho_isAfterCut3",&pho_isAfterCut3,"pho_isAfterCut3/I");
                myTree_->Branch("pho_isAfterCut4",&pho_isAfterCut4,"pho_isAfterCut4/I");
                myTree_->Branch("pho_isAfterCut5",&pho_isAfterCut5,"pho_isAfterCut5/I");
                myTree_->Branch("pho_isAfterCut6",&pho_isAfterCut6,"pho_isAfterCut6/I");
                myTree_->Branch("pho_isAfterCut7",&pho_isAfterCut7,"pho_isAfterCut7/I");
                myTree_->Branch("pho_isAfterCut8",&pho_isAfterCut8,"pho_isAfterCut8/I");
                myTree_->Branch("pho_jetEMfraction",&pho_jetEMfraction,"pho_jetEMfraction/F");
                myTree_->Branch("pho_DrJetClosest",&pho_DrJetClosest,"pho_DrJetClosest/F");
                myTree_->Branch("pho_isLoose",&pho_isLoose,"pho_isLoose/I");
                myTree_->Branch("pho_isTight",&pho_isTight,"pho_isTight/I");
                myTree_->Branch("pho_isEBGap",&pho_isEBGap,"pho_isEBGap/I");
                myTree_->Branch("pho_isEEGap",&pho_isEEGap,"pho_isEEGap/I");
                myTree_->Branch("pho_GenId",&pho_GenId,"pho_GenId/I");
                myTree_->Branch("pho_MotherId",&pho_MotherId,"pho_MotherId/I");
                myTree_->Branch("pho_isPromptGenPho",&pho_isPromptGenPho,"pho_isPromptGenPho/I");
                myTree_->Branch("pho_isFromQuarkGen",&pho_isFromQuarkGen,"pho_isFromQuarkGen/I");
                myTree_->Branch("pho_isPi0Gen",&pho_isPi0Gen,"pho_isPi0Gen/I");
                myTree_->Branch("pho_isEtaGen",&pho_isEtaGen,"pho_isEtaGen/I");
                myTree_->Branch("pho_isRhoGen",&pho_isRhoGen,"pho_isRhoGen/I");
                myTree_->Branch("pho_isOmegaGen",&pho_isOmegaGen,"pho_isOmegaGen/I");
                myTree_->Branch("pho_isGenElectron",&pho_isGenElectron,"pho_isGenElectron/I");
                myTree_->Branch("pho_eventPassHLT_Photon10_L1R",&pho_eventPassHLT_Photon10_L1R,"pho_eventPassHLT_Photon10_L1R/I");
                myTree_->Branch("pho_eventPassHLT_Photon15_L1R",&pho_eventPassHLT_Photon15_L1R,"pho_eventPassHLT_Photon15_L1R/I");
                myTree_->Branch("pho_eventPassHLT_DoublePhoton10_L1R",&pho_eventPassHLT_DoublePhoton10_L1R,"pho_eventPassHLT_DoublePhoton10_L1R/I");
                myTree_->Branch("pho_eventPtHat",&pho_eventPtHat,"pho_eventPtHat/F");
                myTree_->Branch("pho_nVertex",&pho_nVertex,"pho_nVertex/I");
                myTree_->Branch("pho_PromptGenIsoEnergyStatus1_cone02",&pho_PromptGenIsoEnergyStatus1_cone02,"pho_PromptGenIsoEnergyStatus1_cone02/F");
                myTree_->Branch("pho_PromptGenIsoEnergyStatus2_cone02",&pho_PromptGenIsoEnergyStatus2_cone02,"pho_PromptGenIsoEnergyStatus2_cone02/F");
                myTree_->Branch("pho_PromptGenIsoEnergyStatus1_cone03",&pho_PromptGenIsoEnergyStatus1_cone03,"pho_PromptGenIsoEnergyStatus1_cone03/F");
                myTree_->Branch("pho_PromptGenIsoEnergyStatus2_cone03",&pho_PromptGenIsoEnergyStatus2_cone03,"pho_PromptGenIsoEnergyStatus2_cone03/F");
                myTree_->Branch("pho_PromptGenIsoEnergyStatus1_cone035",&pho_PromptGenIsoEnergyStatus1_cone035,"pho_PromptGenIsoEnergyStatus1_cone035/F");
                myTree_->Branch("pho_PromptGenIsoEnergyStatus2_cone035",&pho_PromptGenIsoEnergyStatus2_cone035,"pho_PromptGenIsoEnergyStatus2_cone035/F");
                myTree_->Branch("pho_PromptGenIsoEnergyStatus1_cone04",&pho_PromptGenIsoEnergyStatus1_cone04,"pho_PromptGenIsoEnergyStatus1_cone04/F");
                myTree_->Branch("pho_PromptGenIsoEnergyStatus2_cone04",&pho_PromptGenIsoEnergyStatus2_cone04,"pho_PromptGenIsoEnergyStatus2_cone04/F");
                myTree_->Branch("pho_seedSeverity",&pho_seedSeverity,"pho_seedSeverity/I");
                myTree_->Branch("pho_recoFlag",&pho_recoFlag,"pho_recoFlag/I");
                myTree_->Branch("pho_HLT30",&pho_HLT30,"pho_HLT30/I");
                myTree_->Branch("pho_HLT50",&pho_HLT50,"pho_HLT50/I");
                myTree_->Branch("pho_HLTdiPho15",&pho_HLTdiPho15,"pho_HLTdiPho15/I");
                myTree_->Branch("pho_HLTdiPho20",&pho_HLTdiPho20,"pho_HLTdiPho20/I");
                myTree_->Branch("pho_SCeta",&pho_SCeta,"pho_SCeta/F");
                myTree_->Branch("pho_SCphi",&pho_SCphi,"pho_SCphi/F");
                myTree_->Branch("pho_SCEtraw",&pho_SCEtraw,"pho_SCEtraw/F");
                myTree_->Branch("pho_SCEt",&pho_SCEt,"pho_SCEt/F");
                myTree_->Branch("pho_SCr9",&pho_SCr9,"pho_SCr9/F");
                myTree_->Branch("pho_SCbr",&pho_SCbr,"pho_SCbr/F");
                myTree_->Branch("pho_SCnbBC",&pho_SCnbBC,"pho_SCnbBC/I");
                myTree_->Branch("pho_SCnXtal",&pho_SCnXtal,"pho_SCnXtal/I");



}

void endMacro(){
	myFile->Write();
	myFile->Close();
}


//miniTreeMaker(){
int main(){
	cout << "coucou" << endl;
	gSystem->Load("/sps/cms/hbrun/CMSSW_3_6_2_new/src/Morgan/IpnTreeProducer/src/libToto.so");
        inputEventTree->Add("/sps/cms/hbrun/dataset_3_6_1_patch7/Run2010A-PromptReco-v4//run144011_data_EG_goodVtx_noscrapping_2_1_nFO.root");
//        inputRunTree->Add("/sps/cms/hbrun/dataset_3_6_1_patch7/Run2010A-PromptReco-v4//run144011_data_EG_goodVtx_noscrapping_2_1_nFO.root");


	myFile=new TFile("theMiniTree.root","RECREATE");

	beginMacro();

	int NbEvents = inputEventTree->GetEntries();	cout << "NbEvents = " << NbEvents << endl;
	//NbEvents = 100;
	int NbHLT20 = 0;
	int NbPhotons = 0;
	int NbPhotonsCleaned = 0;
	int NbPhotonsAccept = 0;
	int NbPhotonsCutEt = 0;
	int NbPhotonsPb = 0;
	
//	cout << "trigger = " << runInfos->nHLTPaths() << endl;
/// Find the trigger bits !
        int *bits;
        int NbRunEntries = inputRunTree->GetEntries();
        if (NbRunEntries > 0) {
                inputRunTree->GetEntry(0);
                bits = InitializeHLTinfo(runInfos, ListWantedHLTnames, nbHlt);
        }
        if (NbRunEntries > 1 ) {
	        int *bits2;
        	for (int j = 1 ; j < NbRunEntries ; j++){
                	inputRunTree->GetEntry(j);
	                bits2 = InitializeHLTinfo(runInfos, ListWantedHLTnames, nbHlt);
	                if (bits2[0]!=bits[0]) cout << "Alert ! ! = les bits ont changes ! " << endl;
        	}
       }

	for (int ievt  = 0 ; ievt < NbEvents ; ievt++){
	      inputEventTree->GetEvent(ievt);
	      	if( (ievt%10==0 && ievt<=100)  || (ievt%100==0 && ievt<=1000)   || (ievt%1000==0 && ievt>1000)  )
		{
		  cout <<"Analyzing "<< ievt << "th event: " << endl;
		}
		if (doHLT){
			if (bits[0] == -1) break;
			if (!(event->trigHLT(bits[0])==1)) continue;
//			if (!((event->trigHLT(bits[0])==1)||(event->trigHLT(bits[5])==1))) continue; // continue if not the bit photon20
			NbHLT20++;
			if (bits[1]>-1) {if (event->trigHLT(bits[1])) pho_HLT30=1; else pho_HLT30=0;} else pho_HLT30=-1;
			if (bits[2]>-1) {if (event->trigHLT(bits[2])) pho_HLT50=1; else pho_HLT50=0;} else pho_HLT50=-1;
			if (bits[3]>-1) {if (event->trigHLT(bits[3])) pho_HLTdiPho15=1; else pho_HLTdiPho15=0;} else pho_HLTdiPho15=-1;
			if (bits[4]>-1) {if (event->trigHLT(bits[4])) pho_HLTdiPho20=1; else pho_HLTdiPho20=0;} else pho_HLTdiPho20=-1;
		}
	      
		pho_nVertex = vertices->GetEntriesFast();		

		for (int iphoton=0; iphoton< photons->GetEntriesFast(); iphoton++){
			NbPhotons++;
			TRootPhoton *myphoton = (TRootPhoton*) photons->At(iphoton);
			if (doMC) {
				pho_eventPtHat = event->ptHat();
				doGenInfo(myphoton, mcParticles, &(pho_GenId), &(pho_MotherId), &(pho_isGenElectron), &(pho_isPromptGenPho), &(pho_isFromQuarkGen), &(pho_isPi0Gen), &(pho_isEtaGen), &(pho_isRhoGen), &(pho_isOmegaGen), &(pho_PromptGenIsoEnergyStatus1_cone02), &(pho_PromptGenIsoEnergyStatus2_cone02), 0.2);
				doGenInfo(myphoton, mcParticles, &(pho_GenId), &(pho_MotherId), &(pho_isGenElectron), &(pho_isPromptGenPho), &(pho_isFromQuarkGen), &(pho_isPi0Gen), &(pho_isEtaGen), &(pho_isRhoGen), &(pho_isOmegaGen), &(pho_PromptGenIsoEnergyStatus1_cone03), &(pho_PromptGenIsoEnergyStatus2_cone03), 0.3);
				doGenInfo(myphoton, mcParticles, &(pho_GenId), &(pho_MotherId), &(pho_isGenElectron), &(pho_isPromptGenPho), &(pho_isFromQuarkGen), &(pho_isPi0Gen), &(pho_isEtaGen), &(pho_isRhoGen), &(pho_isOmegaGen), &(pho_PromptGenIsoEnergyStatus1_cone035), &(pho_PromptGenIsoEnergyStatus2_cone035), 0.35);
				doGenInfo(myphoton, mcParticles, &(pho_GenId), &(pho_MotherId), &(pho_isGenElectron), &(pho_isPromptGenPho), &(pho_isFromQuarkGen), &(pho_isPi0Gen), &(pho_isEtaGen), &(pho_isRhoGen), &(pho_isOmegaGen), &(pho_PromptGenIsoEnergyStatus1_cone04), &(pho_PromptGenIsoEnergyStatus2_cone04), 0.4);


			}
			float abs_eta = fabs(myphoton->superCluster()->Eta());
			float scRawEt = myphoton->superCluster()->rawEnergy() * sin(myphoton->superCluster()->Theta());
			if ( (myphoton->superCluster()->seedSeverity()==3)||(myphoton->superCluster()->seedSeverity()==4)||(myphoton->superCluster()->seedRecoFlag()==2)) continue;
			NbPhotonsCleaned++;
			/*if ( abs_eta>2.5 ) continue;
                        if ( abs_eta>1.4442 && abs_eta<1.566 ) continue;
			NbPhotonsAccept++;
			if ((scRawEt < 10)&&(myphoton->Et()>20)) NbPhotonsPb++;
			if (scRawEt < 10) continue;
			NbPhotonsCutEt++;*/

			pho_iEvent = ievt;
			pho_EventNumber = event->eventId();
			pho_RunNumber = event->runId();
			pho_lumiSection = event->luminosityBlock();
		 	pho_isEB = myphoton->isEBPho();
			pho_isEE = myphoton->isEEPho();
	 		pho_isEEP = (myphoton->isEEPho() && myphoton->Eta()>0);
			pho_isEEM = (myphoton->isEEPho() && myphoton->Eta()<0);
			pho_isAlsoElectron = myphoton->isAlsoElectron();
			pho_isLoose = myphoton->isLoosePhoton();
			pho_isTight = myphoton->isTightPhoton();
			pho_isEBGap = myphoton->isEBGap();
			pho_isEEGap = myphoton->isEEGap();
			pho_hasPixelSeed = myphoton->hasPixelSeed();
			pho_seedSeverity = myphoton->superCluster()->seedSeverity(); 
			pho_recoFlag = myphoton->superCluster()->seedRecoFlag();

	               //Covariances, cluster width...
	               pho_cEE = myphoton->covEtaEta();
	               pho_cEP = myphoton->covEtaPhi();
	               pho_cPP = myphoton->covPhiPhi();


		      //Energy
		      pho_et = myphoton->Et();
		      pho_energy = myphoton->Energy();
		      pho_px = myphoton->Px();
		      pho_py = myphoton->Py();
		      pho_pz = myphoton->Pz();	      
	
		      //Cluster shape
		      pho_eta = myphoton->Eta();
		      pho_phi = myphoton->Phi();
		      pho_s4 = myphoton->superCluster()->s4();
		      pho_s9 = myphoton->e3x3();
		      pho_eMax = myphoton->eMax();
		      pho_e2nd = myphoton->e2nd();
		      pho_r19 = myphoton->r19();
		      pho_r9 = myphoton->r9();
		      if (myphoton->e3x3()-myphoton->eMax()-myphoton->e2nd()!=0) pho_S9overS9minusS1S2 = myphoton->e3x3()/(myphoton->e3x3()-myphoton->eMax()-myphoton->e2nd());
		      else pho_S9overS9minusS1S2 = -999;
	

		      pho_phiwidth = myphoton->superCluster()->phiWidth();
		      pho_etawidth = myphoton->superCluster()->etaWidth();
		      pho_sigmaIetaIeta = myphoton->sigmaIetaIeta();
		      pho_sigmaEtaEta = myphoton->sigmaEtaEta();	
	
		      pho_sigmaphi = -1;
		      pho_sigmaeta = -1;
		      if (myphoton->superCluster()!=0){
			int tagseed = myphoton->superCluster()->seedBasicClusterIndex();
			TRootCluster* myseedcluster = (TRootCluster*) clusters->At(tagseed);
			if (myseedcluster!=NULL){
			  std::vector<Int_t> clustersIndex = myphoton->superCluster()->subBasicClusterIndexVector();
			  pho_sigmaphi = 0;
			  pho_sigmaeta = 0;
			  for (unsigned int idx=0; idx<myphoton->superCluster()->nBasicClusters(); idx++){
			    TRootCluster* mybasiccluster = (TRootCluster*) clusters->At(clustersIndex[idx]);
			    if (mybasiccluster->e5x5()>0){
  			      double dphi = mybasiccluster->Phi()-myseedcluster->Phi();
			      if (dphi>TMath::Pi()) dphi=2*TMath::Pi()-dphi;
			      if (dphi<-TMath::Pi()) dphi=-2*TMath::Pi()-dphi;
			      pho_sigmaphi += sqrt(mybasiccluster->e5x5()/myseedcluster->e5x5()*dphi*dphi);
	
			      double deta = mybasiccluster->Eta()-myseedcluster->Eta();
			      pho_sigmaeta += sqrt(mybasiccluster->e5x5()/myseedcluster->e5x5()*deta*deta);
			    }
			  }
			}
		      }
	
		      //isolation
		      pho_hoe = myphoton->hoe();
		      pho_IsoEcalRechit = myphoton->dR04IsolationEcalRecHit();
		      pho_IsoHcalRechit = myphoton->dR04IsolationHcalRecHit();
		      pho_IsoSolidTrkCone = myphoton->dR04IsolationSolidTrkCone();
		      pho_IsoHollowTrkCone = myphoton->dR04IsolationHollowTrkCone();
		      pho_IsoEcalRechit03 = myphoton->dR03IsolationEcalRecHit();
		      pho_IsoHcalRechit03 = myphoton->dR03IsolationHcalRecHit();
		      pho_IsoSolidTrkCone03 = myphoton->dR03IsolationSolidTrkCone();
		      pho_IsoHollowTrkCone03 = myphoton->dR03IsolationHollowTrkCone();
		      pho_esRatio = myphoton->superCluster()->esRatio();
	
		      if (pho_isAlsoElectron==true){
			pho_ElectronClassification = -1;
			for (int ielec=0; ielec<electrons->GetEntries(); ielec++){
			  TRootElectron* myelectron = (TRootElectron*) electrons->At(ielec);
			  if (myphoton->superClusterRawEnergy()==myelectron->superClusterRawEnergy())
			    pho_ElectronClassification = myelectron->classification();
	
			}
		      }
	
		      //conversions
		      pho_convNTracks = myphoton->convNTracks();
	
	
		      //Environment variables for gamma/pi0 disc
		      pho_jetEMfraction = GetClosestJetEMFraction(myphoton, jets, 15);
		      pho_DrJetClosest = GetDeltaRClosestJet(myphoton, jets, 15);
	
	 	      pho_ptoverjetpt =  GetPtOverJetPt(myphoton, jets, 0);
	              pho_ptoverjetpt_pt2 = GetPtOverJetPt(myphoton, jets, 2);
	              pho_ptoverjetpt_pt5 = GetPtOverJetPt(myphoton, jets, 5);
	              pho_ptoverjetpt_pt10 = GetPtOverJetPt(myphoton, jets, 10);
		      pho_ptoverjetpt_pt15 = GetPtOverJetPt(myphoton, jets, 15);
	
		      pho_ptoverjetpt_PFlow = GetPtOverJetPt(myphoton, pflowjets, 0);
	              pho_ptoverjetpt_PFlow_pt2 = GetPtOverJetPt(myphoton, pflowjets, 2);
	              pho_ptoverjetpt_PFlow_pt5 = GetPtOverJetPt(myphoton, pflowjets, 5);
        	      pho_ptoverjetpt_PFlow_pt10 = GetPtOverJetPt(myphoton, pflowjets, 10);

		      pho_ptoverjetpt_sisCone = GetPtOverJetPt(myphoton, sisconejets, 0);
        	      pho_ptoverjetpt_sisCone_pt2 = GetPtOverJetPt(myphoton, sisconejets, 2);
	              pho_ptoverjetpt_sisCone_pt5 = GetPtOverJetPt(myphoton, sisconejets, 5);
	              pho_ptoverjetpt_sisCone_pt10 = GetPtOverJetPt(myphoton, sisconejets, 10);

		      pho_transverseMomentumToJetDirection = GetTransverseMomentumToJetDirection(myphoton, jets,0);
		      pho_transverseMomentumToJetDirection_pt2 = GetTransverseMomentumToJetDirection(myphoton, jets,2);
		      pho_transverseMomentumToJetDirection_pt5 = GetTransverseMomentumToJetDirection(myphoton, jets,5);
		      pho_transverseMomentumToJetDirection_pt10 = GetTransverseMomentumToJetDirection(myphoton, jets,10);
		      pho_transverseMomentumToJetDirection_pt20 = GetTransverseMomentumToJetDirection(myphoton, jets,20);
	
		      pho_transverseToJetRatio = pho_transverseMomentumToJetDirection/pho_et;
		      pho_transverseToJetRatio_pt2 = pho_transverseMomentumToJetDirection_pt2/pho_et;
		      pho_transverseToJetRatio_pt5 = pho_transverseMomentumToJetDirection_pt5/pho_et;
		      pho_transverseToJetRatio_pt10 = pho_transverseMomentumToJetDirection_pt10/pho_et;
		      pho_transverseToJetRatio_pt20 = pho_transverseMomentumToJetDirection_pt20/pho_et;


		      //conpho_closestSC_dR
		      double dR=10;
		      for (int isc=0; isc<superClusters->GetEntriesFast(); isc++){
			TRootSuperCluster* mysupercluster = (TRootSuperCluster*) superClusters->At(isc);
		        if (mysupercluster!=NULL){ 
	 		  if (myphoton->superCluster()->rawEnergy()!=mysupercluster->rawEnergy()){
			    double dRtmp = DeltaR(myphoton->Phi(), mysupercluster->Phi(), myphoton->Eta(), mysupercluster->Eta());
			    if (dRtmp<dR) dR=dRtmp;
			  }
		        }
		      }
		      pho_DrSCclosest = dR;
		      
		      //conpho_dR_SCtrkclosest
		      dR=10;
		      double dR_pt2=10, dR_pt5=10, dR_pt10=10;
		      for (unsigned int itk=0; itk<tracks->GetEntriesFast(); itk++){
			TRootTrack* mytrack = (TRootTrack*) tracks->At(itk);
			//if (mytrack->Et()>10){
			double dRtmp = DeltaR(myphoton->Phi(), mytrack->Phi(), myphoton->Eta(), mytrack->Eta());
			if (dRtmp<dR) dR=dRtmp;
			if (dRtmp<dR_pt2 && mytrack->Et()>2) dR_pt2=dRtmp;
	                if (dRtmp<dR_pt5 && mytrack->Et()>5) dR_pt5=dRtmp;
	                if (dRtmp<dR_pt10 && mytrack->Et()>10) dR_pt10=dRtmp;
			//}
		      }
		      pho_DrTrkclosest = dR;
		      pho_DrTrkclosest_pt2 = dR_pt2;
	              pho_DrTrkclosest_pt5 = dR_pt5;
	              pho_DrTrkclosest_pt10 = dR_pt10;

			// now SC variables 
			pho_SCphi = myphoton->superCluster()->Phi();
			pho_SCeta = myphoton->superCluster()->Eta();
			pho_SCEtraw = scRawEt;
			pho_SCEt = myphoton->superCluster()->Pt();
			pho_SCr9 = myphoton->r9();
			if (myphoton->superCluster()->etaWidth()!=0) pho_SCbr = myphoton->superCluster()->phiWidth()/myphoton->superCluster()->etaWidth(); else pho_SCbr = -1;
			pho_SCnbBC = myphoton->superCluster()->nBasicClusters();
			pho_SCnXtal = myphoton->superCluster()->nXtals();
				
				myTree_->Fill();

		}

	}


	cout << "NbEvents=" << NbEvents << " NbEventsHLT20= " << NbHLT20 << endl;
	cout << "NbPhotons=" << NbPhotons << " NbPhotonsCleaned=" << NbPhotonsCleaned <<  " NbPhotonsCutEt=" << NbPhotonsCutEt << " NbPhotonsPb=" << NbPhotonsPb << endl;
	endMacro();
}