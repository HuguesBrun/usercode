#include "miniTreeMaker.h"

void beginMacro(){

	doHLT                    = true;
	doHLTobject		 = true;
  	doMC                     = false;
  	doJetMC                  = false;
  	doMETMC                  = false;
  	doPDFInfo                = true;
  	doSignalMuMuGamma        = false;
	doLeadingPhoton		 = true;
  	doSignalTopTop           = false;
  	doPhotonConversionMC     = false;
  	doElectronConversionMC   = false;
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
	doDiphotons		 = true;


	mcParticles = new TClonesArray("TRootMCParticle", 0);
	genJets = new TClonesArray("TRootParticle", 0);
	genMETs = new TClonesArray("TRootParticle", 0);
	mcPhotons = new TClonesArray("TRootMCPhoton", 0);
	mcElectrons = new TClonesArray("TRootMCElectron", 0);
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
	HLTObjects = new TClonesArray("TRootHLTObject", 0);


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

  if (doElectronConversionMC)
    {
      inputEventTree->SetBranchAddress("MCElectrons",&mcElectrons, &mcElectrons_br);
      inputEventTree->SetBranchStatus("MCElectrons",1);
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
      inputEventTree->SetBranchAddress("muons", &muons, &muons_br);
      inputEventTree->SetBranchStatus("muons", 1);
    }
  
  if(doElectron)
    {
      inputEventTree->SetBranchAddress("gsfElectrons", &electrons, &electrons_br);
      inputEventTree->SetBranchStatus("gsfElectrons", 1);
    }
  
  if(doPhoton)
    {
      inputEventTree->SetBranchAddress("photons", &photons, &photons_br);
      inputEventTree->SetBranchStatus("photons", 1);
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
      inputEventTree->SetBranchAddress("met", &met, &met_br);
      inputEventTree->SetBranchStatus("met", 1);
    }
  
  if(doBardak)
    {
      inputEventTree->SetBranchAddress("bardak", &bardak, &bardak_br);
      inputEventTree->SetBranchStatus("bardak", 1);
    }
  if (doHLTobject)
    {
      inputEventTree->SetBranchAddress("HLTObjects",&HLTObjects, &HLTObjects_br);
      inputEventTree->SetBranchStatus("HLTObjects", 1);
    }



	myTree_ = new TTree("photons","");
               myTree_->Branch("pho_iEvent",&pho_iEvent,"pho_iEvent/I");
                myTree_->Branch("pho_EventNumber",&pho_EventNumber,"pho_EventNumber/I");
                myTree_->Branch("pho_RunNumber",&pho_RunNumber,"pho_RunNumber/I");
		myTree_->Branch("pho_lumiSection",&pho_lumiSection,"pho_lumiSection/I");
		myTree_->Branch("pho_eventProcessId",&pho_eventProcessId,"pho_eventProcessId/I");
		myTree_->Branch("pho_rho",&pho_rho,"pho_rho/F");
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
		myTree_->Branch("pho_e2x2",&pho_e2x2,"pho_e2x2/F");
		myTree_->Branch("pho_e5x5",&pho_e5x5,"pho_e5x5/F");
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
		myTree_->Branch("pho_IsoSolidNtrackCone",&pho_IsoSolidNtrackCone,"pho_IsoSolidNtrackCone/I");
		myTree_->Branch("pho_IsoHollowNtrackCone",&pho_IsoHollowNtrackCone,"pho_IsoHollowNtrackCone/I");
		myTree_->Branch("pho_IsoNNiceTracks",&pho_IsoNNiceTracks,"pho_IsoNNiceTracks/I");
                myTree_->Branch("pho_IsoEcalRechit03",&pho_IsoEcalRechit03,"pho_IsoEcalRechit03/F");
                myTree_->Branch("pho_IsoHcalRechit03",&pho_IsoHcalRechit03,"pho_IsoHcalRechit03/F");
                myTree_->Branch("pho_IsoSolidTrkCone03",&pho_IsoSolidTrkCone03,"pho_IsoSolidTrkCone03/F");
                myTree_->Branch("pho_IsoHollowTrkCone03",&pho_IsoHollowTrkCone03,"pho_IsoHollowTrkCone03/F");
                myTree_->Branch("pho_esRatio",&pho_esRatio,"pho_esRatio/F");
		myTree_->Branch("pho_nbOther",&pho_nbOther,"pho_nbOther/I");
		myTree_->Branch("pho_nbOtherIso",&pho_nbOtherIso,"pho_nbOtherIso/I");
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
				myTree_->Branch("pho_isAlsoRecoAsElectron",&pho_isAlsoRecoAsElectron,"pho_isAlsoRecoAsElectron/I");
				myTree_->Branch("pho_fBrem",&pho_fBrem,"pho_fBrem/F");
				myTree_->Branch("pho_momentumCorrected",&pho_momentumCorrected,"pho_momentumCorrected/F");
				myTree_->Branch("pho_d0",&pho_d0,"pho_d0/F");
				myTree_->Branch("pho_isAlsoRecoAsJet",&pho_isAlsoRecoAsJet,"pho_isAlsoRecoAsJet/I");
                 myTree_->Branch("pho_matchWithEle",&pho_matchWithEle,"pho_matchWithEle/F");
                 myTree_->Branch("pho_tightEleId",&pho_tightEleId,"pho_tightEleId/F");
                 myTree_->Branch("pho_eleTrkIso",&pho_eleTrkIso,"pho_eleTrkIso/F");
                 myTree_->Branch("pho_eleEcalIso",&pho_eleEcalIso,"pho_eleEcalIso/F");
                 myTree_->Branch("pho_eleHcalIso",&pho_eleHcalIso,"pho_eleHcalIso/F");
                 myTree_->Branch("pho_eleDeltaPhiIn",&pho_eleDeltaPhiIn,"pho_eleDeltaPhiIn/F");
                 myTree_->Branch("pho_eleDeltaEtaIn",&pho_eleDeltaEtaIn,"pho_eleDeltaEtaIn/F");
                 myTree_->Branch("pho_eleHoE",&pho_eleHoE,"pho_eleHoE/F");
                 myTree_->Branch("pho_eleSigmaIetaIeta",&pho_eleSigmaIetaIeta,"pho_eleSigmaIetaIeta/F");
                 myTree_->Branch("pho_eleMissHits",&pho_eleMissHits,"pho_eleMissHits/I");
	myTree_->Branch("pho_eleDistConvPartner",&pho_eleDistConvPartner,"pho_eleDistConvPartner/F");
	myTree_->Branch("pho_eleDcotConvPartner",&pho_eleDcotConvPartner,"pho_eleDcotConvPartner/F");
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
                myTree_->Branch("pho_nGenVertex",&pho_nGenVertex,"pho_nGenVertex/I");
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
		myTree_->Branch("pho_seedEnergy",&pho_seedEnergy,"pho_seedEnergy/F");
		myTree_->Branch("pho_seedTime",&pho_seedTime,"pho_seedTime/F");
                myTree_->Branch("pho_HLT_bit0",&pho_HLT_bit0,"pho_HLT_bit0/I");
                myTree_->Branch("pho_HLT_bit1",&pho_HLT_bit1,"pho_HLT_bit1/I");
                myTree_->Branch("pho_HLT_bit2",&pho_HLT_bit2,"pho_HLT_bit2/I");
                myTree_->Branch("pho_HLT_bit3",&pho_HLT_bit3,"pho_HLT_bit3/I");
                myTree_->Branch("pho_HLT_bit4",&pho_HLT_bit4,"pho_HLT_bit4/I");
                myTree_->Branch("pho_HLT_bit5",&pho_HLT_bit5,"pho_HLT_bit5/I");
                myTree_->Branch("pho_HLT_bit6",&pho_HLT_bit6,"pho_HLT_bit6/I");
                myTree_->Branch("pho_HLT_bit7",&pho_HLT_bit7,"pho_HLT_bit7/I");
                myTree_->Branch("pho_HLT_bit8",&pho_HLT_bit8,"pho_HLT_bit8/I");
                myTree_->Branch("pho_HLT_bit9",&pho_HLT_bit9,"pho_HLT_bit9/I");
                myTree_->Branch("pho_HLT_bit10",&pho_HLT_bit10,"pho_HLT_bit10/I");
                myTree_->Branch("pho_HLT_bit11",&pho_HLT_bit11,"pho_HLT_bit11/I");
                myTree_->Branch("pho_HLT_bit12",&pho_HLT_bit12,"pho_HLT_bit12/I");
                myTree_->Branch("pho_HLT_bit13",&pho_HLT_bit13,"pho_HLT_bit13/I");
                myTree_->Branch("pho_HLT_bit14",&pho_HLT_bit14,"pho_HLT_bit14/I");
                myTree_->Branch("pho_HLT_bit15",&pho_HLT_bit15,"pho_HLT_bit15/I");
                myTree_->Branch("pho_HLT_bit16",&pho_HLT_bit16,"pho_HLT_bit16/I");
                myTree_->Branch("pho_HLT_bit17",&pho_HLT_bit17,"pho_HLT_bit17/I");
                myTree_->Branch("pho_HLT_bit18",&pho_HLT_bit18,"pho_HLT_bit18/I");
                myTree_->Branch("pho_HLT_bit19",&pho_HLT_bit19,"pho_HLT_bit19/I");
                myTree_->Branch("pho_HLT_bit20",&pho_HLT_bit20,"pho_HLT_bit20/I");
                myTree_->Branch("pho_HLT_bit21",&pho_HLT_bit21,"pho_HLT_bit21/I");
                myTree_->Branch("pho_HLT_bit22",&pho_HLT_bit22,"pho_HLT_bit22/I");
                myTree_->Branch("pho_HLT_bit23",&pho_HLT_bit23,"pho_HLT_bit23/I");
                myTree_->Branch("pho_HLT_bit24",&pho_HLT_bit24,"pho_HLT_bit24/I");
                myTree_->Branch("pho_HLT_bit25",&pho_HLT_bit25,"pho_HLT_bit25/I");
                myTree_->Branch("pho_HLT_bit26",&pho_HLT_bit26,"pho_HLT_bit26/I");
                myTree_->Branch("pho_HLT_bit27",&pho_HLT_bit27,"pho_HLT_bit27/I");
                myTree_->Branch("pho_HLT_bit28",&pho_HLT_bit28,"pho_HLT_bit28/I");
                myTree_->Branch("pho_HLT_bit29",&pho_HLT_bit29,"pho_HLT_bit29/I");
                myTree_->Branch("pho_HLT_bit30",&pho_HLT_bit30,"pho_HLT_bit30/I");
                myTree_->Branch("pho_HLT_bit31",&pho_HLT_bit31,"pho_HLT_bit31/I");
                myTree_->Branch("pho_HLT_bit32",&pho_HLT_bit32,"pho_HLT_bit32/I");
                myTree_->Branch("pho_HLT_bit33",&pho_HLT_bit33,"pho_HLT_bit33/I");
                myTree_->Branch("pho_HLT_bit34",&pho_HLT_bit34,"pho_HLT_bit34/I");
                myTree_->Branch("pho_HLT_bit35",&pho_HLT_bit35,"pho_HLT_bit35/I");
                myTree_->Branch("pho_HLT_bit36",&pho_HLT_bit36,"pho_HLT_bit36/I");
                myTree_->Branch("pho_HLT_bit37",&pho_HLT_bit37,"pho_HLT_bit37/I");
                myTree_->Branch("pho_HLT_bit38",&pho_HLT_bit38,"pho_HLT_bit38/I");
                myTree_->Branch("pho_HLT_bit39",&pho_HLT_bit39,"pho_HLT_bit39/I");
                myTree_->Branch("pho_HLT_bit40",&pho_HLT_bit40,"pho_HLT_bit40/I");
                myTree_->Branch("pho_HLT_bit41",&pho_HLT_bit41,"pho_HLT_bit41/I");
                myTree_->Branch("pho_HLT_bit42",&pho_HLT_bit42,"pho_HLT_bit42/I");
                myTree_->Branch("pho_HLT_bit43",&pho_HLT_bit43,"pho_HLT_bit43/I");
                myTree_->Branch("pho_HLT_bit44",&pho_HLT_bit44,"pho_HLT_bit44/I");
                myTree_->Branch("pho_HLT_bit45",&pho_HLT_bit45,"pho_HLT_bit45/I");
                myTree_->Branch("pho_SCeta",&pho_SCeta,"pho_SCeta/F");
                myTree_->Branch("pho_SCphi",&pho_SCphi,"pho_SCphi/F");
                myTree_->Branch("pho_SCEtraw",&pho_SCEtraw,"pho_SCEtraw/F");
		myTree_->Branch("pho_SCEraw",&pho_SCEraw,"pho_SCEraw/F");
                myTree_->Branch("pho_SCEt",&pho_SCEt,"pho_SCEt/F");
                myTree_->Branch("pho_SCr9",&pho_SCr9,"pho_SCr9/F");
                myTree_->Branch("pho_SCbr",&pho_SCbr,"pho_SCbr/F");
                myTree_->Branch("pho_SCnbBC",&pho_SCnbBC,"pho_SCnbBC/I");
                myTree_->Branch("pho_SCnXtal",&pho_SCnXtal,"pho_SCnXtal/I");
		myTree_->Branch("isAspike",&isAspike,"isAspike/I");
		myTree_->Branch("pho_etaLAT",&pho_etaLAT,"pho_etaLAT/F");
		myTree_->Branch("pho_phiLAT",&pho_phiLAT,"pho_phiLAT/F");
		myTree_->Branch("pho_LAT",&pho_LAT,"pho_LAT/F");
		myTree_->Branch("pho_Zernike20",&pho_Zernike20,"pho_Zernike20/F");
		myTree_->Branch("pho_Zernike42",&pho_Zernike42,"pho_Zernike42/F");
		myTree_->Branch("pho_secondMomentMaj",&pho_secondMomentMaj,"pho_secondMomentMaj/F");
		myTree_->Branch("pho_secondMomentMin",&pho_secondMomentMin,"pho_secondMomentMin/F");
		myTree_->Branch("pho_secondMomentAlpha",&pho_secondMomentAlpha,"pho_secondMomentAlpha/F");
		myTree_->Branch("pho_ESratio",&pho_ESratio,"pho_ESratio/F");
		myTree_->Branch("pho_trueE",&pho_trueE,"pho_trueE/F");
		myTree_->Branch("pho_truePx",&pho_truePx,"pho_truePx/F");
		myTree_->Branch("pho_truePy",&pho_truePy,"pho_truePy/F");
		myTree_->Branch("pho_truePz",&pho_truePz,"pho_truePz/F");
		myTree_->Branch("pho_trueEta",&pho_trueEta,"pho_trueEta/F");
		myTree_->Branch("pho_truePhi",&pho_truePhi,"pho_truePhi/F");
		myTree_->Branch("pho_isMatchingWithHLTObject",&pho_isMatchingWithHLTObject,"pho_isMatchingWithHLTObject/I");
		myTree_->Branch("pho_isConverted",&pho_isConverted,"pho_isConverted/I");
		myTree_->Branch("pho_NtrackConv",&pho_NtrackConv,"pho_NtrackConv/I");
		myTree_->Branch("pho_convEoverP",&pho_convEoverP,"pho_convEoverP/F");
		myTree_->Branch("pho_convMass",&pho_convMass,"pho_convMass/F");
		myTree_->Branch("pho_convCotanTheta",&pho_convCotanTheta,"pho_convCotanTheta/F");
		myTree_->Branch("pho_convLikely",&pho_convLikely,"pho_convLikely/F");
		myTree_->Branch("pho_convVertexX",&pho_convVertexX,"pho_convVertexX/F");
		myTree_->Branch("pho_convVertexY",&pho_convVertexY,"pho_convVertexY/F");
		myTree_->Branch("pho_convVertexZ",&pho_convVertexZ,"pho_convVertexZ/F");
		myTree_->Branch("pho_MCisConverted",&pho_MCisConverted,"pho_MCisConverted/I");
		myTree_->Branch("pho_MCconvEoverP",&pho_MCconvEoverP,"pho_MCconvEoverP/F");
		myTree_->Branch("pho_MCconvCotanTheta",&pho_MCconvCotanTheta,"pho_MCconvCotanTheta/F");
		myTree_->Branch("pho_MCconvVertexX",&pho_MCconvVertexX,"pho_MCconvVertexX/F");
		myTree_->Branch("pho_MCconvVertexY",&pho_MCconvVertexY,"pho_MCconvVertexY/F");
		myTree_->Branch("pho_MCconvVertexZ",&pho_MCconvVertexZ,"pho_MCconvVertexZ/F");
		myTree_->Branch("pho_xVertex",&pho_xVertex,"pho_xVertex/F");
		myTree_->Branch("pho_yVertex",&pho_yVertex,"pho_yVertex/F");
		myTree_->Branch("pho_zVertex",&pho_zVertex,"pho_zVertex/F");
		myTree_->Branch("pho_eleMCtruthBrem",&pho_eleMCtruthBrem,"pho_eleMCtruthBrem/F");
		myTree_->Branch("pho_eleMCtruthNBrem",&pho_eleMCtruthNBrem,"pho_eleMCtruthNBrem/I");
		myTree_->Branch("pho_isLeadingPhoton",&pho_isLeadingPhoton,"pho_isLeadingPhoton/I");
		myTree_->Branch("pho_isMatchWithMuon",&pho_isMatchWithMuon,"pho_isMatchWithMuon/I");
}

void endMacro(){
	myFile->Write();
	myFile->Close();
}


//miniTreeMaker(){
int main(){
	cout << "coucou" << endl;
	gSystem->Load("/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/src/libToto.so");
        inputEventTree->Add("../test/MC_EG_goodVtx_noscrapping.root");

	myFile=new TFile("theMiniTree.root","RECREATE");

	beginMacro();

	int NbEvents = inputEventTree->GetEntries();	cout << "NbEvents = " << NbEvents << endl;
//	NbEvents = 500;
	int NbHLT20 = 0;
	int NbPhotons = 0;
	int NbPhotonsCleaned = 0;
	int NbPhotonsAccept = 0;
	int NbPhotonsCutEt = 0;
	int NbPhotonsPb = 0;
	

	cout << "nomFichier=" << inputEventTree->GetFile()->GetName() << " nbEntries=" << NbEvents << endl;
	for (int ievt  = 0 ; ievt < NbEvents ; ievt++){
	      inputEventTree->GetEvent(ievt);
	      	if( (ievt%10==0 && ievt<=100)  || (ievt%100==0 && ievt<=1000)   || (ievt%1000==0 && ievt>1000)  )
		{
		  cout <<"Analyzing "<< ievt << "th event: " << endl;
		}
		if (doHLT){
			
			//if (!((event->hltAccept(ListWantedHLTnames[0]))||(event->hltAccept(ListWantedHLTnames[1]))||(event->hltAccept(ListWantedHLTnames[2]))||(event->hltAccept(ListWantedHLTnames[3])))) continue;
			if (nbHlt > 0) {if (event->hltAccept(ListWantedHLTnames[0])) pho_HLT_bit0 = 1; else pho_HLT_bit0 = 0;}
			if (nbHlt > 1) {if (event->hltAccept(ListWantedHLTnames[1])) pho_HLT_bit1 = 1; else pho_HLT_bit1 = 0;}
			if (nbHlt > 2) {if (event->hltAccept(ListWantedHLTnames[2])) pho_HLT_bit2 = 1; else pho_HLT_bit2 = 0;}
			if (nbHlt > 3) {if (event->hltAccept(ListWantedHLTnames[3])) pho_HLT_bit3 = 1; else pho_HLT_bit3 = 0;}
			if (nbHlt > 4) {if (event->hltAccept(ListWantedHLTnames[4])) pho_HLT_bit4 = 1; else pho_HLT_bit4 = 0;}
			if (nbHlt > 5) {if (event->hltAccept(ListWantedHLTnames[5])) pho_HLT_bit5 = 1; else pho_HLT_bit5 = 0;}
			if (nbHlt > 6) {if (event->hltAccept(ListWantedHLTnames[6])) pho_HLT_bit6 = 1; else pho_HLT_bit6 = 0;}
			if (nbHlt > 7) {if (event->hltAccept(ListWantedHLTnames[7])) pho_HLT_bit7 = 1; else pho_HLT_bit7 = 0;}
			if (nbHlt > 8) {if (event->hltAccept(ListWantedHLTnames[8])) pho_HLT_bit8 = 1; else pho_HLT_bit8 = 0;}
			if (nbHlt > 9) {if (event->hltAccept(ListWantedHLTnames[9])) pho_HLT_bit9 = 1; else pho_HLT_bit9 = 0;}
			if (nbHlt > 10) {if (event->hltAccept(ListWantedHLTnames[10])) pho_HLT_bit10 = 1; else pho_HLT_bit10 = 0;}
			if (nbHlt > 11) {if (event->hltAccept(ListWantedHLTnames[11])) pho_HLT_bit11 = 1; else pho_HLT_bit11 = 0;}
			if (nbHlt > 12) {if (event->hltAccept(ListWantedHLTnames[12])) pho_HLT_bit12 = 1; else pho_HLT_bit12 = 0;}
			if (nbHlt > 13) {if (event->hltAccept(ListWantedHLTnames[13])) pho_HLT_bit13 = 1; else pho_HLT_bit13 = 0;}
			if (nbHlt > 14) {if (event->hltAccept(ListWantedHLTnames[14])) pho_HLT_bit14 = 1; else pho_HLT_bit14 = 0;}
			if (nbHlt > 15) {if (event->hltAccept(ListWantedHLTnames[15])) pho_HLT_bit15 = 1; else pho_HLT_bit15 = 0;}
			if (nbHlt > 16) {if (event->hltAccept(ListWantedHLTnames[16])) pho_HLT_bit16 = 1; else pho_HLT_bit16 = 0;}
			if (nbHlt > 17) {if (event->hltAccept(ListWantedHLTnames[17])) pho_HLT_bit17 = 1; else pho_HLT_bit17 = 0;}
			if (nbHlt > 18) {if (event->hltAccept(ListWantedHLTnames[18])) pho_HLT_bit18 = 1; else pho_HLT_bit18 = 0;}
			if (nbHlt > 19) {if (event->hltAccept(ListWantedHLTnames[19])) pho_HLT_bit19 = 1; else pho_HLT_bit19 = 0;}
			if (nbHlt > 20) {if (event->hltAccept(ListWantedHLTnames[20])) pho_HLT_bit20 = 1; else pho_HLT_bit20 = 0;}
			if (nbHlt > 21) {if (event->hltAccept(ListWantedHLTnames[21])) pho_HLT_bit21 = 1; else pho_HLT_bit21 = 0;}
			if (nbHlt > 22) {if (event->hltAccept(ListWantedHLTnames[22])) pho_HLT_bit22 = 1; else pho_HLT_bit22 = 0;}
			if (nbHlt > 23) {if (event->hltAccept(ListWantedHLTnames[23])) pho_HLT_bit23 = 1; else pho_HLT_bit23 = 0;}
			if (nbHlt > 24) {if (event->hltAccept(ListWantedHLTnames[24])) pho_HLT_bit24 = 1; else pho_HLT_bit24 = 0;}
			if (nbHlt > 25) {if (event->hltAccept(ListWantedHLTnames[25])) pho_HLT_bit25 = 1; else pho_HLT_bit25 = 0;}
			if (nbHlt > 26) {if (event->hltAccept(ListWantedHLTnames[26])) pho_HLT_bit26 = 1; else pho_HLT_bit26 = 0;}
			if (nbHlt > 27) {if (event->hltAccept(ListWantedHLTnames[27])) pho_HLT_bit27 = 1; else pho_HLT_bit27 = 0;}
			if (nbHlt > 28) {if (event->hltAccept(ListWantedHLTnames[28])) pho_HLT_bit28 = 1; else pho_HLT_bit28 = 0;}
			if (nbHlt > 29) {if (event->hltAccept(ListWantedHLTnames[29])) pho_HLT_bit29 = 1; else pho_HLT_bit29 = 0;}
			if (nbHlt > 30) {if (event->hltAccept(ListWantedHLTnames[30])) pho_HLT_bit30 = 1; else pho_HLT_bit30 = 0;}
			if (nbHlt > 31) {if (event->hltAccept(ListWantedHLTnames[31])) pho_HLT_bit31 = 1; else pho_HLT_bit31 = 0;}
			if (nbHlt > 32) {if (event->hltAccept(ListWantedHLTnames[32])) pho_HLT_bit32 = 1; else pho_HLT_bit32 = 0;}
			if (nbHlt > 33) {if (event->hltAccept(ListWantedHLTnames[33])) pho_HLT_bit33 = 1; else pho_HLT_bit33 = 0;}
			if (nbHlt > 34) {if (event->hltAccept(ListWantedHLTnames[34])) pho_HLT_bit34 = 1; else pho_HLT_bit34 = 0;}
			if (nbHlt > 35) {if (event->hltAccept(ListWantedHLTnames[35])) pho_HLT_bit35 = 1; else pho_HLT_bit35 = 0;}
			if (nbHlt > 36) {if (event->hltAccept(ListWantedHLTnames[36])) pho_HLT_bit36 = 1; else pho_HLT_bit36 = 0;}
			if (nbHlt > 37) {if (event->hltAccept(ListWantedHLTnames[37])) pho_HLT_bit37 = 1; else pho_HLT_bit37 = 0;}
			if (nbHlt > 38) {if (event->hltAccept(ListWantedHLTnames[38])) pho_HLT_bit38 = 1; else pho_HLT_bit38 = 0;}
			if (nbHlt > 39) {if (event->hltAccept(ListWantedHLTnames[39])) pho_HLT_bit39 = 1; else pho_HLT_bit39 = 0;}
			if (nbHlt > 40) {if (event->hltAccept(ListWantedHLTnames[40])) pho_HLT_bit40 = 1; else pho_HLT_bit40 = 0;}
			if (nbHlt > 41) {if (event->hltAccept(ListWantedHLTnames[41])) pho_HLT_bit41 = 1; else pho_HLT_bit41 = 0;}
			if (nbHlt > 42) {if (event->hltAccept(ListWantedHLTnames[42])) pho_HLT_bit42 = 1; else pho_HLT_bit42 = 0;}
			if (nbHlt > 43) {if (event->hltAccept(ListWantedHLTnames[43])) pho_HLT_bit43 = 1; else pho_HLT_bit43 = 0;}
			if (nbHlt > 44) {if (event->hltAccept(ListWantedHLTnames[44])) pho_HLT_bit44 = 1; else pho_HLT_bit44 = 0;}
			if (nbHlt > 45) {if (event->hltAccept(ListWantedHLTnames[45])) pho_HLT_bit45 = 1; else pho_HLT_bit45 = 0;}
		}
	    	NbHLT20++;
		pho_nVertex = vertices->GetEntriesFast();	
		int nbDipho = 0;
		int nbDiphoIso = 0;
		float maxEnergy = 0;
		if ((doDiphotons)||(doLeadingPhoton)) {
			for (int iphoton=0; iphoton< photons->GetEntriesFast(); iphoton++){
				TRootPhoton *myphoton = (TRootPhoton*) photons->At(iphoton);
				if (myphoton->Et() > secondPhotonCut) {
					nbDipho++;
					if ((myphoton->dR04IsolationEcalRecHit()<IsoEcal)&&(myphoton->dR04IsolationHcalRecHit()<IsoHcal)&&(myphoton->dR04IsolationHollowTrkCone()<IsoTrk)){
						nbDiphoIso++;
					}
				} 
				if (myphoton->Et() > maxEnergy) maxEnergy = myphoton->Et();
				
			}
		}

		for (int iphoton=0; iphoton< photons->GetEntriesFast(); iphoton++){
			NbPhotons++;
			TRootPhoton *myphoton = (TRootPhoton*) photons->At(iphoton);
			if (doMC) {
				pho_eventPtHat = event->ptHat();
				doGenInfo(myphoton, mcParticles, &(pho_GenId), &(pho_MotherId), &(pho_isGenElectron), &(pho_isPromptGenPho), &(pho_isFromQuarkGen), &(pho_isPi0Gen), &(pho_isEtaGen), &(pho_isRhoGen), &(pho_isOmegaGen), &(pho_PromptGenIsoEnergyStatus1_cone02), &(pho_PromptGenIsoEnergyStatus2_cone02),&(pho_trueE),&(pho_truePx),&(pho_truePy),&(pho_truePz),&(pho_trueEta),&(pho_truePhi), 0.2);
				doGenInfo(myphoton, mcParticles, &(pho_GenId), &(pho_MotherId), &(pho_isGenElectron), &(pho_isPromptGenPho), &(pho_isFromQuarkGen), &(pho_isPi0Gen), &(pho_isEtaGen), &(pho_isRhoGen), &(pho_isOmegaGen), &(pho_PromptGenIsoEnergyStatus1_cone03), &(pho_PromptGenIsoEnergyStatus2_cone03), &(pho_trueE),&(pho_truePx),&(pho_truePy),&(pho_truePz),&(pho_trueEta),&(pho_truePhi),0.3);
				doGenInfo(myphoton, mcParticles, &(pho_GenId), &(pho_MotherId), &(pho_isGenElectron), &(pho_isPromptGenPho), &(pho_isFromQuarkGen), &(pho_isPi0Gen), &(pho_isEtaGen), &(pho_isRhoGen), &(pho_isOmegaGen), &(pho_PromptGenIsoEnergyStatus1_cone035), &(pho_PromptGenIsoEnergyStatus2_cone035), &(pho_trueE),&(pho_truePx),&(pho_truePy),&(pho_truePz),&(pho_trueEta),&(pho_truePhi),0.35);
				doGenInfo(myphoton, mcParticles, &(pho_GenId), &(pho_MotherId), &(pho_isGenElectron), &(pho_isPromptGenPho), &(pho_isFromQuarkGen), &(pho_isPi0Gen), &(pho_isEtaGen), &(pho_isRhoGen), &(pho_isOmegaGen), &(pho_PromptGenIsoEnergyStatus1_cone04), &(pho_PromptGenIsoEnergyStatus2_cone04), &(pho_trueE),&(pho_truePx),&(pho_truePy),&(pho_truePz),&(pho_trueEta),&(pho_truePhi),0.4);
				pho_eventProcessId = event->processID();
				pho_eventPtHat = event->ptHat();
			}
			float abs_eta = fabs(myphoton->superCluster()->Eta());
			float scRawEt = myphoton->superCluster()->rawEnergy() * sin(myphoton->superCluster()->Theta());
			if ( ((myphoton->superCluster()->seedSeverity()==5)||(myphoton->superCluster()->seedSeverity()==4)||((myphoton->superCluster()->seedRecoFlag()==2)&&(myphoton->superCluster()->seedEnergy()<130))||((myphoton->superCluster()->seedEnergy()>=130)&&(myphoton->superCluster()->seedTime()<0)&&(myphoton->superCluster()->seedRecoFlag()==2))||(myphoton->sigmaIetaIeta()<0.001)||(TMath::Sqrt(myphoton->covPhiPhi())<0.001))&&myphoton->isEBPho()==1) isAspike =1; 
			else {
				NbPhotonsCleaned++;
				isAspike = 0;
			}
			pho_seedTime = myphoton->superCluster()->seedTime();
			pho_seedEnergy = myphoton->superCluster()->seedEnergy();

			if (myphoton->Et()== maxEnergy) pho_isLeadingPhoton = 1;
			else pho_isLeadingPhoton = 0; 
			pho_nbOther = nbDipho;
			pho_nbOtherIso = nbDiphoIso;
			pho_iEvent = ievt;
			pho_EventNumber = event->eventId();
			pho_RunNumber = event->runId();
			pho_rho = event->rho();
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

		       // new cluster shape variables
			pho_etaLAT = myphoton->etaLAT();
			pho_phiLAT = myphoton->phiLAT();
			pho_LAT = myphoton->lat();
			pho_Zernike20 = myphoton->zernike20();
			pho_Zernike42 = myphoton->zernike42();
			pho_secondMomentMaj = myphoton->secondMomentMaj();
			pho_secondMomentMin = myphoton->secondMomentMin();
			pho_secondMomentAlpha = myphoton->secondMomentAlpha();
			pho_ESratio = myphoton->superCluster()->esRatio();

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
		      pho_e2x2 = myphoton->e2x2();
		      pho_e5x5 = myphoton->e5x5();
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
		      pho_IsoSolidNtrackCone = myphoton->dR04IsolationNTracksSolidCone();
		      pho_IsoHollowNtrackCone = myphoton->dR04IsolationNTracksHollowCone();
		      pho_IsoNNiceTracks = myphoton->dR04IsolationNNiceTracks();
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

			  pho_isAlsoRecoAsElectron = 0;
			  pho_fBrem = -1;
			  pho_momentumCorrected = -1;
			  pho_d0 = -1;
			pho_tightEleId = -1;
			pho_eleTrkIso = -1;
			pho_eleEcalIso = -1;
			pho_eleHcalIso = -1;
			pho_eleDeltaPhiIn = -1;
			pho_eleDeltaEtaIn = -1;
			pho_eleHoE = -1;
			pho_eleSigmaIetaIeta = -1;
			pho_eleMissHits = -1;
			pho_eleMCtruthBrem = -10000;
			pho_eleMCtruthNBrem = -10000;

			  matchWithAnElectron(myphoton, electrons, &pho_isAlsoRecoAsElectron, &pho_fBrem, &pho_momentumCorrected, &pho_d0,&pho_tightEleId, &pho_eleTrkIso, &pho_eleEcalIso, &pho_eleHcalIso, &pho_eleDeltaPhiIn, &pho_eleDeltaEtaIn, &pho_eleHoE, &pho_eleSigmaIetaIeta, &pho_eleMissHits, &pho_eleDistConvPartner, &pho_eleDcotConvPartner);
	
			if ((pho_isAlsoRecoAsElectron==1)&&(doElectronConversionMC)) {
				findTheMCelectron(myphoton,mcElectrons, pho_eleMCtruthBrem, pho_eleMCtruthNBrem);	
			}			
	
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
			pho_SCEraw = myphoton->superCluster()->rawEnergy();
			
			if (doHLTobject){	
				pho_isMatchingWithHLTObject = findMatchingWithAnHLTObjet(myphoton, HLTObjects, theHTLobject);
			}

			// now fill the conversions variables
			pho_isConverted = 0;
			if (myphoton->convNTracks() > 0 ) pho_isConverted = 1; 
			pho_NtrackConv = myphoton->convNTracks();
			pho_convEoverP = myphoton->convEoverP();
			pho_convMass = myphoton->convMass();
			pho_convCotanTheta = myphoton->convCotanTheta();
			pho_convLikely = myphoton->convLikely();
			pho_convVertexX = myphoton->convVertex().x();
			pho_convVertexY = myphoton->convVertex().y();
			pho_convVertexZ = myphoton->convVertex().z();
			
			if (doPhotonConversionMC){
				findConversionMCtruth(myphoton, mcPhotons, pho_MCisConverted, pho_MCconvEoverP, pho_MCconvMass, pho_MCconvCotanTheta, pho_MCconvVertexX, pho_MCconvVertexY, pho_MCconvVertexZ);
			}
			pho_xVertex = myphoton->vx();
			pho_yVertex = myphoton->vy();
			pho_zVertex = myphoton->vz();

			pho_isMatchWithMuon = isMatchingWithAMuon(myphoton, muons);

			myTree_->Fill();

		}

	}


	cout << "NbEventsTot=" << NbEvents << " NbEventsHLT20=" << NbHLT20 << endl;
	cout << "NbPhotonsTot=" << NbPhotons << " NbPhotonsCleaned=" << NbPhotonsCleaned <<  " NbPhotonsCutEt=" << NbPhotonsCutEt << " NbPhotonsPb=" << NbPhotonsPb << endl;
	endMacro();
}
