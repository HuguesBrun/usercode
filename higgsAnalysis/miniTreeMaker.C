#include "miniTreeMaker.h"

void beginMacro(){

	doHLT                    = true;
	doHLTobject		 = true;
  	doMC                     = true;
  	doJetMC                  = false;
  	doMETMC                  = false;
  	doPDFInfo                = true;
  	doSignalMuMuGamma        = false;
	doLeadingPhoton		 = true;
  	doSignalTopTop           = false;
  	doPhotonConversionMC     = true;
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
	doWorstIsolation 	 = true;


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

    	myTree_ = new TTree("diPhotons","DiPhotonsInfos");
                myTree_->Branch("event_number",&event_number,"event_number/I");
		myTree_->Branch("event_runNumber",&event_runNumber, "event_runNumber/I");
		myTree_->Branch("event_LumiSection",&event_LumiSection,"event_LumiSection/I");
		myTree_->Branch("event_eventPtHat",&event_eventPtHat, "event_eventPtHat/F");
		myTree_->Branch("event_processId",&event_processId, "event_processId/I");
		myTree_->Branch("event_nRecoVertex",&event_nRecoVertex,"event_nRecoVertex/I");
		myTree_->Branch("event_nGenInTimeVertex",&event_nGenInTimeVertex,"event_nGenInTimeVertex/I");
		myTree_->Branch("event_nGenOutOfTimeVertex",&event_nGenOutOfTimeVertex,"event_nGenOutOfTimeVertex/I");
		myTree_->Branch("event_nPhotons",&event_nPhotons,"event_nPhotons/I");
		myTree_->Branch("event_rho",&event_rho,"event_rho/F");
             	myTree_->Branch("dipho_mgg",&dipho_mgg,"dipho_mgg/F");
                myTree_->Branch("dipho_qt",&dipho_qt,"dipho_qt/F");                                 
                myTree_->Branch("dipho_ql",&dipho_ql,"dipho_ql/F");                                 
                myTree_->Branch("dipho_deltaR",&dipho_deltaR,"dipho_deltaR/F");                     
                myTree_->Branch("dipho_costhetastar",&dipho_costhetastar,"dipho_costhetastar/F");   
                myTree_->Branch("dipho_eta",&dipho_eta,"dipho_eta/F");                              
                myTree_->Branch("dipho_etastar",&dipho_etastar,"dipho_etastar/F");         
             	myTree_->Branch("diphoMC_mgg",&diphoMC_mgg,"diphoMC_mgg/F");
                myTree_->Branch("diphoMC_qt",&diphoMC_qt,"diphoMC_qt/F");                                 
                myTree_->Branch("diphoMC_ql",&diphoMC_ql,"diphoMC_ql/F");                                 
                myTree_->Branch("diphoMC_deltaR",&diphoMC_deltaR,"diphoMC_deltaR/F");                     
                myTree_->Branch("diphoMC_costhetastar",&diphoMC_costhetastar,"diphoMC_costhetastar/F");   
                myTree_->Branch("diphoMC_eta",&diphoMC_eta,"diphoMC_eta/F");                              
                myTree_->Branch("diphoMC_etastar",&diphoMC_etastar,"diphoMC_etastar/F");        
		myTree_->Branch("pholead_isMatchingWithMC",&pholead_isMatchingWithMC,"pholead_isMatchingWithMC/I"); 
		myTree_->Branch("photrail_isMatchingWithMC",&photrail_isMatchingWithMC,"photrail_isMatchingWithMC/I");

		myTree_->Branch("pholead_pt",&pholead_pt,"pholead_pt/F"); 
		myTree_->Branch("pholead_eta",&pholead_eta,"pholead_eta/F");   
		myTree_->Branch("pholead_etaSC",&pholead_etaSC,"pholead_etaSC/F");

		myTree_->Branch("pholead_r9",&pholead_r9,"pholead_r9/F");                                             
		myTree_->Branch("pholead_cPP",&pholead_cPP,"pholead_cPP/F");
		myTree_->Branch("pholead_cEP",&pholead_cEP,"pholead_cEP/F");
		myTree_->Branch("pholead_cEE",&pholead_cEE,"pholead_cEE/F");
		myTree_->Branch("pholead_r19",&pholead_r19,"pholead_r19/F");                                             
		myTree_->Branch("pholead_SCEraw",&pholead_SCEraw,"pholead_SCEraw/F");                                             
		myTree_->Branch("pholead_eMax",&pholead_eMax,"pholead_eMax/F");                                             
		myTree_->Branch("pholead_e2x2",&pholead_e2x2,"pholead_e2x2/F");                                             
		myTree_->Branch("pholead_e5x5",&pholead_e5x5,"pholead_e5x5/F");                                             
		myTree_->Branch("pholead_ratioSeed",&pholead_ratioSeed,"pholead_ratioSeed/F");
		myTree_->Branch("pholead_ratioS4",&pholead_ratioS4,"pholead_ratioS4/F");
		myTree_->Branch("pholead_lambdaRatio",&pholead_lambdaRatio,"pholead_lambdaRatio/F");
		myTree_->Branch("pholead_lamdbaDivCov",&pholead_lamdbaDivCov,"pholead_lamdbaDivCov/F");
		myTree_->Branch("pholead_secondMomentMaj",&pholead_secondMomentMaj,"pholead_secondMomentMaj/F");
		myTree_->Branch("pholead_secondMomentMin",&pholead_secondMomentMin,"pholead_secondMomentMin/F");
		myTree_->Branch("pholead_secondMomentAlpha",&pholead_secondMomentAlpha,"pholead_secondMomentAlpha/F");
		myTree_->Branch("pholead_covAngle",&pholead_covAngle,"pholead_covAngle/F");
		myTree_->Branch("pholead_covAngle2",&pholead_covAngle2,"pholead_covAngle2/F");
		myTree_->Branch("pholead_S9overS9minusS1S2",&pholead_S9overS9minusS1S2,"pholead_S9overS9minusS1S2/F");
		myTree_->Branch("pholead_etawidth",&pholead_etawidth,"pholead_etawidth/F");                              
		myTree_->Branch("pholead_phiwidth",&pholead_phiwidth,"pholead_phiwidth/F");                              
		myTree_->Branch("pholead_sigieta",&pholead_sigieta,"pholead_sigieta/F");                                 
		myTree_->Branch("pholead_SCbr",&pholead_SCbr,"pholead_SCbr/F");
		myTree_->Branch("pholead_HcalIso",&pholead_HcalIso,"pholead_HcalIso/F");                                     
		myTree_->Branch("pholead_EcalIso",&pholead_EcalIso,"pholead_EcalIso/F");
		myTree_->Branch("pholead_TrackerIso",&pholead_TrackerIso,"pholead_TrackerIso/F");
		myTree_->Branch("pholead_HcalIsodR03",&pholead_HcalIsodR03,"pholead_HcalIsodR03/F");             
		myTree_->Branch("pholead_EcalIsodR03",&pholead_EcalIsodR03,"pholead_EcalIsodR03/F");
		myTree_->Branch("pholead_TrackerIsodR03",&pholead_TrackerIsodR03,"pholead_TrackerIsodR03/F");
		myTree_->Branch("pholead_hoe",&pholead_hoe,"pholead_hoe/F");

                myTree_->Branch("pholead_HasPixSeed",&pholead_HasPixSeed,"pholead_HasPixSeed/I");
                myTree_->Branch("pholead_seedSeverity",&pholead_seedSeverity,"pholead_seedSeverity/I");
                myTree_->Branch("pholead_recoFlag",&pholead_recoFlag,"pholead_recoFlag/I");
		myTree_->Branch("pholead_isEB",&pholead_isEB,"pholead_isEB/I");
		myTree_->Branch("pholead_isEE",&pholead_isEE,"pholead_isEE/I");
		myTree_->Branch("pholead_NNshapeOutput",&pholead_NNshapeOutput, "pholead_NNshapeOutput/F");

		myTree_->Branch("photrail_pt",&photrail_pt,"photrail_pt/F"); 
		myTree_->Branch("photrail_eta",&photrail_eta,"photrail_eta/F");   
		myTree_->Branch("photrail_etaSC",&photrail_etaSC,"photrail_etaSC/F");

		myTree_->Branch("photrail_r9",&photrail_r9,"photrail_r9/F");                                             
		myTree_->Branch("photrail_cPP",&photrail_cPP,"photrail_cPP/F");
		myTree_->Branch("photrail_cEP",&photrail_cEP,"photrail_cEP/F");
		myTree_->Branch("photrail_cEE",&photrail_cEE,"photrail_cEE/F");
		myTree_->Branch("photrail_r19",&photrail_r19,"photrail_r19/F");                                             
		myTree_->Branch("photrail_SCEraw",&photrail_SCEraw,"photrail_SCEraw/F");                                             
		myTree_->Branch("photrail_eMax",&photrail_eMax,"photrail_eMax/F");                                             
		myTree_->Branch("photrail_e2x2",&photrail_e2x2,"photrail_e2x2/F");                                             
		myTree_->Branch("photrail_e5x5",&photrail_e5x5,"photrail_e5x5/F");                                             
		myTree_->Branch("photrail_ratioSeed",&photrail_ratioSeed,"photrail_ratioSeed/F");
		myTree_->Branch("photrail_ratioS4",&photrail_ratioS4,"photrail_ratioS4/F");
		myTree_->Branch("photrail_lambdaRatio",&photrail_lambdaRatio,"photrail_lambdaRatio/F");
		myTree_->Branch("photrail_lamdbaDivCov",&photrail_lamdbaDivCov,"photrail_lamdbaDivCov/F");
		myTree_->Branch("photrail_secondMomentMaj",&photrail_secondMomentMaj,"photrail_secondMomentMaj/F");
		myTree_->Branch("photrail_secondMomentMin",&photrail_secondMomentMin,"photrail_secondMomentMin/F");
		myTree_->Branch("photrail_secondMomentAlpha",&photrail_secondMomentAlpha,"photrail_secondMomentAlpha/F");
		myTree_->Branch("photrail_covAngle",&photrail_covAngle,"photrail_covAngle/F");
		myTree_->Branch("photrail_covAngle2",&photrail_covAngle2,"photrail_covAngle2/F");
		myTree_->Branch("photrail_S9overS9minusS1S2",&photrail_S9overS9minusS1S2,"photrail_S9overS9minusS1S2/F");
		myTree_->Branch("photrail_etawidth",&photrail_etawidth,"photrail_etawidth/F");                              
		myTree_->Branch("photrail_phiwidth",&photrail_phiwidth,"photrail_phiwidth/F");                              
		myTree_->Branch("photrail_sigieta",&photrail_sigieta,"photrail_sigieta/F");                                 
		myTree_->Branch("photrail_SCbr",&photrail_SCbr,"photrail_SCbr/F");
		myTree_->Branch("photrail_HcalIso",&photrail_HcalIso,"photrail_HcalIso/F");                                     
		myTree_->Branch("photrail_EcalIso",&photrail_EcalIso,"photrail_EcalIso/F");
		myTree_->Branch("photrail_TrackerIso",&photrail_TrackerIso,"photrail_TrackerIso/F");
		myTree_->Branch("photrail_HcalIsodR03",&photrail_HcalIsodR03,"photrail_HcalIsodR03/F");             
		myTree_->Branch("photrail_EcalIsodR03",&photrail_EcalIsodR03,"photrail_EcalIsodR03/F");
		myTree_->Branch("photrail_TrackerIsodR03",&photrail_TrackerIsodR03,"photrail_TrackerIsodR03/F");
		myTree_->Branch("photrail_hoe",&photrail_hoe,"photrail_hoe/F");

                myTree_->Branch("photrail_HasPixSeed",&photrail_HasPixSeed,"photrail_HasPixSeed/I");
                myTree_->Branch("photrail_seedSeverity",&photrail_seedSeverity,"photrail_seedSeverity/I");
                myTree_->Branch("photrail_recoFlag",&photrail_recoFlag,"photrail_recoFlag/I");
		myTree_->Branch("photrail_isEB",&photrail_isEB,"photrail_isEB/I");
		myTree_->Branch("photrail_isEE",&photrail_isEE,"photrail_isEE/I");
		myTree_->Branch("photrail_NNshapeOutput",&photrail_NNshapeOutput, "photrail_NNshapeOutput/F");


}

void saveThisEvent(TRootEvent *theEvent, pair <TRootPhoton*, TRootPhoton*> theDiphotonPair, TClonesArray *thePhotonArray, TClonesArray *theVerticeArray){
		event_number = theEvent->eventId();
		event_runNumber = theEvent->runId();
		event_LumiSection = theEvent->luminosityBlock();
		event_eventPtHat = theEvent->ptHat();
		event_processId = theEvent->processID();
		event_rho = theEvent->rho();
		event_nRecoVertex = theVerticeArray->GetEntriesFast();
		event_nPhotons = thePhotonArray->GetEntriesFast();
		event_nGenInTimeVertex = theEvent->nInTimePUVertices();
		event_nGenOutOfTimeVertex = theEvent->nOOTPUVertices();
		// now calculation of the dipho kine variable 
		TLorentzVector Plead, Ptrail, Psum;
		Plead.SetPxPyPzE((theDiphotonPair.first)->Px(),(theDiphotonPair.first)->Py(),(theDiphotonPair.first)->Pz(),(theDiphotonPair.first)->Energy());
		Ptrail.SetPxPyPzE((theDiphotonPair.second)->Px(),(theDiphotonPair.second)->Py(),(theDiphotonPair.second)->Pz(),(theDiphotonPair.second)->Energy());
		Psum = Plead + Ptrail;
		dipho_mgg = Psum.M();
		dipho_qt = Psum.Pt();
		dipho_ql = Psum.Pz();
		dipho_deltaR = DeltaR((theDiphotonPair.first)->Phi(),(theDiphotonPair.second)->Phi(),(theDiphotonPair.first)->Eta(),(theDiphotonPair.second)->Eta());
		dipho_costhetastar = fabs(CosThetaStar(Plead,Ptrail));
		dipho_eta = Psum.Eta();
		dipho_etastar = 1.0*(Plead.Eta()-Ptrail.Eta())/2;
		int catDipho = findTheDiphoCat(*(theDiphotonPair.first), *(theDiphotonPair.second));
		cout << "passing pair : Run = " << event_runNumber << " LS = " << event_LumiSection << " Event = " << event_number << " SelVtx = 0  CAT4 = " << catDipho << " ggM = " << dipho_mgg << " gg_Pt =  " << dipho_qt << endl;
		TRootMCParticle theLeadMC, theTrailMC;
		TLorentzVector PleadMC, PtrailMC, PsumMC;		
		if (findGenParticle(theDiphotonPair.first, mcParticles, &theLeadMC)) pholead_isMatchingWithMC =1;
		else pholead_isMatchingWithMC = 0;

		if (findGenParticle(theDiphotonPair.second, mcParticles, &theTrailMC)) photrail_isMatchingWithMC =1;
		else photrail_isMatchingWithMC = 0;

		if ((pholead_isMatchingWithMC)&&(photrail_isMatchingWithMC)){
			PleadMC.SetPxPyPzE(theLeadMC.Px(),theLeadMC.Py(),theLeadMC.Pz(),theLeadMC.Energy());
			PtrailMC.SetPxPyPzE(theTrailMC.Px(),theTrailMC.Py(),theTrailMC.Pz(),theTrailMC.Energy());
			PsumMC = PleadMC + PtrailMC;
			diphoMC_mgg = PsumMC.M();
	                diphoMC_qt = PsumMC.Pt();
        	        diphoMC_ql = PsumMC.Pz();
                	diphoMC_deltaR = DeltaR(theLeadMC.Phi(),theTrailMC.Phi(),theLeadMC.Eta(),theTrailMC.Eta());
               		diphoMC_costhetastar = fabs(CosThetaStar(PleadMC,PtrailMC));
                	diphoMC_eta = PsumMC.Eta();
                	diphoMC_etastar = 1.0*(PleadMC.Eta()-PtrailMC.Eta())/2;
		}

		
		pholead_r9 = (theDiphotonPair.first)->r9();
		pholead_cPP = (theDiphotonPair.first)->covPhiPhi();
		pholead_cEE = (theDiphotonPair.first)->covEtaEta();
		pholead_cEP = (theDiphotonPair.first)->covEtaPhi();

		pholead_r19 = (theDiphotonPair.first)->r19();
		pholead_SCEraw = (theDiphotonPair.first)->superCluster()->rawEnergy();
		pholead_eMax = (theDiphotonPair.first)->eMax();
		pholead_e2x2 = (theDiphotonPair.first)->e2x2();
		pholead_e5x5 = (theDiphotonPair.first)->e5x5();
		if (pholead_SCEraw!=0) pholead_ratioSeed = pholead_eMax/pholead_SCEraw; else pholead_ratioSeed = 0;
		if (pholead_e5x5!=0) pholead_ratioS4 = pholead_e2x2 / pholead_e5x5; else pholead_ratioS4 = 0;
		if ( (pholead_cEE+pholead_cPP+sqrt((pholead_cEE-pholead_cPP)*(pholead_cEE-pholead_cPP)+4*pholead_cEP*pholead_cEP))!=0) pholead_lambdaRatio = (pholead_cEE+pholead_cPP-sqrt((pholead_cEE-pholead_cPP)*(pholead_cEE-pholead_cPP)+4*pholead_cEP*pholead_cEP))/(pholead_cEE+pholead_cPP+sqrt((pholead_cEE-pholead_cPP)*(pholead_cEE-pholead_cPP)+4*pholead_cEP*pholead_cEP)); else pholead_lambdaRatio = 0;
		if (pholead_cEE != 0) pholead_lamdbaDivCov = (pholead_cEE+pholead_cPP-sqrt((pholead_cEE-pholead_cPP)*(pholead_cEE-pholead_cPP)+4*pholead_cEP*pholead_cEP))/pholead_cEE; else pholead_lamdbaDivCov = 0;
	
		pholead_secondMomentMaj = (theDiphotonPair.first)->secondMomentMaj();
		pholead_secondMomentMin = (theDiphotonPair.first)->secondMomentMin();
		pholead_secondMomentAlpha = (theDiphotonPair.first)->secondMomentAlpha();
		pholead_etawidth = (theDiphotonPair.first)->superCluster()->etaWidth();
		pholead_sigieta = (theDiphotonPair.first)->superCluster()->phiWidth();
		pholead_SCbr = (theDiphotonPair.first)->superCluster()->phiWidth()/(theDiphotonPair.first)->superCluster()->etaWidth();

		pholead_HcalIso = (theDiphotonPair.first)->dR04IsolationHcalRecHit();
		pholead_EcalIso = (theDiphotonPair.first)->dR04IsolationEcalRecHit();
		pholead_TrackerIso = (theDiphotonPair.first)->dR04IsolationHollowTrkCone();
		pholead_HcalIsodR03 = (theDiphotonPair.first)->dR03IsolationHcalRecHit();
		pholead_EcalIsodR03 = (theDiphotonPair.first)->dR03IsolationEcalRecHit();
		pholead_TrackerIsodR03 = (theDiphotonPair.first)->dR03IsolationHollowTrkCone();
		pholead_hoe = (theDiphotonPair.first)->hoe();

		pholead_HasPixSeed = (theDiphotonPair.first)->hasPixelSeed();
		pholead_seedSeverity = (theDiphotonPair.first)->superCluster()->seedSeverity(); 
		pholead_recoFlag = (theDiphotonPair.first)->superCluster()->seedRecoFlag();
		pholead_isEB = (theDiphotonPair.first)->isEBPho();
		pholead_isEE = (theDiphotonPair.first)->isEEPho();
	//	pholead_NNshapeOutput = (theDiphotonPair.first)

		photrail_r9 = (theDiphotonPair.second)->r9();
		photrail_cPP = (theDiphotonPair.second)->covPhiPhi();
		photrail_cEE = (theDiphotonPair.second)->covEtaEta();
		photrail_cEP = (theDiphotonPair.second)->covEtaPhi();

		photrail_r19 = (theDiphotonPair.second)->r19();
		photrail_SCEraw = (theDiphotonPair.second)->superCluster()->rawEnergy();
		photrail_eMax = (theDiphotonPair.second)->eMax();
		photrail_e2x2 = (theDiphotonPair.second)->e2x2();
		photrail_e5x5 = (theDiphotonPair.second)->e5x5();
		if (photrail_SCEraw!=0) photrail_ratioSeed = photrail_eMax/photrail_SCEraw; else photrail_ratioSeed = 0;
		if (photrail_e5x5!=0) photrail_ratioS4 = photrail_e2x2 / photrail_e5x5; else photrail_ratioS4 = 0;
		if ( (photrail_cEE+photrail_cPP+sqrt((photrail_cEE-photrail_cPP)*(photrail_cEE-photrail_cPP)+4*photrail_cEP*photrail_cEP))!=0) photrail_lambdaRatio = (photrail_cEE+photrail_cPP-sqrt((photrail_cEE-photrail_cPP)*(photrail_cEE-photrail_cPP)+4*photrail_cEP*photrail_cEP))/(photrail_cEE+photrail_cPP+sqrt((photrail_cEE-photrail_cPP)*(photrail_cEE-photrail_cPP)+4*photrail_cEP*photrail_cEP)); else photrail_lambdaRatio = 0;
		if (photrail_cEE != 0) photrail_lamdbaDivCov = (photrail_cEE+photrail_cPP-sqrt((photrail_cEE-photrail_cPP)*(photrail_cEE-photrail_cPP)+4*photrail_cEP*photrail_cEP))/photrail_cEE; else photrail_lamdbaDivCov = 0;
	
		photrail_secondMomentMaj = (theDiphotonPair.second)->secondMomentMaj();
		photrail_secondMomentMin = (theDiphotonPair.second)->secondMomentMin();
		photrail_secondMomentAlpha = (theDiphotonPair.second)->secondMomentAlpha();
		photrail_etawidth = (theDiphotonPair.second)->superCluster()->etaWidth();
		photrail_sigieta = (theDiphotonPair.second)->superCluster()->phiWidth();
		photrail_SCbr = (theDiphotonPair.second)->superCluster()->phiWidth()/(theDiphotonPair.second)->superCluster()->etaWidth();

		photrail_HcalIso = (theDiphotonPair.second)->dR04IsolationHcalRecHit();
		photrail_EcalIso = (theDiphotonPair.second)->dR04IsolationEcalRecHit();
		photrail_TrackerIso = (theDiphotonPair.second)->dR04IsolationHollowTrkCone();
		photrail_HcalIsodR03 = (theDiphotonPair.second)->dR03IsolationHcalRecHit();
		photrail_EcalIsodR03 = (theDiphotonPair.second)->dR03IsolationEcalRecHit();
		photrail_TrackerIsodR03 = (theDiphotonPair.second)->dR03IsolationHollowTrkCone();
		photrail_hoe = (theDiphotonPair.second)->hoe();

		photrail_HasPixSeed = (theDiphotonPair.second)->hasPixelSeed();
		photrail_seedSeverity = (theDiphotonPair.second)->superCluster()->seedSeverity(); 
		photrail_recoFlag = (theDiphotonPair.second)->superCluster()->seedRecoFlag();
		photrail_isEB = (theDiphotonPair.second)->isEBPho();
		photrail_isEE = (theDiphotonPair.second)->isEEPho();
	//	photrail_NNshapeOutput = (theDiphotonPair.second)


		myTree_->Fill();
    }


void endMacro(){
//	myTree_->Write();
	myFile->Write();
	myFile->Close();
}


//miniTreeMaker(){
int main(){
	cout << "coucou" << endl;
	gSystem->Load("/sps/cms/hbrun/CMSSW_4_2_3_testElectron/src/Morgan/IpnTreeProducer/src/libToto.so");
      // inputEventTree->Add("../test/MC_EG_goodVtx_noscrapping.root"); 
	//inputEventTree->Add("root://ccxroot.in2p3.fr:1094//pnfs/in2p3.fr/data/cms/t2data/store/user/obondu/GluGlu_120_PU_S4_v2/GluGlu_120_5_1_e5Q.root");
	//inputEventTree->Add("root://ccxroot.in2p3.fr:1094//pnfs/in2p3.fr/data/cms/t2data/store/user/obondu/GluGlu_120_PU_S4_v1/GluGlu_120_3_1_Gwj.root");
	//inputEventTree->Add("/sps/cms/hbrun/dataset42X/theHiggsFile.root");
	inputEventTree->Add("root://ccxroot.in2p3.fr:1094//pnfs/in2p3.fr/data/cms/t2data/store/user/obondu/GluGlu_120_v2//GluGlu_120_3_1_8BK.root");
	inputEventTree->Add("root://ccxroot.in2p3.fr:1094//pnfs/in2p3.fr/data/cms/t2data/store/user/obondu/GluGlu_120_v2/GluGlu_120_2_1_aHq.root");
	inputEventTree->Add("root://ccxroot.in2p3.fr:1094//pnfs/in2p3.fr/data/cms/t2data/store/user/obondu/GluGlu_120_v2/GluGlu_120_1_1_K6q.root");
	inputEventTree->Add("root://ccxroot.in2p3.fr:1094//pnfs/in2p3.fr/data/cms/t2data/store/user/obondu/GluGlu_120_v2/GluGlu_120_4_1_jnV.root");
	inputEventTree->Add("root://ccxroot.in2p3.fr:1094//pnfs/in2p3.fr/data/cms/t2data/store/user/obondu/GluGlu_120_v2/GluGlu_120_5_1_lfA.root");
//	inputEventTree->Add("../test/MC_EG_goodVtx_noscrapping_good.root");


	myFile=new TFile("theMiniTree.root","RECREATE");

	beginMacro();

	int NbEvents = inputEventTree->GetEntries();	cout << "NbEvents = " << NbEvents << endl;
	//NbEvents = 10000;
	int NbHLT20 = 0;
	int NbPhotons = 0;
	int NbPhotonsCleaned = 0;
	int NbPhotonsAccept = 0;
	int NbPhotonsCutEt = 0;
	int NbPhotonsPb = 0;

	int phoPair = 0;
	int phoPairLeadCiC0 = 0;	
	int phoPairLeadCiC1 = 0;	
	int phoPairLeadCiC2 = 0;	
	int phoPairLeadCiC3 = 0;	
	int phoPairLeadCiC4 = 0;	
	int phoPairLeadCiC5 = 0;	
	int phoPairLeadCiC6 = 0;	
	int phoPairLeadCiC7 = 0;	
	int phoPairTrailCiC0 = 0;	
	int phoPairTrailCiC1 = 0;	
	int phoPairTrailCiC2 = 0;	
	int phoPairTrailCiC3 = 0;	
	int phoPairTrailCiC4 = 0;	
	int phoPairTrailCiC5 = 0;	
	int phoPairTrailCiC6 = 0;	
	int phoPairTrailCiC7 = 0;
	int phoAfterCutEt = 0;
	int phoGoodPair = 0;	

	cout << "nomFichier=" << inputEventTree->GetFile()->GetName() << " nbEntries=" << NbEvents << endl;
	for (int ievt  = 0 ; ievt < NbEvents ; ievt++){
	      inputEventTree->GetEvent(ievt);
	      	if( (ievt%10==0 && ievt<=100)  || (ievt%100==0 && ievt<=1000)   || (ievt%1000==0 && ievt>1000)  )
		{
		  cout <<"Analyzing "<< ievt << "th event: " << endl;
		}
		if (doHLT){
			
			if (nbHlt > 0) {if (event->hltAccept(ListWantedHLTnames[0])) dipho_HLT_bit0 = 1; else dipho_HLT_bit0 = 0;}
			if (nbHlt > 1) {if (event->hltAccept(ListWantedHLTnames[1])) dipho_HLT_bit1 = 1; else dipho_HLT_bit1 = 0;}
			if (nbHlt > 2) {if (event->hltAccept(ListWantedHLTnames[2])) dipho_HLT_bit2 = 1; else dipho_HLT_bit2 = 0;}
			if (nbHlt > 3) {if (event->hltAccept(ListWantedHLTnames[3])) dipho_HLT_bit3 = 1; else dipho_HLT_bit3 = 0;}
			if (nbHlt > 4) {if (event->hltAccept(ListWantedHLTnames[4])) dipho_HLT_bit4 = 1; else dipho_HLT_bit4 = 0;}
			if (nbHlt > 5) {if (event->hltAccept(ListWantedHLTnames[5])) dipho_HLT_bit5 = 1; else dipho_HLT_bit5 = 0;}
			if (nbHlt > 6) {if (event->hltAccept(ListWantedHLTnames[6])) dipho_HLT_bit6 = 1; else dipho_HLT_bit6 = 0;}
			if (nbHlt > 7) {if (event->hltAccept(ListWantedHLTnames[7])) dipho_HLT_bit7 = 1; else dipho_HLT_bit7 = 0;}
			if (nbHlt > 8) {if (event->hltAccept(ListWantedHLTnames[8])) dipho_HLT_bit8 = 1; else dipho_HLT_bit8 = 0;}
			if (nbHlt > 9) {if (event->hltAccept(ListWantedHLTnames[9])) dipho_HLT_bit9 = 1; else dipho_HLT_bit9 = 0;}
			if (nbHlt > 10) {if (event->hltAccept(ListWantedHLTnames[10])) dipho_HLT_bit10 = 1; else dipho_HLT_bit10 = 0;}
			if (nbHlt > 11) {if (event->hltAccept(ListWantedHLTnames[11])) dipho_HLT_bit11 = 1; else dipho_HLT_bit11 = 0;}
			if (nbHlt > 12) {if (event->hltAccept(ListWantedHLTnames[12])) dipho_HLT_bit12 = 1; else dipho_HLT_bit12 = 0;}
			if (nbHlt > 13) {if (event->hltAccept(ListWantedHLTnames[13])) dipho_HLT_bit13 = 1; else dipho_HLT_bit13 = 0;}
			if (nbHlt > 14) {if (event->hltAccept(ListWantedHLTnames[14])) dipho_HLT_bit14 = 1; else dipho_HLT_bit14 = 0;}
			if (nbHlt > 15) {if (event->hltAccept(ListWantedHLTnames[15])) dipho_HLT_bit15 = 1; else dipho_HLT_bit15 = 0;}
			if (nbHlt > 16) {if (event->hltAccept(ListWantedHLTnames[16])) dipho_HLT_bit16 = 1; else dipho_HLT_bit16 = 0;}
			if (nbHlt > 17) {if (event->hltAccept(ListWantedHLTnames[17])) dipho_HLT_bit17 = 1; else dipho_HLT_bit17 = 0;}
			if (nbHlt > 18) {if (event->hltAccept(ListWantedHLTnames[18])) dipho_HLT_bit18 = 1; else dipho_HLT_bit18 = 0;}
			if (nbHlt > 19) {if (event->hltAccept(ListWantedHLTnames[19])) dipho_HLT_bit19 = 1; else dipho_HLT_bit19 = 0;}
			if (nbHlt > 20) {if (event->hltAccept(ListWantedHLTnames[20])) dipho_HLT_bit20 = 1; else dipho_HLT_bit20 = 0;}
			if (nbHlt > 21) {if (event->hltAccept(ListWantedHLTnames[21])) dipho_HLT_bit21 = 1; else dipho_HLT_bit21 = 0;}
			if (nbHlt > 22) {if (event->hltAccept(ListWantedHLTnames[22])) dipho_HLT_bit22 = 1; else dipho_HLT_bit22 = 0;}
			if (nbHlt > 23) {if (event->hltAccept(ListWantedHLTnames[23])) dipho_HLT_bit23 = 1; else dipho_HLT_bit23 = 0;}
			if (nbHlt > 24) {if (event->hltAccept(ListWantedHLTnames[24])) dipho_HLT_bit24 = 1; else dipho_HLT_bit24 = 0;}
			if (nbHlt > 25) {if (event->hltAccept(ListWantedHLTnames[25])) dipho_HLT_bit25 = 1; else dipho_HLT_bit25 = 0;}
			if (nbHlt > 26) {if (event->hltAccept(ListWantedHLTnames[26])) dipho_HLT_bit26 = 1; else dipho_HLT_bit26 = 0;}
			if (nbHlt > 27) {if (event->hltAccept(ListWantedHLTnames[27])) dipho_HLT_bit27 = 1; else dipho_HLT_bit27 = 0;}
			if (nbHlt > 28) {if (event->hltAccept(ListWantedHLTnames[28])) dipho_HLT_bit28 = 1; else dipho_HLT_bit28 = 0;}
			if (nbHlt > 29) {if (event->hltAccept(ListWantedHLTnames[29])) dipho_HLT_bit29 = 1; else dipho_HLT_bit29 = 0;}
			if (nbHlt > 30) {if (event->hltAccept(ListWantedHLTnames[30])) dipho_HLT_bit30 = 1; else dipho_HLT_bit30 = 0;}
			if (nbHlt > 31) {if (event->hltAccept(ListWantedHLTnames[31])) dipho_HLT_bit31 = 1; else dipho_HLT_bit31 = 0;}
			if (nbHlt > 32) {if (event->hltAccept(ListWantedHLTnames[32])) dipho_HLT_bit32 = 1; else dipho_HLT_bit32 = 0;}
			if (nbHlt > 33) {if (event->hltAccept(ListWantedHLTnames[33])) dipho_HLT_bit33 = 1; else dipho_HLT_bit33 = 0;}
			if (nbHlt > 34) {if (event->hltAccept(ListWantedHLTnames[34])) dipho_HLT_bit34 = 1; else dipho_HLT_bit34 = 0;}
			if (nbHlt > 35) {if (event->hltAccept(ListWantedHLTnames[35])) dipho_HLT_bit35 = 1; else dipho_HLT_bit35 = 0;}
			if (nbHlt > 36) {if (event->hltAccept(ListWantedHLTnames[36])) dipho_HLT_bit36 = 1; else dipho_HLT_bit36 = 0;}
			if (nbHlt > 37) {if (event->hltAccept(ListWantedHLTnames[37])) dipho_HLT_bit37 = 1; else dipho_HLT_bit37 = 0;}
			if (nbHlt > 38) {if (event->hltAccept(ListWantedHLTnames[38])) dipho_HLT_bit38 = 1; else dipho_HLT_bit38 = 0;}
			if (nbHlt > 39) {if (event->hltAccept(ListWantedHLTnames[39])) dipho_HLT_bit39 = 1; else dipho_HLT_bit39 = 0;}
			if (nbHlt > 40) {if (event->hltAccept(ListWantedHLTnames[40])) dipho_HLT_bit40 = 1; else dipho_HLT_bit40 = 0;}
			if (nbHlt > 41) {if (event->hltAccept(ListWantedHLTnames[41])) dipho_HLT_bit41 = 1; else dipho_HLT_bit41 = 0;}
			if (nbHlt > 42) {if (event->hltAccept(ListWantedHLTnames[42])) dipho_HLT_bit42 = 1; else dipho_HLT_bit42 = 0;}
			if (nbHlt > 43) {if (event->hltAccept(ListWantedHLTnames[43])) dipho_HLT_bit43 = 1; else dipho_HLT_bit43 = 0;}
			if (nbHlt > 44) {if (event->hltAccept(ListWantedHLTnames[44])) dipho_HLT_bit44 = 1; else dipho_HLT_bit44 = 0;}
			if (nbHlt > 45) {if (event->hltAccept(ListWantedHLTnames[45])) dipho_HLT_bit45 = 1; else dipho_HLT_bit45 = 0;}
		}

		vector<pair <TRootPhoton*, TRootPhoton*> > theDiphotonPairs; //lead, trail

		int NbPhoton = photons->GetEntriesFast();
//		cout << "nb of photon " << NbPhoton << endl;
		for (int iphoton=0; iphoton< NbPhoton ; iphoton++){
			NbPhotons++;
			TRootPhoton *myphoton = (TRootPhoton*) photons->At(iphoton);
			if (!(photonPassingPreselection(myphoton))) continue;
			for (int iphoton2=(iphoton+1) ; iphoton2 < NbPhoton ; iphoton2++){
				TRootPhoton *myphoton2 = (TRootPhoton*) photons->At(iphoton2);
				if (!(photonPassingPreselection(myphoton2))) continue;
				pair <TRootPhoton*, TRootPhoton*> theDiphotonPair;
				if (myphoton->Et()>myphoton2->Et()) theDiphotonPair = make_pair(myphoton, myphoton2);
				else theDiphotonPair = make_pair(myphoton2, myphoton);
				theDiphotonPairs.push_back(theDiphotonPair);

			}
		}
		// now run over diphoton to see if pass CiC	
		float theMaxSumEt = 0; int  theMaxSumEtInd = -10; int nbOfGood = 0;
		for (int i = 0 ; i < theDiphotonPairs.size() ; i++){
			int leadCiCStop, trailCiCStop;
			phoPair++;
		TLorentzVector Plead, Ptrail, Psum;
		Plead.SetPxPyPzE((theDiphotonPairs[i].first)->Px(),(theDiphotonPairs[i].first)->Py(),(theDiphotonPairs[i].first)->Pz(),(theDiphotonPairs[i].first)->Energy());
		Ptrail.SetPxPyPzE((theDiphotonPairs[i].second)->Px(),(theDiphotonPairs[i].second)->Py(),(theDiphotonPairs[i].second)->Pz(),(theDiphotonPairs[i].second)->Energy());
		Psum = Plead + Ptrail;
			int theDiphoCat = findTheDiphoCat(*(theDiphotonPairs[i].first), *(theDiphotonPairs[i].second));
			int phoLeadCat = findThePhoCat(*(theDiphotonPairs[i].first)); 
			cout << "///////////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
			cout << "Pair event Number=" << event->eventId() << "   runNumber=" << event->runId() << "   LS=" << event->luminosityBlock() <<"    cat of the pair=" << theDiphoCat << "   selectedVertex=0 "<< endl;
			bool leadPassingCiC = photonIsPassingCIC(*(theDiphotonPairs[i].first), vertices, tracks, *beamSpot, electrons, &leadCiCStop, phoLeadCat);
		//	cout << leadCiCStop << endl;
			switch (leadCiCStop){
				case 7 : phoPairLeadCiC7++; 
				case 6 : phoPairLeadCiC6++; 
				case 5 : phoPairLeadCiC5++;
				case 4 : phoPairLeadCiC4++;
				case 3 : phoPairLeadCiC3++;
				case 2 : phoPairLeadCiC2++;
				case 1 : phoPairLeadCiC1++;
				case 0 : phoPairLeadCiC0++;
			}
			int phoTrailCat = findThePhoCat(*(theDiphotonPairs[i].second)); 
			bool trailPassingCiC = photonIsPassingCIC(*(theDiphotonPairs[i].second), vertices, tracks, *beamSpot, electrons, &trailCiCStop, phoTrailCat);	
			switch (trailCiCStop){
				case 7 : phoPairTrailCiC7++;
				case 6 : phoPairTrailCiC6++;
				case 5 : phoPairTrailCiC5++;
				case 4 : phoPairTrailCiC4++;
				case 3 : phoPairTrailCiC3++;
				case 2 : phoPairTrailCiC2++;
				case 1 : phoPairTrailCiC1++;
				case 0 : phoPairTrailCiC0++; 
			}
			if (!(leadPassingCiC&&((theDiphotonPairs[i].first)->Et()>40))) continue;//test if lead pass CiC
			if (!(trailPassingCiC&&((theDiphotonPairs[i].second)->Et()>30))) continue; // tets if trail pass CiC
			phoAfterCutEt++;
	//		cout << "le second passe " << endl;
			nbOfGood++;
			float sumEt = (theDiphotonPairs[i].first)->Et() + (theDiphotonPairs[i].second)->Et();
			if (sumEt > theMaxSumEt) {
			//	cout << "Et lead " << (theDiphotonPairs[i].first)->Et() << " photon trail " << (theDiphotonPairs[i].second)->Et() << endl;
				theMaxSumEt = sumEt;
				theMaxSumEtInd = i;
			}
			//cout << "good pair in the event  " << nbOfGood << endl;
		}
		if ( theMaxSumEtInd >= 0 ) {
			phoGoodPair++;
			saveThisEvent(event, theDiphotonPairs[theMaxSumEtInd], photons, vertices);
		}
	}
	cout << "phopair " << phoPair << endl;
	cout << "phoPairLeadCiC0 = " << phoPairLeadCiC0 << endl;
	cout << "phoPairLeadCiC1 = " << phoPairLeadCiC1 << endl;
	cout << "phoPairLeadCiC2 = " << phoPairLeadCiC2 << endl;
	cout << "phoPairLeadCiC3 = " << phoPairLeadCiC3 << endl;
	cout << "phoPairLeadCiC4 = " << phoPairLeadCiC4 << endl;
	cout << "phoPairLeadCiC5 = " << phoPairLeadCiC5 << endl;
	cout << "phoPairLeadCiC6 = " << phoPairLeadCiC6 << endl;
	cout << "phoPairLeadCiC7 = " << phoPairLeadCiC7 << endl;
	cout << "phoPairTrailCiC0 = " << phoPairTrailCiC0 << endl;
	cout << "phoPairTrailCiC1 = " << phoPairTrailCiC1 << endl;
	cout << "phoPairTrailCiC2 = " << phoPairTrailCiC2 << endl;
	cout << "phoPairTrailCiC3 = " << phoPairTrailCiC3 << endl;
	cout << "phoPairTrailCiC4 = " << phoPairTrailCiC4 << endl;
	cout << "phoPairTrailCiC5 = " << phoPairTrailCiC5 << endl;
	cout << "phoPairTrailCiC6 = " << phoPairTrailCiC6 << endl;
	cout << "phoPairTrailCiC7 = " << phoPairTrailCiC7 << endl;
	cout << "phoAfterCutEt = " << phoAfterCutEt << endl;
	cout << "phoGoodPair = " << phoGoodPair << endl;


	endMacro();
}
