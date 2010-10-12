#include "makeDiphotonTree.h"

int beforeMacro(){
	chain->SetBranchAddress("pho_iEvent",&pho_iEvent);         
	chain->SetBranchAddress("pho_EventNumber",&pho_EventNumber);
	chain->SetBranchAddress("pho_RunNumber",&pho_RunNumber);   
	chain->SetBranchAddress("pho_lumiSection",&pho_lumiSection); 
	chain->SetBranchAddress("pho_isEB",&pho_isEB);              
	chain->SetBranchAddress("pho_isEE",&pho_isEE);              
	chain->SetBranchAddress("pho_isEEP",&pho_isEEP);            
	chain->SetBranchAddress("pho_isEEM",&pho_isEEM);            
	chain->SetBranchAddress("pho_isAlsoElectron",&pho_isAlsoElectron);
	chain->SetBranchAddress("pho_ElectronClassification",&pho_ElectronClassification);
	chain->SetBranchAddress("pho_hasPixelSeed",&pho_hasPixelSeed);                    
	chain->SetBranchAddress("pho_et",&pho_et);                                        
	chain->SetBranchAddress("pho_energy",&pho_energy);                                
	chain->SetBranchAddress("pho_px",&pho_px);                                        
	chain->SetBranchAddress("pho_py",&pho_py);                                        
	chain->SetBranchAddress("pho_pz",&pho_pz);                                        
	chain->SetBranchAddress("pho_eta",&pho_eta);                                      
	chain->SetBranchAddress("pho_phi",&pho_phi);                                      
	chain->SetBranchAddress("pho_s4",&pho_s4);                                        
	chain->SetBranchAddress("pho_s9",&pho_s9);                                        
	chain->SetBranchAddress("pho_eMax",&pho_eMax);                                    
	chain->SetBranchAddress("pho_e2nd",&pho_e2nd);                                    
	chain->SetBranchAddress("pho_r19",&pho_r19);                                      
	chain->SetBranchAddress("pho_r9",&pho_r9);                                        
	chain->SetBranchAddress("pho_S9overS9minusS1S2",&pho_S9overS9minusS1S2);          
	chain->SetBranchAddress("pho_cEE",&pho_cEE);                                      
	chain->SetBranchAddress("pho_cEP",&pho_cEP);                                      
	chain->SetBranchAddress("pho_cPP",&pho_cPP);                                      
	chain->SetBranchAddress("pho_phiwidth",&pho_phiwidth);                            
	chain->SetBranchAddress("pho_etawidth",&pho_etawidth);                            
	chain->SetBranchAddress("pho_sigmaphi",&pho_sigmaphi);                            
	chain->SetBranchAddress("pho_sigmaeta",&pho_sigmaeta);                            
	chain->SetBranchAddress("pho_sigmaIetaIeta",&pho_sigmaIetaIeta);                  
	chain->SetBranchAddress("pho_sigmaEtaEta",&pho_sigmaEtaEta);                      
	chain->SetBranchAddress("pho_hoe",&pho_hoe);                                      
	chain->SetBranchAddress("pho_IsoEcalRechit",&pho_IsoEcalRechit);                  
	chain->SetBranchAddress("pho_IsoHcalRechit",&pho_IsoHcalRechit);                  
	chain->SetBranchAddress("pho_IsoSolidTrkCone",&pho_IsoSolidTrkCone);              
	chain->SetBranchAddress("pho_IsoHollowTrkCone",&pho_IsoHollowTrkCone);            
	chain->SetBranchAddress("pho_IsoEcalRechit03",&pho_IsoEcalRechit03);              
	chain->SetBranchAddress("pho_IsoHcalRechit03",&pho_IsoHcalRechit03);              
	chain->SetBranchAddress("pho_IsoSolidTrkCone03",&pho_IsoSolidTrkCone03);          
	chain->SetBranchAddress("pho_IsoHollowTrkCone03",&pho_IsoHollowTrkCone03);        
	chain->SetBranchAddress("pho_esRatio",&pho_esRatio);                              
	chain->SetBranchAddress("pho_convNTracks",&pho_convNTracks);                      
	chain->SetBranchAddress("pho_ptoverjetpt",&pho_ptoverjetpt);                      
	chain->SetBranchAddress("pho_DrSCclosest",&pho_DrSCclosest);                      
	chain->SetBranchAddress("pho_DrTrkclosest",&pho_DrTrkclosest);                    
	chain->SetBranchAddress("pho_ptoverjetpt_pt2",&pho_ptoverjetpt_pt2);              
	chain->SetBranchAddress("pho_DrTrkclosest_pt2",&pho_DrTrkclosest_pt2);            
	chain->SetBranchAddress("pho_ptoverjetpt_pt5",&pho_ptoverjetpt_pt5);              
	chain->SetBranchAddress("pho_DrTrkclosest_pt5",&pho_DrTrkclosest_pt5);            
	chain->SetBranchAddress("pho_ptoverjetpt_pt10",&pho_ptoverjetpt_pt10);            
	chain->SetBranchAddress("pho_DrTrkclosest_pt10",&pho_DrTrkclosest_pt10);          
	chain->SetBranchAddress("pho_ptoverjetpt_pt15",&pho_ptoverjetpt_pt15);            
	chain->SetBranchAddress("pho_ptoverjetpt_PFlow",&pho_ptoverjetpt_PFlow);          
	chain->SetBranchAddress("pho_ptoverjetpt_PFlow_pt2",&pho_ptoverjetpt_PFlow_pt2);  
	chain->SetBranchAddress("pho_ptoverjetpt_PFlow_pt5",&pho_ptoverjetpt_PFlow_pt5);  
	chain->SetBranchAddress("pho_ptoverjetpt_PFlow_pt10",&pho_ptoverjetpt_PFlow_pt10);
	chain->SetBranchAddress("pho_ptoverjetpt_sisCone",&pho_ptoverjetpt_sisCone);      
	chain->SetBranchAddress("pho_ptoverjetpt_sisCone_pt2",&pho_ptoverjetpt_sisCone_pt2);
	chain->SetBranchAddress("pho_ptoverjetpt_sisCone_pt5",&pho_ptoverjetpt_sisCone_pt5);
	chain->SetBranchAddress("pho_ptoverjetpt_sisCone_pt10",&pho_ptoverjetpt_sisCone_pt10);
	chain->SetBranchAddress("pho_isAfterCut1",&pho_isAfterCut1);                          
	chain->SetBranchAddress("pho_isAfterCut2",&pho_isAfterCut2);                          
	chain->SetBranchAddress("pho_isAfterCut3",&pho_isAfterCut3);                          
	chain->SetBranchAddress("pho_isAfterCut4",&pho_isAfterCut4);                          
	chain->SetBranchAddress("pho_isAfterCut5",&pho_isAfterCut5);                          
	chain->SetBranchAddress("pho_isAfterCut6",&pho_isAfterCut6);                          
	chain->SetBranchAddress("pho_isAfterCut7",&pho_isAfterCut7);                          
	chain->SetBranchAddress("pho_isAfterCut8",&pho_isAfterCut8);                          
	chain->SetBranchAddress("pho_jetEMfraction",&pho_jetEMfraction);                      
	chain->SetBranchAddress("pho_DrJetClosest",&pho_DrJetClosest);                        
	chain->SetBranchAddress("pho_isLoose",&pho_isLoose);                                  
	chain->SetBranchAddress("pho_isTight",&pho_isTight);                                  
	chain->SetBranchAddress("pho_isEBGap",&pho_isEBGap);                                  
	chain->SetBranchAddress("pho_isEEGap",&pho_isEEGap);                                  
	chain->SetBranchAddress("pho_GenId",&pho_GenId);                                      
	chain->SetBranchAddress("pho_MotherId",&pho_MotherId);                                
	chain->SetBranchAddress("pho_isPromptGenPho",&pho_isPromptGenPho);                    
	chain->SetBranchAddress("pho_isFromQuarkGen",&pho_isFromQuarkGen);
	chain->SetBranchAddress("pho_isPi0Gen",&pho_isPi0Gen);
	chain->SetBranchAddress("pho_isEtaGen",&pho_isEtaGen);
	chain->SetBranchAddress("pho_isRhoGen",&pho_isRhoGen);
	chain->SetBranchAddress("pho_isOmegaGen",&pho_isOmegaGen);
	chain->SetBranchAddress("pho_isGenElectron",&pho_isGenElectron);
	chain->SetBranchAddress("pho_eventPassHLT_Photon10_L1R",&pho_eventPassHLT_Photon10_L1R);
	chain->SetBranchAddress("pho_eventPassHLT_Photon15_L1R",&pho_eventPassHLT_Photon15_L1R);
	chain->SetBranchAddress("pho_eventPassHLT_DoublePhoton10_L1R",&pho_eventPassHLT_DoublePhoton10_L1R);
	chain->SetBranchAddress("pho_eventPtHat",&pho_eventPtHat);
	chain->SetBranchAddress("pho_nVertex",&pho_nVertex);
	chain->SetBranchAddress("pho_PromptGenIsoEnergyStatus1_cone02",&pho_PromptGenIsoEnergyStatus1_cone02);
	chain->SetBranchAddress("pho_PromptGenIsoEnergyStatus2_cone02",&pho_PromptGenIsoEnergyStatus2_cone02);
	chain->SetBranchAddress("pho_PromptGenIsoEnergyStatus1_cone03",&pho_PromptGenIsoEnergyStatus1_cone03);
	chain->SetBranchAddress("pho_PromptGenIsoEnergyStatus2_cone03",&pho_PromptGenIsoEnergyStatus2_cone03);
	chain->SetBranchAddress("pho_PromptGenIsoEnergyStatus1_cone035",&pho_PromptGenIsoEnergyStatus1_cone035);
	chain->SetBranchAddress("pho_PromptGenIsoEnergyStatus2_cone035",&pho_PromptGenIsoEnergyStatus2_cone035);
	chain->SetBranchAddress("pho_PromptGenIsoEnergyStatus1_cone04",&pho_PromptGenIsoEnergyStatus1_cone04);
	chain->SetBranchAddress("pho_PromptGenIsoEnergyStatus2_cone04",&pho_PromptGenIsoEnergyStatus2_cone04);
	chain->SetBranchAddress("pho_seedSeverity",&pho_seedSeverity);
	chain->SetBranchAddress("pho_recoFlag",&pho_recoFlag);
	chain->SetBranchAddress("pho_HLT30",&pho_HLT30);
	chain->SetBranchAddress("pho_HLT50",&pho_HLT50);
	chain->SetBranchAddress("pho_HLTdiPho15",&pho_HLTdiPho15);
	chain->SetBranchAddress("pho_HLTdiPho20",&pho_HLTdiPho20);
	chain->SetBranchAddress("pho_SCeta",&pho_SCeta);
	chain->SetBranchAddress("pho_SCphi",&pho_SCphi);
	chain->SetBranchAddress("pho_SCEtraw",&pho_SCEtraw);
	chain->SetBranchAddress("pho_SCEt",&pho_SCEt);
	chain->SetBranchAddress("pho_SCr9",&pho_SCr9);
	chain->SetBranchAddress("pho_SCbr",&pho_SCbr);
	chain->SetBranchAddress("pho_SCnbBC",&pho_SCnbBC);
	chain->SetBranchAddress("pho_SCnXtal",&pho_SCnXtal);
	
		myTree_ = new TTree("diPhotons","DiPhotonsInfos");
                myTree_->Branch("event_number",&event_number,"event_number/I");
		myTree_->Branch("event_realNumber",&event_realNumber,"event_realNumber/I");
		myTree_->Branch("event_runNumber",&event_runNumber, "event_runNumber/I");
		myTree_->Branch("event_LumiSection",&event_LumiSection,"event_LumiSection/I");
		myTree_->Branch("event_eventPtHat",&event_eventPtHat, "event_eventPtHat/F");
		myTree_->Branch("event_nVertex",&event_nVertex,"event_nVertex/I");
                myTree_->Branch("event_weightLO",&event_weightLO,"event_weightLO/I");
                myTree_->Branch("event_weightNLO",&event_weightNLO,"event_weightNLO/I");
                myTree_->Branch("dipho_cat4",&dipho_cat4,"dipho_cat4/I");               
                myTree_->Branch("dipho_cat4bis",&dipho_cat4bis,"dipho_cat4bis/I");      
                myTree_->Branch("dipho_cat6",&dipho_cat6,"dipho_cat6/I");               
                myTree_->Branch("dipho_cat12",&dipho_cat12,"dipho_cat12/I");            
                myTree_->Branch("dipho_mgg",&dipho_mgg,"dipho_mgg/F");
                myTree_->Branch("dipho_qt",&dipho_qt,"dipho_qt/F");                                 
                myTree_->Branch("dipho_ql",&dipho_ql,"dipho_ql/F");                                 
                myTree_->Branch("dipho_deltaR",&dipho_deltaR,"dipho_deltaR/F");                     
                myTree_->Branch("dipho_costhetastar",&dipho_costhetastar,"dipho_costhetastar/F");   
                myTree_->Branch("dipho_eta",&dipho_eta,"dipho_eta/F");                              
                myTree_->Branch("dipho_etastar",&dipho_etastar,"dipho_etastar/F");                  
                myTree_->Branch("dipho_rapidity",&dipho_rapidity,"dipho_rapidity/F");               
                myTree_->Branch("dipho_deltaphi",&dipho_deltaphi,"dipho_deltaphi/F");               
                myTree_->Branch("dipho_costhetastar_parton",&dipho_costhetastar_parton,"dipho_costhetastar_parton/F");
                myTree_->Branch("dipho_tanhrapiditydiff",&dipho_tanhrapiditydiff,"dipho_tanhrapiditydiff/F");         
                myTree_->Branch("pholead_pt",&pholead_pt,"pholead_pt/F");                                             
                myTree_->Branch("photrail_pt",&photrail_pt,"photrail_pt/F");                                          
                myTree_->Branch("pho1_pt",&pho1_pt,"pho1_pt/F");                                                      
                myTree_->Branch("pho2_pt",&pho2_pt,"pho2_pt/F");                                                      
                myTree_->Branch("pholead_eta",&pholead_eta,"pholead_eta/F");                                          
                myTree_->Branch("photrail_eta",&photrail_eta,"photrail_eta/F");                                       
                myTree_->Branch("pholead_r9",&pholead_r9,"pholead_r9/F");                                             
                myTree_->Branch("photrail_r9",&photrail_r9,"photrail_r9/F");                                          
                myTree_->Branch("pholead_cPP",&pholead_cPP,"pholead_cPP/F");                                          
                myTree_->Branch("photrail_cPP",&photrail_cPP,"photrail_cPP/F");                                       
                myTree_->Branch("pholead_sigmaphi",&pholead_sigmaphi,"pholead_sigmaphi/F");                           
                myTree_->Branch("photrail_sigmaphi",&photrail_sigmaphi,"photrail_sigmaphi/F");                        
                myTree_->Branch("pholead_S9overS9minusS1S2",&pholead_S9overS9minusS1S2,"pholead_S9overS9minusS1S2/F");
                myTree_->Branch("photrail_S9overS9minusS1S2",&photrail_S9overS9minusS1S2,"photrail_S9overS9minusS1S2/F");
                myTree_->Branch("pholead_etawidth",&pholead_etawidth,"pholead_etawidth/F");                              
                myTree_->Branch("photrail_etawidth",&photrail_etawidth,"photrail_etawidth/F");                           
                myTree_->Branch("pholead_sigieta",&pholead_sigieta,"pholead_sigieta/F");                                 
                myTree_->Branch("photrail_sigieta",&photrail_sigieta,"photrail_sigieta/F");                              
                myTree_->Branch("pholead_ptoverjetpt",&pholead_ptoverjetpt,"pholead_ptoverjetpt/F");                     
                myTree_->Branch("photrail_ptoverjetpt",&photrail_ptoverjetpt,"photrail_ptoverjetpt/F");                  
                myTree_->Branch("pholead_DrScClosest",&pholead_DrScClosest,"pholead_DrScClosest/F");                     
                myTree_->Branch("photrail_DrScClosest",&photrail_DrScClosest,"photrail_DrScClosest/F");                  
                myTree_->Branch("pholead_DrTrkClosest",&pholead_DrTrkClosest,"pholead_DrTrkClosest/F");                  
                myTree_->Branch("photrail_DrTrkClosest",&photrail_DrTrkClosest,"photrail_DrTrkClosest/F");               
                myTree_->Branch("pholead_MLPoutputClusterShape",&pholead_MLPoutputClusterShape,"pholead_MLPoutputClusterShape/F");                                                                                                                        
                myTree_->Branch("photrail_MLPoutputClusterShape",&photrail_MLPoutputClusterShape,"photrail_MLPoutputClusterShape/F");                                                                                                                     
                myTree_->Branch("pholead_MLPoutputEnv",&pholead_MLPoutputEnv,"pholead_MLPoutputEnv/F");                      
                myTree_->Branch("photrail_MLPoutputEnv",&photrail_MLPoutputEnv,"photrail_MLPoutputEnv/F");                   
                myTree_->Branch("pholead_MLPoutputClusterShapeEnv",&pholead_MLPoutputClusterShapeEnv,"pholead_MLPoutputClusterShapeEnv/F");                                                                                                               
                myTree_->Branch("photrail_MLPoutputClusterShapeEnv",&photrail_MLPoutputClusterShapeEnv,"photrail_MLPoutputClusterShapeEnv/F");                                                                                                            
                myTree_->Branch("phoNNmin_MLPoutputClusterShape",&phoNNmin_MLPoutputClusterShape,"phoNNmin_MLPoutputClusterShape/F");                                                                                                                     
                myTree_->Branch("phoNNmax_MLPoutputClusterShape",&phoNNmax_MLPoutputClusterShape,"phoNNmax_MLPoutputClusterShape/F");                                                                                                                     
                myTree_->Branch("phoNNmin_MLPoutputEnv",&phoNNmin_MLPoutputEnv,"phoNNmin_MLPoutputEnv/F");                   
                myTree_->Branch("phoNNmax_MLPoutputEnv",&phoNNmax_MLPoutputEnv,"phoNNmax_MLPoutputEnv/F");                   
                myTree_->Branch("phoNNmin_MLPoutputClusterShapeEnv",&phoNNmin_MLPoutputClusterShapeEnv,"phoNNmin_MLPoutputClusterShapeEnv/F");                                                                                                            
                myTree_->Branch("phoNNmax_MLPoutputClusterShapeEnv",&phoNNmax_MLPoutputClusterShapeEnv,"phoNNmax_MLPoutputClusterShapeEnv/F");                                                                                                            
                myTree_->Branch("pholead_MLPoutputNNisol",&pholead_MLPoutputNNisol,"pholead_MLPoutputNNisol/F");             
                myTree_->Branch("photrail_MLPoutputNNisol",&photrail_MLPoutputNNisol,"photrail_MLPoutputNNisol/F");          
                myTree_->Branch("pho1_MLPoutputNNisol",&pho1_MLPoutputNNisol,"pho1_MLPoutputNNisol/F");                      
                myTree_->Branch("pho2_MLPoutputNNisol",&pho2_MLPoutputNNisol,"pho2_MLPoutputNNisol/F");                      
                myTree_->Branch("phoNNmin_MLPoutputNNisol",&phoNNmin_MLPoutputNNisol,"phoNNmin_MLPoutputNNisol/F");          
                myTree_->Branch("phoNNmax_MLPoutputNNisol",&phoNNmax_MLPoutputNNisol,"phoNNmax_MLPoutputNNisol/F");          
                myTree_->Branch("pholead_HcalIso",&pholead_HcalIso,"pholead_HcalIso/F");                                     
                myTree_->Branch("photrail_HcalIso",&photrail_HcalIso,"photrail_HcalIso/F");                                  
                myTree_->Branch("pholead_EcalIso",&pholead_EcalIso,"pholead_EcalIso/F");
                myTree_->Branch("photrail_EcalIso",&photrail_EcalIso,"photrail_EcalIso/F");
                myTree_->Branch("pholead_TrackerIso",&pholead_TrackerIso,"pholead_TrackerIso/F");
                myTree_->Branch("photrail_TrackerIso",&photrail_TrackerIso,"photrail_TrackerIso/F");
                myTree_->Branch("pholead_HcalIsoPerso",&pholead_HcalIsoPerso,"pholead_HcalIsoPerso/F");
                myTree_->Branch("photrail_HcalIsoPerso",&photrail_HcalIsoPerso,"photrail_HcalIsoPerso/F");
                myTree_->Branch("pholead_EcalIsoPerso",&pholead_EcalIsoPerso,"pholead_EcalIsoPerso/F");
                myTree_->Branch("photrail_EcalIsoPerso",&photrail_EcalIsoPerso,"photrail_EcalIsoPerso/F");
                myTree_->Branch("pholead_TrackerIsoPerso",&pholead_TrackerIsoPerso,"pholead_TrackerIsoPerso/F");
                myTree_->Branch("photrail_TrackerIsoPerso",&photrail_TrackerIsoPerso,"photrail_TrackerIsoPerso/F");
                myTree_->Branch("pholead_hoe",&pholead_hoe,"pholead_hoe/F");
                myTree_->Branch("photrail_hoe",&photrail_hoe,"photrail_hoe/F");
                myTree_->Branch("pholead_DrJetClosest",&pholead_DrJetClosest,"pholead_DrJetClosest/F");
                myTree_->Branch("photrail_DrJetClosest",&photrail_DrJetClosest,"photrail_DrJetClosest/F");
                myTree_->Branch("pholead_jetEMfraction",&pholead_jetEMfraction,"pholead_jetEMfraction/F");
                myTree_->Branch("photrail_jetEMfraction",&photrail_jetEMfraction,"photrail_jetEMfraction/F");
                myTree_->Branch("pholead_isPromptGenPho",&pholead_isPromptGenPho,"pholead_isPromptGenPho/I");
                myTree_->Branch("photrail_isPromptGenPho",&photrail_isPromptGenPho,"photrail_isPromptGenPho/I");
                myTree_->Branch("pholead_isFromQuarkGen",&pholead_isFromQuarkGen,"pholead_isFromQuarkGen/I");
                myTree_->Branch("photrail_isFromQuarkGen",&photrail_isFromQuarkGen,"photrail_isFromQuarkGen/I");
                myTree_->Branch("pholead_HasPixSeed",&pholead_HasPixSeed,"pholead_HasPixSeed/I");
                myTree_->Branch("photrail_HasPixSeed",&photrail_HasPixSeed,"photrail_HasPixSeed/I");
                myTree_->Branch("pholead_seedSeverity",&pholead_seedSeverity,"pholead_seedSeverity/I");
                myTree_->Branch("photrail_seedSeverity",&photrail_seedSeverity,"photrail_seedSeverity/I");
                myTree_->Branch("pholead_recoFlag",&pholead_recoFlag,"pholead_recoFlag/I");
                myTree_->Branch("photrail_recoFlag",&photrail_recoFlag,"photrail_recoFlag/I");
		myTree_->Branch("pholead_isEB",&pholead_isEB,"pholead_isEB/I");
		myTree_->Branch("photrail_isEB",&photrail_isEB,"photrail_isEB/I");
		myTree_->Branch("pholead_isEE",&pholead_isEE,"pholead_isEE/I");
		myTree_->Branch("photrail_isEE",&photrail_isEE,"photrail_isEE/I");
                myTree_->Branch("dipho_eventPassHLT_Photon30_L1R",&dipho_eventPassHLT_Photon30_L1R,"dipho_eventPassHLT_Photon30_L1R/I");
                myTree_->Branch("dipho_eventPassHLT_Photon50_L1R",&dipho_eventPassHLT_Photon50_L1R,"dipho_eventPassHLT_Photon50_L1R/I");
                myTree_->Branch("dipho_eventPassHLT_DoublePhoton10_L1R",&dipho_eventPassHLT_DoublePhoton10_L1R,"dipho_eventPassHLT_DoublePhoton10_L1R/I");
                myTree_->Branch("dipho_eventPassHLT_DoublePhoton15_L1R",&dipho_eventPassHLT_DoublePhoton15_L1R,"dipho_eventPassHLT_DoublePhoton15_L1R/I");
                myTree_->Branch("dipho_Kfactor",&dipho_Kfactor,"dipho_Kfactor/F");


}

int afterMacro(){
	myFile->Write();
	myFile->Close();

}

int fillThisEvent(int iteLead, int iteTrail){
	
	pholead_pt = P_loc[iteLead].Pt();
	pholead_eta = Eta_loc[iteLead];
	pholead_r9 = R9_loc[iteLead];
	pholead_cPP = cPP_loc[iteLead];
	pholead_sigmaphi = sigmaphi_loc[iteLead];
	pholead_etawidth = EtaWidth_loc[iteLead];
	pholead_sigieta = SigIeta_loc[iteLead];
 	pholead_HcalIso = HcalIso_loc[iteLead];
	pholead_EcalIso = EcalIso_loc[iteLead];
	pholead_TrackerIso = TrackerIso_loc[iteLead];
	pholead_hoe = HoE_loc[iteLead];
	pholead_isPromptGenPho =  isPromptGenPho_loc[iteLead];
	pholead_isFromQuarkGen = isFromQuarkGen_loc[iteLead];
	pholead_HasPixSeed = HasPixSeed_loc[iteLead];
	pholead_isEB = pho_isEB_loc[iteLead];
	pholead_isEE = pho_isEE_loc[iteLead];

	photrail_pt = P_loc[iteTrail].Pt();
	photrail_eta = Eta_loc[iteTrail];
	photrail_r9 = R9_loc[iteTrail];
	photrail_cPP = cPP_loc[iteTrail];
	photrail_sigmaphi = sigmaphi_loc[iteTrail];
	photrail_etawidth = EtaWidth_loc[iteTrail];
	photrail_sigieta = pho_sigmaIetaIeta_loc[iteTrail];
 	photrail_HcalIso = HcalIso_loc[iteTrail];
	photrail_EcalIso = EcalIso_loc[iteTrail];
	photrail_TrackerIso = TrackerIso_loc[iteTrail];
	photrail_hoe = HoE_loc[iteTrail];
	photrail_isPromptGenPho =  isPromptGenPho_loc[iteTrail];
	photrail_isFromQuarkGen = isFromQuarkGen_loc[iteTrail];
	photrail_HasPixSeed = HasPixSeed_loc[iteTrail];
	photrail_isEB = pho_isEB_loc[iteTrail];
	photrail_isEE = pho_isEE_loc[iteTrail];

	TLorentzVector Psum = P_loc[iteLead] + P_loc[iteTrail];
	dipho_mgg = Psum.M();

	NbEntries++;

	myTree_->Fill();
}



//makeDiphotonTree(){
	int main(){
	doMC = true; 
	NbEntries = 0;
	myFile = new TFile("diphoton_1.root","RECREATE");
	float MggMin = 20;
	beforeMacro();
	chain->Add("output_1.root");
	int Nevents = chain->GetEntries();
//	Nevents = 50;
	std::cout << "Nevents = " << Nevents << std::endl;
	int ite = 0;
	int numEventBefore, numEventHere, theNumbEvent, theNumbRun;
	int* leadTrail;
	for (int ievt = 0 ; ievt < Nevents ; ievt++){
		if (ievt != 0 ) { // fill the Event level infos
			event_runNumber = pho_RunNumber; 
			event_realNumber = pho_EventNumber;
			event_LumiSection = pho_lumiSection;
			event_nVertex = pho_nVertex;
			dipho_eventPassHLT_Photon30_L1R = pho_HLT30;
			dipho_eventPassHLT_Photon50_L1R = pho_HLT50;
			dipho_eventPassHLT_DoublePhoton10_L1R = pho_HLTdiPho15;
			dipho_eventPassHLT_DoublePhoton15_L1R = pho_HLTdiPho20;
			if (doMC) event_eventPtHat = pho_eventPtHat; 
		}
		chain->GetEntry(ievt);
		if ((ievt % 100000) == 0 ) std::cout << "on est au " << ievt << "evenement " << std::endl;
		if (ievt == 0) {numEventBefore = pho_iEvent;}
		numEventHere = pho_iEvent;
		if (!(numEventHere == numEventBefore)) {
			if (ite >= 2) {
				leadTrail = findLeadAndTrail(P_loc,ite);
//				std::cout << "lead = " << P_loc[leadTrail[0]].Et() << " trail = " << P_loc[leadTrail[1]].Et() << std::endl;
				TLorentzVector Psum = P_loc[leadTrail[0]] + P_loc[leadTrail[1]];
				float Mgg = Psum.M();
				if (Mgg < 0 ) std::cout << "alert Mgg est negatif !! " << std::endl;
				if (Mgg >= MggMin) fillThisEvent(leadTrail[0], leadTrail[1]); 
			}
			ite = 0;

		}
		P_loc[ite].SetPxPyPzE(pho_px,pho_py,pho_pz,pho_energy);
//		std::cout << "P_loc = " << P_loc[ite].Et() << std::endl;
		Eta_loc[ite] = pho_eta;
		R9_loc[ite] = pho_r9;  
		cPP_loc[ite] = pho_cPP;
		EtaWidth_loc[ite] = pho_etawidth;
		ptoverjetpt_loc[ite] = pho_ptoverjetpt;
		dR_SCtrkclosest_loc[ite] = pho_DrTrkclosest;
		closestSC_dR_loc[ite] = pho_DrSCclosest;
		HoE_loc[ite] = pho_hoe;
		TrackerIso_loc[ite] = pho_IsoHollowTrkCone;
		EcalIso_loc[ite] = pho_IsoEcalRechit;
		HcalIso_loc[ite] = pho_IsoHcalRechit;
		HasPixSeed_loc[ite] = pho_hasPixelSeed;
		isPromptGenPho_loc[ite] = pho_isPromptGenPho;
		isFromQuarkGen_loc[ite] = pho_isFromQuarkGen;
		pho_HLT30_loc[ite] = pho_HLT30;
		pho_HLT50_loc[ite] = pho_HLT50;
		pho_HLTdiPho15_loc[ite] = pho_HLTdiPho15;
		pho_HLTdiPho20_loc[ite] = pho_HLTdiPho20;
		pho_sigmaIetaIeta_loc[ite] = pho_sigmaIetaIeta;
		pho_isEB_loc[ite] = pho_isEB;
		pho_isEE_loc[ite] = pho_isEE;
		ite++;

		numEventBefore = pho_iEvent;


		
	}
	std::cout << "NbEventBefore=" << Nevents << " NbEventAfter=" << NbEntries << std::endl;

	afterMacro();
}
