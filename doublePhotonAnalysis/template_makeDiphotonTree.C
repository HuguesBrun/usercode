#include "makeDiphotonTree.h"
#include "TMath.h"

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
	chain->SetBranchAddress("pho_transverseToJetRatio_pt20",&pho_transverseToJetRatio); 
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
	chain->SetBranchAddress("pho_eventPtHat",&pho_eventPtHat);
	if (doMC) {chain->SetBranchAddress("pho_eventProcessId",&pho_eventProcessId);
	chain->SetBranchAddress("pho_trueE",&pho_trueE);
	chain->SetBranchAddress("pho_truePx",&pho_truePx);
	chain->SetBranchAddress("pho_truePy",&pho_truePy);
	chain->SetBranchAddress("pho_truePz",&pho_truePz);
	chain->SetBranchAddress("pho_truePhi",&pho_truePhi);
	chain->SetBranchAddress("pho_trueEta",&pho_trueEta);
	}
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
        chain->SetBranchAddress("pho_HLT_bit0",&pho_HLT_bit0);
        chain->SetBranchAddress("pho_HLT_bit1",&pho_HLT_bit1);
        chain->SetBranchAddress("pho_HLT_bit2",&pho_HLT_bit2);
        chain->SetBranchAddress("pho_HLT_bit3",&pho_HLT_bit3);
        chain->SetBranchAddress("pho_HLT_bit4",&pho_HLT_bit4);
	if (doAllHlt){
        	chain->SetBranchAddress("pho_HLT_bit5",&pho_HLT_bit5);
        	chain->SetBranchAddress("pho_HLT_bit6",&pho_HLT_bit6);
        	chain->SetBranchAddress("pho_HLT_bit7",&pho_HLT_bit7);
        	chain->SetBranchAddress("pho_HLT_bit8",&pho_HLT_bit8);
        	chain->SetBranchAddress("pho_HLT_bit9",&pho_HLT_bit9);
        	chain->SetBranchAddress("pho_HLT_bit10",&pho_HLT_bit10);
        	chain->SetBranchAddress("pho_HLT_bit11",&pho_HLT_bit11);
	}
	chain->SetBranchAddress("pho_SCeta",&pho_SCeta);
	chain->SetBranchAddress("pho_SCphi",&pho_SCphi);
	chain->SetBranchAddress("pho_SCEtraw",&pho_SCEtraw);
	chain->SetBranchAddress("pho_SCEt",&pho_SCEt);
	chain->SetBranchAddress("pho_SCr9",&pho_SCr9);
	chain->SetBranchAddress("pho_SCbr",&pho_SCbr);
	chain->SetBranchAddress("pho_SCnbBC",&pho_SCnbBC);
	chain->SetBranchAddress("pho_SCnXtal",&pho_SCnXtal);
	chain->SetBranchAddress("isAspike",&isAspike);
	if (doNN) {
		chain->SetBranchAddress("pho_NNshapeOutput",&pho_NNshapeOutput);
		chain->SetBranchAddress("pho_NNenvOutput",&pho_NNenvOutput);
		chain->SetBranchAddress("pho_NNcombOutput",&pho_NNcombOutput);
	}
        if (doHLTobject) {
		chain->SetBranchAddress("pho_isMatchingWithHLTObject",&pho_isMatchingWithHLTObject);
	}
	
		myTree_ = new TTree("diPhotons","DiPhotonsInfos");
                myTree_->Branch("event_number",&event_number,"event_number/I");
		myTree_->Branch("event_realNumber",&event_realNumber,"event_realNumber/I");
		myTree_->Branch("event_runNumber",&event_runNumber, "event_runNumber/I");
		myTree_->Branch("event_LumiSection",&event_LumiSection,"event_LumiSection/I");
		myTree_->Branch("event_eventPtHat",&event_eventPtHat, "event_eventPtHat/F");
		myTree_->Branch("event_processId",&event_processId, "event_processId/I");
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
				myTree_->Branch("dipho_mggTrue",&dipho_mggTrue,"dipho_mggTrue/F");
				myTree_->Branch("dipho_qtTrue",&dipho_qtTrue,"dipho_qtTrue/F");
				myTree_->Branch("dipho_deltaphiTrue",&dipho_deltaphiTrue,"dipho_deltaphiTrue/F");
	myTree_->Branch("dipho_costhetastarTrue",&dipho_costhetastarTrue,"dipho_costhetastarTrue/F");
                myTree_->Branch("dipho_costhetastar_parton",&dipho_costhetastar_parton,"dipho_costhetastar_parton/F");
                myTree_->Branch("dipho_tanhrapiditydiff",&dipho_tanhrapiditydiff,"dipho_tanhrapiditydiff/F");        
		myTree_->Branch("dipho_minNNshape",&dipho_minNNshape,"dipho_minNNshape/F");
		myTree_->Branch("dipho_minNNcomb",&dipho_minNNcomb,"dipho_minNNcomb/F");
		myTree_->Branch("dipho_minNNenv",&dipho_minNNenv,"dipho_minNNenv/F");
		myTree_->Branch("dipho_minIsoTrack",&dipho_minIsoTrack,"dipho_minIsoTrack/F"); 
		myTree_->Branch("dipho_maxIsoTrack",&dipho_maxIsoTrack,"dipho_maxIsoTrack/F"); 
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
                myTree_->Branch("pholead_DrSCclosest",&pholead_DrScClosest,"pholead_DrScClosest/F");                     
                myTree_->Branch("photrail_DrSCclosest",&photrail_DrScClosest,"photrail_DrScClosest/F");                  
                myTree_->Branch("pholead_DrTrkClosest",&pholead_DrTrkClosest,"pholead_DrTrkClosest/F");                  
                myTree_->Branch("photrail_DrTrkClosest",&photrail_DrTrkClosest,"photrail_DrTrkClosest/F");              
		myTree_->Branch("pholead_transverseToJetRatio",&pholead_transverseToJetRatio,"pholead_transverseToJetRatio/F");
		myTree_->Branch("photrail_transverseToJetRatio",&photrail_transverseToJetRatio,"photrail_transverseToJetRatio/F");
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
                myTree_->Branch("pholead_HcalIsodR03",&pholead_HcalIsodR03,"pholead_HcalIsodR03/F");                                     
                myTree_->Branch("photrail_HcalIsodR03",&photrail_HcalIsodR03,"photrail_HcalIsodR03/F");                                  
                myTree_->Branch("pholead_EcalIsodR03",&pholead_EcalIsodR03,"pholead_EcalIsodR03/F");
                myTree_->Branch("photrail_EcalIsodR03",&photrail_EcalIsodR03,"photrail_EcalIsodR03/F");
                myTree_->Branch("pholead_TrackerIsodR03",&pholead_TrackerIsodR03,"pholead_TrackerIsodR03/F");
                myTree_->Branch("photrail_TrackerIsodR03",&photrail_TrackerIsodR03,"photrail_TrackerIsodR03/F");
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
				myTree_->Branch("pholead_NNshapeOutput",&pholead_NNshapeOutput, "pholead_NNshapeOutput/F");
				myTree_->Branch("photrail_NNshapeOutput",&photrail_NNshapeOutput, "photrail_NNshapeOutput/F");
				myTree_->Branch("pholead_NNenvOutput",&pholead_NNenvOutput, "pholead_NNenvOutput/F");
				myTree_->Branch("photrail_NNenvOutput",&photrail_NNenvOutput, "photrail_NNenvOutput/F");
				myTree_->Branch("pholead_NNcombOutput",&pholead_NNcombOutput, "pholead_NNcombOutput/F");
				myTree_->Branch("photrail_NNcombOutput",&photrail_NNcombOutput, "photrail_NNcombOutput/F");
				myTree_->Branch("pholead_etaSC",&pholead_etaSC,"pholead_etaSC/F");
				myTree_->Branch("photrail_etaSC",&photrail_etaSC,"photrail_etaSC/F");
				myTree_->Branch("pholead_SCbr",&pholead_SCbr,"pholead_SCbr/F");
				myTree_->Branch("photrail_SCbr",&photrail_SCbr,"photrail_SCbr/F");
				myTree_->Branch("pholead_isMin",&pholead_isMin,"pholead_isMin/I");
				myTree_->Branch("photrail_isMin",&photrail_isMin,"photrail_isMin/I");
				myTree_->Branch("dipho_HLT_bit0",&dipho_HLT_bit0, "dipho_HLT_bit0/I");
				myTree_->Branch("dipho_HLT_bit1",&dipho_HLT_bit1, "dipho_HLT_bit1/I");
				myTree_->Branch("dipho_HLT_bit2",&dipho_HLT_bit2, "dipho_HLT_bit2/I");
				myTree_->Branch("dipho_HLT_bit3",&dipho_HLT_bit3, "dipho_HLT_bit3/I");
				myTree_->Branch("dipho_HLT_bit4",&dipho_HLT_bit4, "dipho_HLT_bit4/I");
				myTree_->Branch("dipho_HLT_bit5",&dipho_HLT_bit5, "dipho_HLT_bit5/I");
				myTree_->Branch("dipho_HLT_bit6",&dipho_HLT_bit6, "dipho_HLT_bit6/I");
				myTree_->Branch("dipho_HLT_bit7",&dipho_HLT_bit7, "dipho_HLT_bit7/I");
				myTree_->Branch("dipho_HLT_bit8",&dipho_HLT_bit8, "dipho_HLT_bit8/I");
				myTree_->Branch("dipho_HLT_bit9",&dipho_HLT_bit9, "dipho_HLT_bit9/I");
				myTree_->Branch("dipho_HLT_bit10",&dipho_HLT_bit10, "dipho_HLT_bit10/I");
				myTree_->Branch("dipho_HLT_bit11",&dipho_HLT_bit11, "dipho_HLT_bit11/I");
				myTree_->Branch("dipho_Kfactor",&dipho_Kfactor,"dipho_Kfactor/F");
				myTree_->Branch("pholead_trueE",&pholead_trueE,"pholead_trueE/F");
				myTree_->Branch("pholead_truePx",&pholead_truePx,"pholead_truePx/F");
				myTree_->Branch("pholead_truePy",&pholead_truePy,"pholead_truePy/F");				
				myTree_->Branch("pholead_truePz",&pholead_truePz,"pholead_truePz/F");
				myTree_->Branch("pholead_truePhi",&pholead_truePhi,"pholead_truePhi/F");
				myTree_->Branch("pholead_trueEta",&pholead_trueEta,"pholead_trueEta/F");
				myTree_->Branch("photrail_trueE",&photrail_trueE,"photrail_trueE/F");
				myTree_->Branch("photrail_truePx",&photrail_truePx,"photrail_truePx/F");
				myTree_->Branch("photrail_truePy",&photrail_truePy,"photrail_truePy/F");				
				myTree_->Branch("photrail_truePz",&photrail_truePz,"photrail_truePz/F");
				myTree_->Branch("photrail_truePhi",&photrail_truePhi,"photrail_truePhi/F");
				myTree_->Branch("photrail_trueEta",&photrail_trueEta,"photrail_trueEta/F");
				myTree_->Branch("pholead_isMatchingWithHLTObject",&pholead_isMatchingWithHLTObject,"pholead_isMatchingWithHLTObject/I");
				myTree_->Branch("photrail_isMatchingWithHLTObject",&photrail_isMatchingWithHLTObject,"photrail_isMatchingWithHLTObject/I");
				myTree_->Branch("event_nPhotons",&event_nPhotons,"event_nPhotons/I");
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
	pholead_sigieta = pho_sigmaIetaIeta_loc[iteLead];
 	pholead_HcalIso = HcalIso_loc[iteLead];
	pholead_EcalIso = EcalIso_loc[iteLead];
	pholead_TrackerIso = TrackerIso_loc[iteLead];
 	pholead_HcalIsodR03 = HcalIsodR03_loc[iteLead];
	pholead_EcalIsodR03 = EcalIsodR03_loc[iteLead];
	pholead_TrackerIsodR03 = TrackerIsodR03_loc[iteLead];
	pholead_hoe = HoE_loc[iteLead];
	pholead_isPromptGenPho =  isPromptGenPho_loc[iteLead];
	pholead_isFromQuarkGen = isFromQuarkGen_loc[iteLead];
	pholead_HasPixSeed = HasPixSeed_loc[iteLead];
	pholead_isEB = pho_isEB_loc[iteLead];
	pholead_isEE = pho_isEE_loc[iteLead];
	pholead_etaSC = etaSC_loc[iteLead];
	pholead_SCbr = SCbr_loc[iteLead];
	pholead_DrScClosest = closestSC_dR_loc[iteLead];
	pholead_DrTrkClosest = dR_SCtrkclosest_loc[iteLead];
	pholead_transverseToJetRatio = transverseToJetRatio_loc[iteLead];
	pholead_trueE = pho_trueE_loc[iteLead];
	pholead_truePx = pho_truePx_loc[iteLead];
	pholead_truePy = pho_truePy_loc[iteLead];
	pholead_truePz = pho_truePz_loc[iteLead];
	pholead_trueEta = pho_trueEta_loc[iteLead];
	pholead_truePhi = pho_truePhi_loc[iteLead];
	if (doNN) {
		pholead_NNshapeOutput = NNshapeOutput_loc[iteLead];
		pholead_NNenvOutput = NNenvOutput_loc[iteLead];
		pholead_NNcombOutput = NNcombOutput_loc[iteLead];
	}
	if (doHLTobject){
		pholead_isMatchingWithHLTObject = pho_isMatchingWithHLTObject_loc[iteLead];
	}
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
 	photrail_HcalIsodR03 = HcalIsodR03_loc[iteTrail];
	photrail_EcalIsodR03 = EcalIsodR03_loc[iteTrail];
	photrail_TrackerIsodR03 = TrackerIsodR03_loc[iteTrail];
	photrail_hoe = HoE_loc[iteTrail];
	photrail_isPromptGenPho =  isPromptGenPho_loc[iteTrail];
	photrail_isFromQuarkGen = isFromQuarkGen_loc[iteTrail];
	photrail_HasPixSeed = HasPixSeed_loc[iteTrail];
	photrail_isEB = pho_isEB_loc[iteTrail];
	photrail_isEE = pho_isEE_loc[iteTrail];
	photrail_etaSC = etaSC_loc[iteTrail];
	photrail_SCbr = SCbr_loc[iteTrail];
	photrail_DrScClosest = closestSC_dR_loc[iteTrail];
	photrail_DrTrkClosest = dR_SCtrkclosest_loc[iteTrail];
	photrail_transverseToJetRatio = transverseToJetRatio_loc[iteTrail];
	photrail_trueE = pho_trueE_loc[iteTrail];
	photrail_truePx = pho_truePx_loc[iteTrail];
	photrail_truePy = pho_truePy_loc[iteTrail];
	photrail_truePz = pho_truePz_loc[iteTrail];
	photrail_trueEta = pho_trueEta_loc[iteTrail];
	photrail_truePhi = pho_truePhi_loc[iteTrail];
	if (doNN) {
		photrail_NNshapeOutput = NNshapeOutput_loc[iteTrail];
		photrail_NNenvOutput = NNenvOutput_loc[iteTrail];
		photrail_NNcombOutput = NNcombOutput_loc[iteTrail];
	}
	if (doHLTobject){
		photrail_isMatchingWithHLTObject = pho_isMatchingWithHLTObject_loc[iteTrail];
	}
	TLorentzVector Psum = P_loc[iteLead] + P_loc[iteTrail];
	TLorentzVector Ptruesum = Ptrue_loc[iteLead]+Ptrue_loc[iteTrail];
	dipho_mgg = Psum.M();
	dipho_mggTrue = Ptruesum.M();
	dipho_costhetastar = fabs(CosThetaStar(P_loc[iteLead],P_loc[iteTrail]));
	dipho_costhetastarTrue = fabs(CosThetaStar(Ptrue_loc[iteLead],Ptrue_loc[iteTrail]));
	dipho_etastar = 1.0*(P_loc[iteLead].Eta()-P_loc[iteTrail].Eta())/2;
	dipho_qt = Psum.Pt();
	dipho_qtTrue = Ptruesum.Pt();
	dipho_ql = Psum.Pz();
	dipho_deltaphi = DeltaR(P_loc[iteLead].Phi(),P_loc[iteTrail].Phi(),0,0);
	dipho_deltaphiTrue = DeltaR(Ptrue_loc[iteLead].Phi(),Ptrue_loc[iteTrail].Phi(),0,0);
	dipho_deltaR = DeltaR(P_loc[iteLead].Phi(),P_loc[iteTrail].Phi(),P_loc[iteLead].Eta(),P_loc[iteTrail].Eta());
	dipho_minIsoTrack = findTheMini(pholead_TrackerIso,photrail_TrackerIso);
	dipho_maxIsoTrack = findTheMaxi(pholead_TrackerIso,photrail_TrackerIso);
	if (findTheMini(pholead_TrackerIso,photrail_TrackerIso)==pholead_TrackerIso) {
		pholead_isMin = 1;
		photrail_isMin = 0;
	}
	else {
		pholead_isMin = 0;
		photrail_isMin = 1;
	}		
	if (doNN){
		dipho_minNNshape = findTheMini(photrail_NNshapeOutput,pholead_NNshapeOutput);
		dipho_minNNcomb = findTheMini(photrail_NNcombOutput,pholead_NNcombOutput);
		dipho_minNNenv = findTheMini(photrail_NNenvOutput,pholead_NNenvOutput);
	}
        if ((photrail_pt>15)&&(pholead_pt>15)&&(pholead_etaSC<2.5)&&(photrail_etaSC<2.5)&&(pholead_hoe<0.1)&&(photrail_hoe<0.1)&&(pholead_TrackerIsodR03<(7+0.002*pholead_pt))&&(photrail_TrackerIsodR03<(7+0.002*photrail_pt))&&(pholead_EcalIsodR03<(8.4+0.012*pholead_pt))&&(photrail_EcalIsodR03<(8.4+0.012*photrail_pt)&&(pholead_HcalIsodR03<(4.4+0.005*pholead_pt))&&(photrail_HcalIsodR03<(4.4+0.005*photrail_pt)))){
          NbEntries++;
         myTree_->Fill();
          }

}



//makeDiphotonTree(){
	int main(){
	doMC = false;
	doAllHlt = true;
	doHLTobject = true;
	doNN = false; 
	NbEntries = 0;
	myFile = new TFile("diphoton_part1.root","RECREATE");
	float MggMin = 0;
	beforeMacro();
	chain->Add("/sps/cms/hbrun/miniTree39X/run149291/output_*.root");
	int Nevents = chain->GetEntries();
//	Nevents = 10;
	std::cout << "Nevents = " << Nevents << std::endl;
	int ite = 0;
	int nbPhotons = 0;
	int numEventBefore, numEventHere, theNumbEvent, theNumbRun;
	int* leadTrail;
	for (int ievt = 0 ; ievt < Nevents ; ievt++){
		nbPhotons++;
		if (ievt != 0 ) { // fill the Event level infos
			event_runNumber = pho_RunNumber; 
			event_realNumber = pho_EventNumber;
			event_LumiSection = pho_lumiSection;
			event_nVertex = pho_nVertex;
			dipho_HLT_bit0 = pho_HLT_bit0;
			dipho_HLT_bit1 = pho_HLT_bit1;
			dipho_HLT_bit2 = pho_HLT_bit2;
			dipho_HLT_bit3 = pho_HLT_bit3;
			dipho_HLT_bit4 = pho_HLT_bit4;
			dipho_HLT_bit5 = pho_HLT_bit5;
			dipho_HLT_bit6 = pho_HLT_bit6;
			dipho_HLT_bit7 = pho_HLT_bit7;
			dipho_HLT_bit8 = pho_HLT_bit8;
			dipho_HLT_bit9 = pho_HLT_bit9;
			dipho_HLT_bit10 = pho_HLT_bit10;
			dipho_HLT_bit11 = pho_HLT_bit11;
			if (doMC){ 
				event_eventPtHat = pho_eventPtHat;
				event_processId = pho_eventProcessId;
			} 
		}
		chain->GetEntry(ievt);
		//if (!((pho_HLT_bit1==1))) continue;
		if ((ievt % 100000) == 0 ) std::cout << "on est au " << ievt << "evenement " << std::endl;
		if (ievt == 0) {numEventBefore = pho_iEvent;}
		numEventHere = pho_iEvent;
		if (!(numEventHere == numEventBefore)) {
			event_nPhotons = (nbPhotons-1);
			nbPhotons = 1; 
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
		//if (!((pho_HLT_bit1==1))) continue;
		if (isAspike == 1) continue;
		P_loc[ite].SetPxPyPzE(pho_px,pho_py,pho_pz,pho_energy);
		Ptrue_loc[ite].SetPxPyPzE(pho_truePx,pho_truePy,pho_truePz,sqrt(pho_truePx*pho_truePx+pho_truePy*pho_truePy+pho_truePz*pho_truePz));
//		std::cout << "P_loc = " << P_loc[ite].Et() << std::endl;
		Eta_loc[ite] = pho_eta;
		R9_loc[ite] = pho_r9;  
		cPP_loc[ite] = pho_cPP;
		EtaWidth_loc[ite] = pho_etawidth;
		ptoverjetpt_loc[ite] = pho_ptoverjetpt;
		dR_SCtrkclosest_loc[ite] = pho_DrTrkclosest_pt2;
		closestSC_dR_loc[ite] = pho_DrSCclosest;
		transverseToJetRatio_loc[ite] = pho_transverseToJetRatio;
		HoE_loc[ite] = pho_hoe;
		TrackerIso_loc[ite] = pho_IsoHollowTrkCone;
		EcalIso_loc[ite] = pho_IsoEcalRechit;
		HcalIso_loc[ite] = pho_IsoHcalRechit;
		TrackerIsodR03_loc[ite] = pho_IsoHollowTrkCone03;
		EcalIsodR03_loc[ite] = pho_IsoEcalRechit03;
		HcalIsodR03_loc[ite] = pho_IsoHcalRechit03;
		HasPixSeed_loc[ite] = pho_hasPixelSeed;
		isPromptGenPho_loc[ite] = pho_isPromptGenPho;
		isFromQuarkGen_loc[ite] = pho_isFromQuarkGen;
		pho_sigmaIetaIeta_loc[ite] = pho_sigmaIetaIeta;
		pho_isEB_loc[ite] = pho_isEB;
		pho_isEE_loc[ite] = pho_isEE;
		NNshapeOutput_loc[ite] = pho_NNshapeOutput;
		NNenvOutput_loc[ite] = pho_NNenvOutput;
		NNcombOutput_loc[ite] = pho_NNcombOutput;
		etaSC_loc[ite] = pho_SCeta;
		SCbr_loc[ite] = pho_SCbr;
		pho_trueE_loc[ite] = pho_trueE;
		pho_truePx_loc[ite] = pho_truePx;
		pho_truePy_loc[ite] = pho_truePy;
		pho_truePz_loc[ite] = pho_truePz;
		pho_truePhi_loc[ite] = pho_truePhi;
		pho_trueEta_loc[ite] = pho_trueEta;
		pho_isMatchingWithHLTObject_loc[ite] = pho_isMatchingWithHLTObject;
		ite++;

		numEventBefore = pho_iEvent;


		
	}
	std::cout << "NbEventBefore=" << Nevents << " NbEventAfter=" << NbEntries << std::endl;

	afterMacro();
}
