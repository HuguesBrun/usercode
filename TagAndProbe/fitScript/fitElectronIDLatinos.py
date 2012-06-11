import FWCore.ParameterSet.Config as cms
### USAGE:
###    cmsRun fitElecID.py <scenario> [ <id> [ <binning1> ... <binningN> ] ]
###
### scenarios:
###   - data_all (default)  
###   - signal_mc

import sys
args = sys.argv[1:]
if (sys.argv[0] == "cmsRun"): args =sys.argv[2:]
scenario = "data_all"
if len(args) > 0: scenario = args[0]
print "Will run scenario ", scenario 

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

Template = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    NumCPU = cms.uint32(1),
    SaveWorkspace = cms.bool(False),

    Variables = cms.PSet(
        mass = cms.vstring("Tag-elec Mass", "70", "130", "GeV/c^{2}"),
        pt = cms.vstring("elec p_{T}", "0", "1000", "GeV/c"),
        eta    = cms.vstring("elec #eta", "-2.5", "2.5", ""),
        abseta = cms.vstring("elec |#eta|", "0", "2.5", ""),
        absSCeta = cms.vstring("elec  |SC #eta|", "0","2.5",""),
        SCeta = cms.vstring("elec SC #eta", "-2.5","2.5",""),
        nVtx   = cms.vstring("Number of vertices", "0", "999", ""),
        weight = cms.vstring("weight","0","1000",""),
       ),

    Categories = cms.PSet(
                          passID   = cms.vstring("passID", "dummy[pass=1,fail=0]"),
                          passID2   = cms.vstring("passID2", "dummy[pass=1,fail=0]"),
                          passIso   = cms.vstring("passIso", "dummy[pass=1,fail=0]"),
                          passIP   = cms.vstring("passIP", "dummy[pass=1,fail=0]"),
                          passConvR   = cms.vstring("passConvR", "dummy[pass=1,fail=0]"),
                          passID2012   = cms.vstring("passID2012", "dummy[pass=1,fail=0]"),
                          passISO2012   = cms.vstring("passISO2012", "dummy[pass=1,fail=0]"),
                          passPogLoose   = cms.vstring("passPogLoose", "dummy[pass=1,fail=0]"),
                          passPogTight   = cms.vstring("passPogTight", "dummy[pass=1,fail=0]"),
                          passTrig   = cms.vstring("passTrig", "dummy[pass=1,fail=0]"),
                          passTrigPlusID   = cms.vstring("passTrigPlusID", "dummy[pass=1,fail=0]"),
                          passIdTrigIso   = cms.vstring("passIdTrigIso", "dummy[pass=1,fail=0]"),
                          tag_HLT_Ele27_WP80 = cms.vstring("tag_HLT_Ele27_WP80","dummy[pass=1,fail=0]"),
                          passIsoTrig = cms.vstring("passIsoTrig","dummy[pass=1,fail=0]"),
                          passElec_FO = cms.vstring("passElec_FO","dummy[pass=1,fail=0]"),
                          passElec_FO_ID = cms.vstring("passElec_FO_ID","dummy[pass=1,fail=0]"),
                          passElec_FO_ISO = cms.vstring("passElec_FO_ISO","dummy[pass=1,fail=0]"),
                          passElec_FO_ISO_ID = cms.vstring("passElec_FO_ISO_ID","dummy[pass=1,fail=0]"),
                          tag_AnyE1 = cms.vstring("tag_AnyE1","dummy[pass=1,fail=0]"),

                          
	
        mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
      
    ),
    PDFs = cms.PSet(
        voigtPlusExpo = cms.vstring(
            "Voigtian::signal(mass, mean[90,80,100], width[2.495], sigma[3,1,20])",
            "Exponential::backgroundPass(mass, lp[0,-5,5])",
            "Exponential::backgroundFail(mass, lf[0,-5,5])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        ),
        vpvPlusExpo = cms.vstring(
            "Voigtian::signal1(mass, mean1[90,80,100], width[2.495], sigma1[2,1,3])",
            "Voigtian::signal2(mass, mean2[90,80,100], width,        sigma2[4,2,10])",
            "SUM::signal(vFrac[0.8,0,1]*signal1, signal2)",
            "Exponential::backgroundPass(mass, lp[-0.1,-1,0.1])",
            "Exponential::backgroundFail(mass, lf[-0.1,-1,0.1])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        ),
        vpvPlusExpoMin70 = cms.vstring(
            "Voigtian::signal1(mass, mean1[90,80,100], width[2.495], sigma1[2,1,3])",
            "Voigtian::signal2(mass, mean2[90,80,100], width,        sigma2[4,3,10])",
            "SUM::signal(vFrac[0.8,0.5,1]*signal1, signal2)",
            "Exponential::backgroundPass(mass, lp[-0.1,-1,0.1])",
            "Exponential::backgroundFail(mass, lf[-0.1,-1,0.1])",
            "efficiency[0.9,0.7,1]",
            "signalFractionInPassing[0.9]"
        )
    ),

    binnedFit = cms.bool(True),
    binsForFit = cms.uint32(40),

    Efficiencies = cms.PSet(), # will be filled later
)

TRIGGER = cms.PSet(tag_Mu24 = cms.vstring("pass"))
if "mc" in scenario or "39X" in scenario or "38X" in scenario:
    TRIGGER = cms.PSet(tag_Mu15 = cms.vstring("pass"), tag_pt = cms.vdouble(24.,9999.))

PT_ETA_BINS = cms.PSet(
    pt     = cms.vdouble(  10, 15, 20, 25, 50, 150 ),
                       #pt     = cms.vdouble(  10, 20, 30, 40, 60, 100 ),
    abseta = cms.vdouble(  0.0, 1.4442, 1.556, 2.5),
    tag_HLT_Ele27_WP80 = cms.vstring("pass"),
                       #abseta = cms.vdouble(  0.0, 1.2, 2.4)
)

PT_ETA_BINS_XCHECK = cms.PSet(
                              pt  = cms.vdouble(20, 30, 40, 60, 100),
                              abseta = cms.vdouble(0, 1.479, 2.5),
                              tag_HLT_Ele27_WP80 = cms.vstring("pass"),
                              )


VTX_BINS  = cms.PSet(
    pt     = cms.vdouble(  10, 150 ),
    abseta = cms.vdouble(  0.0, 2.4),
	nVtx = cms.vdouble(0.5,2.5,4.5,6.5,8.5,10.5,12.5,14.5,16.5,18.5, 20.5, 22.5, 24.5, 26.5, 28.5, 30.5, 32.5, 34.5)
)

PT_ETA_BINS_SMURFS = cms.PSet(
                              pt = cms.vdouble(10,15,20,30,40,50,7000),
                              absSCeta = cms.vdouble(0, 0.8, 1.479, 2, 2.5),
                              tag_HLT_Ele27_WP80 = cms.vstring("pass"),
                              )








#PREFIX="/data/gpetrucc/7TeV/tnp/2011.02.17/"
PREFIX="file:/afs/cern.ch/work/h/hbrun/latinoTnP/newProdElec/"
process.TnP_ElecID = Template.clone(
    InputFileNames = cms.vstring(PREFIX+"tnpZ_withMVAIsoNew.root"),
    InputTreeName = cms.string("fitter_tree"),
    InputDirectoryName = cms.string("tpTreeElEl"),
    OutputFileName = cms.string("TnP_ElecID_%s.root" % scenario),
    Efficiencies = cms.PSet(),
)
if "mc" in scenario:
   process.TnP_ElecID.WeightVariable = cms.string("weight")

if "data" in scenario:
    if   "v1" in scenario: process.TnP_ElecID.InputFileNames = [ PREFIX+"tnpZ_HWWid2012.root" ]
    elif "v2" in scenario: process.TnP_ElecID.InputFileNames = [ PREFIX+"tnpZ_2011A_v2_GOLDEN.root" ]
    elif "huguesTest" in scenario: process.TnP_ElecID.InputFileNames = [ PREFIX+"tnpElectronID_runA.root", PREFIX+"tnpElectronID_runB.root" ]
    else: print "coucou c encore moi"
if "mc" in scenario:
    process.TnP_ElecID.InputFileNames = [ PREFIX+"tnpMC_weightsEle.root" ]

if "tag35" in scenario:
    process.TnP_ElecID.Variables.tag_pt[1]='35'

print "les fichiers que l'on va utiliser = ", process.TnP_ElecID.InputFileNames

#IDS = ["TOGCPFTIPMVA_from_TrackerOrGlobal"]
#IDS = [ "passPogLoose","passPogTight","passTrig","passTrigPlusID","passIdTrigIso","passIsoTrig","passID2012_from_passTrig","passISO2012_from_passTrig","passISO2012_from_passID2012","passID2012_from_passISO2012"]
IDS = ["passElec_FO_ISO_ID_from_passElec_FO_ISO","passElec_FO_ISO_ID_from_passElec_FO_ID","passElec_FO_ISO_ID_from_passISO2012"]

#ALLBINS = [("ptXcheck",PT_ETA_BINS_XCHECK)]
#ALLBINS = [("pt",PT_ETA_BINS)]
ALLBINS = [("ptEta_smurfs", PT_ETA_BINS_SMURFS)]
#ALLBINS += [("vtx",VTX_BINS)]




if len(args) > 1 and args[1] not in IDS: IDS += [ args[1] ]
for ID in IDS:
    if len(args) > 1 and ID != args[1]: continue
    for X,B in ALLBINS:
        if len(args) > 2 and X not in args[2:]: continue
        module = process.TnP_ElecID.clone(OutputFileName = cms.string("TnP_ElecID_%s_%s_%s.root" % (scenario, ID, X)))
        shape = "vpvPlusExpo"
        if "eta" in X and not "abseta" in X: shape = "voigtPlusExpo"
        if "pt_abseta" in X: shape = "voigtPlusExpo"
        if X.find("pt_abseta") != -1: module.Variables.mass[1]="77";
        if X.find("overall") != -1: module.binsForFit = 120
        DEN = B.clone(); num = ID;
        if "24" in ID and hasattr(DEN,'pt') and "pt" not in X: DEN.pt[0] = 25
        if "_from_" in ID:
            parts = ID.split("_from_")
            num = parts[0]
            setattr(DEN, parts[1], cms.vstring("pass"))
        if scenario.find("tagiso") != -1:  
            DEN.tag_combRelIso = cms.vdouble(-1, 0.1)
        if scenario.find("loosetagiso") != -1:  
            DEN.tag_combRelIso = cms.vdouble(-1, 0.2)
        if scenario.find("probeiso") != -1:
            DEN.isoTrk03Abs = cms.vdouble(-1, 3)
        #if scenario.find("calo") != -1: DEN.caloCompatibility = cms.vdouble(0.9,1.1)  # same as above, I think.
        if "mc" in scenario:
            if num == "Mu24": num = "Mu15"
            if num == "IsoMu17": num = "IsoMu15"
            if num == "DoubleMu7": num = "DoubleMu3"
            if num == "Mu8_forEMu": num = "DoubleMu3"
            if num == "Mu17_forEMu": num = "DoubleMu3"
        if "EG5" in scenario: DEN.pair_nL1EG5 = cms.vdouble(0.5,999)
	if "data" in scenario:
	    setattr(module.Efficiencies, ID+"_"+X, cms.PSet(
		EfficiencyCategoryAndState = cms.vstring(num,"pass"),
	 	UnbinnedVariables = cms.vstring("mass"),
		BinnedVariables = DEN,
		BinToPDFmap = cms.vstring(shape)
	    ))
        if scenario.find("mc") != -1:
	    setattr(module.Efficiencies, ID+"_"+X, cms.PSet(
		EfficiencyCategoryAndState = cms.vstring(num,"pass"),
		UnbinnedVariables = cms.vstring("mass","weight"),
		BinnedVariables = DEN,
		BinToPDFmap = cms.vstring(shape)
	    ))
        setattr(process, "TnP_ElecID_"+ID+"_"+X, module)        
        setattr(process, "run_"+ID+"_"+X, cms.Path(module))

