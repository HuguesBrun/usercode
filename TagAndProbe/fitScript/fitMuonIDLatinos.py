import FWCore.ParameterSet.Config as cms
### USAGE:
###    cmsRun fitMuonID.py <scenario> [ <id> [ <binning1> ... <binningN> ] ]
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
        mass = cms.vstring("Tag-muon Mass", "70", "130", "GeV/c^{2}"),
        pt = cms.vstring("muon p_{T}", "0", "1000", "GeV/c"),
        eta    = cms.vstring("muon #eta", "-2.5", "2.5", ""),
        abseta = cms.vstring("muon |#eta|", "0", "2.5", ""),
        phi    = cms.vstring("muon #phi at vertex", "-3.1416", "3.1416", ""),
        nvtx   = cms.vstring("Number of vertices", "0", "999", ""),
		weight = cms.vstring("weight", "0.0", "30.0", ""),

    ),

    Categories = cms.PSet(
                          passID = cms.vstring("passID", "dummy[pass=1,fail=0]"),
                          passID_2010 = cms.vstring("passID_2010", "dummy[pass=1,fail=0]"),
                          passIP = cms.vstring("passIP", "dummy[pass=1,fail=0]"),
                          passIso = cms.vstring("passIso", "dummy[pass=1,fail=0]"),
                          passGlb = cms.vstring("passGlb", "dummy[pass=1,fail=0]"),
                          passTM = cms.vstring("passTM", "dummy[pass=1,fail=0]"),
                          passGlbOrTM = cms.vstring("passGlbOrTM", "dummy[pass=1,fail=0]"),
                          passGlobalOrTracker = cms.vstring("passGlobalOrTracker", "dummy[pass=1,fail=0]"),
                          pass2012Tight = cms.vstring("pass2012Tight", "dummy[pass=1,fail=0]"),
                          pass2012ID_ICHEP = cms.vstring("pass2012ID_ICHEP", "dummy[pass=1,fail=0]"),
                          pass2012ISO_ICHEP = cms.vstring("pass2012ISO_ICHEP", "dummy[pass=1,fail=0]"),
                          pass2012ISOloose_ICHEP = cms.vstring("pass2012ISOloose_ICHEP", "dummy[pass=1,fail=0]"),
                          pass2012ICHEP = cms.vstring("pass2012ICHEP", "dummy[pass=1,fail=0]"),
                          pass2012ICHEPloose = cms.vstring("pass2012ICHEPloose", "dummy[pass=1,fail=0]"),
                          tag_AnyM1 = cms.vstring("tag_AnyM1","dummy[pass=1,fail=0]"),

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

    Efficiencies = cms.PSet(
    ), # will be filled later
)



TRIGGER = cms.PSet(tag_Mu24 = cms.vstring("pass"))
if "mc" in scenario or "39X" in scenario or "38X" in scenario:
    TRIGGER = cms.PSet(tag_Mu15 = cms.vstring("pass"), tag_pt = cms.vdouble(24.,9999.))



PT_ETA_BINS_REAL = cms.PSet(
                       pt   = cms.vdouble( 10, 15, 20, 50, 150),
                       abseta = cms.vdouble(0, 1.479, 2.5),
)

ETA_PT_POG = cms.PSet(
                      eta = cms.vdouble( -2.4, -2.1, -1.6, -1.1, -0.6, 0, 0.6, 1.1, 1.6, 2.1, 2.4),
                      pt = cms.vdouble( 20, 100),
                      
                      )

PT_ETA_POG = cms.PSet(
                      pt = cms.vdouble(20,30,40,60,100),
                      abseta = cms.vdouble(0,1.2,2.4),
                      
                      )
PT_ETA_POG_LARGE = cms.PSet(
                            abseta = cms.vdouble(0,1.2,2.4),
                            pt = cms.vdouble(20, 100),
                            )
                    
                      
PT_ETA_BINS_SMURFS = cms.PSet(
                              pt = cms.vdouble(10,15,20,30,40,50,7000),
                              abseta = cms.vdouble(0, 0.8, 1.2, 2.5),
    )




#PREFIX="/data/gpetrucc/7TeV/tnp/2011.02.17/"
PREFIX="file:/afs/cern.ch/work/h/hbrun/latinoTnP/"
process.TnP_MuonID = Template.clone(
    InputFileNames = cms.vstring(PREFIX+"tnpZ_withMVAIsoNew.root"),
    InputTreeName = cms.string("fitter_tree"),
    InputDirectoryName = cms.string("tpTreeMuMu"),
    OutputFileName = cms.string("TnP_MuonID_%s.root" % scenario),
    Efficiencies = cms.PSet(),
)


if "mc" in scenario:
    process.TnP_MuonID.WeightVariable = cms.string("weight")

if "data" in scenario:
    if   "v1" in scenario: process.TnP_MuonID.InputFileNames = [ PREFIX+"tnpZ_HWWid2012.root" ]
    elif "v2" in scenario: process.TnP_MuonID.InputFileNames = [ PREFIX+"tnpZ_2011A_v2_GOLDEN.root" ]
    elif "huguesTest" in scenario: process.TnP_MuonID.InputFileNames = [ PREFIX+"tnpMuonID.root" ]

if "mc" in scenario:
    process.TnP_MuonID.InputFileNames = [PREFIX + "tnpMC_weightsMu.root"]


if "tag35" in scenario:
    process.TnP_MuonID.Variables.tag_pt[1]='35'

print "les fichiers que l'on va utiliser = ", process.TnP_MuonID.InputFileNames

#IDS = ["TOGCPFTIPMVA_from_TrackerOrGlobal"]
IDS = [ "pass2012Tight", "pass2012ISO_ICHEP_from_pass2012ID_ICHEP","pass2012ID_ICHEP_from_pass2012ISO_ICHEP", "pass2012ICHEP","pass2012ICHEP_from_passGlobalOrTracker"]
#, "TOGCPF_from_TrackerOrGlobal", "TOGCPFT_from_TrackerOrGlobal", "TOGCPFTIP_from_TrackerOrGlobal", "TOGCPFTIPMVA_from_TrackerOrGlobal"]
#IDS += [ "TOGCPF_from_TOGclean", "TOGCPFT_from_TOGCPF","TOGCPFTIP_from_TOGCPFT","TOGCPFTIPMVA_from_TOGCPFTIP"]
ALLBINS = [("eta",ETA_PT_POG)]
ALLBINS += [("pt",PT_ETA_POG)]
ALLBINS += [("absEta",PT_ETA_POG_LARGE)]
ALLBINS += [("ptReal",PT_ETA_BINS_REAL)]

ALLBINS2 = [("smurfs",PT_ETA_BINS_SMURFS)]


if len(args) > 1 and args[1] not in IDS: IDS += [ args[1] ]
for ID in IDS:
    if len(args) > 1 and ID != args[1]: continue
    if "from" in ID:
        ALLBINS = ALLBINS2
    
    for X,B in ALLBINS:
        if len(args) > 2 and X not in args[2:]: continue
        module = process.TnP_MuonID.clone(OutputFileName = cms.string("TnP_MuonID_%s_%s_%s.root" % (scenario, ID, X)))
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
        setattr(module.Efficiencies, ID+"_"+X, cms.PSet(
            EfficiencyCategoryAndState = cms.vstring(num,"pass"),
            UnbinnedVariables = cms.vstring("mass","weight"),
            BinnedVariables = DEN,
            BinToPDFmap = cms.vstring(shape)
        ))
        if scenario.find("mc") != -1:
            setattr(module.Efficiencies, ID+"_"+X+"_mcTrue", cms.PSet(
                EfficiencyCategoryAndState = cms.vstring(num,"pass"),
                UnbinnedVariables = cms.vstring("mass","weight"),
                BinnedVariables = DEN.clone(mcTrue = cms.vstring("true"))
            ))
        setattr(process, "TnP_MuonID_"+ID+"_"+X, module)        
        setattr(process, "run_"+ID+"_"+X, cms.Path(module))

