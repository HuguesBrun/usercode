import FWCore.ParameterSet.Config as cms

process = cms.Process("READTRIGGER")

process.load('Configuration/EventContent/EventContent_cff')

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#"file:/tmp/hbrun/outputFULL.root"
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_1_1_ojk.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_2_1_szQ.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_3_1_323.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_4_1_XvA.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_5_1_IA1.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_6_1_D3P.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_7_1_HBZ.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_8_1_nQz.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_9_1_nkU.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_10_1_3NC.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_11_1_BJB.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_12_1_n5C.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_13_1_HTe.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_14_1_YEl.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_15_1_dQF.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_16_1_k15.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_17_1_NIF.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_18_1_r8q.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_19_1_6wt.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_20_1_tMn.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_21_1_0GH.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_22_1_8Mv.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_23_1_kaW.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_24_1_iUZ.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_25_1_lMI.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_26_1_GTw.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_27_1_lEf.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_28_1_cm6.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_29_1_l32.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_30_1_dST.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_31_1_cMy.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_32_1_Ypz.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_33_1_JyK.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_34_1_fjD.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_35_1_prv.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_36_1_8xI.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_37_1_4zG.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_38_1_002.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_39_1_oVU.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_40_1_M3S.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_41_1_EpL.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_42_1_2sz.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_43_1_iQo.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_44_1_5PT.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_45_1_6fy.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_46_1_vhZ.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_47_1_TM0.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_48_1_zj5.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_49_1_Hjf.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_50_1_bCF.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_51_1_LEK.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_52_1_QME.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_53_1_SGm.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_54_1_buW.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_55_1_XpP.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_56_1_PSO.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_57_1_QtE.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_58_1_8F8.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_59_1_VaN.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_60_1_QVc.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_61_1_Kfu.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_62_1_czi.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_63_1_eYp.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_64_1_GfY.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_65_1_zYx.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_66_1_OZL.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_67_1_qDh.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_68_1_MRN.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_69_1_QEf.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_70_1_Lmf.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_71_1_QWP.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_72_1_2Nd.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_73_1_elB.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_74_1_UDl.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_75_1_xxt.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_76_1_3V1.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_77_1_uZN.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_78_1_cZR.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_79_1_Nac.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_80_1_rLU.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_81_1_IiR.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_82_1_JX1.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_83_1_0Mk.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_84_1_ZzJ.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_85_1_Dwk.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_86_1_HMO.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_87_1_hvF.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_88_1_FOW.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_89_1_hw4.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_90_1_DP7.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_91_1_2hs.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_92_1_CUi.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_93_1_Mms.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_94_1_P1m.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_95_1_Y0I.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_96_1_zZY.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_97_1_HQx.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_98_1_wup.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_99_1_N34.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_100_1_fwu.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_101_1_Di9.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_102_1_kMz.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_103_1_tbS.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_104_1_X7h.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_105_1_B6R.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_106_1_qF7.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_107_1_4UZ.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_108_1_k6y.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_109_1_1mv.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_110_1_Ly1.root",
"rfio:///castor/cern.ch/user/j/jfan/GluGluCopyNew/outputFULLnew_111_1_vnL.root"
)
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)



process.readTheTrigger = cms.EDAnalyzer('ReadTrigger',
	hltProducer = cms.InputTag('TriggerResults','','TEST'),
        outputFile = cms.string('/tmp/hbrun/theTriggerResults.root')
)

process.p = cms.Path(process.readTheTrigger)

