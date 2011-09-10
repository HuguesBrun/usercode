#!/bin/csh
#echo "ca commence" | mail -s "debut job" hbrun@cern.ch
# myName
set LOCALDIR=`pwd`
set EXECDIR="/sps/cms/hbrun/CMSSW_4_1_4_patch2/"
set PYTHONDIR="../DiPhotons42X/diPhotonMC"
set NOMCHAINE="theEXE"
set SORTIEDIR="/sps/cms/hbrun/DiPhotons42X/diPhotonMC/theDiphoMiniTreesOld"
source /afs/cern.ch/cms/sw/cmsset_default.csh
setenv SCRAM_ARCH slc5_amd64_gcc434 
if !(-d "$EXECDIR/$PYTHONDIR") then 
	echo "$theName" | mail -s "Encore un job qui s est degonfle " hbrun@cern.ch
endif
cd $EXECDIR 
cmsenv
cd $EXECDIR/$PYTHONDIR 
cp $NOMCHAINE $LOCALDIR
cd $LOCALDIR
./$NOMCHAINE 
cp cpsortie.root $SORTIEDIR 
rm rmsortie.root
#echo "Le job est termine" | mail -s "Fin job $PYTHONDIR" hbrun@cern.ch
