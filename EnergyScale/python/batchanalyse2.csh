#!/bin/csh
#echo "ca commence" | mail -s "debut job" hugues.ochiaz@orange.fr
set LOCALDIR=`pwd`
set EXECDIR="/afs/cern.ch/user/h/hbrun/scratch0/CMSSW_3_1_0_pre1_Cluster"

source /afs/cern.ch/cms/sw/cmsset_default.csh
cd $EXECDIR/src
eval `scramv1 runtime -csh`
cd Maravin/EnergyScale/python/python_gg2
set list = (`ls | grep py`)
foreach i ($list)
	cmsRun $i
end

