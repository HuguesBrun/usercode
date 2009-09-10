#!/bin/csh
setenv chemin /afs/cern.ch/user/h/hbrun/scratch0/CMSSW_3_1_2_Cluster/src/Maravin/EnergyScale/python/python_gg 
#setenv chemin /castor/cern.ch/user/h/hbrun/outESAnalyser/runElectron

set list = (`ls $chemin | grep Di `)
set j = 0
foreach i ($list)
echo 'chain->Add("'$chemin'/'$i'");'
end
