#!/bin/csh

setenv castorpath /castor/cern.ch/user/h/hbrun/cluster_RECO_CMSSW_3_1_2

set list = (`rfdir $castorpath | grep  allcorrnewcut | awk '{print $9}'`)
foreach i ($list)
echo "$i"
stager_get -M $castorpath/$i
end
