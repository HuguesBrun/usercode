#!/bin/csh

setenv castorpath /castor/cern.ch/user/h/hbrun/cluster_RECO_CMSSW_3_1_2	

set list = (`rfdir $castorpath | grep  fetaet | awk '{print $9}'`)
foreach i ($list)
set j = (`echo "$i" | awk -F "." '{print $1}' | awk -F "_" '{print $2}'`)
rfrm $castorpath/$i 
end
