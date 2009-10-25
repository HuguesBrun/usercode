#!/bin/csh
setenv castorpath /castor/cern.ch/user/h/hbrun/cluster_RECO_CMSSW_3_1_2

set list = (`rfdir $castorpath | grep allcorrnewcut | awk '{print $9}'`)
set old_list = (`cat old_list`)
foreach i ($list)
set j = 0
foreach k ($old_list)
#		echo "salut $i"
	if ($k == $i) then
#			echo "t es dans la black list !"
#		set j = 1
	endif
end
if ($j == 1) then
continue
endif
echo "rfio:$castorpath/$i"
end
