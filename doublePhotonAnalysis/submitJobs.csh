#!/bin/csh

set list = (`ls batch*.csh`) 

foreach i ($list)
	set name = (`echo $i | awk -F"_" '{print $2}' | awk -F"." '{print $1}'`)
	echo $name
	qsub $i -eo -o /sps/cms/hbrun/DiPhotons42X/diPhotonMC/$name.log -q T -N $name -l  M=4096MB,scratch=30720MB,platform=LINUX,u_sps_cmsf
#	qdel $name
end 
