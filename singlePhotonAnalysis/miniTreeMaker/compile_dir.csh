#!/bin/csh
set listLocal = "list"
set theList = (`cat $listLocal `)
foreach j ($theList)
	echo $j
	set castorpath = ${j}_miniTree
	set list = (`ls $castorpath | grep .C | awk -F "." '{ print $1}'`)
	foreach i ($list)  
		echo $castorpath/$i
		g++ $castorpath/$i.C -L/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/src -lToto `root-config --libs --cflags` -o $castorpath/$i.exe
	end
end 
