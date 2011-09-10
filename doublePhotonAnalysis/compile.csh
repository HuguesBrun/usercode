#!/bin/csh

set list = (`ls makeDiphotonTree_*`)

foreach i ($list)
	set nom = (`echo $i | awk -F "_" '{print $2}' | awk -F "." '{print $1}'`)
	if ( $nom == "template" ) then 
		continue
	endif
	echo  "on compile $i"
	#g++ $i `root-config --libs --cflags` -o theEXE_$nom
	g++ $i  -L/afs/in2p3.fr/grid/toolkit/cms2/slc5_amd64_gcc434/lcg/root/5.27.06b-cms14/lib -lTMVA -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMLP -lTreePlayer -lMinuit -pthread -lm -ldl -rdynamic -pthread -m64 -I/afs/in2p3.fr/grid/toolkit/cms2/slc5_amd64_gcc434/lcg/root/5.27.06b-cms14/include -o theEXE_$nom
end 
