#!/bin/csh
		#g++ miniTreeMaker.C `root-config --libs --cflags` -o theEXE.exe
		g++ miniTreeMaker.C -L/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/src -lToto `root-config --libs --cflags` -o theEXE.exe
