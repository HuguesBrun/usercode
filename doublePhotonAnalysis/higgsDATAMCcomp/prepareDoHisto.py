import re, os

fileName = "do_hist_2photons"

file = open("list","r")
files = file.readlines()
file.close()

file = open(fileName+"MC.C","r")
scriptLine = file.readlines()
file.close()

fileNumber = 1
counter = 0
while counter < len(files):
        counter += 1
	dossier = files[counter-1][:-1];
	counter2 = 0
	while counter2 < 1:
	        fileToCopy = "do_hist_2photons_"+dossier+str(counter2)+".C"
		print dossier
      	        file = open(fileToCopy,"w")
        	for line in scriptLine:
			if len(re.split("name = \"mc\"", line)) > 1:
				file.write("TString name = \""+dossier+str(counter2)+"\";\n")
				continue
                	if len(re.split("Add\(",line))> 1:
                        	file.write("chain->Add(\"/sps/cms/hbrun/DiPhotons42X/diPhotonMC/theDiphoMiniTreesHiggs/diphoton_"+ dossier +".root\");\n")
                        	continue
			if len(re.split("myFile =",line))> 1:
				file.write("myFile = new TFile(\"diphoFile_"+dossier+str(counter2)+".root\",\"RECREATE\");\n")
				continue
			if len(re.split("do_hist_2photonsMC",line))> 1:
				file.write("do_hist_2photons_"+dossier+str(counter2)+"(){\n")
				continue
                        if len(re.split("theSelecType=1",line))> 1:
                                file.write("int theSelecType = "+str(counter2)+";")
                                continue
			if len(re.split("theCutAvant=",line)) > 1:
				if (dossier=="GJetPt20"):
					file.write("TString theCutAvant= \"&&(!(event_processId==18))\";\n")
					continue
               		file.write(line)
   		fileNumber +=1
        	file.close()
		counter2+=1

