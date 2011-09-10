import re, os

fileName = "makeDiphotonTree"

file = open("list","r")
files = file.readlines()
file.close()

file = open("template_"+fileName+".C","r")
scriptLine = file.readlines()
file.close()

fileNumber = 1
counter = 0
while counter < len(files):
        counter += 1
	dossier = files[counter-1][:-1];
        fileToCopy = fileName+"_"+dossier+".C"
	print dossier
        file = open(fileToCopy,"w")
        for line in scriptLine:
                if len(re.split("Add\(",line))> 1:
                        file.write("chain->Add(\"/sps/cms/hbrun/miniTree42XMC/"+ dossier +"/output_*.root\");\n")
# 			file.write("chain->Add(\"/sps/cms/hbrun/miniTree41XMC/miniTreeEnriched/"+ dossier + ".root\");\n")
                        continue
		if len(re.split("TFile",line))> 1:
			file.write("myFile = new TFile(\"diphoton_"+dossier+".root\",\"RECREATE\");\n")
			continue
                file.write(line)
        fileNumber +=1
        file.close()

