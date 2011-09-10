import re, os

fileName = "batch"

file = open("list","r")
files = file.readlines()
file.close()

file = open("template_"+fileName+".csh","r")
scriptLine = file.readlines()
file.close()

fileNumber = 1
counter = 0
while counter < len(files):
        counter += 1
	dossier = files[counter-1][:-1];
        fileToCopy = fileName+"_"+dossier+".csh"
	print dossier
        file = open(fileToCopy,"w")
        for line in scriptLine:
                if len(re.split("theEXE",line))> 1:
                        file.write("set NOMCHAINE=\"theEXE_"+dossier+"\"\n")
                        continue
		if len(re.split("cpsortie",line))> 1:
			file.write("cp diphoton_"+dossier+".root $SORTIEDIR\n")
			continue
		if len(re.split("rmsortie",line))> 1:
			file.write("rm diphoton_"+dossier+".root\n")
			continue
                file.write(line)
        fileNumber +=1
        file.close()

