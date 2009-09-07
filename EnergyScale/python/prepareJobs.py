import re, os


file = open("analyser.template","r")
scriptLines =  file.readlines()
file.close()

counter = 1

while counter < 50:
    file = open("tourneanalyser_"+str(counter)+".csh","w")
    for line in scriptLines:
	if len(re.split("increment <",line)) > 1:
	    file.write("	if ($increment < "+str((counter-1)*10)+" ) then\n")
	    continue
	if len(re.split("increment >=",line)) > 1:
 	    file.write("	if ($increment >= "+str(counter*10)+" ) then\n")
	    continue
	file.write(line)
    file.close()
    counter += 1

