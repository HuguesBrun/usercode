import re, random, os

#name of the event list
#theEventFile ="listEvent.txt"
theEventFile ="/sps/cms/sgandurr/CMSSW_4_2_8/src/PhysicsTools/Utilities/scripts/miniTree_v12_2011Bf_40_80_partALL.txt"

#load the Event List File
file = open(theEventFile,"r")
list = file.readlines()
file.close()

li = []

for theEventEntry in list:
    print theEventEntry[:-1]
    run = re.split(":",theEventEntry[:-1])[0]
    lumiSection = re.split(":",theEventEntry[:-1])[1]
    eventNumber = re.split(":",theEventEntry[:-1])[2]
    print "run=",run," lumiSection=",lumiSection," eventNumber=",eventNumber
    theDetEvent = (int(run), int(lumiSection), int(eventNumber))
    li.append(theDetEvent)

li.sort()
print len(li)
print li

localRun = 0
localLumi = 0
theCounter = 0
prevLumi = 0
theNextLi = (0,0,0)

file = open("theJasonFile.json","w")
file.write("{")
while theCounter < len(li):
    localLi = li[theCounter]
    if theCounter < (len(li)-1): 
        theNextLi = li[theCounter+1]
    if localRun != localLi[0]:
        localRun = localLi[0]
        if theCounter == 0:
            file.write('"'+str(localRun)+'":[')
        else:
            file.write('], "'+str(localRun)+'":[')
    else:
        file.write(", ")
    if localLi[1]==theNextLi[1]:
        theCounter += 1
    file.write('['+str(localLi[1])+','+str(localLi[1])+']')
    theCounter += 1
    print theCounter," == ",len(li)
    if theCounter >= len(li):
        file.write("]}")
file.close

file = open("template_python.py","r")
lines = file.readlines()
file.close()


file = open("../interface/theEvents.h","w")

file.write("int nEventToRead = "+str(len(li))+";\n")
file.write("int runNumbers["+str(len(li)+1)+"] = {")
for ite in li:
    file.write(str(ite[0])+",")
file.write("0};\n")

file.write("int eventNumbers["+str(len(li)+1)+"] = {")
for ite in li:
    file.write(str(ite[2])+",")
file.write("0};\n")








