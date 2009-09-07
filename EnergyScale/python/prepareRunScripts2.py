import re, random, os

#============================================================================
# Initial setting

numberOfJobs = 5  

fileName = "hybrid.py"

#============================================================================
# Run the cmsDriver job

file = open("list_ee","r")
files = file.readlines()
file.close()

file = open(fileName, "r")
scriptLines = file.readlines()
file.close()

fileNumber = 1
counter = 1
while counter < len(files):
    fileToCopy = "hybrid2_" + str(fileNumber) + ".py"
    fileNumber += 1

    file = open("python/" + fileToCopy, "w")
    for line in scriptLines:
        if len(re.split("initialSeed ", line)) > 1:
            file.write("        initialSeed = cms.untracked.uint32(" + str(int(random.random()*10000000)) + "),\n")
            continue
        
        # process.source
        if len(re.split("fileNames = cms.untracked.vstring", line)) > 1:
            list = "    fileNames = cms.untracked.vstring("
            print "before the loop", counter, counter%numberOfJobs, len(files)
            while counter%numberOfJobs <> 0 and len(files) > counter:
                print counter
                if counter <> len(files) and counter%numberOfJobs <> 0:
                    list += "'file:" + files[counter-1][:-1] + "',"
                    counter += 1
                if (counter <> len(files) and counter%numberOfJobs == 0) or (counter == len(files)):
                    list += "'file:" + files[counter-1][:-1] + "')\n"
                    counter += 1
                    break
            file.write(list)
            continue

        # process only the ECAL-related part
        if len(re.split("process.reconstruction\_step = cms.Path", line)) > 1:
            file.write("process.ecalLocalRecoSequence_nopreshower_makeClusters = cms.Sequence(process.ecalLocalRecoSequence_nopreshower*process.ecalClusters)\n")
            file.write("process.reconstruction_step = cms.Path(process.ecalLocalRecoSequence_nopreshower_makeClusters)\n")
            continue

        # rename output file
        if len(re.split("    outputFile = cms.string", line)) > 1:
            file.write("    outputFile = cms.string('DiElectronPt2to350_Hybrid_310_" + str(fileNumber-1) + ".root'),\n")
            continue
        file.write(line)
    file.close()
