echo "import FWCore.ParameterSet.Config as cms

readFiles = cms.untracked.vstring()

readFiles.extend([" > tmp_text.txt
nsls /folder/subfolder/subsubfolder |grep 'tnp' | awk -F'.root' '{print "\"rfio:/folder/subfolder/subsubfolder/"$1".root\","}' >> tmp_text.txt
echo "])

readFiles.extend([" >> tmp_text.txt
echo "])

source = cms.Source(\"PoolSource\",
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    fileNames = readFiles,
)" >> tmp_text.txt
mv tmp_text.txt MC_JPsi_XX_cfi.py
