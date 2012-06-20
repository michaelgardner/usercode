#!/bin/csh
cd _in_dir_
#cms
eval `scramv1 runtime -csh`
cd -
cmsRun _in_dir_/_runcfg_.py
cp *.root _in_dir_
#cp *.root /tmp/mgardner/MC_69

