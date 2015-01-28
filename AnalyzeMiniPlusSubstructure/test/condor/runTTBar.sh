#!/bin/bash

#source /uscmst1/prod/sw/cms/shrc uaf
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /uscms_data/d2/skhalil/YPrime/CMSSW_7_3_0/src
#cmsenv
eval `scram runtime -sh`

cd ${_CONDOR_SCRATCH_DIR}

let "sample=${1}+1"
cp /uscms_data/d2/skhalil/YPrime/CMSSW_7_3_0/src/MiniValidation/AnalyzeMiniPlusSubstructure/test/ttbar/ttbar_${sample}.py .
cmsRun ttbar_${sample}.py
mv *.root /eos/uscms/store/user/skhalil/YPrime/Hist/ttbar
rm ttbar_${sample}.py
ls
echo "DONE!"

