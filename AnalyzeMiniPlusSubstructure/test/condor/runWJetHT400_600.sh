#!/bin/bash

#source /uscmst1/prod/sw/cms/shrc uaf
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /uscms_data/d2/skhalil/YPrime/CMSSW_7_3_0/src
#cmsenv
eval `scram runtime -sh`
cd ${_CONDOR_SCRATCH_DIR}

let "sample=${1}+1"
cp /uscms_data/d2/skhalil/YPrime/CMSSW_7_3_0/src/MiniValidation/AnalyzeMiniPlusSubstructure/test/wjets_ht_400to600/wjets_ht_400to600_${sample}.py .
cmsRun wjets_ht_400to600_${sample}.py
mv *.root /eos/uscms/store/user/skhalil/YPrime/Hist/wjets_ht_400to600
rm wjets_ht_400to600_${sample}.py
ls
echo "DONE!"

