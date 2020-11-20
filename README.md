Installation:

```

cmsrel CMSSW_10_6_10
cd CMSSW_10_6_10/src/
cmsenv
git cms-init
git cms-merge-topic isobelojalvo:run3-dev-boosted-CMSSW_10_6_0_pre4
cd L1Trigger
git clone -b 2020_May_21-boosted_tauVeto git@github.com:pallabidas/L1TCaloSummary.git
git clone -b 2020_Apr_5-boosted_tauVeto git@github.com:pallabidas/Run3Ntuplizer.git
cd ..
USER_CXXFLAGS="-Wno-error=reorder -Wno-unused-variable" scram b -j 8

cd L1Trigger/Run3Ntuplizer/test

## to run zero bias data
cmsRun testL1TCaloSummary-ZeroBias.py

## to run susy ggHbb events
cmsRun testL1TCaloSummary-ggHbb.py


```
