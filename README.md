Installation:

```
cmsrel CMSSW_10_6_0_pre4
cd CMSSW_10_6_0_pre4/src
cmsenv
git cms-init
git cms-merge-topic isobelojalvo:run3-dev-$CMSSW_VERSION
cd L1Trigger
git clone git@github.com:isobelojalvo/L1TCaloSummary.git
git clone --branch jul4 git@github.com:isobelojalvo/Run3Ntuplizer.git
cd ..
scram b -j 8
```
