#Installation

## Step 1
The first part of setup is from [this](https://github.com/pallabidas/Run3Ntuplizer/tree/2020_Apr_5-boosted_tauVeto) + email instructions:

```
cmsrel CMSSW_10_6_20
cd CMSSW_10_6_20/src/
cmsenv
scram b -j 8
git cms-init
git cms-merge-topic isobelojalvo:run3-dev-boosted-CMSSW_10_6_0_pre4
```

```
cd L1Trigger/
git clone -b 2020_May_21-boosted_tauVeto git@github.com:pallabidas/L1TCaloSummary.git
git clone -b 2020_Apr_5-boosted_boostedTau git@github.com:pallabidas/Run3Ntuplizer.git
emacs L1TCaloLayer1/src/UCTRegion.hh
```

Move the following lines from `protected` to `public`:
```
std::vector<UCTTower*> towers;
```

Continue:

```
cd ..
USER_CXXFLAGS="-Wno-error=reorder -Wno-unused-variable" scram b -j 8
```

## Step 2
The second part of setup is from [this](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePFTauID#Running_of_the_bug_fix_version_o).
In a new, separate directory:

```
# in a new, separate directory

cmsrel CMSSW_10_2_22

cd CMSSW_10_2_22/src 

cmsenv

git cms-merge-topic -u cms-tau-pog:CMSSW_10_2_X_tau-pog_boostedTausMiniFix
```

## Step 3
The third step is to check-out PhysicsTools/PatAlgos and RecoTauTag/RecoTau in `10_6_20` (where L1Trigger is already installed)

```
# cd back to CMSSW_10_6_20/src

cmsenv

git cms-addpkg PhysicsTools/PatAlgos

git cms-addpkg RecoTauTag/RecoTau
```

## Step 4
Copy and paste the changed files from Step 2 into these packages of Step 3. These are the files:

```
# (from Pallabi's second email)

PhysicsTools/PatAlgos/plugins/PATBoostedTauCleaner.cc

RecoTauTag/RecoTau/plugins/DeepBoostedTauId.cc

RecoTauTag/RecoTau/python/tools/runBoostedTauIdMVA.py

RecoTauTag/RecoTau/test/BoostedTauAnalyzer.cc

RecoTauTag/RecoTau/test/BuildFile.xml

RecoTauTag/RecoTau/test/boostedTauAnalyzer_cfg.py
```

For my convenience, these are the commands I ran from `CMSSW_10_6_20/src` to do this step:

```
# cd to CMSSW_10_6_20/src

cp ../../CMSSW_10_2_22/src/PhysicsTools/PatAlgos/plugins/PATBoostedTauCleaner.cc PhysicsTools/PatAlgos/plugins/ 

cp ../../CMSSW_10_2_22/src/RecoTauTag/RecoTau/plugins/DeepBoostedTauId.cc RecoTauTag/RecoTau/plugins/

cp ../../CMSSW_10_2_22/src/RecoTauTag/RecoTau/python/tools/runBoostedTauIdMVA.py RecoTauTag/RecoTau/python/tools/ 

cp ../../CMSSW_10_2_22/src/RecoTauTag/RecoTau/test/BoostedTauAnalyzer.cc RecoTauTag/RecoTau/test/   

cp ../../CMSSW_10_2_22/src/RecoTauTag/RecoTau/test/BuildFile.xml RecoTauTag/RecoTau/test/BuildFile.xml 

cp ../../CMSSW_10_2_22/src/RecoTauTag/RecoTau/test/boostedTauAnalyzer_cfg.py RecoTauTag/RecoTau/test/
```

## Step 5
Create a new header file `test_TauMVAId.h` within `RecoTauTag/RecoTau/interface` since the original header file `PFRecoTauClusterVariables.h` has major changes between `10_6_x` and `10_2_x`.

```
# still in CMSSW_10_6_20/src from Step 4

cp /afs/cern.ch/work/p/pdas/public/forStephanie/boostedtaus/RecoTau/interface/test_TauMVAId.h RecoTauTag/RecoTau/interface/
```

## Step 6
Copy and paste these files from Pallabi's area to `CMSSW_10_6_20` (some of these have an additional include statement of the `test_TauMVAId.h` header file)

```
# still in CMSSW_10_6_20/src

cp /afs/cern.ch/work/p/pdas/public/forStephanie/boostedtaus/RecoTau/plugins/DeepTauId.cc RecoTauTag/RecoTau/plugins/

cp /afs/cern.ch/work/p/pdas/public/forStephanie/boostedtaus/RecoTau/plugins/DeepBoostedTauId.cc RecoTauTag/RecoTau/plugins

cp /afs/cern.ch/work/p/pdas/public/forStephanie/boostedtaus/RecoTau/plugins/PATTauDiscriminationByMVAIsolationRun2.cc RecoTauTag/RecoTau/plugins/
  
cp /afs/cern.ch/work/p/pdas/public/forStephanie/boostedtaus/RecoTau/plugins/PFRecoTauDiscriminationByMVAIsolationRun2.cc RecoTauTag/RecoTau/plugins/

cp /afs/cern.ch/work/p/pdas/public/forStephanie/boostedtaus/RecoTau/test/BuildFile.xml RecoTauTag/RecoTau/test/

cp /afs/cern.ch/work/p/pdas/public/forStephanie/boostedtaus/RecoTau/test/rerunMVAIsolationOnMiniAOD.cc RecoTauTag/RecoTau/test/ 
```

## Step 7
Compile everything in `CMSSW_10_6_20`:

```
cmsenv
USER_CXXFLAGS="-Wno-error=reorder -Wno-unused-variable" scram b -j 8
```

## Step 8
After succesful compilation, test it:

```
# again in CMSSW_10_6_20
cmsenv
cd L1Trigger/Run3Ntuplizer/test
cmsRun testL1TCaloSummary-ggHtautau.py
```

This creates a file (`l1TNtuple-ggHtautau.root`) in the same directory.



## Original instructions
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
