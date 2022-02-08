# To check out
```
cmsrel CMSSW_11_1_0 #You can use 10_X_X too
cd CMSSW_11_1_0/src
mkdir PhysicsTools
cd PhysicsTools
git clone -b mods https://github.com/SUEPPhysics/SUEPScouting.git 
```

# To setup and compile
do this each time you login
```
cd $CMSSW_BASE/src
cmsenv
scram b
```

# To run 
```
cmsRun SUEPScouting/test/ScoutingNanoAOD_cfg.py inputFiles=file:aod.root outputFile=flatscouting.root maxEvents=1000000
```
