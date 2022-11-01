# To check out
```
cmsrel CMSSW_10_6_0 #you can use CMSSW_11_1_0
cd CMSSW_10_6_0/src
git clone -b mods https://github.com/SUEPPhysics/SUEPScouting.git PhysicsTools/SUEPScouting
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
cmsRun SUEPScouting/test/ScoutingNanoAOD_cfg.py inputFiles=file:aod.root outputFile=flatscouting.root maxEvents=-1
```
