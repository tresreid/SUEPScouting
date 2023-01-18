# To check out
```
cmsrel CMSSW_10_6_26 #you can use CMSSW_11_1_0
cd CMSSW_10_6_26/src
git clone -b mods https://github.com/SUEPPhysics/SUEPScouting.git PhysicsTools/SUEPScouting
git clone -b master https://github.com/tresreid/PatUtils.git PhysicsTools/PatUtils
```

# To setup
do this each time you login
```
cd src
cmsenv
voms-proxy-init
```

# To compile
```
cd src
scram b
```

# To run 
```
cmsRun ScoutingNanoAOD_cfg.py fileList=sample.txt outputFile=sample.root maxEvents=10
```

# To do
* jet energy corrections?
* job submission
