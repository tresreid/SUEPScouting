# To check out
```
cmsrel CMSSW_10_6_0 #you can use CMSSW_11_1_0
cd CMSSW_10_6_0/src
git clone -b mods https://github.com/SUEPPhysics/SUEPScouting.git PhysicsTools/SUEPScouting
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
