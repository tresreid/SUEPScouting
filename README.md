# To check out
```
cmsrel CMSSW_10_2_26 #You can use 10_X_X too
cd CMSSW_10_2_26/src
mkdir PhysicsTools
cd PhysicsTools
git clone https://github.com/SUEPPhysics/SUEPScouting.git 
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
