To check out
```
cmsrel CMSSW_11_1_0 #You can use 10_X_X too
cd CMSSW_11_1_0/src
mkdir PhysicsTools
cd PhysicsTools
git clone https://github.com/SUEPPhysics/ScoutingNanoAOD.git 
```

To run 
```
cmsRun ScoutingNanoAOD_cfg.py fileList=sample.txt outputFile=sample.root maxEvents=10
```
