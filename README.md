# HGCtupleAna
PyROOT analyzer for HGCAL flat tuples

1. Install CMSSW `CMSSW_9_0_0_pre2`
2. Execute `cmsenv` everytime to get python and required libraries
3. `git clone git@github.com:beaudett/reco-ntuples.git` (or your fork) into src
4. Compile with `scram b`
5. `git clone git@github.com:artlbv/HGCtupleAna.git` (or your fork) into anywhere (src)
6. Follow next steps

To get all options, run:

* For simple rechit/particle deltaXYZ
```
python analyzeHGCtree.py --help
```
* For simple cluster2d/particle deltaXYZ
```
python analyze2dClusters.py --help
```
* For simple rechit density/distance demonstration
```
python analyzeClusters.py --help
```
