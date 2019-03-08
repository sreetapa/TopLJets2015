#Pre-compile macro

Check where fastjet is using `scram tool info fastjet`.
Load it in ROOT and then compile the macro

```
root -l
gSystem->Load("/cvmfs/cms.cern.ch/slc6_amd64_gcc700/external/fastjet/3.3.0-omkpbe/lib/libfastjet.so")
.L runTTto2Lselection.C++
```