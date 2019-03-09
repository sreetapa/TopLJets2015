## Search for PbPb->ttbar with 2018 data

The selection/plot filling is implemented in bin/runTTto2Lselection.cc.
It uses directly as inputs the HiForest contents.
The executable can be compiled with `scram b`.
To run on a single file for testing one  can give the command
```
runTTto2Lselection --in /eos/cms/store/group/phys_top/PbPbTTbar_2018/SkimElectrons_PromptRecov1/Chunk_0_ext0.root --out test.root
```

To loop over all the available forest trees better to use condor and then merge the outputs
```
dir=/eos/cms/store/group/phys_top/PbPbTTbar_2018/
out=`pwd`/plots/Chunks
tags=(SkimElectrons  SkimElectrons_PromptRecov1  SkimMuons_PromptRecov1  SkimMuons_PromptRecov2)
for tag in ${tags[@]}; do
    python scripts/launchAnalysis.py $dir/$tag $out $tag;
done   
```

Merge the outputs (hadds all the chunks according to the tag defined)
```
python scripts/mergeOutputs.py plots
```

Plot the results
```
FIXME
```