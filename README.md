#Running the analysis
To run on a single file can give the command
```
runTTto2Lselection --in /eos/cms/store/group/phys_top/PbPbTTbar_2018/SkimElectrons_PromptRecov1/Chunk_0_ext0.root --out test.root
```
To process a full directory use condor
```
dir=/eos/cms/store/group/phys_top/PbPbTTbar_2018/
out=`pwd`/plots/Chunks
tag=SkimElectrons_PromptRecov1
python scripts/launchAnalysis.py $dir/$tag $out $tag  
```
Massive processing example
```
tags=(SkimElectrons  SkimElectrons_PromptRecov1  SkimMuons_PromptRecov1  SkimMuons_PromptRecov2)
for tag in ${tags[@]}; do
    python scripts/launchAnalysis.py $dir/$tag $out $tag;
done   
```
