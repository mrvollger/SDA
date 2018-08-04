#!/usr/bin/env bash
g=$1
b=${g%.*}
ABP=/net/eichler/vol5/home/mchaisso/projects/AssemblyByPhasing/scripts/abp/
$ABP/AddRepulsionEdges.py --graph $g --radius 10000 --out $b.repl.gml --load 10
echo "Done adding repulsion edges"
#$ABP/SimplifyGraph.py --graph $b.repl.gml  --max-nodes 1000 --out $b.repl.simpl.gml --snv assembly.consensus.fragments.snv --pos assembly.consensus.fragments.snv.pos
#echo "Done simplifying graph"
$ABP/RepulsionGraphToLP.py --graph $b.repl.gml --lp $b.repl.gml.lp  --wp 1 --wn 1  --relaxed
echo "Done writing LP"
# mipopt or optimize
$ABP/RunCplex.sh $b.repl.gml.lp optimize
echo "Done with CPLEX"
$ABP/RetainSolutionEdges.py --graph $b.repl.gml --sol $b.repl.gml.lp.sol --out $b.repl.sol.gml  --no-repulsion --min-value 0.5

$ABP/ComponentsToPhasedVCF.py $b.repl.sol.gml assembly.consensus.fragments.snv.pos assembly.consensus.nucfreq.vcf  --base cplex  --minComponent 4

$ABP/../analysis/simulation/CompareSNVS.py --vcf c*.vcf --snvs ../duplications.snvs --asmVcf assembly.consensus.nucfreq.vcf  --writeMat  --out CompareToSNVS.txt
