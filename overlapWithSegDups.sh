
cat */ref.fasta.bed |  sort -k1,1 -k2,2n | bedtools merge -i - > locationsInHumanReference.bed

mkdir -p overlap 

for i in $(seq 0.1 0.1 1.0); do
	echo $i
	bedtools intersect -a locationsInHumanReference.bed -b overlap/merged.segdup38.bed -wa -f $i \
		| bedtools merge -i - > overlap/"$i".bed 
done 



