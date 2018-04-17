

mkdir -p trackHub


if [ $1 == "remake" ]; then
	for region in $(cat regions.txt); do
		echo $region
		cd $region 
		~mvollger/projects/abp/bedForABP.py 
		cd ..
	done
fi



echo 'track name=ABP_Collapses type=bedDetail description="ABP" visibility=2 itemRgb="On"' \
	> trackHub/ABP.details.bed
cat */asm.bed >> trackHub/ABP.details.bed

echo ""
cp $PWD/trackHub/ABP.details.bed /net/eichler/vol2/home/mvollger/public_html/CHM1/ABP.details.bed 
head /net/eichler/vol2/home/mvollger/public_html/CHM1/ABP.details.bed
