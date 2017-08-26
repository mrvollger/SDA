#!/bin/bash
module load ucsc/20140617

path=$(readlink -f $1)
echo $path
region=$(basename $path)
dups=$path/duplications.fixed.fasta


echo $region

project=$region
mkdir -p $project 
pushd $project
mkdir -p $region


URL="https://eee:muller41@eichlerlab.gs.washington.edu/help/mvollger/tracks/$project"
genome=genomes.txt
hub=hub.txt
track=$region/trackDb.txt
me=mvollger@uw.edu
group=$region/groups.txt
html=$region/description.html
twobit="$region/$region".2bit
twoinfo="$region/chromInfo.txt"


touch $genome $hub $track $group

# make 2bit
# fixed dups so I think the sed is no longer nessisary 
sed -e 's/:/_/g; s/-/_/g' $dups | faToTwoBit stdin $twobit
twoBitInfo $twobit $twoinfo 


# define default position
chr=$(head -n 1 $twoinfo | awk '{print $1;}')
end=$(head -n 1 $twoinfo | awk '{print $2;}')
pos="$chr":1-"$end"


#
# populate files required for a hub
#

forgenome="genome $region
trackDb $track
twoBitPath $twobit
groups $group
description $region CHM1 contig
organism Human
defaultPos $pos
orderKey 4800
scientificName Homo sapiens
htmlPath $html 
"
printf "$forgenome" > $genome

forhub="hub myhub
shortLabel my hub
longLabel My Hub
genomesFile $genome
email $me

"
printf "$forhub" > $hub

#
# generate bam tracks 
#
for g in $path/group.*/; do
    groupidx=$(echo $g | grep -o  'group.[0-9]\+')
    idx=$(echo $groupidx | grep -o  '[0-9]\+')
    urlBam=$URL/$region/WH."$idx".bam
    bam=$g/WH."$idx".bam
    lbam=$PWD/$region/WH."$idx".bam
   

    ln -s $bam $lbam
    ln -s "$bam".bai "$lbam".bai

    fortrack="track bam$idx
    shortLabel group $idx
    longLabel CC Group $idx
    bigDataUrl $urlBam
    type bam
    visibility full
    bamColorMode off
    group psvbam

    "
    printf "$fortrack" >> $track

done
#cat $track
echo "Tracks are done"

forgroup="
name psvbam
label PsvBam
priority 2
defaultIsClosed 0

"
printf "$forgroup" > $group
echo "Groups are done "



#
#
#
forhtml="
<h2>$region contig</h2>
 
<p>Assembly of the human $region region from CHM1 Mitchell Vollger.</p>

" 
printf "$forhtml" > $html
echo "HTML is done"




#
# FINISHING UP
#

# make a symbolic link to my home public direcotry 
rm -f $HOME/public_html/tracks/$project
ln -sf $(readlink -f $PWD) $HOME/public_html/tracks/$project
chmod 777 * 
chmod 777 $region/*
popd 
chmod 777 $project
# show the link
echo "$URL/hub.txt"


