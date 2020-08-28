export DATA=Banner.ADvsND.STG.robust.corrected
export ENDPOINT=STG.CC
bedtools sort -i $DATA.txt >$DATA.sorted.bed
mkdir $DATA
for chr in `cut -f 1 ${DATA}.sorted.bed | sort | uniq`; do grep -w $chr $DATA.sorted.bed > $DATA/$chr.output.bed; done
../../scripts/batch_comb-p.pl $DATA
cat $DATA/chr*.regions.sig.bed | awk -F "\t" '{if ($7 <0.05) {print $0}}' | awk -F "\t" '{if ($5 >= 3) {print $0}}' | sort -k1,1V -k2,2n -k3,3n > DMR.${ENDPOINT}.regions.sig.all.bed
/Janssen/bix/pkgs/homer/bin/annotatePeaks.pl DMR.${ENDPOINT}.regions.sig.all.bed hg19 -annStats DMR.${ENDPOINT} -go go.${ENDPOINT} -genomeOntology genomeOntology.${ENDPOINT} > DMR.${ENDPOINT}.regions.sig.all.anno 2>homer_${ENDPOINT}_output
cut -f 2-5,8,15-16,18-20 DMR.${ENDPOINT}.regions.sig.all.anno | awk '{if(NR==1){print $1"\t"$2"\t"$0} else {print $1"\t"$2-1"\t"$0}}'| cut -f 1-2,5- > DMR.${ENDPOINT}.regions.sig.all.anno.m
awk -F "\t" 'NR==FNR{a[$1"\t"$2"\t"$3]=$5"\t"$7}NR>FNR{if($1"\t"$2"\t"$3 in a){print $0"\t"a[$1"\t"$2"\t"$3]}}' DMR.${ENDPOINT}.regions.sig.all.bed DMR.${ENDPOINT}.regions.sig.all.anno.m > DMR.${ENDPOINT}.regions.sig.all.anno.m.final
cat DMR.${ENDPOINT}.regions.sig.all.anno.m.final | cut -f 7 | sort | uniq | wc -l
