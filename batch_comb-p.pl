#!/usr/bin/perl

$exec = "/Janssen/bix/pkgs/combined-pvalues-0.48/cpv";
$projDIR = shift;
$dataDIR = "/Isilon/qli2/projects/Brain/results/DMR/$projDIR";

for ($i = 1; $i <= 22; $i++) {
	open (SH, ">$dataDIR/chr$i.sh") || die "cannot open";
	#$cmd = "/Janssen/bix/pkgs/miniconda2/bin/python $exec/pipeline.py -c 4 --seed 1e-4 --dist 500 -p chr$i --region-filter-p 0.1 --anno hg19 $dataDIR/chr$i.output.bed"; 
	$cmd = "/Janssen/bix/pkgs/miniconda2/bin/python $exec/acf.py -d 1:500:50 -c 4 $dataDIR/chr$i.output.bed > $dataDIR/chr$i.acf.txt";
	print SH "$cmd\n";
	$cmd = "/Janssen/bix/pkgs/miniconda2/bin/python $exec/slk.py --acf $dataDIR/chr$i.acf.txt -c 4 $dataDIR/chr$i.output.bed > $dataDIR/chr$i.pvals.acf.bed\n";
	$cmd .= "/Janssen/bix/pkgs/miniconda2/bin/python $exec/peaks.py --dist 500 --seed 1e-4 $dataDIR/chr$i.pvals.acf.bed > $dataDIR/chr$i.pvals.regions.bed\n";
	$cmd .= "/Janssen/bix/pkgs/miniconda2/bin/python $exec/region_p.py -p $dataDIR/chr$i.output.bed -c 4 -r $dataDIR/chr$i.pvals.regions.bed -s 50 >$dataDIR/chr$i.regions.sig.bed\n";
	print SH "$cmd\n";
	close SH;
	system("qsub $dataDIR/chr$i.sh");	
}
	
