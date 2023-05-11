OUTDIR=/move_bioinfo/anopheles_cnrfp/angsd_output/by_year_by_tx
REF=/move_bioinfo/anopheles_03_21/basicwgs/ref/stash/vb_agamp3.fna
TNUM=20

for f in `ls /move_bioinfo/anopheles_cnrfp/poplists/by_year_by_tx`;
do

	numinds=`eval cat /move_bioinfo/anopheles_cnrfp/poplists/by_year_by_tx/$f | wc -l`
	mininds=`expr $numinds / 2`

	/home/dennistpw/packages/angsd/angsd  -rf /move_bioinfo/anopheles_cnrfp/regions.txt -b /move_bioinfo/anopheles_cnrfp/poplists/by_year_by_tx/$f -gl 1 -doGlf 2 -anc $REF -domajorminor 4 -doCounts 1 -trim 0 -domaf 1 -remove_bads 1 -only_proper_pairs 1 -baq 1 -nThreads $TNUM -doSaf 1 -uniqueonly 1 -out $OUTDIR/$f.wg.nomaffilter -minQ 20 -ref $REF -fai $REF.fai -C 50 -minMapQ 15
done
