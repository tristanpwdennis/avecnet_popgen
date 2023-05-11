OUTDIR=/move_bioinfo/anopheles_cnrfp/angsd_output/by_cluster_by_tp
REF=/move_bioinfo/anopheles_03_21/basicwgs/ref/stash/vb_agamp3.fna
TNUM=20

for f in `ls poplists/by_cluster_by_tp`;
do
	/home/dennistpw/packages/angsd/angsd  -rf /move_bioinfo/anopheles_cnrfp/regions.txt -b /move_bioinfo/anopheles_cnrfp/poplists/by_cluster_by_tp/$f -gl 1 -anc $REF -domajorminor 4 -doCounts 1 -trim 0 -domaf 1 -remove_bads 1 -only_proper_pairs 1 -baq 1 -nThreads $TNUM -doSaf 1 -uniqueonly 1 -minmapQ 15 -out $OUTDIR/$f.wg.nomaffilter -minQ 20 -ref $REF -fai $REF.fai -C 50
done
