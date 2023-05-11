#realSFS script to be run in directory containing SAF file output of angsd
#that has been run on by cluster by yearmonth samples
#outputs summary and windowed thetas

for f in `ls *.nomaffilter.saf.idx | sed s/.nomaffilter.saf.idx//g`
do
~/packages/angsd/misc/realSFS $f.nomaffilter.saf.idx -fold 1 -cores 10 > $f.nomaffilter.sfs
~/packages/angsd/misc/realSFS saf2theta $f.nomaffilter.saf.idx -outname $f.nomaffilter -sfs $f.nomaffilter.sfs -fold 1
~/packages/angsd/misc/thetaStat do_stat $f.nomaffilter.thetas.idx
~/packages/angsd/misc/thetaStat do_stat $f.thetas.idx -win 5000 -step 1000  -outnames $f.nomaffilter.thetasWindow.gz
done
