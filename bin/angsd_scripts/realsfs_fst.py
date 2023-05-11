#run from directory containing all your SAF files
#will run realSFS to make 2D-SFS for each pair of populations
#run Hudson's Fst estimator (see Bhatia 2013)

import glob
import itertools as it
import os

files =  glob.glob('*.saf.idx')

for a, b in it.combinations(files, 2):
			prefa = a.removesuffix('.wg.saf.idx')
			prefb = b.removesuffix('.wg.saf.idx')
			outsfs = prefa + '_' + prefb + '.2dsfs'
			fstfile = prefa + '_' + prefb
			#call realsfs
			#os.system('~/packages/angsd/misc/realSFS ' + a + ' ' + b + ' -fold 1 -P 20 '+ ' > fst/' + outsfs)
			os.system('~/packages/angsd/misc/realSFS ' + 'fst index -whichFst 1 ' + a + ' ' + b  + ' -sfs fst/' + outsfs + ' -fstout fst/' + fstfile + ' -P 20')
			os.system('~/packages/angsd/misc/realSFS ' + 'fst stats2 fst/' + fstfile + '.fst.idx ' + '-win 1000 -step 200 > fst/' + fstfile + '.window1step0p2.txt' + ' -P 20')
			os.system('~/packages/angsd/misc/realSFS ' + 'fst stats fst/' + fstfile + '.fst.idx ' + '> fst/' + fstfile + '.globalstat.txt')
