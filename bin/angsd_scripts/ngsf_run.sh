#for four seeds
for f in 12345 34856 39547 83749 
	do
		#initial estimation
		zcat cnrfp_allsamples_wg_AgamP4_3L.hd.glf.gz | ~/packages/ngsF/ngsF --n_ind 893 --n_sites 2294843 --glf -  -min_epsilon 1e-9 -init_values r --out cnrfp_3l.$f_seed.approxindF --n_threads 15 --seed $f --approx_EM
		#final estimation
		zcat cnrfp_allsamples_wg_AgamP4_3L.hd.glf.gz | ~/packages/ngsF/ngsF --n_ind 893 --n_sites 2294843 --glf -  -min_epsilon 1e-9 -init_values cnrfp_3l.$f_seed.approxindF.pars --out cnrfp_3l.$f__seed.final --n_threads 15 --seed $f
	done
