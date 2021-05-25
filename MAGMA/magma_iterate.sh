#!/bin/bash

while getopts p: option
do
case "${option}"
in
p) pathway=${OPTARG};;
esac
done


magma_starter() {
	mkdir snp_${1}_percentage_${2}_h2_${3}
	cd snp_${1}_percentage_${2}_h2_${3}
     	tar -xzf /media/array/phenome_proj/LSEA/benchmarking/multiple_pathway_test_2.1_2/$4/snp_${1}_percentage_${2}_h2_${3}/Associations.tar.gz
      	$5/filter_script.sh
        for P in `seq 1 30`; do
        	$5/magma/magma --annotate --snp-loc A_target.P_${P}.snploc \
              	--gene-loc $5/data_magma/NCBI37.3.gene.loc --out annotation_P${P}

                $5/magma/magma --bfile $5/data_magma/g1000_eur \
                --pval A_target.P_${P}.qassoc.filtered_magma N=1000 --gene-annot annotation_P${P}.genes.annot --out genes_P${P}

                $5/magma/magma --gene-results genes_P${P}.genes.raw \
                --set-annot $5/data_magma/c2.all.v7.1.entrez.gmt --out set_P${P}
    	done
        rm -r ./Associations
        mkdir Associations_filtered_magma
        mv A_target* ./Associations_filtered_magma/
        tar -zcvf Associations_filtered_magma
        rm -r Associations_filtered_magma
        cd ../
}


startdir=/media/array/phenome_proj/LSEA/benchmarking/magma_multiple_pathway_2.1_2
mkdir $pathway
cd $pathway
for snp_number in 5 10 20; do
	for percentage in 0.2 0.4 0.6 0.8; do
		for h2 in 0.1 0.3 0.45; do
			magma_starter $snp_number $percentage $h2 $pathway $startdir &
			while [ $( ps -AF | grep 'magma' | wc -l ) -ge 40 ]; do sleep 1; done
		done
	done
wait
done

