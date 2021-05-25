# Benchmarking
The process of benchmarking implies several steps, including genotype and phenotype data simulation, GWAS association
and lastly gene set analysis with LSEA or MAGMA.

### 1. Haplotype simulation with HapGen2 based on 1000GP_Phase3 genotype data

```{r, engine=bash}
hapdir=~/1000GP_Phase3
for chr in `seq 1 22`; do
	gunzip $hapdir/1000GP_Phase3_chr${chr}.legend.gz
	gunzip $hapdir/1000GP_Phase3_chr${chr}.hap.gz
	dummyDL=`sed -n '2'p $hapdir/1000GP_Phase3_chr${chr}.legend | cut -d ' ' -f 2`
	hapgen2 -m $hapdir/genetic_map_chr${chr}_combined_b37.txt \
        -l $hapdir/1000GP_Phase3_chr${chr}.legend \
        -h $hapdir/1000GP_Phase3_chr${chr}.hap -o ~/hapgen_results/genotypes_chr${chr}_hapgen \
        -dl $dummyDL 0 0 0 -n 1000 0 -no_haps_output 
done
```

### 2. Conversion of Oxford to PLINK format

```{r, engine=bash}
for chr in `seq 1 22`; do
        plink --data genotypes_chr${chr}_hapgen.controls \
        --oxford-single-chr $chr \
        --make-bed \
        --out genotypes_chr${chr}_hapgen.controls
        
	echo -e "genotypes_chr${chr}_hapgen.controls" >> file_list
done

plink --merge-list file_list --maf 0.05 --make-bed --out genotypes_genome_hapgen.MAF_higher_0.05
```

### 3. Variant selection for causal SNP and genetic background simulation

```{r, engine=bash}
PATHWAY=${OPTARG}
NUMBER=${OPTARG}
H2=${OPTARG}
PERCENTAGE=${OPTARG}


gene_list=$(cat /media/array/phenome_proj/LSEA/data/c2.all.v7.0.symbols.gmt | grep "$PATHWAY" | cut -f3- )
for i in $gene_list; do
	SNPS+="$(grep -w "$i" snp_gene_list_matched | awk -v FS=' ' '{print $1}') "
done


genomedir=/media/array/phenome_proj/LSEA/benchmarking/dalexeev/pathway
currentdir=$(pwd)

let "shuf_number = 1000 - $NUMBER"

snp_number_target_causal=$(echo "scale=0; $NUMBER * $PERCENTAGE" | bc)
snp_number_random_causal=$(echo "scale=0; $NUMBER - $snp_number_target_causal" | bc)

cat $SNPS | tr " " "\n" | shuf -n $snp_number_target_causal > target_snp_list_${NUMBER}_${PERCENTAGE}_${H2}
comm -23 <(cut -f2 $genomedir/genotypes_genome_hapgen.MAF_higher_0.05.bim | sort) <(cat $SNPS | tr " " "\n" | sort) | \
shuf -n $snp_number_random_causal >> target_snp_list_${NUMBER}_${PERCENTAGE}_${H2}

/media/array/phenome_proj/LSEA/benchmarking/dalexeev/plink/plink --bfile $genomedir/genotypes_genome_hapgen.MAF_higher_0.05\
--extract target_snp_list_${NUMBER}_${PERCENTAGE}_${H2} --make-bed --out genotypes_hapgen.controls.target

cut -f2 $genomedir/genotypes_genome_hapgen.MAF_higher_0.05.bim | shuf -n $shuf_number > snps.subset.map
cat target_snp_list_${NUMBER}_${PERCENTAGE}_${H2} >> snps.subset.map

/media/array/phenome_proj/LSEA/benchmarking/dalexeev/plink/plink --bfile $genomedir/genotypes_genome_hapgen.MAF_higher_0.05 --extract snps.subset.map \
--make-bed --out genotypes_subset_hapgen.controls

```

### 4. Phenotype construction via PhenotypeSimulator (look into Phenotype_script.R)
```{r, engine=bash}
Rscript --vanilla --verbose $STARTDIR/Phenotype_script.R $currentdir $NUMBER $H2
```

### 5. GWAS association via PLINK

```{r, engine=bash}
cat Y_caspase | awk 'BEGIN { FS="\t"; OFS="\t" } { $1=$1 "\t" $1 } 1' > Y_caspase1 
cat Y_caspase1 | awk '{FS="\t";OFS="\t"} {sub("id1_","id2_",$2)}1' > Y_caspase2
cat Y_caspase2 | tr -d '"' > Y_caspase3

/media/array/phenome_proj/LSEA/benchmarking/dalexeev/plink/plink --bfile $genomedir/genotypes_genome_hapgen.MAF_higher_0.05 --maf 0.05 --hwe 1e-10 --allow-no-sex --pheno Y_target3 --all-pheno --assoc --out A_target
```
### 6. LSEA gene set analysis

```{r, engine=bash}
awk '{F="\t";OFS="\t"} NR>1{split($2,a,":");print $1,$2,$3,a[3],a[4],$4,$5,$6,$7,$8,$9}' \
A_target.qassoc | sed -e '1i\CHR\tVARIANT\tBP\tA1\tA2\tNMISS\tBETA\tSE\tR2\tT\tP' > A_target.qassoc.filtered_ldsc

python3 /media/array/phenome_proj/LSEA/benchmarking/uk_bio/LSEA/LSEA_2.1.py \
	-input $savename \
	-plink_dir /media/array/phenome_proj/LSEA/benchmarking/plink \
	-bfile /media/array/phenome_proj/EUR_005_nodups \
	-ldsc_dir /media/array/phenome_proj/LSEA/benchmarking/ldsc -n $N \
	-universe /media/array/phenome_proj/LSEA/benchmarking/uk_bio/LSEA/data/universe.json \
	-column_names chr bp variant pval -qval_threshold 1 -out result

```
### 7. Results
You can find results on gene set enrichment analysis conducted by LSEA and MAGMA for different sets of GWAS data. This data differs by number of causal variants (NUMBER), heritability (H2), source of causal variants (PATHWAY) and percentage of causal variants sampled from the taget gene set.
