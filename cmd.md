
```sh
iqtree -v -s core_gene_alignment_filtered.aln -pre core_tree -nt 8 -fast -m GTR
```


```sh
panfeed -v -g ../../gff_antibiotic -p gene_presence_absence.csv -o panfeed1 --upstream 100 --downstream 100 --compress --cores 3
```


```sh
../../../scripts/phylogeny_distance.py --lmm core_tree.treefile > pyseer/phylogeny_K.tsv
```

Panfeed1
```sh
zcat hashes_to_patterns.tsv.gz > hashes_to_patterns.tsv
zcat kmers_to_hashes.tsv.gz > kmers_to_hashes.tsv
zcat kmers.tsv.gz > kmers.tsv
```


```sh
pyseer --lmm \
--phenotypes ../../traits_antibiotic.tsv \
--pres panfeed1/hashes_to_patterns.tsv \
--similarity pyseer/phylogeny_K.tsv \
--output-patterns pyseer/hashed_patterns.txt \
--cpu 4 > pyseer/pyseer.tsv
```


```sh
ls gff_antibiotic | sed 's/.gff//g' > gff_antibiotic_panfeed2.txt 
```

```sh
panfeed-get-clusters -t THRES -a pyseer/pyseer.tsv -p panfeed1/kmers_to_hashes.tsv  > pyseer/gene_clusters.txt
```


```sh
panfeed -v \
-g ../../gff_antibiotic \
-p gene_presence_absence.csv \
-o panfeed2 \
--targets ../../gff_antibiotic_panfeed2.txt \
--genes pyseer/gene_clusters.txt \
--upstream 100 \
--downstream 100 \
--compress \
--cores 4
```

Panfeed2
```bash
zcat hashes_to_patterns.tsv.gz > hashes_to_patterns.tsv
zcat kmers_to_hashes.tsv.gz > kmers_to_hashes.tsv
zcat kmers.tsv.gz > kmers.tsv
sed -i '$ d' kmers.tsv
```

```sh
panfeed-get-kmers -v -t THRES -a pyseer/pyseer.tsv -p panfeed2/kmers_to_hashes.tsv  -k panfeed2/kmers.tsv  > pyseer/kmer_hits.tsv
```
