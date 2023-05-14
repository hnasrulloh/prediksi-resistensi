```
iqtree -v -s core_gene_alignment_filtered.aln -pre core_tree -nt 8 -fast -m GTR
```


```
panfeed -v -g ../../gff_antibitic -p gene_presence_absence.csv -o panfeed1 --upstream 100 --downstream 100 --compress --cores 3
```


```
../../../scripts/phylogeny_distance.py --lmm core_tree.treefile > pyseer/phylogeny_K.tsv
```

Panfeed1
```
zcat hashes_to_patterns.tsv.gz > hashes_to_patterns.tsv
zcat kmers_to_hashes.tsv.gz > kmers_to_hashes.tsv
zcat kmers.tsv.gz > kmers.tsv
```


```sh
pyseer --lmm \
--phenotypes ../../traits_antibiotic.tsv \
--pres panfeed1/hashes_to_patterns.tsv \
--similarity pyseer/phylogeny_K.tsv \
--output-patterns pyseer/kmer_patterns.txt \
--cpu 4 > pyseer/pyseer.tsv
```


```sh
ls gff_antibiotic | sed 's/.gff//g' > gff_antibiotic_panfeed2.txt 
```

```
panfeed-get-clusters -t 1E-10 -a pyseer/pyseer.tsv -p panfeed1/kmers_to_hashes.tsv  > pyseer/gene_clusters.txt
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
Manually remove a broken last line
```bash
zcat hashes_to_patterns.tsv.gz > hashes_to_patterns.tsv
zcat kmers_to_hashes.tsv.gz > kmers_to_hashes.tsv
zcat kmers.tsv.gz > kmers.tsv
```
