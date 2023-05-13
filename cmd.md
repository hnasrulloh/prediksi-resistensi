iqtree -v -s core_gene_alignment_filtered.aln -pre core_tree -nt 8 -fast -m GTR

panfeed -v -g ../../gff_ethambutol -p gene_presence_absence.csv -o panfeed1 --upstream 100 --downstream 100 --compress --cores 3

pyseer --lmm \ 
--phenotypes ../../traits_isoniazid.tsv \
--pres panfeed/hashes_to_patterns.tsv \
--similarity pyseer/phylogeny_K.tsv \
--output-patterns pyseer/kmer_patterns.txt \
--cpu 6 > pyseer/pyseer_kmer.tsv

ls gff_isoniazid/ | sed 's/.gff//g' > gff_rifampin_panfeed2.txt 

panfeed -v \
 -g ../../gff_isoniazid \  
-p gene_presence_absence.csv \ 
-o panfeed2 \
--targets ../../gff_isoniazid_panfeed2.txt \
--genes pyseer/gene_clusters.txt \
--upstream 100 \
--downstream 100 \
--compress \
--cores 4

zcat kmers.tsv.gz > kmers.tsv


