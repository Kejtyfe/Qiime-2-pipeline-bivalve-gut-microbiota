# Qiita output and Qiime 2 pipeline for Bivalve gut microbiota
This repository contains my first analysis workflow using Qiita (https://qiita.ucsd.edu/) and Qiime 2 (https://qiime2.org/).
It processes and analyses 16s rRNA data from my raw reads for diversity analysis and taxonomy

Notes 🧙‍♂️ : 
-Taxonomy was assigned using the SILVA 138 classifier
-Alpha/beta diversity used filtered OTU tables
-Sample metadata is stored in metadata.tsv

## 📁 Folder structure

Qiita workflow summary

├──1. Raw Data Upload

├──2.Split libraries

├──3. Demultiplexed

├──4.Pick closed-reference OTUs (QIIMEq2 1.9.1)

├──5.OTU Table Summary (Feature Table)

└──6. Output Files

**1. Raw Data Upload**
	- Platform: Illumina
**2. Split libraries**
	- Tool: Split libraries FASTQ (QIIME 1.9.1) - To demultiplex raw FASTQ reads and apply initial quality filtering.
	- Key settings:
	Removed reads with ambiguous bases (Ns)
	Kept reads with at least 75% of the original length
	Allowed up to 3 consecutive low-quality bases
	Quality threshold: Q3
	Barcode type: Not-barcoded (reads were already separated per sample)
**3. Demultiplexed**
    - Total: 1455219  
    - Max: 301  
    - Mean: 301  
    - Standard deviation: 301
**4. Pick closed-reference OTUs (QIIMEq2 1.9.1)**
    - Classifier: silva_119_taxonomy_97_7
    - Similarity threshold: 97% (only sequences matching reference OTUs at ≥97% similarity were retained)
    - Other settings:
    sortmerna_e_value: 1
    sortmerna_max_pos: 10,000
    sortmerna_coverage: 0.97
    Threads used: 5
**5.OTU Table Summary (Feature Table)**
**6. Output Files**
    - Feature table (BIOM)
    - Feature table (QZA)
    - Sortmerna picked Otus
    - support files
    - seqs.fasta.gz

Qiime 2 pipeline\

├── scripts/ # Ubuntu WSL

├── metadata.tsv # Sample metadata

└── README.md

Steps in Pipeline 🔧
1. **Import raw data**
2. **OTU picking**
3. **Filtering chloroplast/mitochondria**
4. **Assign taxonomy (SILVA)**
5. **Collapse feature table**
6. **Alpha diversity** (Shannon)
7. **Beta diversity** (UniFrac, Bray-Curtis, pcoa)
8. **Visualisations** (barplots, diversity plots, pcoa)

```bash
activation of the environment
```

`conda activate qiime2-amplicon-2024.5` 

importing the biom and creating a qiime 2 artifact - qza
```bash
qiime tools import \
  --input-path your_depository/216124_otu_table.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path your_depository/216124_otu-feature-table.qza
```

Extracts the contents of my SortMeRNA OTU results

Check what is in sortmerna_picked_otus.tgz

```bash
tar -xvzf your_depository/216124_sortmerna_picked_otus.tgz -C your_depository_test/`
```

unzip seq.fasta.gz

```bash
gunzip -c your_depository/216111_seqs.fasta.gz > your_depository/216111_seqs.fasta
```

Generating a representative sequence Fasta file

A full set of all sequences (in FASTA format)

An OTU clustering file (text file from SortMeRNA)

seq.fasta has to match with seqs_otus

This is a code using a Python package to create a rep_set.fna, which contains one representative sequence per OTU.

```bash
pip install biopython  
```

biopython package (needed to read/write FASTA files using SeqIO)
Now, create a script for Python

```bash
nano generate_rep_set.py

insert into nano

from Bio import SeqIO
```

File paths

```bash
seqs_fasta = "your_depository/216111_seqs.fasta"
otus_txt = "your_depository/sortmerna_picked_otus/seqs_otus.txt"
output_fasta = "your_depository/rep_set.fna"
```

Load sequence records

```bash
seq_records = SeqIO.to_dict(SeqIO.parse(seqs_fasta, "fasta"))
```

#Parse the OTU file and write the representative sequences
```bash
with open(otus_txt, "r") as otus_file, open(output_fasta, "w") as out_fasta:
    for line in otus_file:
        columns = line.strip().split()
        otu_id = columns[0]  # OTU ID
        sequence_ids = columns[1:]  # Sequences in this OTU
        for seq_id in sequence_ids:
            if seq_id in seq_records:
                # Use the first sequence as the representative
                rep_seq = seq_records[seq_id]
                rep_seq.id = otu_id  # Rename header to OTU ID
                rep_seq.description = ""  # Clear description
                SeqIO.write(rep_seq, out_fasta, "fasta")
                break
```

Run the Python script

```bash
python3 generate_rep_set.py
```

converting rep_set.fna into an artifact so I can work with it in qiime

```bash
qiime tools import \
  --input-path your_depository/rep_set.fna \
  --output-path your_depository/otu-rep-seqs.qza \
  --type 'FeatureData[Sequence]'
```

now assigning taxonomy using SILVA -> silva-138-99-nb-classifier.qza

```bash
qiime feature-classifier classify-sklearn \
  --i-classifier your_path/silva-138-99-nb-classifier.qza \    # representative sequences for each OTUs!
  --i-reads your_depository/otu-rep-seqs.qza \
  --o-classification your_depository/otu-taxonomy.qza
```

filtering

filtering Mitochondria, Chloroplast, Samples with less than 5 reads (The reads depend on your choice and the data)

```bash
qiime taxa filter-table \
  --i-table your_depository/216124_otu-feature-table.qza \
  --i-taxonomy your_depository/otu-taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table your_depository/otu-table-no-mitochondria-chloroplast.qza
```

less than 5 reads

```bash
qiime feature-table filter-features \
  --i-table your_depository/otu-table-no-mitochondria-chloroplast.qza \
  --p-min-frequency 5 \
  --o-filtered-table your_depository/otu-table-filtered.qza
```

 Visualisation of the taxonomy
 
```bash
qiime metadata tabulate \
  --m-input-file your_depository/otu-table-filtered.qza \
  --o-visualization your_depository/otu-taxonomy.qzv
```

summary

```bash
qiime feature-table summarize \
  --i-table your_depository/otu-table-filtered.qza \
  --o-visualization your_depository/otu-table-summary.qzv
```

still checks after the filtering and converts to tsv, it is just for checking!
less than 5 reads

```bash
qiime tools export \
  --input-path your_depository/otu-table-filtered.qza \
  --output-path your_depository/otu-table-filtered-exported

biom convert \
  -i your_depository_test/otu-table-filtered-exported/feature-table.biom \
  -o your_depository_test/otu-table-filtered.tsv \
 --to-tsv
```


Extract OTU table (from qza to Bioam)
```bash
qiime tools export \
  --input-path your_depository/otu-table-filtered.qza \
  --output-path your_depository/exported_taxonomy
```
  
Now I need to assign taxonomy

```bash
qiime tools export \
  --input-path your_depository/otu-taxonomy.qza \
  --output-path your_depository/exported_taxonomy
```
Adding taxonomy to each OTU in the biom table

```bash
biom add-metadata \
  -i your_depository_test/exported_otu_genus_table/feature-table.biom \
  -o your_depository_test/feature-table-with-taxonomy.biom \
  --observation-metadata-fp your_depository_test/exported_taxonomy/taxonomy.tsv \
  --sc-separated taxonomy
```

The header must follow the rules

```bash
sed -i '1s/.*/#OTUID\ttaxonomy\tconfidence/' your_depository_test/exported_taxonomy/taxonomy.tsv
```

Converting the Biom table into TSV

```bash
biom convert \
  -i your_depository_test/feature-table-with-taxonomy.biom \
  -o your_depository_test/feature-table-with-taxonomy.tsv \
  --to-tsv \
  --header-key taxonomy
```

Collapsing - Even though I’ve added taxonomy into the feature table file, it’s still better to let QIIME 2 handle the taxonomy-based grouping using Qiime taxa collapse, because it understands the taxonomic hierarchy and does the grouping properly; Genus-level insights are more meaningful than raw OTUs.

```bash
qiime taxa collapse \
  --i-table your_depository_test/otu-table-filtered.qza \
  --i-taxonomy your_depository_test/otu-taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table your_depository_test/otu-genus-feature-table.qza
```

Export to table

```bash
qiime tools export \
  --input-path your_depository_test/otu-genus-feature-table.qza \
  --output-path your_depository_test/exported_otu_genus_table
```
  
And to tsv

```bash
biom convert \
  -i your_depository_test/exported_otu_genus_table/feature-table.biom \
  -o your_depository_test/otu-genus-feature-table.tsv \
  --to-tsv \
  --table-type="OTU table"
```

Blast and long formatting of the table is done in R Studio

# Alpha diversity

Get rid of the prefix

```bash
sed -i 's/15931\.//g' your_depository_test/otu-table-filtered.tsv
```

Now the qza needs to be fixed by the biom

```bash
biom convert \
  -i your_depository_test/otu-table-filtered.tsv \
  -o your_depository_test/otu-table-filtered.biom \
  --to-hdf5 \
  --table-type="OTU table"
```
  
and back to qza

```bash
qiime tools import \
  --input-path your_depository_test/otu-table-filtered.biom \
  --type 'FeatureTable[Frequency]' \
  --output-path your_depository/otu-table-filtered.qza
```

```bash
qiime diversity alpha \
  --i-table your_depository_test/otu-table-filtered.qza \
  --p-metric shannon \
  --o-alpha-diversity your_depository_test/shannon_alpha_diversity.qza
```

Visualisation

```bash
qiime metadata tabulate \
  --m-input-file your_depository_test/shannon_alpha_diversity.qza \
  --o-visualization your_depository_test/shannon_alpha_diversity.qzv
```

Alpha diversity with metadata containing duration_exposure (or whatever you need)
My IDs have a prefix from Qiita 
It is easier to add a prefix to the metadata

Adding the prefix to metadata

```bash
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 {$1=$1} NR>1 {$1="15931."$1} 1' \
```
Creating TSV

```bash
your_depository/metadata_withcombined.tsv \
```

```bash
> your_depository/metadata_prefixed.tsv
```

BIOM to QZA

```bash
qiime tools import \
  --input-path your_depository_test/otu-table-filtred.biom \
  --type 'FeatureTable[Frequency]' \
  --output-path your_depository_test/otu-table-filtered-fixed.qza
```

```bash
qiime diversity alpha \
  --i-table your_depository_test/otu-table-filtered.qza \
  --p-metric shannon \
  --o-alpha-diversity your_depository_test/shannon_alpha_diversity.qza
```
Visualisation

```bash
qiime diversity alpha-group-significance \
  --i-alpha-diversity your_depository/shannon_alpha_diversity.qza \
  --m-metadata-file your_depository/metadata_prefixed.tsv \
  --o-visualization your_depository/shannon-group-significance.qzv
```


# Beta diversity

Aligns OTU sequences, removes noisy alignment regions, builds a tree, and roots it — use it in QIIME to calculate phylogenetic diversity metrics (UniFrac)

```bash
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences your_depository_test/otu-rep-seqs.qza \
  --o-alignment your_depository_test/aligned-rep-seqs.qza \
  --o-masked-alignment your_depository_test/masked-alignment.qza \
  --o-tree your_depository_test/unrooted-tree.qza \
  --o-rooted-tree your_depository_test/rooted-tree.qza
```

 Unweighted Unifrac -> I want presence and absence information
 
```bash
qiime diversity beta-phylogenetic \
  --i-table your_depository/otu-table-filtered.qza \
  --i-phylogeny your_depository/rooted-tree.qza \
  --p-metric unweighted_unifrac \
  --o-distance-matrix your_depository/unweighted_unifrac_distance_matrix.qza
```

Measuring how different microbial communities are between samples

Based on the presence or absence of taxa, and how evolutionarily distinct those taxa are

```bash
qiime diversity beta-phylogenetic \
  --i-table your_depository_test/otu-table-filtered.qza \
  --i-phylogeny your_depository_test/rooted-tree.qza \
  --p-metric unweighted_unifrac \
  --o-distance-matrix your_depository_test/unweighted_unifrac_distance_matrix.qza
```

Extract sample IDs from the OTU table

```bash
head -n 1 your_depository_test/otu-table-filtered.tsv | cut -f2- > your_depository_test/sample-ids.txt
```

Reformat the sample ID

grep to filter another file to keep only lines that match those sample names

```bash
(head -n 1 your_depository/metadata.tsv && \
grep -F -w -f your_depository_test/sample-ids-fixed.txt \
your_depository/metadata.tsv) \
> your_depository_test/metadata-matched.tsv
```


needs to extract the real headers

```bash
sed -n '2p' your_depository/otu-table-filtered.tsv | cut -f2- \
| tr '\t' '\n' > your_depository_test/sample-ids-fixed.txt
```

and check it

```bash
cat your_depository_test/sample-ids-fixed.txt
```

Now convert to biom

```bash
biom convert \
  -i your_depository_test/otu-table-filtered.tsv \
  -o your_depository_test/otu-table-fixed.biom \
  --table-type="OTU table" \
  --to-hdf5
```

and convert biom to qza so I can work with it in qiime

```bash
qiime tools import \
  --input-path your_depository_test/otu-table-fixed.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path your_depository_test/otu-table-fixed.qza
```

# Visualisation

```bash
qiime taxa barplot \
  --i-table your_depository/otu-table-fixed.qza \
  --i-taxonomy your_depository/otu-taxonomy.qza \
  --m-metadata-file your_depository/metadata.tsv \
  --o-visualization your_depository/taxa-barplot.qzv
```

# PCA

```bash
qiime diversity beta \
  --i-table your_depository/otu-table-fixed.qza \
  --p-metric braycurtis \
  --o-distance-matrix your_depository/braycurtis-distance.qza
```

# braycurtis-distance

```bash
qiime diversity pcoa \
  --i-distance-matrix your_depository_test/braycurtis-distance.qza \
  --o-pcoa your_depository_test/braycurtis-pcoa.qza
```
 
# visualisation

```bash
qiime emperor plot \
  --i-pcoa your_depository/braycurtis-pcoa.qza \
  --m-metadata-file your_depository/metadata.tsv \
  --o-visualization your_depository/braycurtis-pcoa.qzv
```

💻😃🌼


Metadata
sample_name	description	sample_type	condition	sex	day	feeding	exposure	condition_duration
15931.GD7A1ChlControl7	mussel	experimental	Chl Control	f	day7	chlorella	no	Chl Control 7
15931.GD7A2ChlControl7	mussel	experimental	Chl Control	m	day7	chlorella	no	Chl Control 7
15931.GD7A3ChlControl7	mussel	experimental	Chl Control	m	day7	chlorella	no	Chl Control 7

٩(^ᗜ^)و Happy coding
