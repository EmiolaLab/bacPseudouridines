# Profiling of pseudouridines in microbiome-derived bacterial RNA
Pseudouridine (Ψ) is a widespread RNA modification that influences RNA stability, structure, and translation. However, its role in bacterial mRNA, particularly within complex microbiomes, remains poorly defined. We describe a computational pipeline for base-resolution, quantitative mapping of pseudouridine in microbiome transcriptomes. This pipeline is applicable to bisulfite-based, microbial sequencing workflow.  More information can be found in: 

Sharma et al (2026) "Quantitative mapping of pseudouridines in bacterial RNA" *Nature Communications* 

For a detailed protocol, please refer to Sharma et al (2026) "High-throughput profiling of pseudouridines in microbiome-derived bacterial RNA" *Current Protocols*

### Installing dependencies
`clone the github page (git clone https://github.com/EmiolaLab/bacPseudouridines.git)`\
`cd bacPseudouridines/`\
`conda env create -f environment.yml`\
`conda activate Pseudouridine_env`\
`wget https://raw.githubusercontent.com/sridhar0605/brc-parser/refs/heads/main/brc-parser.py`\
`wget https://github.com/mozack/abra2/releases/download/v2.23/abra2-2.23.jar`\
`chmod +x brc-parser.py`

### Example run
This example involves a very small dataset that has been optimized for a quick run. Note that running real samples can be memory intensive. Therefore, we recommend allocating a large memory (> 100G) for real-world microbiome samples.

`cd example`

*microbes.fa* is a FASTA file containing genome sequences of selected microbial taxa.

`bwa index microbes.fa`

#### 1. Generate genome annotation file
We use Prokka to generate genome annotations, including strand information.

The commands below generate annotations for the concatenated genomes. Ideally, each genome should be annotated individually, and the resulting GFF files should then be merged. However, for the purposes of this example, we will generate annotations directly from the concatenated genomes.

`prokka --quiet --kingdom Bacteria --outdir out.dir --locustag genomes --prefix genomes microbes.fa`\
`sed '/##FASTA/,$d' out.dir/genomes.gff | grep -v "##" > all_gffs.txt`

#### 2. Read alignment to reference database

The example below is for BS-treated samples. For untreated samples, replace  `treated_subset.R1.fastq.gz`  with  `untreated_subset.R1.fastq.gz`

`bwa mem -t 16 microbes.fa treated_subset.R1.fastq.gz treated_subset.R2.fastq.gz > sample.sam`\
`samtools view -@ 16 -b sample.sam > sample.bam`\
`samtools sort -@ 16 sample.bam > sample.sorted.bam`\
`samtools index sample.sorted.bam`\
`rm sample.sam sample.bam`

#### 3. Conversion, sorting, and indexing of alignments
Convert the SAM file to BAM format, then sort and index using samtools: 

`samtools view -@ 16 -b sample.sam > sample.bam`\
`samtools sort -@ 16 sample.bam > sample.sorted.bam`\
`samtools index sample.sorted.bam`

Optionally, remove intermediate files to conserve disk space:\
`rm sample.sam sample.bam`

#### 4. Local realignment of reads
Perform local realignment to improve detection of insertion and deletion events using ABRA2: 

`mkdir TMP`\
`java -Xmx32G -jar abra2-2.23.jar --in sample.sorted.bam --out sample.sorted.realign.bam --ref microbes.fa --threads 16 --tmpdir TMP --sa > sample.abra.log`

The “--sa" (skip assembly) flag reduces runtime by bypassing assembly-based realignment. 
Optionally, remove intermediate folder\
`rm -rf TMP`

#### 5. Retrieval of nucleotide coverage and variant information
Index the realigned BAM file\
`samtools index sample.sorted.realign.bam`

Generate per-base coverage and nucleotide variant information using bam-readcount:\
`bam-readcount -w1 -f microbes.fa sample.sorted.realign.bam > sample.brc.tsv`

This step produces a tab-delimited file containing per-base read counts and nucleotide composition across the reference genomes.
Parse the output into a structured format using brc-parser.py:\
`python brc-parser.py sample.brc.tsv`

Convert the resulting CSV file to a tab-delimited format for downstream analysis:\
`sed 's/,/\t/g' sample.brc_parsed.csv > sample_treated.txt`

#### 6. Post-processing filteration
Here, we are retrieving only 'A' and 'T' sites with depth >= 20.

For BS-treated samples, run:\
`awk '$3~/^(T|A)$/ {print $0}' sample_treated.txt | awk '$4 >= 20' | grep -w "del" > sample_treated_filtered.txt`

For untreated samples, run:\
`awk '$3~/^(T|A)$/ {print $0}' sample_untreated.txt | awk '$4 >= 20' > sample_untreated_filtered.txt`

#### 7. Retrieve final output
We take both the BS-treated and untreated outputs for comparison to identify putative pseudouridine sites.\
`Rscript PseudoU.R -u sample_untreated_filtered.txt -t sample_treated_filtered.txt -g all_gffs.txt -p 0.01 -o results.txt`

Here, -u untreated output file\
-t BS-treated output file\
-g genome annotation file (generated in the 1st step)\
-p P-value threshold by Fisher’s exact test, comparing deletion counts and total read depth between treated and untreated samples \
-o final output file


