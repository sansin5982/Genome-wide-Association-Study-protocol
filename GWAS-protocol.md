-   <a href="#introduction" id="toc-introduction">Introduction</a>
-   <a href="#material-and-methods" id="toc-material-and-methods">Material
    and Methods</a>
    -   <a href="#basic-input-files" id="toc-basic-input-files">Basic input
        files</a>
        -   <a href="#ped-file" id="toc-ped-file">PED File</a>
        -   <a href="#map-file" id="toc-map-file">MAP file</a>
    -   <a href="#binary-ped-files" id="toc-binary-ped-files">Binary PED
        files</a>
        -   <a href="#bed" id="toc-bed">BED</a>
        -   <a href="#bim" id="toc-bim">BIM</a>
        -   <a href="#fam" id="toc-fam">FAM</a>
    -   <a
        href="#steps-to-perform-quality-control-in-a-genome-wide-association-study"
        id="toc-steps-to-perform-quality-control-in-a-genome-wide-association-study">Steps
        to perform quality control in a Genome-wide association study:</a>
        -   <a href="#per-individual-quality-control"
            id="toc-per-individual-quality-control">Per Individual Quality
            Control</a>
        -   <a href="#step-1-converting-ped-and-map-into-binary-format"
            id="toc-step-1-converting-ped-and-map-into-binary-format">Step 1:
            Converting ped and map into binary format</a>
        -   <a
            href="#step-2-identification-of-individuals-with-discordant-sex-information"
            id="toc-step-2-identification-of-individuals-with-discordant-sex-information">Step
            2: Identification of individuals with discordant sex information</a>
        -   <a
            href="#step-3-identification-of-individuals-with-elevated-missing-data-rates"
            id="toc-step-3-identification-of-individuals-with-elevated-missing-data-rates">Step
            3: Identification of individuals with elevated missing data rates</a>
        -   <a
            href="#step-4-identification-of-individuals-with-outlying-heterozygosity-rate"
            id="toc-step-4-identification-of-individuals-with-outlying-heterozygosity-rate">Step
            4: Identification of individuals with outlying heterozygosity rate</a>
        -   <a href="#step-5-identification-of-duplicate-samples"
            id="toc-step-5-identification-of-duplicate-samples">Step 5:
            Identification of duplicate samples</a>
        -   <a href="#step-6-identification-of-individuals-of-divergent-ancestry"
            id="toc-step-6-identification-of-individuals-of-divergent-ancestry">Step
            6: Identification of Individuals of divergent ancestry</a>
    -   <a href="#per-snp-quality-control" id="toc-per-snp-quality-control">Per
        SNP quality control</a>
-   <a href="#references" id="toc-references">References</a>

# Introduction

A genome-wide association study (GWAS) is a type of study that uses
genomic data to identify genetic variations that are associated with a
particular trait or disease (Marees et al, 2018).The National Institute
of Health (NIH) defines GWAS as a study of common genetic variations
across the entire genome designed to identify genetic associations with
observable traits. A GWAS uses high-throughput genotyping technology. It
examines thousands of polymorphisms to relate them to a clinical
condition or measurable trait.It is a hypothesis-free approach to scan
markers across the whole genome.

It involves genotyping a large number of single nucleotide polymorphisms
(SNPs) across the entire genome of many individuals and comparing the
frequency of these SNPs in individuals with the trait or disease of
interest compared to those without the trait or disease.A traditional
**case-control study design** is the most common approach in GWASs. The
goal of a GWAS is to identify specific genetic loci that are associated
with the trait or disease, which can provide insight into the underlying
biological mechanisms and potentially lead to the development of new
diagnostic tools and treatments.

GWAS has been used to identify genetic risk factors for a wide range of
complex diseases and traits, including cardiovascular disease, type 2
diabetes, cancer, and many others. By utilizing large sample sizes and
high-throughput genotyping technologies, GWAS has revolutionized our
understanding of the genetic basis of complex traits and has become a
powerful tool in the field of genetics and genomics.

# Material and Methods

**Tools, Software and Files Required**

-   **plink**
    -   [link1](https://www.cog-genomics.org/plink2/)
    -   [link2](https://www.cog-genomics.org/plink/2.0/)
-   **bcftools and samtools**
    -   [link](http://www.htslib.org/download/)
-   **vcftools**
    -   [link](https://sourceforge.net/projects/vcftools/files/)
-   **Check for compatibility with the input required by the Sanger
    server.**
    -   [link](http://qbrc.swmed.edu/zhanxw/software/checkVCF/checkVCF-20140116.tar.gz)
-   **Haplotype Reference Consortium v1.1 panel for HRC site list**
    -   [link](HRC.r1-1.GRCh37.wgs.mac5.sites.tab)[link](ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/)
-   **HRC preparation checking tool**
    -   [link](https://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.2.9.zipa)
-   Set the path for these tools so you directly call them from a
    different folder. To set the path first you have enter in that
    specific tools

**export PATH=*P**A**T**H*:(pwd)**

## Basic input files

We expect two basic input files, Raw.ped and Raw.map, for association
analysis. Here, raw is the file name. PED (Plink Pedigree File) and MAP
(Plink Map File) are file formats commonly used in genome-wide
association studies (GWAS). These files can also be converted into
binary format.

### PED File

The PED file is a tab-delimited text file that contains information
about individuals in a study, including their family relationships and
genotype data. Each row in the PED file represents a single individual,
and each column represents a different type of information, such as the
individual’s family ID, individual ID, and genotype data. The first six
columns in a PED file are required and include the following
information:

<table>
<thead>
<tr class="header">
<th>FID</th>
<th>IID</th>
<th>PID</th>
<th>MID</th>
<th>Sex</th>
<th>P</th>
<th>rs1</th>
<th>rs2</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>1</td>
<td>1</td>
<td>0</td>
<td>0</td>
<td>1</td>
<td>1</td>
<td>CT</td>
<td>AG</td>
</tr>
<tr class="even">
<td>2</td>
<td>2</td>
<td>0</td>
<td>0</td>
<td>0</td>
<td>1</td>
<td>CG</td>
<td>AA</td>
</tr>
<tr class="odd">
<td>3</td>
<td>3</td>
<td>0</td>
<td>0</td>
<td>1</td>
<td>0</td>
<td>CC</td>
<td>TT</td>
</tr>
</tbody>
</table>

-   FID = Family ID,a unique identifier for the family to which the
    individual belongs
-   IID = Individual ID, a unique identifier for the individual within
    the family
-   PID = Paternal ID,the individual ID of the father (-9 if unknown)
-   MID = Maternal ID, the individual ID of the mother (-9 if unknown)
-   Sex (1 = male; 2 = female; other = -9)
-   P = Phenotype, the phenotype of the individual, usually coded as a
    binary trait (e.g. affected=1, unaffected=2, missing=0)

The remaining columns in the PED file contain genotype data for each
individual.

-   rs1 = SNP1
-   rs2 = SNP2

Sixth column represents phenotype. The phenotype can be either an
affection status or quantitative trait. Plink can automatically detect
the phenotype based on the code.

### MAP file

Each line of the MAP file explains a single marker and contains four
columns. The MAP file is a tab-delimited text file that contains
information about each SNP (single nucleotide polymorphism) that was
genotyped in the study. Each row in the MAP file represents a single
SNP, and each column represents a different type of information, such as
the chromosome number, SNP identifier, and base-pair position. The
columns in a MAP file include the following information:

<table>
<thead>
<tr class="header">
<th>Chr</th>
<th>SNP</th>
<th>GD</th>
<th>BPP</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>1</td>
<td>rs1</td>
<td>0</td>
<td>870050</td>
</tr>
<tr class="even">
<td>1</td>
<td>rs2</td>
<td>0</td>
<td>870150</td>
</tr>
<tr class="odd">
<td>1</td>
<td>rs3</td>
<td>0</td>
<td>870322</td>
</tr>
<tr class="even">
<td>1</td>
<td>rs4</td>
<td>0</td>
<td>870878</td>
</tr>
</tbody>
</table>

-   Chr = Chromosome (1-22, X, Y or 0 if unplaced), the chromosome
    number on which the SNP is located
-   rs = a unique identifier for the SNP
-   GD = the distance between this SNP and the previous SNP (in Morgans)
-   BPP = the position of the SNP on the chromosome (in base-pairs)

The autosomes should be coded 1 through 22. The following other codes
can be used to specify other chromosome types: <br>

<table>
<thead>
<tr class="header">
<th>Region</th>
<th>Genomic region</th>
<th>Number</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>X</td>
<td>X chromosome</td>
<td>23</td>
</tr>
<tr class="even">
<td>Y</td>
<td>Y chromosome</td>
<td>24</td>
</tr>
<tr class="odd">
<td>XY</td>
<td>Pseudo-autosomal region of X</td>
<td>25</td>
</tr>
<tr class="even">
<td>MT</td>
<td>Mitochondrial</td>
<td>26</td>
</tr>
</tbody>
</table>

Together, the PED and MAP files provide all the information needed to
perform a GWAS. The genotype data in the PED file can be combined with
the information about each SNP in the MAP file to test for associations
between genotype and phenotype. This is typically done using a
statistical software package, such as PLINK.

## Binary PED files

Binary file formats are an alternative to the text-based PED/MAP file
formats commonly used in GWAS. Binary formats offer several advantages
over text-based file formats, including faster loading times, reduced
storage space requirements, and improved data compression. Binary files
contain extended map file (a.bim) having information about allele names,
binary ped file (a\*bed) and pedigree/phenotype information in separate
file (a.fam).

### BED

The BED file contains information on the presence or absence of a
specific allele for each marker at each sample, represented as binary
values (0 or 1).

### BIM

A BIM file is a file format used to store and organize data on genetic
variants that are included in a study. The BIM file contains information
on the chromosome, variant ID, position, and allele information for each
genetic variant.

The BIM file is typically a tab-delimited text file with six columns:
chromosome number, variant ID, genetic position, base-pair position, the
reference allele, and the alternate allele. The chromosome number column
specifies which chromosome the genetic variant is located on, while the
variant ID column provides a unique identifier for the variant. The
genetic position and base-pair position columns provide information on
the location of the variant on the chromosome. The reference allele and
alternate allele columns provide information on the two possible alleles
for the variant.

GWAS typically involve the analysis of large numbers of genetic
variants, and the BIM file helps to ensure that data is organized and
analyzed correctly. The BIM file is used in conjunction with other
files, such as genotype files (e.g., PLINK binary files), to identify
associations between genetic variants and traits or diseases of
interest.

<table>
<thead>
<tr class="header">
<th>Chr</th>
<th>SNP</th>
<th>GD</th>
<th>BPP</th>
<th>Allele1</th>
<th>Allele2</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>1</td>
<td>rs1</td>
<td>0</td>
<td>870050</td>
<td>C</td>
<td>T</td>
</tr>
<tr class="even">
<td>1</td>
<td>rs2</td>
<td>0</td>
<td>870150</td>
<td>A</td>
<td>G</td>
</tr>
<tr class="odd">
<td>1</td>
<td>rs3</td>
<td>0</td>
<td>870322</td>
<td>G</td>
<td>G</td>
</tr>
<tr class="even">
<td>1</td>
<td>rs4</td>
<td>0</td>
<td>870878</td>
<td>C</td>
<td>A</td>
</tr>
</tbody>
</table>

### FAM

In a GWAS, a FAM file is a text file that contains information about the
individuals being studied. The FAM file is used to specify the familial
relationships and phenotype information for each individual in the
study.

The FAM file typically contains three columns of data separated by
whitespace or tab characters. The first column contains a family ID,
which is a unique identifier for a particular family or group of related
individuals. The second column contains an individual ID, which is a
unique identifier for each individual within a family. The third column
contains the phenotype information for each individual, such as disease
status or a quantitative trait.

For example, a FAM file might look like this:

<table>
<thead>
<tr class="header">
<th>familyID</th>
<th>individualID</th>
<th>fatherID</th>
<th>motherID</th>
<th>sex</th>
<th>phenotype</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>1</td>
<td>1</td>
<td>0</td>
<td>0</td>
<td>1</td>
<td>1</td>
</tr>
<tr class="even">
<td>1</td>
<td>2</td>
<td>0</td>
<td>0</td>
<td>2</td>
<td>0</td>
</tr>
<tr class="odd">
<td>1</td>
<td>3</td>
<td>1</td>
<td>2</td>
<td>1</td>
<td>1</td>
</tr>
<tr class="even">
<td>1</td>
<td>4</td>
<td>1</td>
<td>2</td>
<td>2</td>
<td>0</td>
</tr>
</tbody>
</table>

In this example, there is one family with four individuals. The first
column specifies the family ID, which is 1 in this case. The second
column specifies the individual ID, ranging from 1 to 4. The third and
fourth columns specify the IDs of the father and mother for each
individual, respectively. The fifth column specifies the sex of each
individual (1 for male and 2 for female). The sixth column specifies the
phenotype, with 1 indicating affected by the disease or trait of
interest, and 0 indicating unaffected.The phenotype can be either a
quantitative trait or an affection status column: PLINK will
automatically detect which type. \# Methods

## Steps to perform quality control in a Genome-wide association study:

A GWAS use thousands of markers to test a genetic association. Even a
low rate of error or bias will introduce false positive associations. It
has been suggested that removal of handful individuals or small
percentage of markers should not markedly decrease the overall power of
the study. Removal of one marker has greater impact than removal of one
individual. Removal of one marker is potentially an overlooked disease
association. Thus, implementing QC per-individual before QC per-marker
will maximize the number of markers remaining in the study.

### Per Individual Quality Control

Per-individual QC screens genotype to identify subjects that may
introduce bias, if not removed. Poorly genotyped (low call rate)
individuals will increase error in the study and may significantly
affect the results. There are several steps of per-individual QC for a
GWAS data set.

### Step 1: Converting ped and map into binary format

PLINK command: <br> <br> ./plink –file raw\_GWAS\_data –make-bed <br>
<br>

### Step 2: Identification of individuals with discordant sex information

This was the first step of QC, performed to identify subjects that have
inconclusive/contradictory gender information, as it can lead to
spurious associations. In a case-control study, samples that show the
wrong gender information are suggested to be excluded from further QC
and statistical analysis. Subjects were appropriately recoded or
removed, if information was inconclusive, for further analyses.

<br> **PLINK command to check the sex of all individuals** <br> <br>
**./plink –bfile raw\_GWAS\_data –check-sex –out GWAS\_Sex\_Check** <br>

-   Command Create a list of individuals with discordant sex data in
    file “GWAS\_Sex\_Check.sexcheck”. Ascertained sex is indicated in
    column 3 while genotype-based sex is indicated in column 4. The
    genotype data are noted in column 4 with a 0 when the homozygosity
    rate is more than 0.2 but less than 0.8, indicating that the data
    are inconclusive about the sex of an individual.

<table>
<thead>
<tr class="header">
<th>FID</th>
<th>IID</th>
<th>PEDSEX</th>
<th>SNPSEX</th>
<th>STATUS</th>
<th>F</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>1</td>
<td>201320</td>
<td>2</td>
<td>2</td>
<td>OK</td>
<td>0.027</td>
</tr>
<tr class="even">
<td>1</td>
<td>201320</td>
<td>1</td>
<td>1</td>
<td>OK</td>
<td>0.976</td>
</tr>
<tr class="odd">
<td>1</td>
<td>201327</td>
<td>2</td>
<td>1</td>
<td>PROBLEM</td>
<td>0.974</td>
</tr>
<tr class="even">
<td>1</td>
<td>201335</td>
<td>2</td>
<td>0</td>
<td>PROBLEM</td>
<td>0.456</td>
</tr>
<tr class="odd">
<td>1</td>
<td>201342</td>
<td>1</td>
<td>2</td>
<td>PROBLEM</td>
<td>0.632</td>
</tr>
<tr class="even">
<td>1</td>
<td>201359</td>
<td>1</td>
<td>0</td>
<td>PROBLEM</td>
<td>0.321</td>
</tr>
</tbody>
</table>

-   The above command also provide a log file, here
    GWAS\_Sex\_Check.log, that provides many information including
    number of cases and controls, males and females count, individuals
    with ambiguous code, etc

-   Extract the IDs of individuals with discordant sex information. In
    situations where discrepancies cannot be resolved, remove the
    individuals through following command. <br> <br>

**PLINK command to remove the individuals based on sex information**
<br> <br> **plink –bfile raw\_GWAS\_data –remove
discordant-sex-individuals-file.txt –make-bed –out
1\_QC\_Raw\_GWAS\_data** <br> <br> (File
“discordant-sex-individuals-file.txt”, should contain only FID and IID
of the individuals that have to be removed)

<table>
<thead>
<tr class="header">
<th>FID</th>
<th>IID</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>1</td>
<td>201320</td>
</tr>
<tr class="even">
<td>1</td>
<td>201327</td>
</tr>
<tr class="odd">
<td>1</td>
<td>201335</td>
</tr>
<tr class="even">
<td>1</td>
<td>201342</td>
</tr>
<tr class="odd">
<td>1</td>
<td>201359</td>
</tr>
</tbody>
</table>

    Gender <- read.table("Sex_check_1.sexcheck", header = T, as.is = T) %>%
      na.omit()
    png("Gender_check.png")

    ggplot(Gender, aes(x=F, y= PEDSEX, col = STATUS))+
      geom_point()+
      labs(y="Gender", x = "F score")+
      theme_classic()+
      scale_x_continuous(breaks = round(seq(min(Gender$F), 1, by = 0.01), 1))+
      theme(strip.text.x = element_blank())+
      geom_circle(aes(x0 = 0.61,  y0 = 2, r = 0.42),
                  inherit.aes = FALSE,
                  col = "Red")+
      geom_circle(aes(x0 = 0.05,  y0 = 1, r = 0.46),
                  inherit.aes = FALSE,
                  col = "Red")+
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = "none")+
      annotate(geom="text", x=-0.2, y=2.1, label="Females (< 0.2)",
               color="#e75480")+
      annotate(geom="text", x=0.9, y=1.1, label="Males (> 0.8)",
               color="#00bfff")+
      annotate(geom="text", x=0.61, y=2.1, label="Females failed (> 0.2)",
               color="#e75480")+
      annotate(geom="text", x=0.05, y=1.1, label="Males failed (< 0.8)",
               color="#00bfff")
    dev.off()

The above R script will create a png file containing information based
on STATUS column and differentiate males, females and ambiguous data
separately.

<img src="Gender_check.png" alt="Discordant Sex information"  />
<p class="caption">
Discordant Sex information
</p>

### Step 3: Identification of individuals with elevated missing data rates

This QC step is used to identify and exclude individuals with too much
missing genotype data. Genotype accuracy and genotype call rate can be
significantly affected by variations in DNA quality. A high genotype
failure rate suggest poor DNA sample quality. <br> **PLINK command to
calculate missing rate.** <br> <br> **./plink2 –bfile
1\_QC\_Raw\_GWAS\_data –missing –out missing\_data\_rate** <br> <br> -
Command creates the files “missing\_data\_rate.imiss” and
“missing\_data\_rate.lmiss”. - The “imiss” file (individual missingness)
reports the proportion of missing genotypes per individual in the
dataset. It lists each individual ID and the number and proportion of
missing genotypes for that individual. This can be useful for
identifying samples with high levels of missing data that may need to be
removed from downstream analyses.

<table>
<thead>
<tr class="header">
<th>FID</th>
<th>IID</th>
<th>MISS_PHENO</th>
<th>N_MISS</th>
<th>N_GENO</th>
<th>F_MISS</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>1</td>
<td>201320</td>
<td>N</td>
<td>631</td>
<td>536323</td>
<td>0.00177</td>
</tr>
<tr class="even">
<td>1</td>
<td>201327</td>
<td>N</td>
<td>566</td>
<td>536323</td>
<td>0.00105</td>
</tr>
<tr class="odd">
<td>1</td>
<td>201419</td>
<td>N</td>
<td>94784</td>
<td>536323</td>
<td>0.1767</td>
</tr>
<tr class="even">
<td>1</td>
<td>201567</td>
<td>N</td>
<td>16860</td>
<td>538448</td>
<td>0.03131</td>
</tr>
<tr class="odd">
<td>1</td>
<td>201359</td>
<td>N1103</td>
<td>538448</td>
<td>0.00204</td>
<td></td>
</tr>
</tbody>
</table>

-   The fourth column in the .imiss file (N\_MISS) denotes the number of
    missing SNPs and the sixth column (F\_MISS) denotes the proportion
    of missing SNPs per individual.

### Step 4: Identification of individuals with outlying heterozygosity rate

<br> PLINK command <br> **./plink –bfile 1\_QC\_Raw\_GWAS\_data –het
–out outlying\_heterozygosity\_rate** <br> <br>

<table>
<thead>
<tr class="header">
<th>FID</th>
<th>IID</th>
<th>O(HOM)</th>
<th>E(HOM)</th>
<th>N(NM)</th>
<th>F</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>1</td>
<td>201320</td>
<td>228640</td>
<td>2.259e+05</td>
<td>325867</td>
<td>0.0272</td>
</tr>
<tr class="even">
<td>1</td>
<td>201337</td>
<td>227550</td>
<td>2.26e+05</td>
<td>326067</td>
<td>0.01506</td>
</tr>
</tbody>
</table>

**NOTE: plink2 format will give results in a different way** <br> -
Command creates the file “outlying\_heterozygosity\_rate.het”, in which
the third column denotes the observed number of homozygous genotypes
\[O(Hom)\] and the fifth column denotes the number of non-missing
genotypes \[N(NM)\] per individual.

    # Missing individual & Heterozygosity rate
    miss <- fread("SEX_data/Missing_sample/missing_data_rate.imiss")
    hetro <- fread("SEX_data/Missing_sample/Heterozygosity_rate/outlying_heterozygosity_rate.het")
    head(miss, 2)
    head(hetro, 2)

Next we calculate observed heterozygosity rate

    # Calculate the observed heterozyosity rate
    hetro$obs_hetero_rate <- ((hetro$`N(NM)`)-hetro$`E(HOM)`)/hetro$`N(NM)`

-   Merge the miss and hetro dataframe above created

<!-- -->

    # Merge missing file and heterozygoisty file
    hetro_miss <- miss %>% 
      left_join(hetro, by = "IID")

    # Creating plot
    png("Missing_hetero_check.png")
    ggplot(hetro_miss, aes(x = F_MISS, y = obs_hetero_rate))+
      geom_point(alpha = 0.5, col = "#00bfff")+
      labs(x ="Proportion of Missing genotypes(log scale)", y = "Heterozygosity rate")+
      scale_x_log10(limits = c(0.0001, 1))+
      theme_classic()+
      scale_y_continuous(limits = c((min(hetro_miss$obs_hetero_rate) - 0.02), (max(hetro_miss$obs_hetero_rate) + 0.02)))+
      geom_hline(yintercept = mean(hetro_miss$obs_hetero_rate) + 3*sd(hetro_miss$obs_hetero_rate), col = "Grey")+
      geom_hline(yintercept = mean(hetro_miss$obs_hetero_rate) - 3*sd(hetro_miss$obs_hetero_rate), col = "Grey")+
      geom_vline(xintercept = 0.006, col = "Grey")+
      geom_point(data=hetro_miss %>%
                   filter(F_MISS > 0.01),
                 pch = 19,
                 size=1.6,
                 colour = "#e75480")+
      annotate(geom="text", x=0.035, y=0.303, label="Missing genotypes (>1%)",
               color="#003300")+
      annotate(geom="text", x=0.1, y=0.285, label="Excess heterozygosity rate (Â± 3 sd from mean)",
               color="#003300")+
      annotate(geom="text", x=0.00017, y=0.315, label="+3 sd from mean",
               color="#003300")+
      annotate(geom="text", x=0.00017, y=0.301, label="-3 sd from mean",
               color="#003300")+
      annotate(geom="text", x=0.0015, y=0.31, label="Genotyping > 99%",
               color="#003300")

    dev.off()

The above R script will create a png file containing information based
on missing rate and heterozygosity rate. For missing rate we set a
threshold value &lt; 0.01 and for heterzygosity rate between + and - 3
standard deviation.

<img src="Missing_hetero_check.png" alt="Individual missingness and heterozygoisty rate"  />
<p class="caption">
Individual missingness and heterozygoisty rate
</p>

### Step 5: Identification of duplicate samples

In a population based study, it is important that all samples should be
unrelated. Presence of duplicate, first- or second-degree relatives will
introduce bias in the study as their genotypes will be overrepresented.
This step was used to identify all related and duplicate individuals for
removal. A metric (identity by state, IBS) for each pair of individuals
was calculated to identify duplicate samples. IBS is defined as, at a
locus, two individuals who have an identical nucleotide sequence or the
same allele.

<table>
<thead>
<tr class="header">
<th>PIHAT values</th>
<th>Relation</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>1</td>
<td>Identical or Duplicate</td>
</tr>
<tr class="even">
<td>0.8 and above</td>
<td>Highly related</td>
</tr>
<tr class="odd">
<td>0.5</td>
<td>First Degree</td>
</tr>
<tr class="even">
<td>0.25</td>
<td>Highly related</td>
</tr>
<tr class="odd">
<td>0.125</td>
<td>First Cousin</td>
</tr>
<tr class="even">
<td>0.0625</td>
<td>Second Cousin</td>
</tr>
</tbody>
</table>

-   Check the relatedness. Use the independent SNPs (pruning) for this
    analysis and limit to autosomal chromosome only <br> PLINK command
    <br> **plink –bfile 2\_QC\_Raw\_GWAS\_data –chr 1-22 –make-bed –out
    Autosomal**
-   Create independent SNPs through pruning <br> PLINK command <br>
    **plink –bfile Autosomal –indep-pairwise 50 5 0.2 –out
    raw-GWAS-data** <br> it will generate raw-GWAS-data.prune.in file.
    This file use in next step
-   Check relatedness <br> PLINK command <br> **plink –bfile
    2\_QC\_Raw\_GWAS\_data –extract raw-GWAS-data.prune.in –genome –out
    related\_check** <br>

<img src="Related_samples.png" alt="Relatedness"  />
<p class="caption">
Relatedness
</p>

### Step 6: Identification of Individuals of divergent ancestry

#### Approach 1: Multidimensional scaling

<br> PLINK Command <br> **plink –bfile 3\_QC\_Raw\_GWAS\_data –extract
raw-GWAS-data.prune.in –genome –cluster –mds-plot 10**

-   Visualizing population structure using MDS is useful for identifying
    subpopulations, population stratification and systematic genotyping
    or sequencing errors, and can also be used to detect individual
    outliers that may need to be removed, e.g. European-Americans
    included in a study of African-Americans.

#### Approach 2: Principal component analysis

<br> PLINK Command <br> **plink –bfile 3\_QC\_Raw\_GWAS\_data –genome
–cluster –pca 10**

-   We can also utilize hapmap data to perform PCA

## Per SNP quality control

#### Hardy Weinberg equilibrium

-   If the frequency of observed genotypes of a variant in a population
    can be derived from the observed allele frequencies, the genetic
    variant is said to be in Hardy–Weinberg equilibrium.
-   Genotype frequencies remain stable from one generation to another in
    the absence of any evolutionary pressure (selection, mutation and
    migration).
-   The goodness-of-fit test is used to test HWE is not a powerful test
    and reliability of results depends on the sample size.

<table>
<thead>
<tr class="header">
<th>Sample size</th>
<th>50</th>
<th>100</th>
<th>200</th>
<th>300</th>
<th>400</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>AA count</td>
<td>21</td>
<td>42</td>
<td>84</td>
<td>126</td>
<td>168</td>
</tr>
<tr class="even">
<td>AB count</td>
<td>25</td>
<td>50</td>
<td>100</td>
<td>150</td>
<td>200</td>
</tr>
<tr class="odd">
<td>BB count</td>
<td>4</td>
<td>8</td>
<td>16</td>
<td>24</td>
<td>32</td>
</tr>
<tr class="even">
<td>HWE p value</td>
<td>0.52</td>
<td>0.26</td>
<td>0.08</td>
<td><strong>0.02</strong></td>
<td><strong>0.009</strong></td>
</tr>
</tbody>
</table>

Source: Genetic Epidemiology: Mehmet T Dorak <br> - Rule of thumb: the
heterozygote frequency can only reach a maximum of 50%. If heterozygote
frequencies are more than 50%, it is a clear sign of HWD, regardless of
statistical test result. - Most common reason is not biological,
**genotyping error** is most plausible exploration.

<br> PLINK Command <br> **./plink2 –bfile 4\_QC\_Raw\_GWAS\_data –geno
0.01 –hwe 0.00000001 –make-bed –out 5\_QC\_Raw\_GWAS\_data**

# References

1- Marees, A.T., et al, 2018. A tutorial on conducting genome‐wide
association studies: Quality control and statistical analysis. *Int J
Methods Psychiatr Res*, Jun; 27(2): e1608.

2- Anderson, C.A. et al, 2010. Data quality control in genetic
case-control association studies. *Nat Protoc*, Sep:5(9):1564-73

3- Singh, Sandeep Kumar, “A Case-Only Genome-wide Association Study of
Gender- and Age-specific Risk Markers for Childhood Leukemia” (2015).
FIU Electronic Theses and Dissertations. 1832
