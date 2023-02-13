# Introduction

**Tools, Software and Files Required**

-   **plink** <https://www.cog-genomics.org/plink2/>
    <https://www.cog-genomics.org/plink/2.0/>

-   **bcftools and samtools** <http://www.htslib.org/download/>

-   **vcftools** <https://sourceforge.net/projects/vcftools/files/>

-   **Check for compatibility with the input required by the Sanger
    server.**
    <http://qbrc.swmed.edu/zhanxw/software/checkVCF/checkVCF-20140116.tar.gz>

-   **Haplotype Reference Consortium v1.1 panel for HRC site list**
    (HRC.r1-1.GRCh37.wgs.mac5.sites.tab)
    <ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/>

-   **HRC preparation checking tool**
    <https://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.2.9.zipa>

-   Set the path for these tools so you directly call them from a
    different folder. To set the path first you have enter in that
    specific tools **export PATH=*P**A**T**H*:(pwd)**

## Basic input files

We expect two basic input files, Raw.ped and Raw.map, for association
analysis. Here Raw is the file name. These files can also be converted
into binary format. Description of ped and map files is givne below:

### PED File

The PED file is a white-space (space or tab) delimited file: the first
six mandatory columns are given below.

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

Legend:

-   FID = Family ID
-   IID = Individual ID
-   PID = Paternal ID
-   MID = Maternal ID
-   Sex (1 = male; 2 = female; other = -9)
-   P = Phenotype

Sixth column represents phenotype. The phenotype can be either an
affection status or quantitative trait. Plink can automatically detect
the phenotype based on the code.

### MAP file

Each line of the MAP file explains a single marker and contains four
columns. MAP file for gender- and age-specific analyses were similar.

-   Chromosome (1-22, X, Y or 0 if unplaced)
-   rs# or snp identifier
-   Genetic distance (morgans)
-   Base-pair position (bp units)

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

### Binary file format

Binary files (fam, bed and bin) are smaller in size and run faster than
basic files (ped and map). Quality control, statistical analysis and
data entry (phenotype) are much easier and takes lesser time.

# Methods

## Steps to perform Genome-wide association study:

1.  If we get ped and map file, first should convert it into binary
    format (bed, bim and fam) <br> <br> **./plink –file raw\_GWAS\_data
    –make-bed** <br> <br> Note: Earlier version of plink reads X
    chromosome as 23, Y chromosome as 24, pseudo-autosomal region of X
    as 25 and mitochondrial as 26. However, plink2 reads them with X, Y
    etc. <br> <br>

**Quality controls steps:** <br> GWASs use thousands of markers to test
association. Even a low rate of error or bias will introduce false
positive associations. It has been suggested that removal of handful
individuals or small percentage of markers should not markedly decrease
the overall power of the study. Removal of one marker has greater impact
than removal of one individual. Removal of one marker is potentially an
overlooked disease association. Thus, implementing QC per-individual
before QC per-marker will maximize the number of markers remaining in
the study.

### Per Individual Quality Control

Per-individual QC screens genotype to identify subjects that may
introduce bias, if not removed. There are several steps of
per-individual QC for a GWAS data set.

#### Identification of individuals with discordant sex information

<br> PLINK command <br> **plink –bfile raw\_GWAS\_data –check-sex –out
GWAS\_Sex\_Check** <br>

-   Command Create a list of individuals with discordant sex data in
    file “GWAS\_Sex\_Check.sexcheck”. Column 3 denotes ascertained sex
    and column 4 denotes sex according to genotype data. When the
    homozygosity rate is more than 0.2 but less than 0.8, the genotype
    data are inconclusive regarding the sex of an individual and these
    are marked in column 4 with a 0.
-   Extract the IDs of individuals with discordant sex information. In
    situations in which discrepancy cannot be resolved, remove the
    individuals through following command. <br> <br> PLINK command <br>
    **plink –bfile raw\_GWAS\_data –remove
    discordant-sex-individuals-file.txt –make-bed –out
    1\_QC\_Raw\_GWAS\_data** <br> <br> (File
    “discordant-sex-individuals-file.txt”, should contain only FID and
    IID of the individuals that have to be removed)

<!-- -->

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

<img src="Gender_check.png" alt="Discordant Sex information"  />
<p class="caption">
Discordant Sex information
</p>

#### Identification of individuals with elevated missing data rates

<br> PLINK command <br> **./plink2 –bfile 1\_QC\_Raw\_GWAS\_data
–missing –out missing\_data\_rate** <br> - Command creates the files
“missing\_data\_rate.imiss” and “missing\_data\_rate.lmiss”. - The
fourth column in the .imiss file (N\_MISS) denotes the number of missing
SNPs and the sixth column (F\_MISS) denotes the proportion of missing
SNPs per individual.

#### Identification of individuals with outlying heterozygosity rate

<br> PLINK command <br> **./plink –bfile 1\_QC\_Raw\_GWAS\_data –het
–out outlying\_heterozygosity\_rate** <br> <br> **NOTE: plink2 format
will give results in a different way** <br> - Command creates the file
“outlying\_heterozygosity\_rate.het”, in which the third column denotes
the observed number of homozygous genotypes \[O(Hom)\] and the fifth
column denotes the number of non-missing genotypes \[N(NM)\] per
individual.

    # Missing individual & Heterozygosity rate
    miss <- fread("SEX_data/Missing_sample/missing_data_rate.imiss")
    hetro <- fread("SEX_data/Missing_sample/Heterozygosity_rate/outlying_heterozygosity_rate.het")
    head(miss, 2)
    head(hetro, 2)

    # Calculate the observed heterozyosity rate
    hetro$obs_hetero_rate <- ((hetro$`N(NM)`)-hetro$`E(HOM)`)/hetro$`N(NM)`

-   Merge the “missing\_data\_rate.smiss” and
    “outlying\_heterozygosity\_rate.het” hetroandmiss&lt;- merge(hetro,
    miss, by=“IID”)

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

<img src="Missing_hetero_check.png" alt="Individual missingness and heterozygoisty rate"  />
<p class="caption">
Individual missingness and heterozygoisty rate
</p>

#### Identification of duplicated or related individuals

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

<img src="Related_samples.png" alt="Relatedness"  />
<p class="caption">
Relatedness
</p>

#### Identification of Individuals of divergent ancestry

**Population Stratification**

<img src="Population_stratification.jpg" alt="Population Stratification" width="60%" />
<p class="caption">
Population Stratification
</p>

-   A threat to validity of a genetic association studies and may
    introduce bias. Population stratification may confound results, if
    not adjusted. A population in a genetic association study may arise
    from two distinct sub-populations which may differ in their genetic
    ancestry. The presence of systematic differences in allele
    frequencies between these two sub populations may confound the
    result. There are thousands of neutral loci, known as genomic
    controls used to correct for population stratification in GWASs.

1.  **Multidimensional scaling** <br> PLINK Command <br> **plink –bfile
    3\_QC\_Raw\_GWAS\_data –extract raw-GWAS-data.prune.in –genome
    –cluster –mds-plot 10**

-   Visualizing population structure using MDS is useful for identifying
    subpopulations, population stratification and systematic genotyping
    or sequencing errors, and can also be used to detect individual
    outliers that may need to be removed, e.g. European-Americans
    included in a study of African-Americans.

1.  **Principal component analysis** <br> PLINK Command <br> **plink
    –bfile 3\_QC\_Raw\_GWAS\_data –genome –cluster –pca 10**

-   We can also utilize hapmap data to perform PCA

### Per SNP quality control

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
