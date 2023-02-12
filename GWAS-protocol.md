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
six mandatory columns are given below. - Family ID - Individual ID -
Paternal ID - Maternal ID - Sex (1 = male; 2 = female; other = -9) -
Phenotype

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
