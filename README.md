# Paracottus_knerii_transcriptome
Reproducible code for the analysis of a *Paracottus knerii* ("stone sculpin" from Baikal) fry transcriptome.

# Methods

## Raw data

Total RNA was isolated from a dissected *Paracottus knerii* fry (30–60 days old) fixed in RNAlater with a MagMax kit (Thermo Fisher). After quality control (RNA concentration >100 ng/uL, RIN 7.6) 100 ng of RNA was processsed with a TruSeq Stranded mRNA library preparation kit (Illumina), and the library was sequenced with a NovaSeq 6000 device (2 x 101 bp). Q30 value of sequencing: 92.06 %. 
Demultiplexing of the sequencing reads was performed with Illumina bcl2fastq (version 2.20). Adapters were trimmed with Skewer (version 0.2.2) (Jiang et al. 2014). Up to this point, the analysis was performed by the CeGaT company.

Read quality was analyzed with FastQC v0.11.9 (https://github.com/s-andrews/FastQC) and was fairly good to proceed with assembly.

```{bash}
fastqc *
```

The reads were submitted to NCBI: [BioProject PRJNA1200955](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1200955) and [SRA sample SRR31806589](https://trace.ncbi.nlm.nih.gov/Traces/sra?run=SRR31806589).

## Transcriptome assembly

The main transcriptome assembly used in downstream analyses was performed with rnaSPAdes (Bushmanova et al., 2019, doi:10.1093/gigascience/giz100) v3.13.1 using the `--ss-fr` option. In addition, the Oyster River Protocol ([MacManes, 2018](https://peerj.com/articles/5428/); [readthedocs](https://oyster-river-protocol.readthedocs.io/en/latest/)) was utilized to compare assemblers.
Assembly quality was controlled with BUSCO v5.4.5 (Manni et al., 2021, doi: 10.1002/cpz1.323) using the Actinopterygii database (`actinopterygii_odb10`).

```{bash}
cd ~/S12139/S12139_1/02_assembly
rnaspades.py -t 6 -1 ../01_fastq/RNA_S12139Nr1.1.fastq.gz -2 ../01_fastq/RNA_S12139Nr1.2.fastq.gz --ss-fr -o Pkn_rnaspades_ssfr
busco -i ../02_assembly/Pkn_rnaspades_ssfr/transcripts.fasta -l ./actinopterygii_odb10/ -o Pkn_rnaspades_busco -m transcriptome --offline
```

All analyses except indicated otherwise were performed using a small computing cluster (64 Gb RAM, 6 physical cores, 12 virtual cores).

# Results

## Read QC

51M paired reads (10G bases) were obtained. Quality control with FastQC showed good quality (Q>28) and absense of remaining sequence adapters.

## Assembly QC

The rnaSPAdes assembly had the following BUSCO score:
```{}
       --------------------------------------------------
        |Results from dataset                             |
        --------------------------------------------------
        |C:80.5%[S:64.2%,D:16.3%],F:5.0%,M:14.5%,n:3640   |
        |2929   Complete BUSCOs (C)                       |
        |2336   Complete and single-copy BUSCOs (S)       |
        |593    Complete and duplicated BUSCOs (D)        |
        |181    Fragmented BUSCOs (F)                     |
        |530    Missing BUSCOs (M)                        |
        |3640   Total BUSCO groups searched               |
        --------------------------------------------------
```
It is not close to the desired 100 % but comparable with other fish transcriptome assemblies (e.g., [Zhou et al., 2020](https://doi.org/10.1038/s41597-020-0361-6) or [Kokkonen et al., 2024](https://doi.org/10.1111/eva.13735)) and is suitable for the planned downstream applications, primer design and as a database for MS-MS proteome analysis.

## Biological assembly QC

### Is it the right species?

## Annotation

## Gene finding
