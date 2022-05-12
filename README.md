# CpGSelectionPipeline

This repository comprises scripts to select CpGs from bulk DNA methylation data for subsequent analysis with scTAM-seq. 

In general, three classes of CpGs are selected:
- Differentially Methylated CpGs (DMCs, in [RnBeads](RnBeads))
- Intermediately Methylated CpGs (IMCs, in [IMS](IMS))
- CpGs harboring Within-Sample Heterogeneity (WSH, in [WSH](WSH))

# Snakemake pipeline

The full pipeline can be executed using the Snakemake file located in [snakemake/full/Snakefile]. For a successful execution, [snakemake](https://snakemake.readthedocs.io/en/stable/) has to be installed as well as [sra-tools](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software), [TrimGalore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/), and [samtools](http://www.htslib.org/).

Additionally, the reference genome ([mouse mm10](http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/)) has to be downloaded and stored into a folder named *references/mm10/*.

Further necessary annotations include enhancer annotations, which can be downloaded as the file 'H3K27Ac\_mm10.bed' from https://www.science.org/doi/10.1126/science.1256271 and has to be put it into the folder 'annotations'. Similarly, download the file 'GSM774291\_FLDN1PU.1I.bb' from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM774291, process it with bedtools and liftover to convert it into a bed file mapped to mm10.

The pipeline can be started using:

```
snakemake --profile sge --use-conda --cluster-config cluster.yaml
```

# Contact

You can contact [Michael Scherer](michael.scherer@crg.eu) for questions and suggestions.
