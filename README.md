# Measuring phasing performance
## Gold standard phased VCF files for NA12878

The files were generated from the GIAB v0.2 VCF file, available from the National Institute of Standards and Technology (NIST) Genome-In-A-Bottle (GIAB) consortium FTP site.
- ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/analysis/GIAB_integration/NIST_RTG_PlatGen_merged_highconfidence_v0.2_Allannotate.vcf.gz

## Running the performance evaluation script

The script measures phasing performance for a chromosome. The following example command assesses the phasing quality of an
input phased VCF file (examples/chr20.vcf.gz) by comparing with the gold standard phased VCF file ("gold_standard/chr20.snp.vcf.gz") and store the results into files with the prefix "examples/test".
```
$ ./scripts/measure_phasing_performance.pl -i examples/chr20.vcf.gz -r gold_standard/chr20.snp.vcf.gz -o examples/test
```
- Output files
  - examples/test.out: Summary of results
  - examples/test.err_pos: Each error position along with the distance to its upstream phased position (in bp) is
stored
  - examples/test.err_pos.more: This file is created to record phasing outcome classes for every position as the following:
    - CORRECT: correctly phased site
    - FIRST_POS: the first site for each phasing block (so not measured for switch error)
    - SWITCH_ERROR_LONG: long switch error
    - SWITCH_ERROR_POINT: point switch error
    - SWITCH_ERROR_UNDEF: undetermined switch error
    - UNPHASED: unphased site

When the optional "-a" parameter is specified with a file name, data for pairwise SNVs yield and accuracy measures are stored in the file, which can be used for plotting. This option usually takes much longer time than without the option.
```
$ ./scripts/measure_phasing_performance.pl -a examples/test.tsv -i examples/chr20.vcf.gz -r gold_standard/chr20.snp.vcf.gz -o examples/test
```

## In-house variant calls

This in-house set of NA12878 SNV variant calls was generated using the reference human genome GRCh37 coordinates with the Illumina
paired-end DNA reads (45X coverage, 148 bp), using BWA-MEM for read mapping and GATK for variant calling.

## Haplotype diversity

The haplotype diversity for each 1kb non-overlapping region was computed from the SHAPEIT 1000GP reference panel.

## Links
### Haplotype reference panels
- [SHAPEIT 1000GP reference panel](https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html)
- [Beagle 1000GP reference panel](http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/)
- [Haplotype Reference Consortium Release 1.1](https://ega-archive.org/datasets/EGAD00001002729)
### Phasing methods and haplotype data for NA12878
- [10X Genomics](https://www.ncbi.nlm.nih.gov/pubmed/26829319), data: ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/analysis/10XGenomics_calls_08142015/
- [Contiguity Preserving Transposon](https://www.ncbi.nlm.nih.gov/pubmed/25326703)
- [Fosmid-pool-based Phasing Strategy](https://www.ncbi.nlm.nih.gov/pubmed/22102577), [data](http://owww.molgen.mpg.de/~genetic-variation/SIH/data/)
- [Moleculo](https://www.ncbi.nlm.nih.gov/pubmed/24561555)
- [Beagle](https://faculty.washington.edu/browning/beagle/beagle.html)
- [Eagle2](https://data.broadinstitute.org/alkesgroup/Eagle/)
- [SHAPEIT](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html)
- [HapCUT](https://github.com/vibansal/hapcut)
