# QDNAseq shallow coverage copy number calling pipeline

This repository contains code to run a modified version of QDNAseq on shallow WGS, on WXS and on targeted sequencing data.

## Running

Whole genome sequencing requires these parameters:
```
Rscript call_copynumber.R [samplename] [tumourbam] [normalname] [normalbam]
```

Whole exome sequencing - where the offtarget reads are used - takes an extra parameter which is a bed file specifying the regions that were targeted, i.e. all exons in the genome
```
Rscript call_copynumber_exome.R [samplename] [tumourbam] [normalname] [normalbam] [targetregions]
```

Targeted sequencing - where offtarget reads are used - takes an extra parameter which is a bed file specifying the regions that were targeted
```
Rscript call_copynumber_offtarget.R [samplename] [tumourbam] [normalname] [normalbam] [targetregions]
```

## Installing

### Dependencies
 * A slightly modified version of QDNAseq from here: [https://github.com/sdentro/QDNAseq/tree/dev](https://github.com/sdentro/QDNAseq/tree/dev)

### Installation
```
R -q -e 'devtools::install_github("sdentro/QDNAseq", ref="dev")'
```

## Reference data


