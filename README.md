# snuupy
[![DOI](https://zenodo.org/badge/302899070.svg)](https://zenodo.org/badge/latestdoi/302899070)

single-nucleus utility in python

# Tutorial
![Schematic_diagram](./Schematic_diagram.png)

Git clone or download this repository then modify the paths in the configuration file before running snakemake. 

# Usage
```text
snuupy.py [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  addGeneName                   parse polished reads mapping results and
                                get...

  addPolyATag                   add polyA length tag for bam files
  addUnmappedBaseTag            get unmapped base tag
  barcodeAssignment             assign barcode for each Nanopore read;
                                based...

  generateH5adFromKb            get adata from kbpython result --tg:...
  generateIlluminaWindow        output illumina reads based on mapping info
  generateIlluminaWindowFromKb  generate illumina windows from kb_python...
  generateMtx                   generate matrices
  generateNanoporeWindow        output nanopore reads based on mapping info
  getMismatch                   calculate mismatch based on blast results
  getSpliceInfo                 get splice information which used for...
  multilayerClustering          cluster nuclei by snf or mofa
  parseIllumina                 parse Illumina bam file and generate...
  polishReads                   polish barcode assigned Nanopore reads
  polyAClusterDetected          detect PolyA Cluster
  windowBlast                   blast find potential UMI/Bc
```
<img src="./snakemake/pipeline.svg" width="500" height="700">

# Packages required
- python3 
  - pandas 
  - scipy 
  - numpy 
  - scanpy 
  - joblib 
  - loguru 
  - portion 
  - ont-fast5-api 
  - more_itertools
  - biopython
  - pyfastx
  - snfpy
  - mofapy2
- R
  - scran
- Snakemake 
- minimap2 
- poaV2 
- seqkit 
- racon 
- picard 


# Acknowledge
This pipeline was inspired by Sicelore

Function multilayerClustering is inspired by scLAPA
