# snuupy
[![DOI](https://zenodo.org/badge/302899070.svg)](https://zenodo.org/badge/latestdoi/302899070)

single-nucleus utility in python

# tutorial
![Schematic_diagram](./Schematic_diagram.png)

Git clone or download this repository then modify the paths in the configuration file before running snakemake. 

# Usage

snuupy.py [OPTIONS] COMMAND [ARGS]

Commands:

  - addGeneName
  - addPolyATag
  - addUnmappedBaseTag
  - barcodeAssignment
  - generateH5adFromKb
  - generateIlluminaWindow
  - generateIlluminaWindowFromKb
  - generateMtx
  - generateNanoporeWindow
  - getMismatch
  - getSpliceInfo
  - parseIllumina
  - polishReads
  - polyAClusterDetected
  - windowBlast


# packages required
Snakemake python3.7 pandas scipy numpy scanpy joblib loguru minimap2 poaV2 seqkit racon biopython picard portaion ont-fast5-api more_itertools

# acknowledge
This pipeline was inspired by Sicelore

# open source license
The source code is released under MIT license. 
