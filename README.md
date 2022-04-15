# Dx_BAF

Software to parse VCF file and create B-allele frequency IGV track.

# Usage
usage: make_BAF_igv.py [-h] [-o OUTPUTFILE] [-c] [--mindepth MINDEPTH]
                       inputfile

positional arguments:
  inputfile             input VCF file (genotyped VCF)

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUTFILE, --outputfile OUTPUTFILE
                        output filename. Without argument, output will be
                        printed in stdout
  -c, --compressed      VCF input is compressed (.gz)
  --mindepth MINDEPTH   Threshold for minimum depth (DP) of SNV (default = 15)

## Installation
To run the BAF scripts we need to create a virtual python environment

### Make virtual python environment (note this could be UMCU-HPC specific)
```bash
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

All scripts are tested using Python 3.6.8
