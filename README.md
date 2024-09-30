# Antibody characterisation pipeline

## Overview
This repository contains code to run the antibody characterisation pipeline described in "Computational design of 
developable therapeutic antibodies: efficient traversal of binder landscapes and rescue of escape mutations" 
(see citation below). 

## Installation
### Environment setup
Use the `environment.yml` file to create a conda environment
```shell
conda env create -f environment.yml -n ab-characterisation
conda activate ab-characterisation
```
### Dependencies
Ensure that you have working installs of the following:

   1) Rosetta: https://www.rosettacommons.org/demos/latest/tutorials/install_build/install_build
   
   2) ChimeraX: https://www.cgl.ucsf.edu/chimerax/download.html 
      
      a) ensure that you can run the basic ChimeraX script in tests/data:
      ```
      ChimeraX --script tests/data/chimera_test_script.py --nogui
      ```
      b) set environment variable `DEBIAN_FRONTEND="noninteractive"`
   
   3) ANARCI: https://github.com/oxpig/ANARCI. Note: On MacOS machines, install the hmmer dependency via brew, otherwise via conda.
      
   4) Ensure you have the correct licences for all linked software.

## Testing your installation
You can test the installation of the environment using `pytest`. 
For this, first set the Rosetta base directory as an environment variable, for example like this:
```shell
export ROSETTA_BASE=/path/to/rosetta/rosetta.binary.linux.release-315
```
Then run pytest
```shell
pytest
```
Which will run an end-to-end example run of the pipeline on a set of 4 antibody sequences (note that depending on your
setup this may take 1h). 

## Running the pipeline
With the conda environment active, the pipeline can be run as follows:
```shell
ab-characterisation --input-file tests/data/test_pipeline.csv --rosetta-base-dir $ROSETTA_BASE
```
(assuming ROSETTA_BASE to have been set as described above).

If you want to multiprocess the pipeline, instead run as
```shell
mpiexec -n N_PROCESSES ab-characterisation --input-file tests/data/test_pipeline.csv --rosetta-base-dir $ROSETTA_BASE 
```

``` 
Usage: ab-characterisation [OPTIONS]

╭─ Options ─────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ *  --input-file                TEXT     Input .csv file, containing sequence_name, heavy_sequence, light_sequence and │
│                                         reference_complex columns. [required]                                         │
│    --chimera-resolution        FLOAT    Resolution of the map used for alignment within ChimeraX. [default: 6.0]      │
│    --output-dir                TEXT     Directory to which output files are written.                                  │
│                                         [default: ./ab_characterisation_output]                                       │
│    --rosetta-replicates        INTEGER  How many replicates to run for Rosetta characterisation steps. [default: 1]   │
│ *  --rosetta-base-dir          TEXT     Base directory for the Roestta software suite,                                │
│                                         e.g. /path/to/rosetta/rosetta.binary.linux.release-315 [required]             │
│    --top-n                     INTEGER  Top N candidate antibodies to provide from the provided .csv file of          │
│                                         antibodies [default: 10]                                                      │
│    --help                               Show this message and exit.                                                   │
╰────────────────────────────────────────────────
```
## Acknowledgements
The antibody characterisation pipeline was developed  by researchers and engineers at Exscientia:

- Frederic Dreyer

- Constantin Schneider

- Aleksandr Kovaltsuk

- Daniel Cutting
 
- Matthew J. Byrne
 
- Daniel A. Nissley
 
- Newton Wahome
 
- Henry Kenlay
 
- Claire Marks
 
- David Errington
 
- Richard J. Gildea
 
- David Damerell
 
- Pedro Tizei
 
- Wilawan Bunjobpol
 
- Sachin Surade
 
- Douglas E. V. Pires

- Charlotte M. Deane

## Citation
If you use this code in your research, please cite the following paper:

[BIBTEX HERE]