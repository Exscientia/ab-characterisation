#!/bin/bash

ROSETTA3=<ROSETTA_BASE_DIR>

$ROSETTA3/main/source/bin/rosetta_scripts.static.linuxgccrelease \
-database $ROSETTA3/main/database \
-in:file:s <INPUT_FILE> \
-in:file:native <INPUT_FILE> \
-parser:protocol ./rosetta_metrics_complex.xml \
-beta \
-include_sugars \
-alternate_3_letter_codes pdb_sugar \
-load_PDB_components false \
-auto_detect_glycan_connections \
-write_glycan_pdb_codes \
-output_alternate_atomids \
-write_pdb_link_records