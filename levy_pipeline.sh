#!/bin/bash

dir=
genelist=
taxlist=
R1_pattern=
R2_pattern=
extra_seqs=
asv_rra=
taxa_rra=
identity_cutoff=

for gene a:

bash step_1_download_refs

bash step_2_make_blastdb

bash step_3_pears

bash step_4_collapse

bash step_5_mk_seqfiles

bash step_6_blast

bash step_7_make_taxatable
