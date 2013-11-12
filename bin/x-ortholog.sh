#!/bin/bash

#
# A generic BLAST command wrapper for reciprocal blast
#
# CLI Parameters:
# - $1: Blast strategy to use
# - $2: Blast DB
# - $3: Query file

set -e
set -u

case "$1" in
'ncbi-blast')
blastn -db $2/db -query $3 -evalue 0.00001 -word_size 7 -max_target_seqs 1 -max_hsps_per_subject 1 -num_threads 1  -outfmt '6 qseqid sseqid evalue score qgi bitscore length nident positive mismatch pident ppos qacc gaps gaopen qaccver qlen qframe qstart qend sstrand sstart send'
;;

'wu-blast')
wu-blastn $2/db $3 -mformat=2 -V=1 -B=1 -spoutmax=1 -warnings -errors -notes W=7 M=5 N=-4 Q=20 R=10 -e 0.00001 -cpus 1 -filter=seg -lcfilter
;;

*) echo "Not a valid BLAST strategy: $1"; exit 1
;;

esac
