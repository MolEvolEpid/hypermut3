python ../hypermut.py data/manuscript_seqs.fasta G A . RD -p output/strict-keepgaps- -k -m strict
python ../hypermut.py data/manuscript_seqs.fasta G A . RD -p output/strict-skipgaps- -m strict
python ../hypermut.py data/manuscript_seqs.fasta G A . RD -p output/partial-keepgaps- -k -m partial
python ../hypermut.py data/manuscript_seqs.fasta G A . RD -p output/partial-skipgaps- -m partial
