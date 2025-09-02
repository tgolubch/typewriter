# typewriter
Aggregate blast results and/or perform MLST for a set of sample sequences vs a given panel of queries.

Searches a given list of blast databases or fasta files for either:
- a panel of one or more sequences (eg resistance genes) from a multifasta query file, or
- the MLST profile for each sample and produces combined output for the dataset.

Similarity is reported as relative coverage, calculated as (sequence identity * query coverage).

With --save_results also outputs the hits for each query sequence in fasta format.

With --do_mutations reports sequence differences (eg. drug resistance mutations) in the form XNY, where X is the reference (subject) residue, N is the position in the reference sequence, and Y is the query residue at this position.

Example: MLST on assembled data from a single sample:
```
  {typewriter.py} --infastapath /path/to/my/sample/contigs.fasta
                  --outprefix my_mlst_results
                  --outdir    /path/to/outut/directory
                  --mlst True
                  --word_size 11
                  --allelefile  /path/to/my/allelefile.tsv
                  --stfile      /path/to/my/STfile.tsv
```
