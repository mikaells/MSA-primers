# MSA-primers

A program to suggest primers from a messy multiple sequence alignments (MSA).

Needs an alignment, primer lenghts, amplicon size limits and a GC% difference cutoff, and will then suggest 3 candidate pairs optimized for fewest wobbly bases and lowest GC difference.

When run on tdaA.aln, the first pair will look like this:

E.g.

    First:

    pos             kmer     divs  GC
    136 atcaaaacgatggaag 0.4155976 0.4

    Consensus:	    atcaaaacgatggaag
    Substitutions:	0000010010000000
    Primer:		    atcaaracratggaag

    Second:

    pos             kmer     divs  GC
    31 atggatatcagtcaac 0.04786485 0.4

    Consensus:	    atggatatcagtcaac
    Substitutions:	0000010000000000
    Primer:		    atggayatcagtcaac
