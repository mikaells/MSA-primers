# MSA-primers

A program to suggest primers (with IUPAC bases) from a messy multiple sequence alignments (MSA). They should absolutely be looked at in CLC or similar. It is also not smart enough to calculate any sort of melting point or primer dimers.

Needs an alignment, primer lenghts, amplicon size limits and a GC% difference cutoff, and will then suggest 3 candidate pairs optimized for fewest wobbly bases and lowest GC difference.

When run on tdaA.aln, the first pair will look like this:

E.g.

    ***pair 1 
    First:

    pos             kmer      divs  GC
    136 136 atcaaaacgatggaag 0.4155976 0.4

    Consensus:	atcaaaacgatggaag
    #wobbles:	0000010010000000
    #substitutions:	0000050080000000
    Primer:		atcaaracratggaag

    Second:

    pos             kmer      divs  GC
    36  36 tatcagtcaactcaag 0.2451616 0.4

    Consensus:	tatcagtcaactcaag
    #wobbles:	1000000000000001
    #substitutions:	1000000000000006
    Primer:		yatcagtcaactcaar

    --------------

One might then consider changing the GC-cutoff or maybe ignoring the first wobble of the second primer, since only one sequence has a divergence.
