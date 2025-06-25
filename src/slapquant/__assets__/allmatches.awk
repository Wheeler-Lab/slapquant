@include "util.awk"

BEGIN {
    c["A"] = "T"; c["C"] = "G"; c["G"] = "C"; c["T"] = "A" ; OFS="\t"

    slsequence_r = revcomp(slsequence);
    polyAsequence_r = revcomp(polyAsequence);
}

/^@/ {}  # Don't care about the header
{
    flags = int($2)
    if (and(flags, 4)) { # unmapped flag set, sequence not aligned
        next
    }

    rname = $3
    pos = int($4)
    cigar = $6

    if (match($6, /^([0-9]+)S([0-9]+)M/, groups_at_start)) {
        nr_clipped = int(groups_at_start[1])
        clipped = substr(sequence, 1, nr_clipped)
        if ((clipped ~ slsequence) || (clipped ~ polyAsequence_r)) {
            next
        }
    }
    if (match($6, /([0-9]+)M([0-9]+)S$/, groups_at_end)) {
        nr_clipped = int(groups_at_end[2])
        seqlen = length(sequence)
        clipped = substr(sequence, seqlen-nr_clipped+1, nr_clipped)
        if ((clipped ~ slsequence_r) || (clipped ~ polyAsequence)) {
            next
        }
    }

    print(rname, pos, reference_consumed(cigar), "*", "*", "*")
}