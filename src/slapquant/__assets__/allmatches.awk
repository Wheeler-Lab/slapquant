function revcomp(arg) {
    o = ""
    for(i = length(arg); i > 0; i--)
        o = o c[substr(arg, i, 1)]
    return(o)
}

BEGIN {
    c["A"] = "T"; c["C"] = "G"; c["G"] = "C"; c["T"] = "A" ; OFS="\t"

    slsequence_r = revcomp(slsequence);
    polyAsequence_r = revcomp(polyAsequence);
}

/^@/ {}  # Don't care about the header
{
    # if ($1 ~ /^@/) { # don't care about the header
    #     next
    # }
    flags = int($2)
    if (and(flags, 4)) { # unmapped flag set, sequence not aligned
        next
    }
    is_forward = and(int($2), 16) == 0;  # SEQ reversed flag unset

    # Skip sequences containing the sl sequence and the poly A motif.
    if ( \
        (($10 ~ slsequence)   || ($10 ~ polyAsequence)) || \
        (($10 ~ slsequence_r) || ($10 ~ polyAsequence_r)) \
    ) {
        next
    }

    rname = $3
    pos = int($4)
    
    # Look for alignments that are soft clipped at the start of the alignment
    if (match($6, /([0-9]+S|)([0-9]+)M/, groups)) {
        nr_matched = int(groups[2])
        if ((nr_matched > 20)) {  # We want at least 20 aligned bases
            nr_clipped = int(groups[1])  # int() ignores characters at the end, handy!
            print(rname, pos + nr_clipped, nr_matched, "start", "*", "*")  # We don't need the clipped sequence in this case.
        }
    }
}