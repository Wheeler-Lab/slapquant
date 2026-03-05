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

        # allowed proportion of mismatches (pass via -v mismatch_prop=0.1); default 0
        if (mismatch_prop < 0) mismatch_prop = 0
        if (mismatch_prop > 1) mismatch_prop = 1

        # if mismatch_prop == 0 use exact/pattern matching logic
        if (mismatch_prop == 0) {
            if ((clipped ~ slsequence) || (clipped ~ polyAsequence_r)) {
                next
            }
        } else {
            lenc = length(clipped)
            allowed = int(lenc * mismatch_prop + 0.999999)

            # compare clipped to the prefix of slsequence and prefix of polyAsequence_r
            pat1 = substr(slsequence, 1, lenc)
            pat2 = substr(polyAsequence_r, 1, lenc)

            mism1 = 0
            lim1 = (lenc < length(pat1) ? lenc : length(pat1))
            for (i = 1; i <= lim1; i++) if (substr(clipped, i, 1) != substr(pat1, i, 1)) mism1++
            mism1 += (lenc > length(pat1) ? lenc - length(pat1) : length(pat1) - lenc)

            mism2 = 0
            lim2 = (lenc < length(pat2) ? lenc : length(pat2))
            for (i = 1; i <= lim2; i++) if (substr(clipped, i, 1) != substr(pat2, i, 1)) mism2++
            mism2 += (lenc > length(pat2) ? lenc - length(pat2) : length(pat2) - lenc)

            if ((mism1 <= allowed) || (mism2 <= allowed)) {
                next
            }
        }
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