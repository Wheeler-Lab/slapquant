{
    if ($1 ~ /^@/) { # don't care about the header
        next
    }
    flags = int($2)
    if (and(flags, 4)) { # unmapped flag set, sequence not aligned
        next
    }
    is_forward = and(int($2), 16) == 0;  # SEQ reversed flag unset

    rname = $3
    pos = int($4)
    
    # Look for alignments that are soft clipped at the start of the alignment
    if (match($6, /^([0-9]+)S([0-9]+)M/, groups_at_start)) {
        nr_matched = int(groups_at_start[2])
        nr_clipped = int(groups_at_start[1])
        if ((nr_matched > 20) || (nr_clipped >= 6)) {  # We want at least 6 soft clipped bases, and 20 aligned bases
            clipped = substr($10, 1, nr_clipped)
            if (is_forward) {
                strand = "+"
            } else {
                strand = "-"
            }
            print(rname, pos + nr_clipped, nr_matched, strand, "start", clipped)
        }
    }
    
    if (match($6, /([0-9]+)M([0-9]+)S$/, groups_at_end)) {
        nr_matched = int(groups_at_end[1])
        nr_clipped = int(groups_at_end[2])
        if ((nr_matched > 20) || (nr_clipped >= 6)) {
            clipped = substr($10, length($10)-nr_clipped+1, nr_clipped)
            if (is_forward) {
                strand = "+"
            } else {
                strand = "-"
            }
            print(rname, pos + nr_matched, nr_matched, strand, "end", clipped)
        }
    }
}