@include "util.awk"

{
    if ($1 ~ /^@/) { # don't care about the header
        next
    }
    flags = int($2)
    if (and(flags, 4)) { # unmapped flag set, sequence not aligned
        next
    }

    rname = $3
    pos = int($4)
    sequence = $10
    cigar = $6
    
    if (sequence == "*") {
        next
    }

    # Look for alignments that are soft clipped at the start of the alignment
    if (match(cigar, /^([0-9]+H)?([0-9]+)S/, groups_at_start)) {
        nr_matched = reference_consumed(cigar)
        nr_clipped = int(groups_at_start[2])
        if ((nr_matched > 20) && (nr_clipped >= 6)) {  # We want at least 6 soft clipped bases, and 20 aligned bases
            clipped = substr(sequence, 1, nr_clipped)
            remainder = substr(sequence, nr_clipped+1)
            print(rname, pos, nr_matched, "start", clipped, remainder)
        }
    }
    
    if (match(cigar, /([0-9]+)S([0-9]+H)?$/, groups_at_end)) {
        nr_matched = reference_consumed(cigar)
        nr_clipped = int(groups_at_end[1])

        seqlen = length(sequence)

        if ((nr_matched > 20) && (nr_clipped >= 6)) {
            clipped = substr(sequence, seqlen-nr_clipped+1, nr_clipped)
            remainder = substr(sequence, 1, seqlen-nr_clipped)
            print(rname, pos + nr_matched, nr_matched, "end", clipped, remainder)
        }
    }
}