BEGIN {
    c["A"] = "T"; c["C"] = "G"; c["G"] = "C"; c["T"] = "A" ; OFS="\t"
}

function revcomp(arg) {
    o = ""
    for(i = length(arg); i > 0; i--)
        o = o c[substr(arg, i, 1)]
    return(o)
}

function reference_consumed(cigar) {
    gsub(/[0-9]+[ISHP]/, "", cigar)
    patsplit(cigar, advances, /[0-9]+/)
    consumed = 0
    for (i in advances) {
        consumed += advances[i]
    }
    return(consumed)
}