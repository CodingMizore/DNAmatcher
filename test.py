def get_reverse_complement_of(s):
    complement_of = {'C': 'G', 'A': 'T', 'G': 'C', 'T': 'A'}
    s_rc = ''
    for ch in reversed(s):
        s_rc += complement_of[ch]
    return s_rc

print(
    get_reverse_complement_of('CAGCAACGCTAGTGGATAATTTTGTTGCTATAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGC')
)
