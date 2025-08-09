'''
In DNA strings, symbols 'A' and 'T' are complements of each other, as are 'C' and 'G'.

The "reverse complement" of a DNA string 's' is the string 'sc'
formed by reversing the symbols of s, then taking the complement
of each symbol (e.g., the reverse complement of "GTCA" is "TGAC").

Given: A DNA string s
 of length at most 1000 bp.

Return: The reverse complement sc
 of s
'''

def reverse_complement(dna_file):
    translate_table = str.maketrans({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'})
    with open(dna_file) as f:
        content = f.read()
    rev = content.translate(translate_table)[::-1]
    return rev

answer = reverse_complement('rosalind_revc.txt')
print(answer)