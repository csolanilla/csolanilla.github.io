def qt(s):
    with open(s) as f:
        content = f.read()
    return content.count("A"), content.count("G"), content.count("C"), content.count("T")


print(qt('rosalind_dna (1).txt'))