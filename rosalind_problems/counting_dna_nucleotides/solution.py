def nucelotide_counter(filename):
    try:
        with open(str(filename)) as f:
            content = str(f.read())
            for i in 'ACGT':
                print(f"{i}: {content.count(i)}")

    except FileNotFoundError:
        print(f"File '{filename}' not found")



nucelotide_counter('rosalind_dna (1).txt')