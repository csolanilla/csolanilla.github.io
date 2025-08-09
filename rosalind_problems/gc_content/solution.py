from collections import Counter

def gc_content(filename):
    combo = ""
    dna_id = []
    seq = []

    with open(filename) as f:
        content = f.read().split()
        for line in content:
            if line.startswith('>'):
                dna_id.append(line)
                if combo:
                    seq.append(combo)
                combo = ""
            else:
                combo += line
        if combo:
            seq.append(combo)

    dna_data = dict(zip(dna_id, seq))
    gc_targets = ['G', 'C']
    ranked_gc = []
    for key, value in dna_data.items():
        count_string = Counter(value)
        gc_toal = (sum(count_string[ele] for ele in gc_targets) / len(value)) * 100
        ranked_gc.append([key.strip('>'), gc_toal])

    return sorted(ranked_gc, key=lambda x: x[1])

dna_data = gc_content('rosalind_gc.txt')

print(dna_data[-1][0])
print(dna_data[-1][1])

