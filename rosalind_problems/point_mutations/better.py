def count_point_mutations(filename):
    with open(filename) as f:
        s1, s2 = f.read().split()
        mutations = sum((a != b for  a,b in zip(s1, s2)))

    return mutations

answer = count_point_mutations('rosalind_hamm.txt')

print(answer)