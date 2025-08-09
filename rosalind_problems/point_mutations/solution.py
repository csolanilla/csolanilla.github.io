def count_mutations(filename):
    strings = []
    with open(filename) as f:
        for line in f.readlines():
            strings.append(line)
    point_mutations = 0
    for char in range(len(strings[0])):
        if str(strings[0][char]) == str(strings[1][char]):
            print(f"{strings[0][char]} | {strings[1][char]}")
        else:
            print(f"***{strings[0][char]} | {strings[1][char]}***")
            point_mutations += 1
    return point_mutations



answer = count_mutations('rosalind_hamm.txt')

print(answer)