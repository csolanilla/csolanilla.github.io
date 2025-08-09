def recurance_rabbits(months, pairs):
    a, b = 1, 1
    for month in range(2,months):
        a,b = b, b + (pairs*a)
    return b

print(recurance_rabbits(34,4))