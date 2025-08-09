'''

Given k, m, n organisms
k: homozygous dominant for a factor
m: heterozygous
n: homozygous recessive

return probability that 2 randomly selected mating organisms
will produce offspring possessing a dominant allele
assuming they can all mate

2,2,2 = 0.78333

'''

def is_dom(k,m,n):
    # k = T[0]
    # m = T[1]
    # n = T[2]
    T = [list(range(x)) for x in [k,m,n]]
    total = sum(len(sublst) for sublst in T)

    final_prob_lst = []
    # for each k, m, and n value
    # get probability of picking parent 1
    for index1, parent1 in enumerate(T):
        # get prob of picking parent 1
        prob1 = len(parent1) / total

        # update the copy with that parent out of population
        testlist = [lst[:-1] if idx == index1 else lst for idx, lst in enumerate(T)]

        # find prob of picking parent 2
        for index2, parent2 in enumerate(testlist):
            total2 = sum(len(sublst) for sublst in testlist)
            prob2 = len(parent2) / total2

            # apply homo + hetero probabilities
            domprob = 1
            if index1 == 1 and index2 == 1:
                domprob = 0.75
            elif index1 == 1 and index2 == 2:
                domprob = 0.5
            elif index1 == 2 and index2 == 1:
                domprob = 0.5
            elif index1 == 2 and index2 == 2:
                domprob = 0

            # apply final domprob:
            final_prob = prob1 * prob2 * domprob
            final_prob_lst.append(final_prob)

    return sum(final_prob_lst)



print(is_dom(24,28,25))