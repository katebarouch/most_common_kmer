# Input:  Two strings p and q
# Output: An integer value representing the Hamming Distance between p and q.
def HammingDistance(p, q):
    hamming_distance = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            hamming_distance += 1
    return hamming_distance

print(HammingDistance("CTACAGCAATACGATCATATGCGGATCCGCAGTGGCCGGTAGACACACGT", "CTACCCCGCTGCTCAATGACCGGGACTAAAGAGGCGAAGATTATGGTGTG"))