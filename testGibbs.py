# Gibbs Sampling

import random
        
def AdjustedRandomizedMotifSearch(Dna, k, t):
    M = ['TGA', 'GTT', 'GAA', 'TGT']
    BestMotifs = M
    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs 


# Input:  A list of strings Dna, and integers k and t
# Output: RandomMotifs(Dna, k, t)
def RandomMotifs(Dna, k, t):
    motifs = []
    for i in range(t):
        start = random.randint(0, len(Dna[i]) - k)
        motifs.append(Dna[i][start:start + k])
    return motifs

def Motifs(profile, dna):
    k = len(profile['A'])
    D = []
    for i in range(0, len(dna)):
        k_mers = []
        pr = []
        for x in range(len(dna[i])-k+1):
            k_mers.append(dna[i][x:x+k])
        for k_mer in k_mers:
            pr.append(Pr(k_mer, profile))
        D.append(k_mers[pr.index(max(pr))])
    return D

# Input:  A set of kmers Motifs
# Output: ProfileWithPseudocounts(Motifs)
def ProfileWithPseudocounts(Motifs):
    t = len(Motifs) + 4
    k = len(Motifs[0])
    profile = {} # output variable
    
    count = CountWithPseudocounts(Motifs)
    
    for symbol in "ACGT":
        profile[symbol] = []
        for j in range(k):
            profile[symbol].append(count[symbol][j] / t)
        
    return profile


# Input:  A set of kmers Motifs
# Output: CountWithPseudocounts(Motifs)
def CountWithPseudocounts(Motifs):
    count = {} # initializing the count dictionary
    
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(1)
                
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
            
    return count

def Count(Motifs):
    count = {} # initializing the count dictionary
    
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(0)
                
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
                
    return count

# Input:  A set of kmers Motifs
# Output: A consensus string of Motifs.
def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)
    
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
        
    return consensus

# Input:  A set of k-mers Motifs
# Output: The score of these k-mers.
def Score(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = Consensus(Motifs)
    
    score = 0 
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            if symbol != consensus[j]:
                score += 1
                
    return score

def ProfileMostProbableKmer(text, k, profile):
    most_probable_kmer = ""
    max_p = -1
    
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        p = Pr(kmer, profile)
        if p > max_p:
            max_p = p
            most_probable_kmer = kmer
    return most_probable_kmer

def Normalize(Probabilities):
    total = sum(Probabilities.values())
    
    for key in Probabilities:
        Probabilities[key] = Probabilities[key]/total
        
    return Probabilities

# Input:  String Text and profile matrix Profile
# Output: Pr(Text, Profile)
def Pr(Text, Profile):
    p = 1
    l = len(Text)
    
    for i in range(l):
        p *= Profile[Text[i]][i]
    return p

# Input:  A dictionary Probabilities whose keys are k-mers and whose values are the probabilities of these kmers
# Output: A randomly chosen k-mer with respect to the values in Probabilities
def WeightedDie(Probabilities):
    total_prob = 0
    rand_num = random.uniform(0, 1)

    for kmer, prob in Probabilities.items():
        total_prob += prob
        if rand_num <= total_prob:
            return kmer
        
# Input:  A string Text, a profile matrix Profile, and an integer k
# Output: ProfileGeneratedString(Text, profile, k)
def ProfileGeneratedString(Text, profile, k):
    n = len(Text)
    probabilities = {} 
    
    for i in range(0,n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
        
    probabilities = Normalize(probabilities)
    
    return WeightedDie(probabilities)

def main():
    Dna = ['TGACGTTC', 'TAAGAGTT', 'GGACGAAA','CTGTTCGC']
    k = 3
    t = len(Dna)
    print(AdjustedRandomizedMotifSearch(Dna, k, t))    

main()