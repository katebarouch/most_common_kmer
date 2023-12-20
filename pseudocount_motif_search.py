#pseudocount motif search

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
# HINT:   You need to use CountWithPseudocounts as a subroutine of ProfileWithPseudocounts
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

# Input:  String Text and profile matrix Profile
# Output: Pr(Text, Profile)
def Pr(Text, Profile):
    p = 1
    l = len(Text)
    
    for i in range(l):
        p *= Profile[Text[i]][i]
    return p

# Write your ProfileMostProbableKmer() function here.
# The profile matrix assumes that the first row corresponds to A, the second corresponds to C,
# the third corresponds to G, and the fourth corresponds to T.
# You should represent the profile matrix as a dictionary whose keys are 'A', 'C', 'G', and 'T' and whose values are lists of floats
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

def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = [] # output variable

    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
            
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs