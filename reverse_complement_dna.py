# Find the reverse complement of a DNA string.
# Input: A DNA string Pattern.
# Output: The reverse complement of Pattern.


def ReverseComplement(Pattern):   

    Pattern = Reverse(Pattern)
    Pattern = Complement(Pattern)
    
    return Pattern


def Reverse(Pattern):
    
    reverse_pattern = ""
    for i in range(len(Pattern)-1, -1, -1):
        reverse_pattern += Pattern[i]
    
    return reverse_pattern

def Complement(Pattern):

    comp_pattern = ""
    
    for nuc in Pattern:
        if nuc == "A":
            comp_pattern += "T"
        if nuc == "T":
            comp_pattern += "A"
        if nuc == "G":
            comp_pattern += "C"
        if nuc == "C":
            comp_pattern += "G"
            
    return comp_pattern
