def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
        else:
            continue
    return words

def FrequencyMap(Text, k):
    freq_dict = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        if Pattern in freq_dict:
            freq_dict[Pattern] += 1
        else:
            freq_dict[Pattern] = 1
    return freq_dict

k_mer = FrequentWords('CGGAGGACTCTAGGTAACGCTTATCAGGTCCATAGGACATTCA', 3)

print(k_mer)