
def one_one(P, G):
    count = 0
    for i in range(len(G)-len(P)):
        if G[i:i+len(P)] == P:
            count+=1
    return count

def one_two(Text, k):
    dict = {}
    for i in range(len(Text)-k):
        if Text[i:i+k] not in dict:
            dict[Text[i:i+k]]=1
        else:
            dict[Text[i:i+k]]+=1
    m = max(dict.values())
    res = [key for key, value in dict.items() if value == m]
    res = ' '.join(res)
    return res


def one_three(P):
    dict = {'A':'T','C':'G','G':'C','T':'A'}
    reverce = list(P)
    for i in range(len(P)):
        p = dict[P[i]]
        reverce[-i-1] = dict[P[i]]
    return reverce

def two_one(P):
    dict = {'AAA':'K','AAC':'N','AAG':'K','AAU':'N','ACA':'T','ACC':'T','ACG':'T','ACU':'T','AGA':'R','AGC':'S','AGG':'R','AGU':'S','AUA':'I','AUC':'I','AUG':'M','AUU':'I','CAA':'Q','CAC':'H','CAG':'Q','CAU':'H','CCA':'P','CCC':'P','CCG':'P','CCU':'P','CGA':'R','CGC':'R','CGG':'R','CGU':'R','CUA':'L','CUC':'L','CUG':'L','CUU':'L','GAA':'E','GAC':'D','GAG':'E','GAU':'D','GCA':'A','GCC':'A','GCG':'A','GCU':'A','GGA':'G','GGC':'G','GGG':'G','GGU':'G','GUA':'V','GUC':'V','GUG':'V','GUU':'V','UAA':'','UAC':'Y','UAG':'','UAU':'Y','UCA':'S','UCC':'S','UCG':'S','UCU':'S','UGA':'','UGC':'C','UGG':'W','UGU':'C','UUA':'L','UUC':'F','UUG':'L','UUU':'F'}
    res = []
    for i in range(len(P)//3-1):
        res.append(dict[P[i*3:i*3+3]])
    return res

def two_two (DNA, P):
    dict = {'AAA':'K','AAC':'N','AAG':'K','AAU':'N','ACA':'T','ACC':'T','ACG':'T','ACU':'T','AGA':'R','AGC':'S','AGG':'R','AGU':'S','AUA':'I','AUC':'I','AUG':'M','AUU':'I','CAA':'Q','CAC':'H','CAG':'Q','CAU':'H','CCA':'P','CCC':'P','CCG':'P','CCU':'P','CGA':'R','CGC':'R','CGG':'R','CGU':'R','CUA':'L','CUC':'L','CUG':'L','CUU':'L','GAA':'E','GAC':'D','GAG':'E','GAU':'D','GCA':'A','GCC':'A','GCG':'A','GCU':'A','GGA':'G','GGC':'G','GGG':'G','GGU':'G','GUA':'V','GUC':'V','GUG':'V','GUU':'V','UAA':'','UAC':'Y','UAG':'','UAU':'Y','UCA':'S','UCC':'S','UCG':'S','UCU':'S','UGA':'','UGC':'C','UGG':'W','UGU':'C','UUA':'L','UUC':'F','UUG':'L','UUU':'F'}
    complem = {'A':'T','C':'G','G':'C','T':'A'}
    res = []
    for i in range(len(DNA)-len(P)*3+1):
        tmp1 = DNA[i:i+len(P)*3]
        tmp2 = ''.join([complem[key] for key in tmp1[::-1]])
        tmp1 = tmp1.replace('T','U')
        tmp2 = tmp2.replace('T','U')
        pat1=[tmp1[i:i+3] for i in range(0, len(tmp1), 3)]
        pat2=[tmp2[i:i+3] for i in range(0, len(tmp2), 3)]
        pat1 = [dict[tripl] for tripl in pat1]
        pat2 = [dict[tripl] for tripl in pat2]
        if ''.join(pat1) == P or ''.join(pat2) == P:
            res.append(''.join(tmp1.replace('U','T')))
    return res

def two_four (P):
    dict = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}
    tmp = []
    tmp.append(P)
    res = []
    n = len(P)-1
    P += P
    for i in range(n+1):
        for j in range(n):
            tmp.append(P[i:i+j+1])

    res.append(0)
    for elem in tmp:
        temp = 0
        for letter in elem:
            temp += dict[letter]
        res.append(temp)

    res.sort()
    return ' '.join(str(x) for x in res)

print(two_two('ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA','MA'))
