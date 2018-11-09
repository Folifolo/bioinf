
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



def expand(Peptides):
    dict = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}
    res = set()
    if len(Peptides) != 0:
        for peptide in Peptides:
            for amino in dict:
                res.add(peptide+amino)
    #else:
    #    for amino in dict:
    #        res.add(amino)
    return res

def consist(pep_spect, Spectrum):
    temp = pep_spect.copy()
    Spec = Spectrum.copy()
    for elem in pep_spect:
        if elem in Spec:
            Spec.remove(elem)
            temp.remove(elem)
    if len(temp) == 0:
        return True
    return False

def three_one (Spectrum):
    dict = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}
    Peptides = {''}
    res = set()
    while len(Peptides) != 0:
        Peptides = expand(Peptides)
        peptides_temp = Peptides.copy()
        for peptide in Peptides:
            pep_spect = two_four(peptide)
            not_cycle = two_four(peptide, is_cycle=False)
            if pep_spect[-1] == Spectrum[-1]:
                if pep_spect == Spectrum:
                    temp = str('')
                    for elem in peptide:
                        temp += str(dict[elem]) +'-'
                    res.add(temp[:-1])
                peptides_temp.remove(peptide)
            elif not consist(not_cycle, Spectrum):
                peptides_temp.remove(peptide)
        Peptides = peptides_temp
    return  res

def three_two (Peptide, Spectrum, is_cycle = False):
    score = 0
    copy_spec = Spectrum.copy()
    if Peptide == '':
        theor_spectrum = ['0']
    else:
        theor_spectrum = two_four(Peptide, is_cycle = is_cycle)
    for element in theor_spectrum:
        if element in copy_spec:
            copy_spec.remove(element)
            score+=1
    return score

def Trim(Leaderboard, Spectrum, N):
    leader = dict()
    result = set()
    if len(Leaderboard) == 0:
        return set()
    for peptide in Leaderboard:
        leader[peptide] = three_two(peptide, Spectrum)
    sort_leader = sorted(leader.items(), key=lambda x: x[1], reverse = True)
    m = sort_leader[N-1][1]
    if sort_leader[N-1][1] != sort_leader[N][1]:
        res = (sort_leader[:N][:])
        for elem in res:
            result.add(elem[0])
        return result

    for i in range(len(Leaderboard)-N+1):
        if sort_leader[N-1][1] != sort_leader[N+i-1][1]:
            break

    res = (sort_leader[:N+i-1][:])
    for elem in res:
        result.add(elem[0])
    return result


def three_three (N, Spectrum):
    dict = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}
    Leaderboard = {''}
    LeaderPeptide = ''
    while len(Leaderboard) !=0:
        Leaderboard = expand(Leaderboard)
        leaderboard_temp = Leaderboard.copy()
        for peptide in Leaderboard:
            pep_spect = two_four(peptide)
            if pep_spect[-1] == Spectrum[-1]:
                if three_two(peptide, Spectrum) > three_two(LeaderPeptide, Spectrum):
                    LeaderPeptide = peptide
            elif int(pep_spect[-1]) > int(Spectrum[-1]):
                leaderboard_temp.remove(peptide)
        leaderboard_temp = Trim(leaderboard_temp, Spectrum, N)
        Leaderboard = leaderboard_temp

    result = []
    for char in LeaderPeptide:
        result.append(dict[char])
    return result


def difference(Pattern1, Pattern2):
    dif = 0
    if len(Pattern1) != len(Pattern2):
        return 10000000
    for i in range(len(Pattern1)):
        if Pattern1[i] != Pattern2[i]:
           dif += 1
    return dif

def pattern_generator(pattern, d):
    peptides = {""}
    for i in range(len(pattern)):
        peptides = expand_atcg(peptides)
    peptides_copy = peptides.copy()
    for elem in peptides:
        if difference(pattern, elem) > d:
            peptides_copy.remove(elem)
    return peptides_copy

def expand_atcg(Peptides):
    dict = {"A", "T", "C", "G"}
    res = set()
    if len(Peptides) != 0:
        for peptide in Peptides:
            for amino in dict:
                res.add(peptide+amino)

    return res

def four_one(Dna, k, d):
    Patterns = set()
    m = len(Dna[0])
    for i in range(len(Dna[0])-k+1):
        pattern = Dna[0][i:i+k]
        pattern_set = pattern_generator(pattern, d)
        for pattern_ in pattern_set:
            pattern_in_dna = False
            for string in Dna:#[1:]
                pattern_in_string = False
                for j in range(len(string)-k+1):
                    if difference(pattern_, string[j:j+k]) <= d:
                        pattern_in_string = True
                if pattern_in_string == False:
                    break
            else:
                Patterns.add(pattern_)

    return Patterns



def d(Pattern, Dna):
    res = 0
    for string in Dna:
        minimum = 10000000
        for i in range(len(string)-len(Pattern)+1):
            dif = difference(string[i:i+len(Pattern)], Pattern)
            if dif < minimum:
                minimum = dif
        res += minimum
    return res


def four_two(k, Dna):
    distance = 100000000
    tmp = "A"*k
    pattern_set = pattern_generator(tmp, k*2)
    for pattern in pattern_set:
        if distance > d(pattern, Dna):
            distance = d
            Median = pattern
    return Median
