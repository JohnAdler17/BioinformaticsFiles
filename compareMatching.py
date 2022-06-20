"""
This program allows a user to run both the naive pattern matching Algorithm
and Boyer-Moore on a specified pattern and text and reports the locations
of any matches and how many character comparisons were done by each algorithm.

Code for Boyer-Moore modified from that provided by Ben Langmead at
https://nbviewer.jupyter.org/github/BenLangmead/comp-genomics-class/blob/master/notebooks/CG_ZAlgorithm.ipynb
by Layla Oesper on September 12, 2019, and Catie Welsh on January 23, 2020
"""
import time

def naive(p, t):
    """
    A function that will do the naive algorithm for exact pattern matching.
    Input: p (string) - the pattern string
           t (string) - the text string to search through

    Output: returns the following two variables (in this order)
            occurences (list of ints) - a list of the indices in t in which p
                        appears.  Each index is the starting index of the match.
            compare (int) - the number of character comparisons computed during
                            the execution of the algorithm.  A comparison is when
                            a character from t is checked for equality with a
                            character from p.
    """
    
    occurrences = []
    compare = 0


    #TODO
    print("Length of p:", len(p))
    print("Length of t:", len(t))
    for i in range(len(t) - len(p) + 1):
        compare += 1
        j = 0

        while j < len(t):
            if t[i + j] != p[j - 1]:
                """compare += 1"""
                
                break
            j += 1
            compare += 1
        """print("J:", j)"""

        if j == len(p):
            occurrences.append(i)
        
    return occurrences, compare

def read_single_fasta(filename):
    """
    Reads in a single fasta file (only contains one sequence) and returns
    the sequence from the file.
    Input: filename (String) - the name of the fasta file to read.  The
            sequence may be on multiple lines in the fasta file.
    Output: sequence (String) - a single string representing the sequence in
            the fasta file (be sure to convert all characters to the same case)
    """
    #TODO
    file = open(filename, "r")

    if file.mode == "r":
        contents = file.read()
        

    return contents #Fix this once you've added your code above

def z_array(s):
    """ Use Z algorithm (Gusfield theorem 1.4.1) to preprocess s """
    assert len(s) > 1
    z = [len(s)] + [0] * (len(s)-1)

    # Initial comparison of s[1:] with prefix
    for i in range(1, len(s)):
        if s[i] == s[i-1]:
            z[1] += 1
        else:
            break

    r, l = 0, 0
    if z[1] > 0:
        r, l = z[1], 1

    for k in range(2, len(s)):
        assert z[k] == 0
        if k > r:
            # Case 1
            for i in range(k, len(s)):
                if s[i] == s[i-k]:
                    z[k] += 1
                else:
                    break
            r, l = k + z[k] - 1, k
        else:
            # Case 2
            # Calculate length of beta
            nbeta = r - k + 1
            zkp = z[k - l]
            if nbeta > zkp:
                # Case 2a: zkp wins
                z[k] = zkp
            else:
                # Case 2b: Compare characters just past r
                nmatch = 0
                for i in range(r+1, len(s)):
                    if s[i] == s[i - k]:
                        nmatch += 1
                    else:
                        break
                l, r = k, r + nmatch
                z[k] = r - k + 1
    return z


def n_array(s):
    """ Compile the N array (Gusfield theorem 2.2.2) from the Z array """
    return z_array(s[::-1])[::-1]


def big_l_prime_array(p, n):
    """ Compile L' array (Gusfield theorem 2.2.2) using p and N array.
        L'[i] = largest index j less than n such that N[j] = |P[i:]| """
    lp = [0] * len(p)
    for j in range(len(p)-1):
        i = len(p) - n[j]
        if i < len(p):
            lp[i] = j + 1
    return lp


def big_l_array(p, lp):
    """ Compile L array (Gusfield theorem 2.2.2) using p and L' array.
        L[i] = largest index j less than n such that N[j] >= |P[i:]| """
    l = [0] * len(p)
    l[1] = lp[1]
    for i in range(2, len(p)):
        l[i] = max(l[i-1], lp[i])
    return l


def small_l_prime_array(n):
    """ Compile lp' array (Gusfield theorem 2.2.4) using N array. """
    small_lp = [0] * len(n)
    for i in range(len(n)):
        if n[i] == i+1:  # prefix matching a suffix
            small_lp[len(n)-i-1] = i+1
    for i in range(len(n)-2, -1, -1):  # "smear" them out to the left
        if small_lp[i] == 0:
            small_lp[i] = small_lp[i+1]
    return small_lp


def good_suffix_table(p):
    """ Return tables needed to apply good suffix rule. """
    n = n_array(p)
    lp = big_l_prime_array(p, n)
    return lp, big_l_array(p, lp), small_l_prime_array(n)


def good_suffix_mismatch(i, big_l_prime, small_l_prime):
    """ Given a mismatch at offset i, and given L/L' and l' arrays,
        return amount to shift as determined by good suffix rule. """
    length = len(big_l_prime)
    assert i < length
    if i == length - 1:
        return 0
    i += 1  # i points to leftmost matching position of P
    if big_l_prime[i] > 0:
        return length - big_l_prime[i]
    return length - small_l_prime[i]


def good_suffix_match(small_l_prime):
    """ Given a full match of P to T, return amount to shift as
        determined by good suffix rule. """
    return len(small_l_prime) - small_l_prime[1]


def dense_bad_char_tab(p, amap):
    """ Given pattern string and list with ordered alphabet characters, create
        and return a dense bad character table.  Table is indexed by offset
        then by character. """
    tab = []
    nxt = [0] * len(amap)
    for i in range(0, len(p)):
        c = p[i]
        assert c in amap
        tab.append(nxt[:])
        nxt[amap[c]] = i+1
    return tab

class BoyerMoore(object):
    """ Encapsulates pattern and associated Boyer-Moore preprocessing. """

    def __init__(self, p, alphabet='ACGT'):
        # Create map from alphabet characters to integers
        self.amap = {alphabet[i]: i for i in range(len(alphabet))}
        # Make bad character rule table
        self.bad_char = dense_bad_char_tab(p, self.amap)
        # Create good suffix rule table
        _, self.big_l, self.small_l_prime = good_suffix_table(p)

    def bad_character_rule(self, i, c):
        """ Return # skips given by bad character rule at offset i """
        assert c in self.amap
        assert i < len(self.bad_char)
        ci = self.amap[c]
        return i - (self.bad_char[i][ci]-1)

    def good_suffix_rule(self, i):
        """ Given a mismatch at offset i, return amount to shift
            as determined by (weak) good suffix rule. """
        length = len(self.big_l)
        assert i < length
        if i == length - 1:
            return 0
        i += 1  # i points to leftmost matching position of P
        if self.big_l[i] > 0:
            return length - self.big_l[i]
        return length - self.small_l_prime[i]

    def match_skip(self):
        """ Return amount to shift in case where P matches T """
        return len(self.small_l_prime) - self.small_l_prime[1]

def boyer_moore(p, p_bm, t):
    """ Do Boyer-Moore matching """
    i = 0
    occurrences = []
    compare = 0
    while i < len(t) - len(p) + 1:
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            compare += 1
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences, compare

def main():
    p = "needle"
    t = "haystack needle haystack"
    occur, compare = naive(p,t)
    print("Pattern: ",p)
    print("Text:",t)
    print("\n----------------------Naive Algorithm--------------")
    print("Occurrence Position(s):",occur)
    print("Naive Comparisons:", compare)

    p_bm = BoyerMoore(p, alphabet='abcdefghijklmnopqrstuvwxyz ')
    bm_occur, bm_compares = boyer_moore(p, p_bm, t)
    print("\n----------------------Boyer-Moore Algorithm--------------")
    print("Occurence Position(s):",bm_occur)
    print("BM Comparisons:", bm_compares)

    print("\n----------------------------------------------\n")


    #Only uncomment the lines below once you have implemented and tested your
    #functions
    
    
    p = read_single_fasta('Alu.fasta')
    print(p)
    
    
    t = read_single_fasta('chr19.fasta')
    print(len(t))
    
    occur, compare = naive(p,t)
    print("\n----------------------Naive Algorithm--------------")
    print("Occurrence Position(s):",occur)
    print("Naive Comparisons:", compare)

    p_bm = BoyerMoore(p, alphabet='acgtACGTnN')
    bm_occur, bm_compares = boyer_moore(p, p_bm, t)
    print("\n----------------------Boyer-Moore Algorithm--------------")
    print("Occurence Position(s):",bm_occur)
    print("BM Comparisons:", bm_compares)
    


if __name__ == '__main__':
    main()
