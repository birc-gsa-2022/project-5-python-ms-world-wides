import argparse
import sys
import ast

def suffixArray(x: str) -> list:
    """Given x return suffix array SA(x). 
       We use Python's sorted function here 
       for simplycity, but we can do better.

    Args:
        x (str): input string.
    """  
    satups = sorted([(x[i:], i) for i in range(len(x))])
    
    return list(map(lambda t: t[1], satups))

def count_to_bucket(count: str) -> dict:
    '''
    >>> count_to_bucket("$iiiimppss")
    {'$': 0, 'i': 1, 'm': 5, 'p': 6, 's': 8}
    '''
    C = {}
    for i,c in enumerate(count):
        if c in C:
            continue
        else:
            C[c] = i

    return C

def calc_O(bwt: str, C: dict) -> dict:
    '''
    >>> calc_O('aaba$', { '$' : 0, 'a' : 1, 'b' : 4})
    {'$': [0, 0, 0, 0, 0, 1], 'a': [0, 1, 2, 2, 3, 3], 'b': [0, 0, 0, 1, 1, 1]}
    '''
    O = C.copy()
    for k in O.keys():
        O[k] = [0]
    
    for char in bwt:
        for k in O.keys():
            if k == char:
                O[k].append(O[k][-1]+1)
            else:
                O[k].append(O[k][-1])
    
    return O  
   
def bwt_C_O(x: str) -> tuple():
    """Calculates SA, C and O.

    Args:
        x (str): input string we want to find pattern in.

    Returns:
        SA, C and O
    """    
    
    last_idx = len(x)-1 # last_idx

    Rsa = suffixArray(reversed(x)+'$')
    Rbwt = ''.join([reversed(x)[(i + last_idx)%len(x)] for i in Rsa])

    sa = suffixArray(x+'$') 
    bwt = ''.join([x[(i + last_idx)%len(x)] for i in sa]) # Calculates bwt(x)

    count = ''.join([x[i] for i in sa]) # Sort x
    C = count_to_bucket(count) # dict with cumulative counts
    
    RO = calc_O(Rbwt, C)
    O = calc_O(bwt, C) # dict (table) with offsets
    return Rsa, sa, C, RO, O

def fasta_func(fastafile: str) -> dict:
    """Function that can take file or list of strings

    Returns:
        dict: return dictionary with fasta sequence coupled with its sequence name
    """    

    sequence = []
    name = ''
    fasta_dict = {}
    for line in fastafile:
        if type(line) == list:
            line = line[0]
        if line.startswith('>'):
            if name != '':
                fasta_dict[name] = ''.join(sequence)
                sequence = []
            name = line[1:].strip()
        else:
            sequence.append(line.strip())

    if name != '':
        fasta_dict[name] = ''.join(sequence)

    return fasta_dict

def fastq_func(fastqfile: str) -> dict: 
    read = []
    name = ''
    fastq_dict = {}
    for line in fastqfile:
        if line.startswith('@'):
            if name != '':
                fastq_dict[name] = ''.join(read)
                read = []
            name = line[1:].strip()
        else:
            read.append(line.strip())

    if name != '':
        fastq_dict[name] = ''.join(read)

    return fastq_dict
    
def process_file(fasta_dict: dict, filename: str):
    """Create a file containing the name of each string 
    with their suffix array, their bucket dict and their O table
    """    
    filename = filename.split('.')[0]
    file = filename + '_prepro.txt'

    with open(file, 'w') as f:
        final = ''
        for k, v in fasta_dict.items():
            final += '>' + k + '\n'
            Rsa, sa, C, , RO, O = bwt_C_O(v)
            final += str(Rsa) + '\n'
            final += str(sa) + '\n'
            final += str(C) + '\n'
            final += str(RO) + '\n'
            final += str(O) + '\n'
        f.writelines(final)

def approximate_matching(prepro_file: str, fastq_dict: dict, edit_limit: int):
    # give room for muliple fasta sequences
    fastanames, Rsa_list, sa_list, C_list, O_list, RO_list = [], [], [], [], []

    with open(prepro_file, 'r') as f:    
        lines = f.readlines()
        linecounter = 0
        for line in lines:
            linecounter += 1
            if line.startswith('>'):
                linecounter = 0
                name = line[1:].strip()
                fastanames.append(name)
            l = ast.literal_eval(line[1:].strip())
            if linecounter == 1:
                Rsa_list.append(l)
            elif linecounter == 2:
                sa_list.append(l)
            elif linecounter == 3:
                C_list.append(l) # convert back to dict
            elif linecounter == 4:
                RO_list.append(l) # convert back to dict
            elif linecounter == 5:
                O_list.append(l)
    
    res = []
    for readname, read in fastq_dict.items():
        for i, sa in enumerate(sa_list):
            D = D_table(Rsa_list[i], C_list[i], RO_list[i], read, edit_limit)
            if D != []:
                simplesam = approx_fm_search(fastanames[i], sa_list[i], C_list[i], O_list[i], read, readname, D)
                yield simplesam


def D_table(Rsa: list, C: dict, RO: dict, fastq: str, edit_limit) -> list:
        
    L, R = 0, len(Rsa)
    D = []
    edits = 0
    
    for char in fastq: # O(m)
        while edits <= edit_limit:
            if L == R or char not in C:
                edits += 1
                L, R = 0, len(Rsa)
            else:
                L = C[char] + RO[char][L] # O(1)
                R = C[char] + RO[char][R] # O(1)
            D.append(edits)
        else:
            D = []
    return D

def approx_fm_search(genomename: str, sa: list, C: dict, O: dict, fastq: str, readname: str, D: list, edit_limit) -> str:

    res = []
    sigma = list(C.keys())
    
    L, R = 0, len(sa)
    p = reversed(fastq)
    j = 0
    res = {} # {edits: (L, R, cigar)}
    queue = [(0, L, R, "", j)] # (edits, L, R, cigar, j)
    while queue.empty() == False:
        edits, L, R, cigar, j = queue.pop()
        if j < len(p) and edits <= edit_limit: # O(m)

            # USE D
            
            if L == R or p[j] not in C: # Do we need this? no, can do deletion instead, just not for the first character
                L = R # for char not in C
                break
            # ADD: Mismatch and match
            

            # Insertion
            queue.append((edits, L, R, cigar, j+1))
            # Deletion
            # ADD: No deletion at the beginning of the pattern!!!
            queue.append((edits, C[p[j]] + O[p[j]][L], C[[j]] + O[[j]][R], cigar, j+1))
        elif j == len(p)-1:
            res[edits] = (L, R, cigar)
        else:
            continue
            
    
    for edits, match in res.items():
        for L,R,cigar in match:
            for a in range(L,R):
                match = sa[a]+1
                res.append('\t'.join([readname, genomename, str(match), f'{cigar}M', fastq]))
    
    return '\n'.join(res)


def main():

    argparser = argparse.ArgumentParser(
        description="FM-index exact pattern matching",
        usage="\n\tfm -p genome\n\tfm genome reads"
    )
    argparser.add_argument(
        "-p", action="store_true",
        help="preprocess the genome."
    )
    argparser.add_argument(
        "genome",
        help="Simple-FASTA file containing the genome.",
        type=argparse.FileType('r')
    )
    argparser.add_argument(
        "reads", nargs="?",
        help="Simple-FASTQ file containing the reads.",
        type=argparse.FileType('r')
    )
    args = argparser.parse_args()

    if args.p:
        print(f"Preprocess {args.genome}")
        fasta_dict = fasta_func(args.genome)

        process_file(fasta_dict, args.genome.name)
    else:
        # here we need the optional argument reads
        if args.reads is None:
            argparser.print_help()
            sys.exit(1)
        prepro_file = args.genome.name.split('.')[0]+'_prepro.txt'
        
        fastq_dict = fastq_func(args.reads)
        print(fm_search(prepro_file, fastq_dict))
        
        


if __name__ == '__main__':
    main()
