import argparse
import sys
import ast
import re

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

    last_idx = len(x) # last_idx
    rev_x = ''.join(reversed(x))+'$'
    Rsa = suffixArray(rev_x)
    Rbwt = ''.join([rev_x[(i + last_idx)%len(rev_x)] for i in Rsa])

    x += '$'
    sa = suffixArray(x) 
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
    
def process_file(fasta_dict: dict, filename: str) -> str:
    """Create a file containing the name of each string 
    with their suffix array, their bucket dict and their O table
    """    
    filename = filename.split('.')[0]
    file = filename + '_prepro.txt'

    with open(file, 'w') as f:
        final = ''
        for k, v in fasta_dict.items():
            final += '>' + k + '\n'
            Rsa, sa, C, RO, O = bwt_C_O(v)
            final += str(Rsa) + '\n'
            final += str(sa) + '\n'
            final += str(C) + '\n'
            final += str(RO) + '\n'
            final += str(O) + '\n'
        f.writelines(final)
    
    return file

def D_table(Rsa: list, C: dict, RO: dict, fastq: str, edit_limit) -> list:
        
    L, R = 0, len(Rsa)
    D = []
    edits = 0
    for char in fastq: # O(m)
        L = C[char] + RO[char][L] # O(1)
        R = C[char] + RO[char][R]
        if edits <= edit_limit:
            if L == R or char not in C:
                edits += 1
                D.append(edits)
                L, R = 0, len(Rsa)
            else:
                D.append(edits)
        else:
            D = []
    return D

def split_blocks(x: str) -> list[str]:
    """Split a string into blocks of equal character.

    Args:
        x (str): A string, but we sorta think it would be edits.

    Returns:
        list[str]: A list of blocks.

    >>> split_blocks('MDMMMMMMIMMMM')
    ['M', 'D', 'MMMMMM', 'I', 'MMMM']

    """
    # In any other language, this would likely not be the most efficient
    # approach to this, but since re.findall calls into C, it is faster
    # than implementing a more reasonable algorithm in pure Python.
    return [m[0] for m in re.findall(r"((.)\2*)", x)]

def edits_to_cigar(edits: str) -> str:
    """Encode a sequence of edits as a CIGAR.

    Args:
        edits (str): A sequence of edit operations

    Returns:
        str: The CIGAR encoding of edits.

    >>> edits_to_cigar('MDMMMMMMIMMMM')
    '1M1D6M1I4M'

    """
    result = split_blocks(edits)
    cigar = ''
    for C in result :
        cigar += str(len(C)) + C[0]
    return cigar

def approx_fm_search(genomename: str, sa: list, C: dict, O: dict, fastq: str, readname: str, D: list, edit_limit: int) -> str:

    
    L, R = 0, len(sa)
    j = 0
    p = ''.join(reversed(fastq))
        
    res = {} # {edits: (L, R, cigar)}
    queue = [(0, L, R, "", j)] # (edits, L, R, cigar, j)
    while queue:

        edits, L, R, cigar, j = queue.pop()
        
        if j == len(p) and edits <= edit_limit: # add result
            if edits in res:
                res[edits].append((L, R, cigar))
            else:
                res[edits] = [(L, R, cigar)]
            j = 0
            continue
        
        d = D[j] # Worth doing D or I or going on with a mismatch?

        if j < len(p) and d <= edit_limit: # O(m)
            # Update 29-11-22: add letters to cigar from the front since we move through p from the back

            if L != R: # Match
                queue.append((edits, C[p[j]] + O[p[j]][L], C[p[j]] + O[p[j]][R], 'M'+cigar, j+1))

            if edits < edit_limit: 
                # Mismatch
                for char in C.keys():
                    if char != '$' and char != p[j]:
                        queue.append((edits+1, C[char] + O[char][L], C[char] + O[char][R], 'M'+cigar, j+1))
            
                # Insertion
                queue.append((edits+1, L, R, 'I'+cigar, j+1))

                # Deletion
                if j != 0 and p[j] in C:
                    queue.append((edits+1, C[p[j]] + O[p[j]][L], C[p[j]] + O[p[j]][R], 'D'+cigar, j))
            
    final = []
    
    for edits, match in res.items():
        
        for L,R,cigar in match:
            for i in range(L,R):
                match = sa[i]+1
                if edit_limit==0:
                    text=str(len(p))+'M'
                else:
                    text=f'{edits_to_cigar(cigar)}'
                final.append('\t'.join([readname, genomename, str(match), text, fastq]))
    if final == []:
        return ''
    return '\n'.join(final)

def approximate_matching(prepro_file: str, fastq_dict: dict, edit_limit: int):
    # give room for muliple fasta sequences
    fastanames, Rsa_list, sa_list, C_list, O_list, RO_list = [], [], [], [], [], []

    with open(prepro_file, 'r') as f: 
        lines = f.readlines()
        linecounter = 0
        for line in lines:
            linecounter += 1
            if line.startswith('>'):
                linecounter = 0
                name = line[1:].strip()
                fastanames.append(name)
                continue
            l = ast.literal_eval(line.strip())
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
    
    for readname, read in fastq_dict.items():
        for i, sa in enumerate(sa_list):
            D = D_table(Rsa_list[i], C_list[i], RO_list[i], read, edit_limit)
            if D != []:
                simplesam = approx_fm_search(fastanames[i], sa_list[i], C_list[i], O_list[i], read, readname, D, edit_limit)
                if simplesam != '':
                    yield simplesam

# def main():

#     genome = 'src/fasta.fa'
#     with open (genome, 'r') as f:
#         lines = f.readlines()
#         fa = fasta_func(lines)

#     read = 'src/fastq.fq'
#     with open (read, 'r') as f:
#         lines = f.readlines()
#         fq = fastq_func(lines)
    
#     editlim = 1

#     filename = process_file(fa, genome)

#     sams = approximate_matching(filename, fq, editlim)
#     for s in sams:
#         print(s)

    

#     """ argparser = argparse.ArgumentParser(
#         description="FM-index exact pattern matching",
#         usage="\n\tfm -p genome\n\tfm genome reads"
#     )
#     argparser.add_argument(
#         "-p", action="store_true",
#         help="preprocess the genome."
#     )
#     argparser.add_argument(
#         "genome",
#         help="Simple-FASTA file containing the genome.",
#         type=argparse.FileType('r')
#     )
#     argparser.add_argument(
#         "reads", nargs="?",
#         help="Simple-FASTQ file containing the reads.",
#         type=argparse.FileType('r')
#     )
#     args = argparser.parse_args()

#     if args.p:
#         print(f"Preprocess {args.genome}")
#         fasta_dict = fasta_func(args.genome)

#         process_file(fasta_dict, args.genome.name)
#     else:
#         # here we need the optional argument reads
#         if args.reads is None:
#             argparser.print_help()
#             sys.exit(1)
#         prepro_file = args.genome.name.split('.')[0]+'_prepro.txt'
        
#         fastq_dict = fastq_func(args.reads)
#         print(fm_search(prepro_file, fastq_dict))
#          """
        


def main():
    argparser = argparse.ArgumentParser(
        description="Readmapper",
        usage="\n\treadmap -p genome\n\treadmap -d dist genome reads"
    )
    argparser.add_argument(
        "-p", action="store_true",
        help="preprocess the genome."
    )
    argparser.add_argument(
        "-d", type=int, metavar="integer",
        default=1, help="max edit distance."
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
        sams = approximate_matching(prepro_file, fastq_dict, args.d)
        for s in sams:
            print(s)


if __name__ == '__main__':
    main()
