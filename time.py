import time
import matplotlib.pyplot as plt
from random import choice
import random
from fm import *
import gc
import math
import numpy as np


def process_file_time(sequence, seqname):
    '''Create a file containing the name of each string 
    with their suffix array, their bucket dict and their O table'''
    
    with open('/home/mathilde/Documents/Kandidat/GSA/Project/Project4_group/project-4-python-ms-world-wides/preprocess_time/{}._prepro.txt'.format(seqname), 'w') as f:
        final = ''
        final += '>' + seqname + '\n'
        sa, C, O = bwt_C_O(sequence)
        final += '#' + str(sa) + '\n'
        final += '+' + str(C) + '\n'
        final += '@' + str(O) + '\n'
        f.writelines(final)

def time_data(length, bases):
    """Function that generates DNA sequences and patterns used 
       to calculate the running time 

    Args:
        length (int): The length of the longest sequence generated.

    Returns:
        Lists, n and t: 
                    n is a list of the sequence length.
                    The t lists are the running times.
    """    
    n = []
    t_sa = []
    t_preprocess = []
    t_search = []
    
    for i in range(10, length):
        sequence = ''.join([choice(bases) for j in range(i)])
        start = random.randint(0, i - round(i/10))
        pattern = sequence[start:(start + round(i/10))]

        # Time to construct SA
        t0_sa = time.time()
        suffixArray(sequence)
        t1_sa = time.time()
        total_sa = t1_sa - t0_sa

        del sequence
        gc.collect()
        sequence = ''.join([choice(bases) for j in range(i)])

        # Time to preprocess file
        t0_preprocess = time.time()
        process_file_time(sequence, 'seq_' + str(i))
        t1_preprocess = time.time()
        total_preprocess = t1_preprocess - t0_preprocess

        prepro_file = 'preprocess_time/{}._prepro.txt'.format('seq_' + str(i))
        reads = {}
        read = 'read' + str(i) + '\n' + pattern
        reads[read] = pattern

        del sequence
        gc.collect()
        sequence = ''.join([choice(bases) for j in range(i)])


        # Time to search
        t0_search = time.time()
        approximate_matching(prepro_file, reads, edit_limit = 1)
        t1_search = time.time()
        total_search = t1_search - t0_search
       
        n.append(i)
        t_sa.append(total_sa)
        t_preprocess.append(total_preprocess)
        t_search.append(total_search)

        time.sleep(0.05)
    
    return n,t_sa, t_preprocess, t_search

def time_data_pattern(length, bases):
    """Function that generates DNA sequences and patterns used 
       to calculate the running time 

    Args:
        length (int): The length of the longest sequence generated.

    Returns:
        Lists, n and t: 
                    n is a list of the sequence length.
                    The t lists are the running times.
    """    
    n = []
    t_search = []
     
    sequence = ''.join([choice(bases) for j in range(1000)])
    process_file_time(sequence, 'multiple_patterns')
    prepro_file = '/home/mathilde/Documents/Kandidat/GSA/Project/Project5/project-5-python-ms-world-wides/preprocess_time/multiple_patterns._prepro.txt'
    
    for i in range(10, length):
        start = random.randint(0, 1000 - round(1000/10))
        pattern = sequence[start:(start + round(i/10))]
        reads = {}
        read = 'read' + str(i) + '\n' + pattern
        reads[read] = pattern

        # Time to fm search
        t0_search = time.time()
        approximate_matching(prepro_file, reads, edit_limit = 1)
        t1_search = time.time()
        total_search = t1_search - t0_search
       
        n.append(i)
        t_search.append(total_search)

        del pattern
        gc.collect()

        time.sleep(0.05)
    
    return n, t_search
bases = ['a','t','g','c']
worst_case = ['a']

n, t1, t2, t3 = time_data(1000, bases)
n1, t11, t21, t31 = time_data(1000, worst_case)


def plot_fig(n, t1, t2, t3, name, number):
    plt.figure(number)
    # Plotting all times simultaneously
    plt.scatter(n, t1, label='Build suffix array')
    plt.scatter(n, t2, label='Preprocess')
    plt.scatter(n, t3, label='FM search')
  
    
    # Naming figure, x-axis and y-axis
    plt.title('Time complexity')
    plt.ylabel('Time (s)')
    plt.xlabel('Size (n)')


    # Adding legends
    plt.legend()


    # Save figure in folder figs
    plt.savefig('/home/mathilde/Documents/Kandidat/GSA/Project/Project5/project-5-python-ms-world-wides/time/{}.png'.format(name))



