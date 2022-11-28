from sequence_generator import (genome_sequence_generator,
                                random_sequence_generator)
from lin import lin_runner
from fm import process_file, fm_search, fasta_func, fastq_func

def test_special():
    empty_fasta = fasta_func(random_sequence_generator(0, 0, 'fasta','empty_fasta'))
    #empty_fastq = fastq_func(random_sequence_generator(0, 0, 'fastq','empty_fastq'))
    a1_fasta = fasta_func(random_sequence_generator(1, 0, 'fasta', 'a1_fasta', True))
    a1_fastq = fastq_func(random_sequence_generator(1, 0, 'fastq', 'a1_fastq', True))
    a10_fasta = fasta_func(random_sequence_generator(10, 0, 'fasta', 'a10_fasta', True))
    process_file(empty_fasta, 'empty_fasta')
    process_file(a1_fasta, 'a1_fasta')
    process_file(a10_fasta, 'a10_fasta')
    assert(sorted(fm_search('a10_fasta_prepro.txt',a1_fastq))==sorted(lin_runner(a10_fasta,a1_fastq)))
    #assert(sorted(fm_search('a10_fasta_prepro.txt',empty_fastq))==sorted(lin_runner(a10_fasta,empty_fastq)))
    assert(sorted(fm_search('a1_fasta_prepro.txt',a1_fastq))==sorted(lin_runner(a1_fasta,a1_fastq)))
    #assert(sorted(fm_search('empty_fasta_prepro.txt',empty_fastq))==sorted(lin_runner(empty_fasta,empty_fastq)))
    assert(sorted(fm_search('empty_fasta_prepro.txt',a1_fastq))==sorted(lin_runner(empty_fasta,a1_fastq)))


def test_random():
    fasta100 = fasta_func(random_sequence_generator(100, 8, 'fasta','fasta100'))
    process_file(fasta100, 'fasta100')
    fastq4 = fastq_func(random_sequence_generator(4, 0, 'fastq','fastq4'))

    assert(sorted(fm_search('fasta100_prepro.txt',fastq4))==sorted(lin_runner(fasta100,fastq4)))
    fastq100 = fastq_func(random_sequence_generator(100, 8, 'fastq','fastq100'))
    assert(sorted(fm_search('fasta100_prepro.txt',fastq100))==sorted(lin_runner(fasta100,fastq100)))

    
def test_genome():
    gen = fasta_func(genome_sequence_generator('src/sample_sequence.gz',
     500, 8, 'fasta', 'genome'))
    process_file(gen, 'genome') 
    
    fastq0 = fasta_func(genome_sequence_generator('src/sample_sequence.gz',
     0, 8, 'fastq', 'fastq0'))
    fastq1 = fasta_func(genome_sequence_generator('src/sample_sequence.gz',
     1, 8, 'fastq', 'fastq1'))
    fastq10 = fasta_func(genome_sequence_generator('src/sample_sequence.gz',
     10, 8, 'fastq', 'fastq10'))
    fastq25 = fasta_func(genome_sequence_generator('src/sample_sequence.gz',
     25, 8, 'fastq', 'fastq25'))
    assert(sorted(fm_search('genome_prepro.txt',fastq0))==sorted(lin_runner(gen,fastq0)))
    assert(sorted(fm_search('genome_prepro.txt',fastq1))==sorted(lin_runner(gen,fastq1)))
    assert(sorted(fm_search('genome_prepro.txt',fastq10))==sorted(lin_runner(gen,fastq10)))
    assert(sorted(fm_search('genome_prepro.txt',fastq25))==sorted(lin_runner(gen,fastq25)))