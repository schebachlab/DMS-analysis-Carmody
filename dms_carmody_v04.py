#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Some functions useful for analyzing next generation sequencing data obtained
in a deep mutational scanning experiment.

@author: Charles Kuntz :: cpkuntz@iu.edu
"""
# Import needed modules
import sys
import os
import re
import textwrap
import string
import random
import numpy as np
import pandas as pd
import scipy.stats as stat
import statsmodels.api as sm
import matplotlib.pyplot as plt
import seaborn as sns
import openpyxl
from datetime import datetime
from collections import Counter
from collections import defaultdict
from configparser import ConfigParser

# Calling this function requires no input, and simply returns
# a dictionary containing three-base codons as keys and 
# the corresponding amino acid or Stop signal as values

def genetic_code_func():
    genetic_code  = """
        TTT F      CTT L      ATT I      GTT V
        TTC F      CTC L      ATC I      GTC V
        TTA L      CTA L      ATA I      GTA V
        TTG L      CTG L      ATG M      GTG V
        TCT S      CCT P      ACT T      GCT A
        TCC S      CCC P      ACC T      GCC A
        TCA S      CCA P      ACA T      GCA A
        TCG S      CCG P      ACG T      GCG A
        TAT Y      CAT H      AAT N      GAT D
        TAC Y      CAC H      AAC N      GAC D
        TAA Stop   CAA Q      AAA K      GAA E
        TAG Stop   CAG Q      AAG K      GAG E
        TGT C      CGT R      AGT S      GGT G
        TGC C      CGC R      AGC S      GGC G
        TGA Stop   CGA R      AGA R      GGA G
        TGG W      CGG R      AGG R      GGG G
        """
    codon_finder = re.compile(r'[ATCG]{3}')
    amino_acid_finder = re.compile(r'\ \w{1}[\ |\n]|\Stop')
    codon_list = codon_finder.findall(genetic_code)
    amino_acid_list = [x.strip() for x in amino_acid_finder.findall(genetic_code)]
    genetic_code_dict = {}
    i = 0
    while i < len(codon_list):
        genetic_code_dict[codon_list[i]] = amino_acid_list[i]
        i += 1
    return genetic_code_dict

# Pass a string containing a path to a fastq file and this
# function reads the data into a list with the following organization:
# index 2 contains the sequence of nucleotides (~300 bases)
# index 4 contains the quality information for each base in ASCII+33 format
# index 5 begins the next read. 

def fastq_parser(path_to_fastq):
    filehandle = open(path_to_fastq, 'r')
    fastq_datastring = filehandle.read()
    delimiter = re.compile(r'[\S]+')
    fastq_datalist = delimiter.findall(fastq_datastring)
    return fastq_datalist


# This function reads a domain sequence (such as wild type H2 or H7)
# and returns a tuple of a clean representation of the domain,
# a list of three-base codons of the domain, and the length of the domain

def domain_processor(domain_sequence_raw):
    domain_sequence = re.sub("[^a-zA-Z]", "", domain_sequence_raw)
    domain_codons = textwrap.wrap(domain_sequence, 3)
    domain_length = len(domain_sequence)
    return domain_sequence, domain_codons, domain_length

def offset_sequence(raw_seq, offset):
    offset_seq = raw_seq[offset:]
    return offset_seq

# This function reads a list of lines from fastq files
# and the wild type domain sequence and returns the 
# index of the first base of the domain within the read

def seqstartfinder(data, domain_sequence):
    data_line = 2
    while data_line < len(data):
        sequence_read = data[data_line]
        if domain_sequence in sequence_read:
            sequence_start = sequence_read.find(domain_sequence)
            break
        else:
            data_line += 5
    return sequence_start


# These functions readsa list of lines from fastq files
# as well as the start position within the read of the
# domain, as well as the length of the domain, and
# a threshold for the average quality score over
# the domain and returns a list of trimmed domain reads
# as well as a list of the quality scores of each domain
# and the percent of reads that passed the quality threshold

# This is not a great way to filter reads. It calculates the average q score
# over the whole read and those that don't surpass a threshold are filtered.
# Problem is, reads can pass this test but still have probable errors,
# depending on the threshold. See:
# https://www.drive5.com/usearch/manual/readqualfiltering.html
def qualityfilter1(data, seqstart, seqlen, threshold=30):
    filtered_data = []
    avg_scores = []
    passed = 0
    failed = 0
    i = 2
    while i < len(data):
        domain_seq = data[i][seqstart:seqstart+seqlen]
        domain_qual = data[i+2][seqstart:seqstart+seqlen]
        if len(domain_seq) == seqlen and len(domain_qual) == seqlen:
            avg_score_domain = np.mean([ord(j) - 33 for j in domain_qual])
            pass_conditions = [avg_score_domain >= threshold,
                            'N' not in domain_seq]
            if all(pass_conditions):
                filtered_data.append(domain_seq)
                avg_scores.append(avg_score_domain)
                passed += 1
                i += 5
            else:
                failed += 1
                i += 5
        else:
            failed += 1
            i += 5
    percent_passed = (passed/(passed + failed))*100
    mean_scores = np.mean(avg_scores)
    return filtered_data, mean_scores, percent_passed

# This is better. It has a maximum allowable number of expected reads. Default = 1 expected error
# If a read has 1 or more expected errors, it fails
def qualityfilter2(data, seqstart, seqlen, errors=1):
    filtered_data = []
    avg_scores = []
    passed = 0
    failed = 0
    i = 2
    while i < len(data):
        domain_seq = data[i][seqstart:seqstart+seqlen]
        domain_qual = data[i+2][seqstart:seqstart+seqlen]
        if len(domain_seq) == seqlen and len(domain_qual) == seqlen:
            domain_qual_num = [ord(q_symbol) - 33 for q_symbol in domain_qual]
            prob_error_read = [10**-(base_q/10) for base_q in domain_qual_num]
            sum_prob_error_read = np.sum(prob_error_read)
            avg_score_domain = np.nanmean(domain_qual_num)
            pass_conditions = [sum_prob_error_read < errors,
                            'N' not in domain_seq]
            if all(pass_conditions):
                filtered_data.append(domain_seq)
                avg_scores.append(avg_score_domain)
                passed += 1
                i += 5
            else:
                failed += 1
                i += 5
        else:
            failed += 1
            i += 5
    percent_passed = (passed/(passed + failed))*100
    mean_scores = np.mean(avg_scores)
    return filtered_data, mean_scores, percent_passed

 # this is kind of redundant, but okay. It has an average q score threshold for the whole read 
 # plus a maximum number of expected errors over the read
def qualityfilter3(data_r1, data_r2, read_seqstart, read_seqlen, len_overlap, threshold=30, errors=1):
    #filtered_coord_r1 = []
    #filtered_coord_r2 = []
    
    filtered_data_r1 = {}
    filtered_data_r2 = {}
    
    avg_scores_r1 = []
    avg_scores_r2 = []
    
    passed_r1 = 0
    passed_r2 = 0
    
    failed_r1 = 0
    failed_r2 = 0
    
    i = 2
    while i < len(data_r1):
        domain_coord = data_r1[i-2]
        domain_seq = data_r1[i][read_seqstart:read_seqstart+read_seqlen]
        domain_qual = data_r1[i+2][read_seqstart:read_seqstart+read_seqlen]
        if len(domain_seq) == read_seqlen and len(domain_qual) == read_seqlen:
            domain_qual_num = [ord(q_symbol) - 33 for q_symbol in domain_qual]
            prob_error_read = [10**-(base_q/10) for base_q in domain_qual_num]
            sum_prob_error_read = np.sum(prob_error_read)
            avg_score_domain = np.mean(domain_qual_num)
            pass_conditions = [sum_prob_error_read < errors,
                            avg_score_domain >= threshold,
                            'N' not in domain_seq]
            if all(pass_conditions):
                #filtered_coord_r1.append(domain_coord)
                filtered_data_r1.update({domain_coord:domain_seq})
                avg_scores_r1.append(avg_score_domain)
                passed_r1 += 1
                i += 5
            else:
                failed_r1 += 1
                i += 5
        else:
            failed_r1 += 1
            i += 5

    len_nonoverlap = read_seqlen - len_overlap

    i = 2
    while i < len(data_r2):
        domain_coord = data_r2[i-2]
        domain_seq = data_r2[i][read_seqstart:read_seqstart+len_nonoverlap]
        domain_qual = data_r2[i+2][read_seqstart:read_seqstart+len_nonoverlap]
        if len(domain_seq) + len_overlap == read_seqlen and len(domain_qual) + len_overlap == read_seqlen:
            domain_qual_num = [ord(q_symbol) - 33 for q_symbol in domain_qual]
            prob_error_read = [10**-(base_q/10) for base_q in domain_qual_num]
            sum_prob_error_read = np.sum(prob_error_read)
            avg_score_domain = np.mean(domain_qual_num)
            pass_conditions = [sum_prob_error_read < errors,
                            avg_score_domain >= threshold,
                            'N' not in domain_seq]
            if all(pass_conditions):
                #filtered_coord_r2.append(domain_coord)
                filtered_data_r2.update({domain_coord:domain_seq})
                avg_scores_r2.append(avg_score_domain)
                passed_r2 += 1
                i += 5
            else:
                failed_r2 += 1
                i += 5
        else:
            failed_r2 += 1
            i += 5

    percent_passed_r1 = (passed_r1/(passed_r1 + failed_r1))*100
    mean_scores_r1 = np.mean(avg_scores_r1)
    
    percent_passed_r2 = (passed_r2/(passed_r2 + failed_r2))*100
    mean_scores_r2 = np.mean(avg_scores_r2)
    
    filtered_data_combined = defaultdict(list)
    for mydict in (filtered_data_r1, filtered_data_r2):
        for key, value in mydict.items():
            filtered_data_combined[key].append(value)
    
    filtered_data_dict = dict(filtered_data_combined)
    
    return filtered_data_dict, mean_scores_r1, mean_scores_r2, percent_passed_r1, percent_passed_r2


# this is similar to FASTX toolkit's quality filter. It filters reads by checking
# the quality score of each nucleotide in the read. If a certain percentage of the
# nucleotides in the read surpass a threshold quality score, the read passes
"""
def qualityfilter4(data, read_seqstart, read_seqlen, threshold=30, required_pass_rate=95):
    qual_test = lambda q_score, threshold: 1 if q_score >= threshold else 0
    filtered_data = []
    avg_scores = []
    passed = 0
    failed = 0
    i = 2
    while i < len(data):
        domain_seq = data[i][read_seqstart:read_seqstart+read_seqlen]
        domain_qual = data[i+2][read_seqstart:read_seqstart+read_seqlen]
        if len(domain_seq) == seqlen and len(domain_qual) == read_seqlen:
            domain_qual_num = [ord(q_symbol) - 33 for q_symbol in domain_qual]
            domain_qual_num_test = [qual_test(q_score, threshold) for q_score in domain_qual_num]
            domain_pass_rate = sum(domain_qual_num_test)/len(domain_qual_num_test) * 100
            avg_score_domain = np.mean(domain_qual_num)
            pass_conditions = [domain_pass_rate >= required_pass_rate,'N' not in domain_seq]
            if all(pass_conditions):
                filtered_data.append(domain_seq)
                avg_scores.append(avg_score_domain)
                passed += 1
                i += 5
            else:
                failed += 1
                i += 5
        else:
            failed += 1
            i += 5
    percent_passed = (passed/(passed + failed))*100
    mean_scores = np.mean(avg_scores)
    return filtered_data, mean_scores, percent_passed
"""
# This function reads a string of bases (a sequence) and
# returns the reverse complementary sequence

def get_overlap(s1_wt, s2_wt):
    n = 0
    while n < len(s1_wt):
        query = s1_wt[n:]
        if query in s2_wt:
            overlap_seq = query
            break
        else:
            n += 1
    s1_slice_idx = [s1_wt.find(overlap_seq), len(s1_wt)]
    s2_slice_idx = [s2_wt.find(overlap_seq), len(overlap_seq)]
    consensus_seq_wt = s1_wt+s2_wt[s2_slice_idx[1]:]
    return overlap_seq, s1_slice_idx, s2_slice_idx, consensus_seq_wt


def make_consensus_seq_check_overlap(filtered_data_dict, s1_slice, s2_slice, offset):
    filtered_data = []
    both_passed_qscore = 0
    for coord, seq_list in filtered_data_dict.items():
        if len(seq_list) == 2:
            both_passed_qscore += 1
            s1 = seq_list[1]
            s2 = reverse_complement_func(seq_list[2])
            if s1[s1_slice[0]:s1_slice[1]] == s2[s2_slice[0]:s2_slice[1]]:
                consensus_seq_raw = s1+s2[s2_slice[1]:]
                consensus_seq = consensus_seq_raw[offset:]
                filtered_data.append(consensus_seq)
    len_dict = len(filtered_data_dict)
    len_filtered = len(filtered_data)
    return filtered_data, both_passed_qscore, len_dict, len_filtered

def make_consensus_seq(filtered_data_dict, s1_slice, s2_slice, offset):
    filtered_data = []
    both_passed_qscore = 0
    for coord, seq_list in filtered_data_dict.items():
        if len(seq_list) == 2:
            both_passed_qscore += 1
            if len(seq_list[0]) > len(seq_list[1]):
                s1 = seq_list[0]
                s2 = reverse_complement_func(seq_list[1])
            if len(seq_list[0]) < len(seq_list[1]):
                s1 = seq_list[1]
                s2 = reverse_complement_func(seq_list[0])
            consensus_seq_raw = s1+s2
            consensus_seq = consensus_seq_raw[offset:]
            filtered_data.append(consensus_seq)
    len_dict = len(filtered_data_dict)
    len_filtered = len(filtered_data)
    return filtered_data, both_passed_qscore, len_dict, len_filtered


def split_consensus_at_slipsite(consensus_seq, slipsite_seq):
    slipsite_idx = consensus_seq.find(slipsite_seq)
    downstream_buffer = slipsite_idx % 3
    upstream_buffer = 3 - downstream_buffer
    upstream_domain_seq = consensus_seq[0:slipsite_idx+upstream_buffer]
    downstream_domain_seq = consensus_seq[slipsite_idx-downstream_buffer:]
    return upstream_domain_seq, downstream_domain_seq


def split_reads_at_slipsite(filtered_reads, consensus_seq, slipsite_seq):
    slipsite_idx = consensus_seq.find(slipsite_seq)
    downstream_buffer = slipsite_idx % 3
    upstream_buffer = 3 - downstream_buffer
    upstream_domain_seq = consensus_seq[0:slipsite_idx+upstream_buffer]
    downstream_domain_seq = consensus_seq[slipsite_idx-downstream_buffer:]
    upstream_reads = []
    downstream_reads = []
    wtcount = 0
    for read in filtered_reads:
        read_codons = textwrap.wrap(read, 3)
        mismatched_codons = 0
        for codon in range(len(read_codons)):
            if read_codons[codon] != domain_codons[codon]:
                mismatched_codons += 1
        if mismatched_codons == 0:
            wtcount += 1
        else:
            current_upstream_read = read[0:slipsite_idx+upstream_buffer]
            upstream_reads.append(current_upstream_read)
            current_downstream_read = read[slipsite_idx-downstream_buffer:]
            downstream_reads.append(current_downstream_read)
    return upstream_reads, downstream_reads, wtcount

def reverse_complement_func(forward_sequence):
    forward_alphabet = 'AGCTagct'
    revcomp_alphabet = 'TCGAtcga'
    reverse_sequence = forward_sequence[::-1]
    reverse_complement_sequence = reverse_sequence.translate({ord(x):
        y for (x, y) in zip(forward_alphabet, revcomp_alphabet)})
    return reverse_complement_sequence

def revcomp_data_func(filtered_data):
    reverse_complement_data = []
    idx = 0
    while idx < len(filtered_data):
        read = filtered_data[idx]
        revcomp_read = reverse_complement_func(read)
        reverse_complement_data.append(revcomp_read)
        idx += 1
    return reverse_complement_data
  

# This function reads the domain in codons and translates it to the amino
# acid sequence of the domain

def ribosome(domain_codons, genetic_code_dict):
    domain_aminoacids = []
    for codon in domain_codons:
        amino_acid = genetic_code_dict[codon]
        domain_aminoacids.append(amino_acid)
    return domain_aminoacids

def filter_multiple_codon_variants(quality_filtered_reads, domain_codons):
    master_variant_pool = []
    reads_with_n_mutated_codons = {0:0, 1:0, 2:0, 3:0, 4:0, 5:0, 6:0,
                                   7:0, 8:0, 9:0, 10:0, 11:0}
    for read in quality_filtered_reads:
        read_codons = textwrap.wrap(read, 3)
        mismatched_codons = 0
        for codon in range(len(read_codons)):
            if read_codons[codon] != domain_codons[codon]:
                mismatched_codons += 1
        if mismatched_codons == 0:
            master_variant_pool.append(read)
            reads_with_n_mutated_codons[0] += 1
        elif mismatched_codons == 1:
            master_variant_pool.append(read)
            reads_with_n_mutated_codons[1] += 1
        elif mismatched_codons == 2:
            reads_with_n_mutated_codons[2] += 1
        elif mismatched_codons == 3:
            reads_with_n_mutated_codons[3] += 1
        elif mismatched_codons == 4:
            reads_with_n_mutated_codons[4] += 1
        elif mismatched_codons == 5:
            reads_with_n_mutated_codons[5] += 1
        elif mismatched_codons == 6:
            reads_with_n_mutated_codons[6] += 1
        elif mismatched_codons == 7:
            reads_with_n_mutated_codons[7] += 1
        elif mismatched_codons == 8:
            reads_with_n_mutated_codons[8] += 1
        elif mismatched_codons == 9:
            reads_with_n_mutated_codons[9] += 1
        elif mismatched_codons == 10:
            reads_with_n_mutated_codons[10] += 1
        else:
            reads_with_n_mutated_codons[11] += 1 
    return master_variant_pool, reads_with_n_mutated_codons


def variant_identifier(rarefied_reads, domain_codons, codon_df_switch):
    df_index, df_columns, row_length, column_length, array_size = df_shape_func(domain_codons, codon_df_switch)
    genetic_code_dict = genetic_code_func()
    domain_translated = ribosome(domain_codons, genetic_code_dict)

    # Initialize array for counting single-codon variants
    single_codon_variant_array = np.zeros((row_length, column_length))
    
    # Initialize counters for wildtype and single-codon variant reads
    wildtype_reads_count = 0
    single_codon_variants = 0
    single_silent_mutations = 0
    single_nonsense_mutations = 0
    single_missense_mutations = 0
    
    wt_stops = domain_translated.count('Stop')

    # Identify the wild-type reads and single-codon variants
    for read in rarefied_reads:
        read_codons = textwrap.wrap(read, 3)
        mismatched_codons = 0
        for codon in range(len(read_codons)):
            if read_codons[codon] != domain_codons[codon]:
                mismatched_codons += 1
        if mismatched_codons == 0:
            wildtype_reads_count += 1
        else: #mismatched_codons == 1:
            single_codon_variants += 1
            translated_read = ribosome(read_codons, genetic_code_dict)
            read_stops = translated_read.count('Stop')
            if translated_read == domain_translated:
                single_silent_mutations += 1
            elif read_stops > wt_stops:
                single_nonsense_mutations += 1
            else: #translated_read != domain_translated and 'Stop' not in translated_read:
                single_missense_mutations += 1
            for codon in range(len(read_codons)):
                if read_codons[codon] != domain_codons[codon]:
                    if codon_df_switch == False:
                        mutation = translated_read[codon]
                    else:
                        mutation = read_codons[codon]
                    for idx in range(len(df_index)):
                        if df_index[idx] == mutation:
                            single_codon_variant_array[codon][idx] += 1
                        else:
                            continue

    # Make dataframe from single codon variant array
    single_codon_variant_array_transpose = single_codon_variant_array.T
    single_codon_variant_df = pd.DataFrame(single_codon_variant_array_transpose, index=df_index)
    single_codon_variant_df.columns = df_columns
    
    return (single_codon_variant_df, 
            wildtype_reads_count, 
            single_codon_variants, 
            single_silent_mutations, 
            single_nonsense_mutations, 
            single_missense_mutations)

# Pass this function a list of dataframes (one df per bin) and it calculates the total counts
# for all four bins before normalizing by wild type counts. It also returns the actual
# weighted intensities for each bin
# This function is also passed a list of wild type counts (scaled or not) and a list
# of bin intensities
# This is the central algorithm of our approach to measuring trafficking enrichment or hinderance
# as a function of amino acid subsitutions over the domain of interest

def weighted_intensity_func(var_df_polarity_bins, wtcount_bins, intensity_bins):
    # Process the variant count dataframes (populated by pandas dataframes in each bin):
    sum_over_bins__varcounts = sum(var_df_polarity_bins)
    weighted_varcounts = [intensity * var_df for intensity, var_df in zip(
            intensity_bins, var_df_polarity_bins)]
    sum_over_bins__weighted_varcounts = sum(weighted_varcounts)
    # Process the wild-type count list (populated by wtcounts in each bin):
    sum_over_bins__wtcounts = sum(wtcount_bins)
    weighted_wtcounts = [intensity * wtcount for intensity, wtcount in zip(
            intensity_bins, wtcount_bins)]
    sum_over_bins__weighted_wtcounts = sum(weighted_wtcounts)
    # Get the variant count weighted intensity and wtcount weighted intensity
    varcount_intensity_weighted = sum_over_bins__weighted_varcounts / sum_over_bins__varcounts
    wtcount_intensity_weighted = sum_over_bins__weighted_wtcounts / sum_over_bins__wtcounts
    # Get the final weighted intensity
    var_intensity_over_wt_intensity = varcount_intensity_weighted / wtcount_intensity_weighted
    # Return the last three objects calculated
    return (varcount_intensity_weighted, wtcount_intensity_weighted, 
            var_intensity_over_wt_intensity, sum_over_bins__varcounts)

# Function to eliminate variants represented in the experiment a number of reads lower than a user-defined cutoff
def import_file(path_to_file):
    file_obj = open(path_to_file)
    data_raw = []
    for line in file_obj:
        data_raw.append(line.split())
    file_obj.close()
    return data_raw


def df_shape_func(domain_codons, codon_df_switch):
    genetic_code_dict = genetic_code_func()
    domain_translated = ribosome(domain_codons, genetic_code_dict)
    # Get shape of dataframes
    polarity_rank = 'ILFVCMAWTYGSNHPQREKD'
    polindex = textwrap.wrap(polarity_rank, 1) + ['Stop']
    old_polindex_codons = ['ATT', 'ATC', 'ATA',                          # I
                       'TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG',     # L
                       'TTT', 'TTC',                                 # F
                       'GTT', 'GTC', 'GTA', 'GTG',                   # V
                       'TGT', 'TGC',                                 # C
                       'ATG',                                        # M
                       'GCT', 'GCC', 'GCA', 'GCG',                   # A
                       'TGG',                                        # W
                       'ACT', 'ACC', 'ACA', 'ACG',                   # T
                       'TAT', 'TAC',                                 # Y
                       'GGT', 'GGC', 'GGA', 'GGG',                   # G
                       'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC',     # S
                       'AAT', 'AAC',                                 # N
                       'CAT', 'CAC',                                 # H
                       'CCT', 'CCC', 'CCA', 'CCG',                   # P
                       'CAA', 'CAG',                                 # Q
                       'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG',     # R
                       'GAA', 'GAG',                                 # E
                       'AAA', 'AAG',                                 # K
                       'GAT', 'GAC',                                 # D
                       'TAA', 'TGA', 'TAG']                          # Stop
    
    polindex_codons = ['GAC',
                       'GAT',
                       'GAA',
                       'AAA',
                       'CCA',
                       'CCC',
                       'CCG',
                       'CCT',
                       'TAC',
                       'TAT',
                       'GAG',
                       'AAC',
                       'AAT',
                       'CAA',
                       'CAG',
                       'CGC',
                       'CGT',
                       'ACG',
                       'CGA',
                       'CGG',
                       'TCA',
                       'TCC',
                       'TCT',
                       'ACC',
                       'ACT',
                       'AGG',
                       'TGC',
                       'TGT',
                       'GGA',
                       'AAG',
                       'ACA',
                       'TTG',
                       'GGC',
                       'GGG',
                       'GGT',
                       'AGC',
                       'AGT',
                       'TCG',
                       'TTC',
                       'TTT',
                       'AGA',
                       'GTC',
                       'GTG',
                       'GTT',
                       'GTA',
                       'ATG',
                       'TGG',
                       'ATC',
                       'ATT',
                       'CTA',
                       'CTC',
                       'CTT',
                       'GCC',
                       'GCT',
                       'TTA',
                       'ATA',
                       'GCA',
                       'GCG',
                       'CTG',
                       'CAC',
                       'CAT',
                       'TAA',
                       'TAG',
                       'TGA']


    if codon_df_switch == True:
        df_index = polindex_codons
        df_columns = domain_codons
    else:
        df_index = polindex
        df_columns = domain_translated

    row_length = len(df_columns)
    column_length = len(df_index)
    array_size = row_length * column_length
    return df_index, df_columns, row_length, column_length, array_size


def filter_rare_variants_func(list_counts_bin_dfs, rarefaction_df, domain_codons, codon_df_switch, cutoff):
    df_index, df_columns, row_length, column_length, array_size = df_shape_func(domain_codons, codon_df_switch)
    num_bins = len(list_counts_bin_dfs)
    sum_over_bins_df = sum(list_counts_bin_dfs)
    sum_over_bins_array = sum_over_bins_df.values
    passing_positions = sum_over_bins_array >= cutoff
    rarefaction_array = rarefaction_df.values
    rarefaction_array_copy = np.copy(rarefaction_array)
    for i in range(column_length):
        for j in range(row_length):
            if passing_positions[i][j] == False:
                rarefaction_array_copy[i][j] = np.nan
            else:
                continue
    filtered_df = pd.DataFrame(rarefaction_array_copy, index=df_index)
    filtered_df.columns = df_columns
    return filtered_df


# Functions for doing the statistical analysis that takes two rarefaction sample dataframes and outputs
# a number of statistical parameters of interest (Pearson's r). Also returns a new dataframe with the 
# average of the values that were kept after doing the 95% confidence interval
# Finally, this function plots a scatterplot correlation between the two rarefactions in matplotlib

def dfs_to_pearson(df1, df2, domain_codons, codon_df_switch):
    df_index, df_columns, row_length, column_length, array_size = df_shape_func(domain_codons, codon_df_switch)
    array1 = df1.values
    array2 = df2.values
    array1_flat_raw = array1.flatten()
    array2_flat_raw = array2.flatten()
    array1_flat = np.copy(array1_flat_raw)
    array2_flat = np.copy(array2_flat_raw)
    not_nan_values = ~np.logical_or(np.isnan(array1_flat), np.isnan(array2_flat))
    for idx in range(array_size):
        if not_nan_values[idx] == True:
            array1_flat[idx] = array1_flat_raw[idx]
            array2_flat[idx] = array2_flat_raw[idx]
        else:
            array1_flat[idx] = -1
            array2_flat[idx] = -1
    array1_flat_copy = np.copy(array1_flat)
    array2_flat_copy = np.copy(array2_flat)
    array1_flat_dropnan = np.compress(not_nan_values, array1_flat_copy)
    array2_flat_dropnan = np.compress(not_nan_values, array2_flat_copy)
    pearson_r, p_value = stat.pearsonr(array1_flat_dropnan, array2_flat_dropnan)
    return array1_flat_raw, array2_flat_raw, array1_flat, array2_flat, pearson_r, p_value


def flat_arrays_to_pearson(array1_flat, array2_flat):
    not_nan_values = ~np.logical_or(np.isnan(array1_flat), np.isnan(array2_flat))
    array1_flat_copy = np.copy(array1_flat)
    array2_flat_copy = np.copy(array2_flat)
    array1_flat_dropnan = np.compress(not_nan_values, array1_flat_copy)
    array2_flat_dropnan = np.compress(not_nan_values, array2_flat_copy)
    pearson_r, p_value = stat.pearsonr(array1_flat_dropnan, array2_flat_dropnan)
    return pearson_r, p_value

def flat_arrays_to_df(array1_flat, array2_flat, domain_codons, codon_df_switch):
    df_index, df_columns, row_length, column_length, array_size = df_shape_func(domain_codons, codon_df_switch)
    array_flat_combined = np.column_stack((array1_flat, array2_flat))
    array_flat_avg = np.average(array_flat_combined, axis=1)
    array_avg = array_flat_avg.reshape(column_length, row_length)
    output_df = pd.DataFrame(array_avg, index=df_index)
    output_df.columns = df_columns
    return output_df
    

def scatterplot_arrays(array1, array2, x_axis, y_axis, filename, title):
    scatterplotname = filename+".png"
    scatterplot_fig = plt.figure()
    plt.plot(array1, array2, 'o', color='black');
    plt.xlabel(x_axis)
    plt.ylabel(y_axis)
    plt.title(title)
    scatterplot_fig.savefig(scatterplotname, dpi=400)


def outlier_elimination_func(array1, array2, domain_codons, codon_df_switch, domain_start):
    df_index, df_columns, row_length, column_length, array_size = df_shape_func(domain_codons, codon_df_switch)
    print("Searching for outliers using Bonferroni correction...")
    OLS_model = sm.OLS(array1, array2)
    OLS_model_results = OLS_model.fit()
    print(OLS_model_results.summary())
    outlier_data = OLS_model_results.outlier_test()
    outlier_indexes = []
    outliers_01 = []
    outliers_02 = []
    for idx in range(array_size):
        if outlier_data[idx][2] < 0.05:
            outlier_index_new = idx
            outlier_indexes.append(outlier_index_new)
            outlier_01_new = array1[idx]
            outlier_02_new = array2[idx]
            outliers_01.append(outlier_01_new)
            outliers_02.append(outlier_02_new)
    outliers = list(zip(outliers_01, outliers_02))
    for idx in range(len(outlier_indexes)):
        outlier_tuple = outliers[idx]
        array_idx = outlier_indexes[idx]
        polindex_position = array_idx // row_length
        domain_position = array_idx % row_length
        domain_aa = df_columns[domain_position]
        polindex_aa = df_index[polindex_position]
        print("Outlier eliminated: ", domain_aa, domain_position + domain_start, "->", polindex_aa, "::", outlier_tuple)
    list_flat_01_filtered = []
    list_flat_02_filtered = []
    for idx in range(array_size):
        if idx in outlier_indexes:
            list_flat_01_filtered.append(-1)
            list_flat_02_filtered.append(-1)
        else:
            list_flat_01_filtered.append(array1[idx])
            list_flat_02_filtered.append(array2[idx])
    array1_filtered = np.array(list_flat_01_filtered)
    array2_filtered = np.array(list_flat_02_filtered)
    array1_output = np.copy(array1_filtered)
    array2_output = np.copy(array2_filtered)
    for idx in range(array_size):
        if array1_filtered[idx] == -1:
            array1_output[idx] = np.nan
            array2_output[idx] = np.nan
        else:
            array1_output[idx] = array1_filtered[idx]
            array2_output[idx] = array2_filtered[idx]
    return array1_output, array2_output

def rarefaction_corr_func(df1, df2, domain_codons, jobname, outlier_switch, codon_df_switch, domain_start=1):
    df_index, df_columns, row_length, column_length, array_size = df_shape_func(domain_codons, codon_df_switch)
    # Get flatted arrays from dataframes and calculate the pearson's r before outlier elimination
    (array1_flat_raw, 
     array2_flat_raw, 
     array1_flat_input, 
     array2_flat_input, 
     pearson_r_input, 
     p_value_input) = dfs_to_pearson(df1, df2, domain_codons, codon_df_switch)
    str_pearson_r1 = "Pearson's r, no outlier removal: {:.4f}".format(pearson_r_input)
    print(str_pearson_r1)
    print("p value, no outlier removal: ", p_value_input)
    # Save scatterplots before outlier removal:
    filename1 = jobname+"_no_outlier_removal"
    scatter_title1 = "Correlation of two rarefaction samplings"+"\n"+str_pearson_r1
    scatterplot_arrays(array1_flat_raw, array2_flat_raw, 'Rarefaction 1', 'Rarefaction 2', filename1, scatter_title1)
    df_avg = flat_arrays_to_df(array1_flat_raw, array2_flat_raw, domain_codons, codon_df_switch)
    output_list = []
    output_list.append(df_avg)
    if outlier_switch == True:
        # Eliminate outliers after correlating two replicates:
        (array1_flat_output, 
        array2_flat_output) = outlier_elimination_func(array1_flat_input, array2_flat_input, domain_codons, codon_df_switch, domain_start)
        # Get Pearson's r after outlier removal
        pearson_r_output, p_value_output = flat_arrays_to_pearson(array1_flat_output, array2_flat_output)
        str_pearson_r2 = "Pearson's r after outlier removal: {:.4f}".format(pearson_r_output)
        print(str_pearson_r2)
        print("p value after outlier removal: ", p_value_output)
        # Save scatterplot after outlier removal:
        filename2 = jobname+"_after_outlier_removal"
        scatter_title2 = "Correlation of two replicates"+"\n"+str_pearson_r2
        scatterplot_arrays(array1_flat_output, array2_flat_output, 'Rarefaction 1', 'Rarefaction 2', filename2, scatter_title2)
        df_avg_outliers_removed = flat_arrays_to_df(array1_flat_output, array2_flat_output, domain_codons, codon_df_switch)
        output_list.append(df_avg_outliers_removed)
    else:
        print("Outliers after rarefaction correlation not removed.")
        empty_array = np.zeros((column_length, row_length))
        empty_df = pd.DataFrame(empty_array, index=df_index)
        empty_df.columns = df_columns
        output_list.append(empty_df)
    return output_list



# This function is passed two dataframes and a list of the amino acid sequence of the domain
# and the function returns the element-wise average of the two dataframes

def average_two_dfs_func(pd_df_01, pd_df_02, domain_codons, codon_df_switch):
    df_index, df_columns, row_length, column_length, array_size = df_shape_func(domain_codons, codon_df_switch)

    np_array_01 = pd_df_01.values
    np_array_02 = pd_df_02.values
    np_array_flat_01 = np_array_01.flatten()
    np_array_flat_02 = np_array_02.flatten()

    np_array_flat_combined = np.column_stack((np_array_flat_01, np_array_flat_02))
    np_array_flat_average = np.average(np_array_flat_combined, axis=1)
    np_array_average = np_array_flat_average.reshape(column_length, row_length)
    output_array = np_array_average
    output_df = pd.DataFrame(output_array, index=df_index)
    output_df.columns = df_columns  
    return output_df

def nucleotide_variant_identifier(reads, orf_codons, orf_seq):
    columns = list(orf_seq)
    rows = ['A', 'C', 'G', 'T']
    column_length = len(rows)
    row_length = len(columns)
    genetic_code_dict = genetic_code_func()
    translated_wildtype = ribosome(orf_codons, genetic_code_dict)

    # Initialize arrays for counting all SNPs, DNPs, and TNPs
    all_count_array = np.zeros((row_length, column_length))
    all_snp_count_array = np.zeros((row_length, column_length))
    all_dnp_count_array = np.zeros((row_length, column_length))
    all_tnp_count_array = np.zeros((row_length, column_length))

    # Initialize arrays for counting silent SNPs, DNPs, and TNPs
    syn_count_array = np.zeros((row_length, column_length))
    syn_snp_count_array = np.zeros((row_length, column_length))
    syn_dnp_count_array = np.zeros((row_length, column_length))
    syn_tnp_count_array = np.zeros((row_length, column_length))

    # Initialize arrays for counting missense SNPs, DNPs, and TNPs
    nonsyn_count_array = np.zeros((row_length, column_length))
    nonsyn_snp_count_array = np.zeros((row_length, column_length))
    nonsyn_dnp_count_array = np.zeros((row_length, column_length))
    nonsyn_tnp_count_array = np.zeros((row_length, column_length))

    # Initialize arrays for counting nonsense SNPs, DNPs, and TNPs
    nonsense_count_array = np.zeros((row_length, column_length))
    nonsense_snp_count_array = np.zeros((row_length, column_length))
    nonsense_dnp_count_array = np.zeros((row_length, column_length))
    nonsense_tnp_count_array = np.zeros((row_length, column_length))

    wildtype_reads = 0
    single_codon_variants = 0
    single_silent_reads = 0
    single_missense_reads = 0
    single_nonsense_reads = 0

    SNPs = 0
    DNPs = 0
    TNPs = 0

    missense_SNPs = 0
    missense_DNPs = 0
    missense_TNPs = 0

    silent_SNPs = 0
    silent_DNPs = 0
    silent_TNPs = 0

    nonsense_SNPs = 0
    nonsense_DNPs = 0
    nonsense_TNPs = 0

    for read in reads:
        read_codons = textwrap.wrap(read, 3)
        mismatched_codons = 0
        for codon in range(len(read_codons)):
            if read_codons[codon] != orf_codons[codon]:
                mismatched_codons += 1
        if mismatched_codons == 0:
            wildtype_reads += 1
        if mismatched_codons == 1:
            single_codon_variants += 1
            translated_read = ribosome(read_codons, genetic_code_dict)
            for base in range(len(read)):
                if read[base] != orf_seq[base]:
                    if read[base] == 'A':
                        all_count_array[base][0] += 1
                    if read[base] == 'C':
                        all_count_array[base][1] += 1
                    if read[base] == 'G':
                        all_count_array[base][2] += 1
                    if read[base] == 'T':
                        all_count_array[base][3] += 1
            for codon in range(len(read_codons)):
                if read_codons[codon] != orf_codons[codon]:
                    nucleotide_substitutions = 0
                    for nucleotide in range(len(read_codons[codon])):
                        if read_codons[codon][nucleotide] != orf_codons[codon][nucleotide]:
                            nucleotide_substitutions += 1
                    if nucleotide_substitutions == 1:
                        SNPs += 1
                        for base in range(len(read)):
                            if read[base] != orf_seq[base]:
                                if read[base] == 'A':
                                    all_snp_count_array[base][0] += 1
                                if read[base] == 'C':
                                    all_snp_count_array[base][1] += 1
                                if read[base] == 'G':
                                    all_snp_count_array[base][2] += 1
                                if read[base] == 'T':
                                    all_snp_count_array[base][3] += 1
                    if nucleotide_substitutions == 2:
                        DNPs += 1
                        for base in range(len(read)):
                            if read[base] != orf_seq[base]:
                                if read[base] == 'A':
                                    all_dnp_count_array[base][0] += 1
                                if read[base] == 'C':
                                    all_dnp_count_array[base][1] += 1
                                if read[base] == 'G':
                                    all_dnp_count_array[base][2] += 1
                                if read[base] == 'T':
                                    all_dnp_count_array[base][3] += 1
                    if nucleotide_substitutions == 3:
                        TNPs += 1
                        for base in range(len(read)):
                            if read[base] != orf_seq[base]:
                                if read[base] == 'A':
                                    all_tnp_count_array[base][0] += 1
                                if read[base] == 'C':
                                    all_tnp_count_array[base][1] += 1
                                if read[base] == 'G':
                                    all_tnp_count_array[base][2] += 1
                                if read[base] == 'T':
                                    all_tnp_count_array[base][3] += 1
            if translated_read == translated_wildtype:
                single_silent_reads += 1
                for base in range(len(read)):
                    if read[base] != orf_seq[base]:
                        if read[base] == 'A':
                            syn_count_array[base][0] += 1
                        if read[base] == 'C':
                            syn_count_array[base][1] += 1
                        if read[base] == 'G':
                            syn_count_array[base][2] += 1
                        if read[base] == 'T':
                            syn_count_array[base][3] += 1
                for codon in range(len(read_codons)):
                    if read_codons[codon] != orf_codons[codon]:
                        nucleotide_substitutions = 0
                        for nucleotide in range(len(read_codons[codon])):
                            if read_codons[codon][nucleotide] != orf_codons[codon][nucleotide]:
                                nucleotide_substitutions += 1
                        if nucleotide_substitutions == 1:
                            silent_SNPs += 1
                            for base in range(len(read)):
                                if read[base] != orf_seq[base]:
                                    if read[base] == 'A':
                                        syn_snp_count_array[base][0] += 1
                                    if read[base] == 'C':
                                        syn_snp_count_array[base][1] += 1
                                    if read[base] == 'G':
                                        syn_snp_count_array[base][2] += 1
                                    if read[base] == 'T':
                                        syn_snp_count_array[base][3] += 1
                        if nucleotide_substitutions == 2:
                            silent_DNPs += 1
                            for base in range(len(read)):
                                if read[base] != orf_seq[base]:
                                    if read[base] == 'A':
                                        syn_dnp_count_array[base][0] += 1
                                    if read[base] == 'C':
                                        syn_dnp_count_array[base][1] += 1
                                    if read[base] == 'G':
                                        syn_dnp_count_array[base][2] += 1
                                    if read[base] == 'T':
                                        syn_dnp_count_array[base][3] += 1
                        if nucleotide_substitutions == 3:
                            silent_TNPs += 1
                            for base in range(len(read)):
                                if read[base] != orf_seq[base]:
                                    if read[base] == 'A':
                                        syn_tnp_count_array[base][0] += 1
                                    if read[base] == 'C':
                                        syn_tnp_count_array[base][1] += 1
                                    if read[base] == 'G':
                                        syn_tnp_count_array[base][2] += 1
                                    if read[base] == 'T':
                                        syn_tnp_count_array[base][3] += 1
            elif 'Stop' in translated_read:
                single_nonsense_reads += 1
                for base in range(len(read)):
                    if read[base] != orf_seq[base]:
                        if read[base] == 'A':
                            nonsense_count_array[base][0] += 1
                        if read[base] == 'C':
                            nonsense_count_array[base][1] += 1
                        if read[base] == 'G':
                            nonsense_count_array[base][2] += 1
                        if read[base] == 'T':
                            nonsense_count_array[base][3] += 1
                for codon in range(len(read_codons)):
                    if read_codons[codon] != orf_codons[codon]:
                        nucleotide_substitutions = 0
                        for nucleotide in range(len(read_codons[codon])):
                            if read_codons[codon][nucleotide] != orf_codons[codon][nucleotide]:
                                nucleotide_substitutions += 1
                        if nucleotide_substitutions == 1:
                            nonsense_SNPs += 1
                            for base in range(len(read)):
                                if read[base] != orf_seq[base]:
                                    if read[base] == 'A':
                                        nonsense_snp_count_array[base][0] += 1
                                    if read[base] == 'C':
                                        nonsense_snp_count_array[base][1] += 1
                                    if read[base] == 'G':
                                        nonsense_snp_count_array[base][2] += 1
                                    if read[base] == 'T':
                                        nonsense_snp_count_array[base][3] += 1
                        if nucleotide_substitutions == 2:
                            nonsense_DNPs += 1
                            for base in range(len(read)):
                                if read[base] != orf_seq[base]:
                                    if read[base] == 'A':
                                        nonsense_dnp_count_array[base][0] += 1
                                    if read[base] == 'C':
                                        nonsense_dnp_count_array[base][1] += 1
                                    if read[base] == 'G':
                                        nonsense_dnp_count_array[base][2] += 1
                                    if read[base] == 'T':
                                        nonsense_dnp_count_array[base][3] += 1
                        if nucleotide_substitutions == 3:
                            nonsense_TNPs += 1
                            for base in range(len(read)):
                                if read[base] != orf_seq[base]:
                                    if read[base] == 'A':
                                        nonsense_tnp_count_array[base][0] += 1
                                    if read[base] == 'C':
                                        nonsense_tnp_count_array[base][1] += 1
                                    if read[base] == 'G':
                                        nonsense_tnp_count_array[base][2] += 1
                                    if read[base] == 'T':
                                        nonsense_tnp_count_array[base][3] += 1
            else:
                single_missense_reads += 1
                for base in range(len(read)):
                    if read[base] != orf_seq[base]:
                        if read[base] == 'A':
                            nonsyn_count_array[base][0] += 1
                        if read[base] == 'C':
                            nonsyn_count_array[base][1] += 1
                        if read[base] == 'G':
                            nonsyn_count_array[base][2] += 1
                        if read[base] == 'T':
                            nonsyn_count_array[base][3] += 1
                for codon in range(len(read_codons)):
                    if read_codons[codon] != orf_codons[codon]:
                        nucleotide_substitutions = 0
                        for nucleotide in range(len(read_codons[codon])):
                            if read_codons[codon][nucleotide] != orf_codons[codon][nucleotide]:
                                nucleotide_substitutions += 1
                        if nucleotide_substitutions == 1:
                            missense_SNPs += 1
                            for base in range(len(read)):
                                if read[base] != orf_seq[base]:
                                    if read[base] == 'A':
                                        nonsyn_snp_count_array[base][0] += 1
                                    if read[base] == 'C':
                                        nonsyn_snp_count_array[base][1] += 1
                                    if read[base] == 'G':
                                        nonsyn_snp_count_array[base][2] += 1
                                    if read[base] == 'T':
                                        nonsyn_snp_count_array[base][3] += 1
                        if nucleotide_substitutions == 2:
                            missense_DNPs += 1
                            for base in range(len(read)):
                                if read[base] != orf_seq[base]:
                                    if read[base] == 'A':
                                        nonsyn_dnp_count_array[base][0] += 1
                                    if read[base] == 'C':
                                        nonsyn_dnp_count_array[base][1] += 1
                                    if read[base] == 'G':
                                        nonsyn_dnp_count_array[base][2] += 1
                                    if read[base] == 'T':
                                        nonsyn_dnp_count_array[base][3] += 1
                        if nucleotide_substitutions == 3:
                            missense_TNPs += 1
                            for base in range(len(read)):
                                if read[base] != orf_seq[base]:
                                    if read[base] == 'A':
                                        nonsyn_tnp_count_array[base][0] += 1
                                    if read[base] == 'C':
                                        nonsyn_tnp_count_array[base][1] += 1
                                    if read[base] == 'G':
                                        nonsyn_tnp_count_array[base][2] += 1
                                    if read[base] == 'T':
                                        nonsyn_tnp_count_array[base][3] += 1

    # Make dataframe from all_count_array
    all_count_array_transpose = all_count_array.T
    all_count_df = pd.DataFrame(all_count_array_transpose, index=rows)
    all_count_df.columns = columns
    # Make dataframe from all_snp_count_array
    all_snp_count_array_transpose = all_snp_count_array.T
    all_snp_count_df = pd.DataFrame(all_snp_count_array_transpose, index=rows)
    all_snp_count_df.columns = columns
    # Make dataframe from all_dnp_count_array
    all_dnp_count_array_transpose = all_dnp_count_array.T
    all_dnp_count_df = pd.DataFrame(all_dnp_count_array_transpose, index=rows)
    all_dnp_count_df.columns = columns
    # Make dataframe from all_tnp_count_array
    all_tnp_count_array_transpose = all_tnp_count_array.T
    all_tnp_count_df = pd.DataFrame(all_tnp_count_array_transpose, index=rows)
    all_tnp_count_df.columns = columns

    # Make dataframe from syn_count_array
    syn_count_array_transpose = syn_count_array.T
    syn_count_df = pd.DataFrame(syn_count_array_transpose, index=rows)
    syn_count_df.columns = columns
    # Make dataframe from syn_snp_count_array
    syn_snp_count_array_transpose = syn_snp_count_array.T
    syn_snp_count_df = pd.DataFrame(syn_snp_count_array_transpose, index=rows)
    syn_snp_count_df.columns = columns
    # Make dataframe from syn_dnp_count_array
    syn_dnp_count_array_transpose = syn_dnp_count_array.T
    syn_dnp_count_df = pd.DataFrame(syn_dnp_count_array_transpose, index=rows)
    syn_dnp_count_df.columns = columns
    # Make dataframe from syn_tnp_count_array
    syn_tnp_count_array_transpose = syn_tnp_count_array.T
    syn_tnp_count_df = pd.DataFrame(syn_tnp_count_array_transpose, index=rows)
    syn_tnp_count_df.columns = columns

    # Make dataframe from nonsyn_count_array
    nonsyn_count_array_transpose = nonsyn_count_array.T
    nonsyn_count_df = pd.DataFrame(nonsyn_count_array_transpose, index=rows)
    nonsyn_count_df.columns = columns
    # Make dataframe from nonsyn_snp_count_array
    nonsyn_snp_count_array_transpose = nonsyn_snp_count_array.T
    nonsyn_snp_count_df = pd.DataFrame(nonsyn_snp_count_array_transpose, index=rows)
    nonsyn_snp_count_df.columns = columns
    # Make dataframe from nonsyn_dnp_count_array
    nonsyn_dnp_count_array_transpose = nonsyn_dnp_count_array.T
    nonsyn_dnp_count_df = pd.DataFrame(nonsyn_dnp_count_array_transpose, index=rows)
    nonsyn_dnp_count_df.columns = columns
    # Make dataframe from nonsyn_tnp_count_array
    nonsyn_tnp_count_array_transpose = nonsyn_tnp_count_array.T
    nonsyn_tnp_count_df = pd.DataFrame(nonsyn_tnp_count_array_transpose, index=rows)
    nonsyn_tnp_count_df.columns = columns

    # Make dataframe from nonsense_count_array
    nonsense_count_array_transpose = nonsense_count_array.T
    nonsense_count_df = pd.DataFrame(nonsense_count_array_transpose, index=rows)
    nonsense_count_df.columns = columns
    # Make dataframe from nonsense_snp_count_array
    nonsense_snp_count_array_transpose = nonsense_snp_count_array.T
    nonsense_snp_count_df = pd.DataFrame(nonsense_snp_count_array_transpose, index=rows)
    nonsense_snp_count_df.columns = columns
    # Make dataframe from nonsense_dnp_count_array
    nonsense_dnp_count_array_transpose = nonsense_dnp_count_array.T
    nonsense_dnp_count_df = pd.DataFrame(nonsense_dnp_count_array_transpose, index=rows)
    nonsense_dnp_count_df.columns = columns
    # Make dataframe from nonsense_tnp_count_array
    nonsense_tnp_count_array_transpose = nonsense_tnp_count_array.T
    nonsense_tnp_count_df = pd.DataFrame(nonsense_tnp_count_array_transpose, index=rows)
    nonsense_tnp_count_df.columns = columns

    return (wildtype_reads,
            single_codon_variants,
            single_silent_reads,
            single_missense_reads,
            single_nonsense_reads,
            SNPs,
            DNPs,
            TNPs,
            missense_SNPs,
            missense_DNPs,
            missense_TNPs,
            silent_SNPs,
            silent_DNPs,
            silent_TNPs,
            nonsense_SNPs,
            nonsense_DNPs,
            nonsense_TNPs,
            all_count_df,
            all_snp_count_df,
            all_dnp_count_df,
            all_tnp_count_df,
            syn_count_df,
            syn_snp_count_df,
            syn_dnp_count_df,
            syn_tnp_count_df,
            nonsyn_count_df,
            nonsyn_snp_count_df,
            nonsyn_dnp_count_df,
            nonsyn_tnp_count_df,
            nonsense_count_df,
            nonsense_snp_count_df,
            nonsense_dnp_count_df,
            nonsense_tnp_count_df)

def interactive_configuration_mode():

    # Get data and analysis parameters from user interactively.
    while True:
        analysis_choice = int(input("Do you want to identify\n[1] amino acid substitutions or\n[2] codon substitutions?\nPlease enter [1] or [2]: "))
        if analysis_choice == 1:
            break
        elif analysis_choice == 2:
            break
        else:
            print("Input error. Please enter [1] or [2]")
    
    if analysis_choice == 1:
        codon_df_switch = False
    else:
        codon_df_switch = True
    """
    while True:
        reverse_read_choice = int(input("Do you want to analyze\n[1] forward reads or \n[2] reverse reads?\nPlease enter [1] or [2]: "))
        if reverse_read_choice == 1:
            break
        elif reverse_read_choice == 2:
            break
        else:
            print("Input error. Please enter [1] for forward reads or [2] for reverse reads.")
    """
    while True:
        num_bins = int(input("Please enter the number of bins in this analysis: "))
        if num_bins >= 2:
            break
        else:
            print("You must enter at least two bins. Please re-enter.")

    num_alpha_dict = dict(zip(range(0, 26), string.ascii_uppercase))    
    path_fastq = []
    intensity = []

    # Ask the user for the path to FASTQ files and bin intensity for each bin
    for facs_bin in range(num_bins):
        path_fastq_bin = str(input(
        "Input path to FASTQ file containing reads for bin {}: ".format(num_alpha_dict[facs_bin])))
        path_fastq.append(path_fastq_bin)
        intensity_bin = int(input(
        "Please enter the intensity weight for bin {}: ".format(num_alpha_dict[facs_bin])))
        intensity.append(intensity_bin)
        
    # Ask user for the wild-type domain sequence:
    domain_seq_raw = str(input("Please enter the wild-type DNA sequence of the domain: "))
    
    domain_start = int(input("Please enter the position within the protein of the first amino acid in the domain: "))

    # Ask user for an output file name:
    jobname = str(input("Please enter a name for this job: "))

    # Some possible answers to yes or no questions
    yes_conditions = ['Yes', 'Y', 'y', 'yes', 'YES', 'yeah', 'uh huh', 'sure', 'alright', 'yea',
                      'yES', 'yeS', 'YEs', 'oui', 'Oui', 'aye', 'Aye', 'Sure', 'yep', 'Yep', 'Yea']
    no_conditions = ['No', 'N', 'n', 'no', 'NO', 'nah', 'nope', 'Nope', 'No way','NOT', 'mmm mmm',
                     'nO', 'Non', 'non', 'nay', 'Nay', 'nada', 'I don\'t think so', 'Nada']

    # Ask if the user wants to provide a seed for debugging the rarefaction
    while True:
        seed_switch_inp = str(input(
        "Do you want to provide a seed for rarefaction sampling?\n"+
        "Please enter yes or no: "))
        if seed_switch_inp in yes_conditions:
            seed_switch = True
            break
        elif seed_switch_inp in no_conditions:
            seed_switch = False
            break
        else:
            print("Answer must be yes or no. Please re-enter: ")
            
    if seed_switch == True:
        seed = int(input("Please enter an integer seed for rarefaction sampling: "))
    else:
        seed = random.randrange(sys.maxsize)

    
    # Ask if the user wants to filter outliers from the averaged rarefaction samplings
    while True:
        outlier_switch_inp = str(input(
        "Do you want to filter outliers from the averaged rarefaction samplings?\nPlease enter yes or no: "))
        if outlier_switch_inp in yes_conditions:
            outlier_switch = True
            break
        elif outlier_switch_inp in no_conditions:
            outlier_switch = False
            break
        else:
            print("Answer must be yes or no. Please re-enter.")


    # Ask if the user wants to filter sparse variants
    while True:
        sparsity_filter_switch_inp = str(input(
        "Do you want to mask rare variants with a total count across bins below some minimum?\n"+
        "Please enter yes or no: "))
        if sparsity_filter_switch_inp in yes_conditions:
            sparsity_filter_switch = True
            break
        elif sparsity_filter_switch_inp in no_conditions:
            sparsity_filter_switch = False
            break
        else:
            print("Answer must be yes or no. Please re-enter.")

    if sparsity_filter_switch == True:
        minimum_total_counts = int(input("Please enter the minimum number of counts over all bins required to pass the rare variant filter: "))
    else:
        minimum_total_counts = 0 # default

    # Ask the user which type of quality filter is preferred: filter using average
    # q score over domain within read passing a threshold or filtering based on
    # likely number of errors within the read (default = 1)
    while True:
        quality_filter_choice = int(input("How do you want to filter reads by quality score?\n"+
        "[1] Filter reads in which the average quality score over the range of the domain is below some threshold\n"+
        "[2] Filter reads in which the expected number of errors in the domain is greater than some number\n"+
        "[3] Filter reads not meeting conditions outlined in both [1] and [2]\n"+
        "[4] Each nucleotide in the read must surpass a minimum required quality score\n"+
        "Please enter [1], [2], [3] or [4]: "))
        if quality_filter_choice == 1:
            threshold = float(input("Enter minimum average quality score over the domain: "))
            errors = 1 # default
            required_passrate = 90 # default
            break
        elif quality_filter_choice == 2:
            threshold = 30 # default
            errors = float(input("Enter maxmium number of expected errors in the domain: "))
            required_passrate = 90 # default
            break
        elif quality_filter_choice == 3:
            threshold = float(input("Enter minimum average quality score over the domain: "))
            errors = float(input("Enter maxmium number of expected errors in the domain: "))
            required_passrate = 90 # default
            break
        elif quality_filter_choice == 4:
            threshold = float(input("Enter the minimum quality score each nucleotide must\n    achieve in order for a read to pass the quality filter: "))
            errors = 1 # default
            required_passrate = float(input("Enter the percentage of nucleotides over the read\n    that must surpass the threshold entered above: "))
            break
        else:
            print("Must choose option [1], [2], [3], or [4]. Please re-enter.")

    return (codon_df_switch,
            reverse_read_choice,
            num_bins,
            num_alpha_dict,
            path_fastq,
            intensity,
            domain_seq_raw,
            domain_start,
            jobname,
            seed,
            outlier_switch,
            sparsity_filter_switch,
            minimum_total_counts,
            quality_filter_choice,
            threshold,
            errors,
            required_passrate)

def parse_config_file(config_file):
    parser = ConfigParser()
    parser.read(config_file)
    codon_df_switch = parser.getboolean('Basic analysis parameters', 'Codon analysis')
    #reverse_read_choice = parser.getint('Basic analysis parameters', 'Forward or reverse reads')
    #if reverse_read_choice not in [1, 2]:
    #    print("Error. Must choose [1] to analyze forward reads or [2] to analyze reverse reads. Check configuration file.")
    #    sys.exit()
    num_bins = parser.getint('Data from DMS experiment', 'Number of bins')
    if num_bins < 2:
        print("Error. You must enter at least two bins. Check configuration file.")
        sys.exit()
    num_alpha_dict = dict(zip(range(0, 26), string.ascii_uppercase))
    path_fastq_raw_r1 = parser.get('Data from DMS experiment', 'R1 FastQ files')
    path_fastq_r1 = path_fastq_raw_r1.split()
    path_fastq_raw_r2 = parser.get('Data from DMS experiment', 'R2 FastQ files')
    path_fastq_r2 = path_fastq_raw_r2.split()
    if len(path_fastq_r1) != num_bins:
        print("Error. Must choose R1 and R2 FastQ files for each bin. Check configuration file.")
        sys.exit()
    if len(path_fastq_r2) != num_bins:
        print("Error. Must choose R1 and R2 FastQ files for each bin. Check configuration file.")
        sys.exit()
    intensity_raw = parser.get('Data from DMS experiment', 'Bin intensities')
    intensity_string = intensity_raw.split()
    intensity = [float(x) for x in intensity_string]
    if len(intensity) != num_bins:
        print("Error. Must enter a bin intensity for each bin. Check configuration file.")
        sys.exit()
    domain_seq_input_r1 = parser.get('Wild-type sequence', 'R1 Domain sequence')
    domain_seq_list_r1 = domain_seq_input_r1.split()
    domain_seq_raw_r1 = domain_seq_list_r1[0]
    domain_seq_input_r2 = parser.get('Wild-type sequence', 'R2 Domain sequence')
    domain_seq_list_r2 = domain_seq_input_r2.split()
    domain_seq_raw_r2 = domain_seq_list_r2[0]
    domain_start = parser.getint('Wild-type sequence', 'Domain start')
    offset = parser.getint('Wild-type sequence', 'Offset')
    slipsite_input = parser.get('Wild-type sequence', 'Slip site sequence')
    slipsite_list = slipsite_input.split()
    slipsite_seq = slipsite_list[0]
    jobname = parser.get('Basic analysis parameters', 'Job name')
    seed_switch = parser.getboolean('Basic analysis parameters', 'Set seed')
    if seed_switch == True:
        seed = parser.getint('Basic analysis parameters', 'Seed')
    else:
        seed = random.randrange(sys.maxsize)
    outlier_switch = parser.getboolean('Basic analysis parameters', 'Eliminate outliers')
    sparsity_filter_switch = parser.getboolean('Basic analysis parameters', 'Filter rare variants')
    if sparsity_filter_switch == True:
        minimum_total_counts = parser.getint('Basic analysis parameters', 'Minimum total counts')
    else:
        minimum_total_counts = 0
 
    threshold = parser.getfloat('Quality filter', 'Quality score threshold')
    errors = parser.getfloat('Quality filter', 'Maximum expected errors in read')


    return (codon_df_switch,
            num_bins,
            num_alpha_dict,
            path_fastq_r1,
            path_fastq_r2,
            intensity,
            domain_seq_raw_r1,
            domain_seq_raw_r2,
            domain_start,
            offset,
            jobname,
            seed,
            outlier_switch,
            sparsity_filter_switch,
            minimum_total_counts,
            threshold,
            errors,
            slipsite_seq)


if __name__ == '__main__':

    if len(sys.argv) == 2:
        (codon_df_switch,
         num_bins,
         num_alpha_dict,
         path_fastq_r1,
         path_fastq_r2,
         intensity,
         domain_seq_raw_r1,
         domain_seq_raw_r2,
         domain_start,
         offset,
         jobname,
         seed,
         outlier_switch,
         sparsity_filter_switch,
         minimum_total_counts,
         threshold,
         errors,
         slipsite_seq) = parse_config_file(sys.argv[1])
    else:
        print("Error! Please include config file.")


    # Set up script work output file
    logfile_string = datetime.now().strftime("%Y-%m-%d-%H-%M_output.txt")
    #sys.stdout = open(jobname+'_'+logfile_string, 'wt')

    # Print the name of the job to the terminal and log file:
    print(datetime.now())
    print("Job name: ", jobname)
    print("Printing codon tables: ", codon_df_switch)
    """
    if reverse_read_choice == 1:
        print("Analyzing forward reads.")
    else:
        print("Analyzing reverse reads.")
    """

    # Set up the genetic code dictionary
    genetic_code_dict = genetic_code_func()
    
    # Print the FASTQ path and mean fluorescence intensities of each bin
    print("Parameters used in this analysis...")
    for facs_bin in range(num_bins):
        print("Path to R1 FASTQ file, bin", num_alpha_dict[facs_bin], ":\n", path_fastq_r1[facs_bin])
        print("Path to R2 FASTQ file, bin", num_alpha_dict[facs_bin], ":\n", path_fastq_r2[facs_bin])
        print("Mean fluoresence intensiy, bin", num_alpha_dict[facs_bin], ":", intensity[facs_bin])
    
    # Process the wild-type domain sequence
    domain_seq_r1, domain_codons_r1, seqlen_r1 = domain_processor(domain_seq_raw_r1)
    domain_seq_r2, domain_codons_r2, seqlen_r2 = domain_processor(domain_seq_raw_r2)
    
    overlap_seq, s1_slice_idx, s2_slice_idx, consensus_seq_wt_raw = get_overlap(domain_seq_raw_r1, domain_seq_raw_r2)

    consensus_seq_wt = offset_sequence(consensus_seq_wt_raw, offset)
    
    len_overlap = len(overlap_seq)

    domain_seq, domain_codons, seqlen = domain_processor(consensus_seq_wt)

    domain_translated = ribosome(domain_codons, genetic_code_dict)

    # Print the domain in terms of a raw string, the codons, and the amino acid sequence
    print("Domain sequence: \n", domain_seq)
    print("Domain sequence in codons: \n", domain_codons)
    print("Amino acid sequence of the wild-type domain: \n", domain_translated)


    # Parse fastq files
    print("Importing reads from FASTQ files...")
    raw_data_r1 = []
    for facs_bin in range(num_bins):
        raw_data_bin = fastq_parser(path_fastq_r1[facs_bin])
        raw_data_r1.append(raw_data_bin)

    raw_data_r2 = []
    for facs_bin in range(num_bins):
        raw_data_bin = fastq_parser(path_fastq_r2[facs_bin])
        raw_data_r2.append(raw_data_bin)  
    print("Data from FASTQ files loaded")

    # Find start position of the wild-type domain within the reads:
    seqstart = seqstartfinder(raw_data_r1[0], domain_seq_r1)
    #seqstart = 0

    print("Domain start position within reads: ", seqstart)
    print("Domain sequence length: ", seqlen)


    # Filter the data using the provided cutoff
    # print("Filtering reads using quality score threshold of: ", threshold)

    # Filter data in each bin
    filtered_data_raw = []
    passed_r1 = []
    passed_r2 = []
    avg_score_r1 = []
    avg_score_r2 = []
    for facs_bin in range(num_bins):
        print("Filtering reads from bin", num_alpha_dict[facs_bin], "...")
        print("Filtering reads by maximum expected errors and average quality score over domain\n"+
                  "Maximum expected errors:", errors, "\n"+
                  "Required mean quality score over domain: ", threshold)
        filtered_dict_bin, avg_score_r1_bin, avg_score_r2_bin, passed_r1_bin, passed_r2_bin = qualityfilter3(
            raw_data_r1[facs_bin], raw_data_r2[facs_bin], seqstart, seqlen_r1, len_overlap, threshold, errors)
        filtered_data_raw.append(filtered_dict_bin)
        passed_r1.append(passed_r1_bin)
        passed_r2.append(passed_r2_bin)
        avg_score_r1.append(avg_score_r1_bin)
        avg_score_r2.append(avg_score_r2_bin)
        print("Percent of passing R1 reads from bin", num_alpha_dict[facs_bin],
              ":", passed_r1[facs_bin])
        print("Percent of passing R2 reads from bin", num_alpha_dict[facs_bin],
              ":", passed_r2[facs_bin])
        print("Average mean score of passing R1 reads from bin",
              num_alpha_dict[facs_bin], ":", avg_score_r1[facs_bin])
        print("Average mean score of passing R2 reads from bin",
              num_alpha_dict[facs_bin], ":", avg_score_r2[facs_bin])

    # Create consensus sequences
    filtered_data = []
    both_passed_qscore = []
    len_dict = []
    len_filtered = []
    for facs_bin in range(num_bins):
        print("Building consensus sequence from bin", num_alpha_dict[facs_bin], "...")
        filtered_data_bin, both_passed_qscore_bin, len_dict_bin, len_filtered_bin = make_consensus_seq(filtered_data_raw[facs_bin], s1_slice_idx, s2_slice_idx, offset)
        filtered_data.append(filtered_data_bin)
        both_passed_qscore.append(both_passed_qscore_bin)
        len_dict.append(len_dict_bin)
        len_filtered.append(len_filtered_bin)
        print("Number of clusters (R1 or R2) passing quality filter from bin", num_alpha_dict[facs_bin], ":", len_dict[facs_bin])
        print("Number of clusters in which R1 and R2 both passed quality filter from bin", num_alpha_dict[facs_bin], ":", both_passed_qscore[facs_bin])
        #print("Number of clusters in which overlapping sequences of R1 = R2 in bin", num_alpha_dict[facs_bin], ":", len_filtered[facs_bin])

    #filtered_data, both_passed_qscore, len_dict, len_filtered = make_consensus_seq(filtered_data_dict, s1_slice, s2_slice)

    """
    if reverse_read_choice == 1:
        filtered_data = []
        for facs_bin in range(num_bins):
            filtered_data_bin = filtered_data_raw[facs_bin]
            filtered_data.append(filtered_data_bin)

    if reverse_read_choice == 2:
        print("Calculating the reverse complement of R2 sequence to get forward reads...")
        filtered_data = []
        for facs_bin in range(num_bins):
            filtered_data_bin = revcomp_data_func(filtered_data_raw[facs_bin])
            filtered_data.append(filtered_data_bin)
    """

    # Get a pool of reads containing only wild-type and single-codon variants
    master_reads_pool = []
    histo = []
    for facs_bin in range(num_bins):
        print("Identifying all reads in bin", num_alpha_dict[facs_bin],
              "with zero and one mutated codons...")
        master_reads_pool_bin, histo_bin = filter_multiple_codon_variants(filtered_data[facs_bin], domain_codons)
        master_reads_pool.append(master_reads_pool_bin)
        histo.append(histo_bin)
        print("Number of reads processed from bin", num_alpha_dict[facs_bin], ":",
              len(filtered_data[facs_bin]))
        print("Number of wild-type and single-codon mutational variants processed from",
              num_alpha_dict[facs_bin], ":", len(master_reads_pool[facs_bin]))
        print("Number of wild-type reads identified in", num_alpha_dict[facs_bin], ":", histo[facs_bin][0])
        print("Number of single-codon mutational variants processed from", num_alpha_dict[facs_bin], ":", histo[facs_bin][1])
        print("Number of reads with 2 mutated codons in bin", num_alpha_dict[facs_bin], ":", histo[facs_bin][2])
        print("Number of reads with 3 mutated codons in bin", num_alpha_dict[facs_bin], ":", histo[facs_bin][3])
        print("Number of reads with 4 mutated codons in bin", num_alpha_dict[facs_bin], ":", histo[facs_bin][4])
        print("Number of reads with 5 mutated codons in bin", num_alpha_dict[facs_bin], ":", histo[facs_bin][5])
        print("Number of reads with 6 mutated codons in bin", num_alpha_dict[facs_bin], ":", histo[facs_bin][6])
        print("Number of reads with 7 mutated codons in bin", num_alpha_dict[facs_bin], ":", histo[facs_bin][7])
        print("Number of reads with 8 mutated codons in bin", num_alpha_dict[facs_bin], ":", histo[facs_bin][8])
        print("Number of reads with 9 mutated codons in bin", num_alpha_dict[facs_bin], ":", histo[facs_bin][9])
        print("Number of reads with 10 mutated codons in bin", num_alpha_dict[facs_bin], ":", histo[facs_bin][10])
        print("Number of reads with more than 10 mutated codons in bin", num_alpha_dict[facs_bin], ":", histo[facs_bin][11])


    # Find bin with smallest number of reads, including WT and one-codon mutants
    len_data = []
    for facs_bin in range(num_bins):
        len_data_bin = len(master_reads_pool[facs_bin])
        len_data.append(len_data_bin)
    minimum = min(len_data)
    index_min = np.argmin(len_data)
    print("Bin with the smallest number of reads: ", num_alpha_dict[index_min])
    print("Number of reads in bin", num_alpha_dict[index_min], ":", minimum)


    # Sample data twice. Test using numpy sampling vs default Python sampling
    # For now, just use default sampling
    print("Rarefaction sampling...")
    random.seed(seed)
    print("Seed for rarefaction sampling: ", seed)
    rarefaction_01_reads = []
    rarefaction_02_reads = []
    for facs_bin in range(num_bins):
        rarefaction_01_reads_bin = random.sample(master_reads_pool[facs_bin], minimum)
        rarefaction_01_reads.append(rarefaction_01_reads_bin)
        rarefaction_02_reads_bin = random.sample(master_reads_pool[facs_bin], minimum)
        rarefaction_02_reads.append(rarefaction_02_reads_bin)

    # Split reads pool into upstream and downstream segments. Upstream segments will be analyzed with amino acid procedure
    # Downstream segments will be analyzed with nucleotide procedure.
    # First, split wild-type consensus sequence
    print("Splitting wild-type consensus sequence at slip site.")
    print("Slip site sequence: ", slipsite_seq)
    upstream_domain_seq_raw, downstream_domain_seq_raw = split_consensus_at_slipsite(consensus_seq_wt, slipsite_seq)
    # Process upstream wild-type domain
    upstream_domain_seq, upstream_domain_codons, upstream_seqlen = domain_processor(upstream_domain_seq_raw)
    upstream_domain_translated = ribosome(upstream_domain_codons, genetic_code_dict)
    print("Wild-type sequence upstream from slip site: \n", upstream_domain_seq)
    print("Wild-type sequence upstream from slip site in codons: \n", upstream_domain_codons)
    print("Amino acid sequence of the wild-type domain, upstream from slip site: \n", upstream_domain_translated)


    # Process downstream wild-type domain
    downstream_domain_seq, downstream_domain_codons, downstream_seqlen = domain_processor(downstream_domain_seq_raw)
    downstream_domain_translated = ribosome(downstream_domain_codons, genetic_code_dict)
    print("Wild-type sequence downstream from slip site: \n", downstream_domain_seq)
    print("Wild-type sequence downstream from slip site in codons: \n", downstream_domain_codons)
    print("Amino acid sequence of the wild-type domain, downstream from slip site: \n", downstream_domain_translated)


    # Now, split rarefaction 1 and rarefaction 2 into upstream and downstream read segments
    rarefaction_01_reads_upstream = []
    rarefaction_01_reads_downstream = []
    wtcount_01 = []
    for facs_bin in range(num_bins):
        print("Splitting reads into segments upstream and downstream of slip site in rarefaction 1 in bin", num_alpha_dict[facs_bin])
        rarefaction_01_reads_upstream_bin, rarefaction_01_reads_downstream_bin, wtcount_01_bin = split_reads_at_slipsite(rarefaction_01_reads[facs_bin], consensus_seq_wt, slipsite_seq)
        rarefaction_01_reads_upstream.append(rarefaction_01_reads_upstream_bin)
        rarefaction_01_reads_downstream.append(rarefaction_01_reads_downstream_bin)
        wtcount_01.append(wtcount_01_bin)
        print("Number of true wild-type reads identified in rarefaction 1, bin", 
              num_alpha_dict[facs_bin], ":", wtcount_01[facs_bin])


    # Identify number of WT reads and reads with each different single-codon mutation in each rarefaction
    # First, rarefaction 1
    rarefaction_01_df = []
    wtcount_01b = []
    single_codon_variants_01 = []
    single_silent_mutations_01 = []
    single_nonsense_mutations_01 = []
    single_missense_mutations_01 = []
    for facs_bin in range(num_bins):
        print("Counting wild-type reads and single-codon variants in rarefaction 1 in bin", num_alpha_dict[facs_bin])
        (rarefaction_01_df_bin,
         wtcount_01_bin,
         single_codon_variants_01_bin,
         single_silent_mutations_01_bin,
         single_nonsense_mutations_01_bin,
         single_missense_mutations_01_bin) = variant_identifier(rarefaction_01_reads_upstream[facs_bin], upstream_domain_codons, codon_df_switch)
        rarefaction_01_df.append(rarefaction_01_df_bin)
        wtcount_01b.append(wtcount_01_bin)
        single_codon_variants_01.append(single_codon_variants_01_bin)
        single_silent_mutations_01.append(single_silent_mutations_01_bin)
        single_nonsense_mutations_01.append(single_nonsense_mutations_01_bin)
        single_missense_mutations_01.append(single_missense_mutations_01_bin)
        print("Number of reads processed from rarefaction 1, bin", num_alpha_dict[facs_bin], ":",
              len(rarefaction_01_reads[facs_bin]))
        print("Number of single-codon mutant variants processed from rarefaction 1, bin",
              num_alpha_dict[facs_bin], ":", single_codon_variants_01[facs_bin])
        print("Reads with single synonymous mutations in rarefaction 1, bin",
              num_alpha_dict[facs_bin], ":", single_silent_mutations_01[facs_bin])
        print("Reads with single nonsense mutations in rarefaction 1, bin",
              num_alpha_dict[facs_bin], ":", single_nonsense_mutations_01[facs_bin])
        print("Reads with single missense mutations in rarefaction 1, bin",
              num_alpha_dict[facs_bin], ":", single_missense_mutations_01[facs_bin])
        print("Variant identification for rarefaction 1, bin", num_alpha_dict[facs_bin],
              "complete.")

    # Now, split rarefaction 2 into upstream and downstream read segments
    rarefaction_02_reads_upstream = []
    rarefaction_02_reads_downstream = []
    wtcount_02 = []
    for facs_bin in range(num_bins):
        print("Splitting reads into segments upstream and downstream of slip site in rarefaction 2 in bin", num_alpha_dict[facs_bin])
        rarefaction_02_reads_upstream_bin, rarefaction_02_reads_downstream_bin, wtcount_02_bin = split_reads_at_slipsite(rarefaction_02_reads[facs_bin], consensus_seq_wt, slipsite_seq)
        rarefaction_02_reads_upstream.append(rarefaction_02_reads_upstream_bin)
        rarefaction_02_reads_downstream.append(rarefaction_02_reads_downstream_bin)
        wtcount_02.append(wtcount_02_bin)
        print("Number of true wild-type reads identified in rarefaction 2, bin", 
              num_alpha_dict[facs_bin], ":", wtcount_02[facs_bin])


    # Next, rarefaction 2
    rarefaction_02_df = []
    wtcount_02b = []
    single_codon_variants_02 = []
    single_silent_mutations_02 = []
    single_nonsense_mutations_02 = []
    single_missense_mutations_02 = []
    for facs_bin in range(num_bins):
        print("Counting wild-type reads and single-codon variants in rarefaction 2 in bin", num_alpha_dict[facs_bin])
        (rarefaction_02_df_bin,
         wtcount_02_bin,
         single_codon_variants_02_bin,
         single_silent_mutations_02_bin,
         single_nonsense_mutations_02_bin,
         single_missense_mutations_02_bin) = variant_identifier(rarefaction_02_reads_upstream[facs_bin], upstream_domain_codons, codon_df_switch)
        rarefaction_02_df.append(rarefaction_02_df_bin)
        wtcount_02b.append(wtcount_02_bin)
        single_codon_variants_02.append(single_codon_variants_02_bin)
        single_silent_mutations_02.append(single_silent_mutations_02_bin)
        single_nonsense_mutations_02.append(single_nonsense_mutations_02_bin)
        single_missense_mutations_02.append(single_missense_mutations_02_bin)
        print("Number of reads processed from rarefaction 2, bin", num_alpha_dict[facs_bin], ":",
              len(rarefaction_02_reads[facs_bin]))
        print("Number of single-codon mutant variants processed from rarefaction 2, bin",
              num_alpha_dict[facs_bin], ":", single_codon_variants_02[facs_bin])
        print("Reads with single synonymous mutations in rarefaction 2, bin",
              num_alpha_dict[facs_bin], ":", single_silent_mutations_02[facs_bin])
        print("Reads with single nonsense mutations in rarefaction 2, bin",
              num_alpha_dict[facs_bin], ":", single_nonsense_mutations_02[facs_bin])
        print("Reads with single missense mutations in rarefaction 2, bin",
              num_alpha_dict[facs_bin], ":", single_missense_mutations_02[facs_bin])
        print("Variant identification for rarefaction 2, bin", num_alpha_dict[facs_bin],
              "complete.")


    # Get weighted intensities
    print("Calculating frameshifting enrichment of amino acid variants before the slip site")
    varcount_ints_weight_01, wtcount_ints_weight_01, var_wt_norm_01, varcount_dfsum_01 = weighted_intensity_func(rarefaction_01_df, wtcount_01, intensity)
    # For rarefaction 02:
    varcount_ints_weight_02, wtcount_ints_weight_02, var_wt_norm_02, varcount_dfsum_02 = weighted_intensity_func(rarefaction_02_df, wtcount_02, intensity)
    
    # Print the sum of the wild-type counts weighted by intensity:
    print("Wild type weighted average (rarefaction 1): ", wtcount_ints_weight_01)
    print("Wild type weighted average (rarefaction 2): ", wtcount_ints_weight_02)
    avg_wtcount_weighted_avg = (wtcount_ints_weight_01 + wtcount_ints_weight_02) / 2
    print("Average over rarefactions, wild-type weighted average: ", avg_wtcount_weighted_avg)

    # Get the dataframe containing the rarefaction-averaged (n=2) weighted-intensity:
    print("Averaging rarefaction samples...")
    varcount_ints_weight_avg = average_two_dfs_func(varcount_ints_weight_01, varcount_ints_weight_02, upstream_domain_codons, codon_df_switch)
    varcount_dfsum_avg = average_two_dfs_func(varcount_dfsum_01, varcount_dfsum_02, upstream_domain_codons, codon_df_switch)

    # Get averages of each bin count:
    rarefaction_counts_avg = []
    for facs_bin in range(num_bins):
        rarefaction_counts_avg_bin = average_two_dfs_func(rarefaction_01_df[facs_bin], rarefaction_02_df[facs_bin], upstream_domain_codons, codon_df_switch)
        rarefaction_counts_avg.append(rarefaction_counts_avg_bin)

    # Get Pearson_r and avg of the two rarefaction-sampled, wt-normalized, bin intensity-weighted counts:
    print("Executing routine to find average of the rarefactions...")
    var_wt_norm_avg_list = rarefaction_corr_func(var_wt_norm_01, var_wt_norm_02, upstream_domain_codons, jobname, outlier_switch, codon_df_switch, domain_start)
    var_wt_norm_avg = var_wt_norm_avg_list[0]
    var_wt_norm_avg_outliers_removed = var_wt_norm_avg_list[1]

    # Save Excel workbooks
    print("Saving Excel workbooks...")
    # Excel workbook for rarefaction 01
    workbook_counts_rarefaction_01 = pd.ExcelWriter(jobname+'_rarefaction_01.xlsx')
    var_wt_norm_01.to_excel(workbook_counts_rarefaction_01, 'Variant count weighted by intensity relative to WT')
    varcount_ints_weight_01.to_excel(workbook_counts_rarefaction_01, 'Variant count weighted by intensity')
    for facs_bin in range(num_bins):
        rarefaction_01_df[facs_bin].to_excel(workbook_counts_rarefaction_01, 'Rarefied counts from bin {}'.format(num_alpha_dict[facs_bin]))
    varcount_dfsum_01.to_excel(workbook_counts_rarefaction_01, 'Total rarefied counts')
    workbook_counts_rarefaction_01.save()

    # Excel workbook for rarefaction 02
    workbook_counts_rarefaction_02 = pd.ExcelWriter(jobname+'_rarefaction_02.xlsx')
    var_wt_norm_02.to_excel(workbook_counts_rarefaction_02, 'Variant count weighted by intensity relative to WT')
    varcount_ints_weight_02.to_excel(workbook_counts_rarefaction_02, 'Variant count weighted by intensity')
    for facs_bin in range(num_bins):
        rarefaction_02_df[facs_bin].to_excel(workbook_counts_rarefaction_02, 'Rarefied counts from bin {}'.format(num_alpha_dict[facs_bin]))
    varcount_dfsum_02.to_excel(workbook_counts_rarefaction_02, 'Total rarefied counts')
    workbook_counts_rarefaction_02.save()

    # Excel workbook for average of rarefactions
    workbook_counts_rarefaction_avg = pd.ExcelWriter(jobname+'_rarefaction_avg.xlsx')
    var_wt_norm_avg.to_excel(workbook_counts_rarefaction_avg, 'Variant count weighted by intensity relative to WT')
    varcount_ints_weight_avg.to_excel(workbook_counts_rarefaction_avg, 'Variant count weighted by intensity')
    for facs_bin in range(num_bins):
        rarefaction_counts_avg[facs_bin].to_excel(workbook_counts_rarefaction_avg, 'Rarefied counts from bin {}'.format(num_alpha_dict[facs_bin]))
    varcount_dfsum_avg.to_excel(workbook_counts_rarefaction_avg, 'Total rarefied counts')
    workbook_counts_rarefaction_avg.save()

    # Excel workbook for average of rarefactions, variant count weighted by intensity only
    workbook_var_ints_rarefaction_avg = pd.ExcelWriter(jobname+'_rarefaction_avg_var_ints.xlsx')
    varcount_ints_weight_avg.to_excel(workbook_var_ints_rarefaction_avg, 'Variant count weighted by intensity')
    workbook_var_ints_rarefaction_avg.save()

    # Excel workbook of WT-normalized weighted surface intensities with outliers removed
    if outlier_switch == True:
        workbook_var_wt_norm_avg_outliers_removed = pd.ExcelWriter(jobname+'_rarefaction_avg_outliers_removed.xlsx')
        var_wt_norm_avg_outliers_removed.to_excel(workbook_var_wt_norm_avg_outliers_removed, 'Variant count weighted by intensity relative to WT')
        workbook_var_wt_norm_avg_outliers_removed.save()

    if (sparsity_filter_switch, outlier_switch) == (True, False):
        print("Filtering sparse variants with cutoff: ", minimum_total_counts)
        rare_variants_filtered_df = filter_rare_variants_func(rarefaction_counts_avg, var_wt_norm_avg, upstream_domain_codons, codon_df_switch, minimum_total_counts)
        workbook_sparsity_filter = pd.ExcelWriter(jobname+'_rare_variants_removed.xlsx')
        rare_variants_filtered_df.to_excel(workbook_sparsity_filter, 'Average variant scores')
        workbook_sparsity_filter.save()
        
    if (sparsity_filter_switch, outlier_switch) == (True, True):
        print("Filtering sparse variants with cutoff: ", minimum_total_counts)
        rare_variants_filtered_df = filter_rare_variants_func(rarefaction_counts_avg, var_wt_norm_avg_outliers_removed, upstream_domain_codons, codon_df_switch, minimum_total_counts)
        workbook_sparsity_filter = pd.ExcelWriter(jobname+'_rare_variants_and_outliers_removed.xlsx')
        rare_variants_filtered_df.to_excel(workbook_sparsity_filter, 'Avg variant scores-outliers')
        workbook_sparsity_filter.save()

    if sparsity_filter_switch == False:
        print("No sparse variant filter applied.")


    # Perform nucleotide analysis on sequence downstream of slip site
    # First, rarefaction 01
    nucleotide_01_df = []
    for facs_bin in range(num_bins):
        print("Performing nucleotide analysis for reads downstream of slip site, rarefaction 1, bin", num_alpha_dict[facs_bin])
        (wildtype_reads,
         single_codon_variants,
         single_silent_reads,
         single_missense_reads,
         single_nonsense_reads,
         SNPs,
         DNPs,
         TNPs,
         missense_SNPs,
         missense_DNPs,
         missense_TNPs,
         silent_SNPs,
         silent_DNPs,
         silent_TNPs,
         nonsense_SNPs,
         nonsense_DNPs,
         nonsense_TNPs,
         all_count_df,
         all_snp_count_df,
         all_dnp_count_df,
         all_tnp_count_df,
         syn_count_df,
         syn_snp_count_df,
         syn_dnp_count_df,
         syn_tnp_count_df,
         nonsyn_count_df,
         nonsyn_snp_count_df,
         nonsyn_dnp_count_df,
         nonsyn_tnp_count_df,
         nonsense_count_df,
         nonsense_snp_count_df,
         nonsense_dnp_count_df,
         nonsense_tnp_count_df) = nucleotide_variant_identifier(rarefaction_01_reads_downstream[facs_bin], downstream_domain_codons, downstream_domain_seq)
        nucleotide_01_df.append(all_snp_count_df)

    # Next, rarefaction 02
    nucleotide_02_df = []
    for facs_bin in range(num_bins):
        print("Performing nucleotide analysis for reads downstream of slip site, rarefaction 2, bin", num_alpha_dict[facs_bin])
        (wildtype_reads,
         single_codon_variants,
         single_silent_reads,
         single_missense_reads,
         single_nonsense_reads,
         SNPs,
         DNPs,
         TNPs,
         missense_SNPs,
         missense_DNPs,
         missense_TNPs,
         silent_SNPs,
         silent_DNPs,
         silent_TNPs,
         nonsense_SNPs,
         nonsense_DNPs,
         nonsense_TNPs,
         all_count_df,
         all_snp_count_df,
         all_dnp_count_df,
         all_tnp_count_df,
         syn_count_df,
         syn_snp_count_df,
         syn_dnp_count_df,
         syn_tnp_count_df,
         nonsyn_count_df,
         nonsyn_snp_count_df,
         nonsyn_dnp_count_df,
         nonsyn_tnp_count_df,
         nonsense_count_df,
         nonsense_snp_count_df,
         nonsense_dnp_count_df,
         nonsense_tnp_count_df) = nucleotide_variant_identifier(rarefaction_02_reads_downstream[facs_bin], downstream_domain_codons, downstream_domain_seq)
        nucleotide_02_df.append(all_snp_count_df)

    # Get weighted intensities
    print("Calculating frameshifting enrichment of SNVs after the slip site...")
    varcount_ints_weight_01_snv, wtcount_ints_weight_01_snv, var_wt_norm_01_snv, varcount_dfsum_01_snv = weighted_intensity_func(nucleotide_01_df, wtcount_01, intensity)
    # For rarefaction 02:
    varcount_ints_weight_02_snv, wtcount_ints_weight_02_snv, var_wt_norm_02_snv, varcount_dfsum_02_snv = weighted_intensity_func(nucleotide_02_df, wtcount_02, intensity)


    # Save Excel workbooks
    print("Saving Excel workbooks, SNVs...")
    # Excel workbook for rarefaction 01
    workbook_counts_rarefaction_01_snv = pd.ExcelWriter(jobname+'_rarefaction_01_snv.xlsx')
    var_wt_norm_01_snv.to_excel(workbook_counts_rarefaction_01_snv, 'vcount relative to WT')
    varcount_ints_weight_01_snv.to_excel(workbook_counts_rarefaction_01_snv, 'vcount weighted by intensity')
    for facs_bin in range(num_bins):
        nucleotide_01_df[facs_bin].to_excel(workbook_counts_rarefaction_01_snv, 'Rarefied counts from bin {}'.format(num_alpha_dict[facs_bin]))
    varcount_dfsum_01_snv.to_excel(workbook_counts_rarefaction_01_snv, 'Total rarefied counts')
    workbook_counts_rarefaction_01_snv.save()


    # Excel workbook for rarefaction 02
    workbook_counts_rarefaction_02_snv = pd.ExcelWriter(jobname+'_rarefaction_02_snv.xlsx')
    var_wt_norm_02_snv.to_excel(workbook_counts_rarefaction_02_snv, 'vcount relative to WT')
    varcount_ints_weight_02_snv.to_excel(workbook_counts_rarefaction_02_snv, 'vcount weighted by intensity')
    for facs_bin in range(num_bins):
        nucleotide_02_df[facs_bin].to_excel(workbook_counts_rarefaction_02_snv, 'Rarefied counts from bin {}'.format(num_alpha_dict[facs_bin]))
    varcount_dfsum_02_snv.to_excel(workbook_counts_rarefaction_02_snv, 'Total rarefied counts')
    workbook_counts_rarefaction_02_snv.save()

    # Inform the user the job is complete
    print(datetime.now())
    print("Job: ", jobname, "... complete.")

    # return 0;
# }



