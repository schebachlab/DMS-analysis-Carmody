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
        TAA *      CAA Q      AAA K      GAA E
        TAG *      CAG Q      AAG K      GAG E
        TGT C      CGT R      AGT S      GGT G
        TGC C      CGC R      AGC S      GGC G
        TGA *      CGA R      AGA R      GGA G
        TGG W      CGG R      AGG R      GGG G
        """
    codon_finder = re.compile(r'[ATCG]{3}')
    amino_acid_finder = re.compile(r'\ \w{1}[\ |\n]|\*')
    codon_list = codon_finder.findall(genetic_code)
    amino_acid_list = [x.strip() for x in amino_acid_finder.findall(genetic_code)]
    genetic_code_dict = {}
    i = 0
    while i < len(codon_list):
        genetic_code_dict[codon_list[i]] = amino_acid_list[i]
        i += 1
    return genetic_code_dict

def domain_processor(domain_sequence_raw):
    domain_sequence = re.sub("[^a-zA-Z]", "", domain_sequence_raw)
    domain_codons = textwrap.wrap(domain_sequence, 3)
    domain_length = len(domain_sequence)
    return domain_sequence, domain_codons, domain_length

def ribosome(domain_codons, genetic_code_dict):
    domain_aminoacids = []
    for codon in domain_codons:
        amino_acid = genetic_code_dict[codon]
        domain_aminoacids.append(amino_acid)
    return domain_aminoacids

def df_shape_func(domain_codons_in, codon_df_switch):
    if codon_df_switch == 3:
        domain_codons_str = ''.join(domain_codons_in)
        domain_codons = textwrap.wrap(domain_codons_str, 3)
    else:
        domain_codons = domain_codons_in
    genetic_code_dict = genetic_code_func()
    domain_translated = ribosome(domain_codons, genetic_code_dict)
    # Get shape of dataframes
    polarity_rank = 'ILFVCMAWTYGSNHPQREKD'
    polindex = textwrap.wrap(polarity_rank, 1) + ['Stop']
    snv_index = ['A', 'C', 'G', 'T']
    nucleotide_seq = list(''.join(domain_codons))
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

    # 1 = peptide, 2 = codons, 3 = snv
    if codon_df_switch == 2:
        df_index = polindex_codons
        df_columns = domain_codons
    elif codon_df_switch == 1:
        df_index = polindex
        df_columns = domain_translated
    else:
        df_index = snv_index
        df_columns = nucleotide_seq

    row_length = len(df_columns)
    column_length = len(df_index)
    array_size = row_length * column_length
    return df_index, df_columns, row_length, column_length, array_size

def flat_arrays_to_df(array1_flat, array2_flat, domain_codons, codon_df_switch):
    df_index, df_columns, row_length, column_length, array_size = df_shape_func(domain_codons, codon_df_switch)
    array_flat_combined = np.column_stack((array1_flat, array2_flat))
    array_flat_avg = np.average(array_flat_combined, axis=1)
    array_avg = array_flat_avg.reshape(column_length, row_length)
    output_df = pd.DataFrame(array_avg, index=df_index)
    output_df.columns = df_columns
    return output_df


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

def excel_parser_func(excel_workbook):
    pd_df = pd.read_excel(excel_workbook, sheet_name=0, header=0, index_col=0)
    return pd_df

def excel_writer_func(dataframe, outfilename):
    workbook = pd.ExcelWriter(outfilename+'.xlsx')
    dataframe.to_excel(workbook, 'Average')
    workbook.save()

#Rarefaction correlation routines
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

def scatterplot_arrays(array1, array2, x_axis, y_axis, filename, title):
    scatterplotname = filename+".png"
    scatterplot_fig = plt.figure()
    plt.plot(array1, array2, 'o', color='black');
    plt.xlabel(x_axis)
    plt.ylabel(y_axis)
    plt.title(title)
    scatterplot_fig.savefig(scatterplotname, dpi=400)


if __name__ == '__main__':

    while True:
        analysis_choice = int(input("Do you want to map\n[1] amino acid substitutions\n[2] codon substitutions or\n[3] nucleotide substitutions?\nPlease enter [1], [2], or [3]: "))
        if analysis_choice == 1:
            break
        elif analysis_choice == 2:
            break
        elif analysis_choice == 3:
            break
        else:
            print("Input error. Please enter [1], [2], or [3]")

    # Ask user for numerator Excel sheet:
    path1 = str(input("Please enter the path to the first excel sheet: "))
    path2 = str(input("Please enter the path to the second excel sheet: "))

    wtseq = str(input("Please enter the sequence of the WT domain: "))

    # Ask user for output filename:
    jobname = str(input("Please enter the output name: "))

    genetic_code_dict = genetic_code_func()
    domain_seq, domain_codons, seqlen = domain_processor(wtseq)
    domain_translated = ribosome(domain_codons, genetic_code_dict)

    if analysis_choice == 1:

        df_columns = domain_translated

        df1 = excel_parser_func(path1)
        df1.columns = df_columns
        df2 = excel_parser_func(path2)
        df2.columns = df_columns

    if analysis_choice == 2:

        df_columns = domain_codons

        df1 = excel_parser_func(path1)
        df1.columns = df_columns
        df2 = excel_parser_func(path2)
        df2.columns = df_columns

    if analysis_choice == 3:

        df_columns = list(domain_seq)

        df1 = excel_parser_func(path1)
        df1.columns = df_columns
        df2 = excel_parser_func(path2)
        df2.columns = df_columns

    #df_avg = average_two_dfs_func(df1, df2, df_columns, analysis_choice)
    #excel_writer_func(df_avg, jobname)


    var_wt_norm_avg_list = rarefaction_corr_func(df1, df2, df_columns, jobname, False, analysis_choice, 1)
    df_avg = var_wt_norm_avg_list[0]
    df_avg_outliers_removed = var_wt_norm_avg_list[1]

    excel_writer_func(df_avg, jobname+'_avg')
    excel_writer_func(df_avg_outliers_removed, jobname+'_avg_bonf')