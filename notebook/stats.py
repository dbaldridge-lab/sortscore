from sequence_parsing import translate_dna
import pandas as pd

############ SUMMARY STATS ############
# calculates summary statistics for activity scores in a given group of variants
def dms_stats_dna(data, act_col):
    # Ensure the column is numeric, coercing errors to NaN
    data[act_col] = pd.to_numeric(data[act_col], errors='coerce')

    # calculate average activity over all sequences
    all_var_avg = data[act_col].mean()
    all_var_min = data[act_col].min()
    all_var_max = data[act_col].max()

    # lookup WT DNA activity value
    wt_dna_act = data.loc[data['annotate_dna'] == 'wt_dna', act_col].iloc[0]

    # lookup spike-in activity value if it exists, else NaN
    spike_in_act = data.loc[data['annotate_dna'] == 'spike-in', act_col].iloc[0]

    # calculate average of missense(DNA level) activity scores
    missense_dna_avg = data.loc[data['annotate_dna'] == 'missense_dna', act_col].mean()
    # calculate range of missense(DNA level) activity scores
    missense_dna_min = data.loc[data['annotate_dna'] == 'missense_dna', act_col].min()
    missense_dna_max = data.loc[data['annotate_dna'] == 'missense_dna', act_col].max()

    # calculate average of synonymous(DNA level) activity scores
    syn_avg = data.loc[data['annotate_dna'] == 'synonymous', act_col].mean()
    syn_median = data.loc[data['annotate_dna'] == 'synonymous', act_col].median()
    syn_min = data.loc[data['annotate_dna'] == 'synonymous', act_col].min()
    syn_max = data.loc[data['annotate_dna'] == 'synonymous', act_col].max()

    return (
    all_var_avg, 
    all_var_min, 
    all_var_max, 
    wt_dna_act, 
    spike_in_act, 
    missense_dna_avg, 
    missense_dna_min, 
    missense_dna_max,
    syn_avg,
    syn_median,
    syn_min,
    syn_max
    )

# calculates summary statistics for activity scores in a given group of variants
def dms_stats_aa(data, act_col):
    # Ensure the column is numeric, coercing errors to NaN
    data[act_col] = pd.to_numeric(data[act_col], errors='coerce')

    nonsense_avg = data.loc[data['annotate_aa'] == 'nonsense', act_col].mean()
    nonsense_median = data.loc[data['annotate_aa'] == 'nonsense', act_col].median()
    nonsense_min = data.loc[data['annotate_aa'] == 'nonsense', act_col].min()
    nonsense_max = data.loc[data['annotate_aa'] == 'nonsense', act_col].max()

    missense_aa_avg = data.loc[data['annotate_aa'] == 'missense_aa', act_col].mean()
    missense_aa_min = data.loc[data['annotate_aa'] == 'missense_aa', act_col].min()
    missense_aa_max = data.loc[data['annotate_aa'] == 'missense_aa', act_col].max()

    return (
        nonsense_avg, 
        nonsense_median,
        nonsense_min, 
        nonsense_max, 
        missense_aa_avg, 
        missense_aa_min, 
        missense_aa_max
    )

############ QUALITY CHECKS ############

def get_dropout(df):

    num_aa = df.shape[1]

    # At each position there are 19 AA variants and 1 stop codon. Include 1 WT sequence.
    possible_sequences = num_aa * 20 + 1
    print("Total number of possible variants (19 AA variants, 1 stop codon)*positions + 1 WT):", possible_sequences)

    dropout_num = df.isna().sum().sum()
    print("Number of variants that were not observed (total dropout):", dropout_num)

    dropout_percent = round((dropout_num / possible_sequences) * 100, 1)
    print("Percent dropout:", dropout_percent)

    return dropout_num, dropout_percent


def highest_freq_seq_dna(counts_df, wt_seq):
    # sort from most abundant to least abundant sequences
    counts_df = counts_df.sort_values(by='counts', ascending=False)
    # Extract reference DNA sequence from the first row of the dataframe.
    abundant_seq = counts_df['bc'].iloc[0]

    if wt_seq != abundant_seq:
        print('The most abundant DNA sequence is not the same as the reference DNA sequence.')
        print('Most abundant DNA sequence:', abundant_seq)
        print('Most abundant AA sequence:', translate_dna(abundant_seq))