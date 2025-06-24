from Bio.Seq import Seq #Translate DNA sequence to protien

############ SEQUENCE PROCESSING ############

def get_reversecomplement(string1):
    '''
    https://stackoverflow.com/questions/25188968/reverse-complement-of-dna-strand-using-python
    '''
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(string1))
    return reverse_complement

# labels subsets of variants at the DNA level
def dms_annotate_dna(data, wt_ref_seq, spike_in_seq):

    # make annotations column with empty strings
    data['annotate_dna'] = ''

    # add WT DNA sequence label
    data.loc[data['bc'] == wt_ref_seq, 'annotate_dna'] = 'wt_dna'
    
    # add a missense label for variants with 1 or more DNA sequence differences
    data.loc[data['dna_seq_diff'] != '', 'annotate_dna'] = 'missense_dna'

    # add spike-in DNA sequence label 
    # (overwrites missense label for opool 3 or wt for other opools)
    data.loc[data['bc'] == spike_in_seq, 'annotate_dna'] = 'spike-in'

    # add synonymous label
    data.loc[(data['bc'] != wt_ref_seq) &
             (data['bc'] != spike_in_seq) & 
             (data['aa_seq_diff'] == ''), 
             'annotate_dna'] = 'synonymous'

    print('DNA level annotation counts: ')
    print(data['annotate_dna'].value_counts())



# labels subset of variants at the AA level
def dms_annotate_aa(data):
    # make annotations column with empty strings
    data['annotate_aa'] = ''

    data.loc[data['annotate_dna'] == 'wt_dna', 'annotate_aa'] = 'wt_dna'
    data.loc[data['annotate_dna'] == 'spike-in', 'annotate_aa'] = 'spike-in'
    data.loc[data['annotate_dna'] == 'synonymous', 'annotate_aa'] = 'synonymous'

    # add a missense label for variants with 1 or more AA sequence differences
    data.loc[data['aa_seq_diff'] != '', 'annotate_aa'] = 'missense_aa'

    # add nonsense label if AA variant sequence contains a stop codon (overwrites missense label for these variants)
    data.loc[data['aa_seq_diff'].str.contains('\*'), 'annotate_aa'] = 'nonsense'

    print('AA level annotation counts:\n')
    print(data['annotate_aa'].value_counts())


def make_codon_table():
    from Bio.Data import CodonTable

    # Standard genetic code
    codon_table = CodonTable.unambiguous_dna_by_id[1]

    # List of codons for amino acids
    amino_acid_codons = list(codon_table.forward_table.keys())

    # List of stop codons
    stop_codons = codon_table.stop_codons

    # Combine both lists to get all 64 codons
    all_codons = amino_acid_codons + stop_codons

    print(f"Total number of codons: {len(all_codons)}")
    return all_codons

# Translate DNA to protein sequence
def translate_dna(dna_sequence):
    dna_seq = Seq(dna_sequence)
    return str(dna_seq.translate())

# Define a function to compare each protein sequence to the reference sequence and report differences
## TODO: some check for runtime errors (e.g. if a * is ever listed first for a protein reference sequence)
def compare_to_reference(ref_seq, sequence):
    differences = []
    min_length = min(len(ref_seq), len(sequence))
    for i in range(min_length):
        if ref_seq[i] != sequence[i]:
            # Format: 'reference_character.position.sequence_character'
            difference = f'{ref_seq[i]}.{i+1}.{sequence[i]}'
            differences.append(difference)
    return ', '.join(differences)

# Function to add starting position along protein to the number in aa_seq_diff
def modify_position(aa_seq_diff, min_pos):
    parts = aa_seq_diff.split('.')
    if len(parts) == 3 and parts[1].isdigit():
        parts[1] = str(int(parts[1]) + min_pos - 1)
    return '.'.join(parts)

# Function to translate DNA sequence to codon list
def get_codons(dna_seq):
    if len(dna_seq) % 3 != 0:
        raise ValueError("The length of the DNA sequence must be divisible by 3 to label codons.")
    return [dna_seq[i:i+3] for i in range(0, len(dna_seq), 3)]

# Function to compare codon lists and output differences
def compare_codon_lists(wt_seq, variant_seq):
    wt_codons = get_codons(wt_seq)
    variant_codons = get_codons(variant_seq)
    differences = []
    for i, (wt_codon, var_codon) in enumerate(zip(wt_codons, variant_codons)):
        if wt_codon != var_codon:
            wt_aa = translate_dna(wt_codon)
            var_aa = translate_dna(var_codon)
            differences.append(f"{wt_aa}({wt_codon}).{i+1}.{var_aa}({var_codon})")
    return '.'.join(differences)