# figures.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec
from sequence_parsing import translate_dna
from stats import get_dropout

############ BEESWARM FUNCTIONS ############
# TODO: add averages for synonymous, nonsense, missense

def dms_swarmplot(data, title, annotate_col, score_col, export=False, output='activity_beeswarm.png'):
    # Create a swarmplot of the average activity scores
    plt.figure(figsize=(15, 10))

    # generate swarmplot
    sns.swarmplot(data, x=annotate_col, y=score_col, hue=annotate_col, legend=False)
    #sns.stripplot(data, x='annotate_aa', y='avgscore', hue='annotate_aa', jitter=True, alpha=0.5)

    # add title and labels
    plt.title(title, fontsize=22)
    plt.xlabel('Annotation', fontsize=18)
    plt.ylabel('Average Activity Score', fontsize=18)

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    if export:
        # Save the beeswarm as a PNG file
        plt.savefig(fname=output_file, dpi=300, format='png')

    plt.show()

############ HEATMAP FUNCTIONS ############

def extract_position(sequence_diff):
    parts = sequence_diff.split('.')
    if len(parts) == 3 and parts[1].isdigit():
        return parts[2], int(parts[1])
    return None, None

def extract_value(cell):
    if cell == 'WT':
        return 'WT'
    elif pd.notna(cell):
        return cell
    else:
        return np.nan

def dms_matrix_template(num_aa, min_pos=1, mutant_type='aa'):
    column_values = list(range(min_pos, min_pos + num_aa))
    if mutant_type == 'aa':
        row_labels = ['W', 'F', 'Y', 'P', 'M', 'I', 'L', 'V', 'A', 'G', 'C', 'S', 'T', 'Q', 'N', 'D', 'E', 'H', 'R', 'K', '*']
    elif mutant_type == 'dna':
        row_labels = [
                      'W(TGG)',                                # Tryptophan
                      'F(TTT)', 'F(TTC)',                      # Phenylalanine
                      'Y(TAT)', 'Y(TAC)',                      # Tyrosine
                      'P(CCT)', 'P(CCC)', 'P(CCA)', 'P(CCG)',  # Proline
                      'M(ATG)',                                # Methionine
                      'I(ATT)', 'I(ATC)', 'I(ATA)',            # Isoleucine
                      'L(TTA)', 'L(TTG)', 'L(CTT)', 'L(CTC)', 'L(CTA)', 'L(CTG)',  # Leucine
                      'V(GTT)', 'V(GTC)', 'V(GTA)', 'V(GTG)',  # Valine
                      'A(GCT)', 'A(GCC)', 'A(GCA)', 'A(GCG)',  # Alanine
                      'G(GGT)', 'G(GGC)', 'G(GGA)', 'G(GGG)',  # Glycine
                      'C(TGT)', 'C(TGC)',                      # Cysteine
                      'S(TCT)', 'S(TCC)', 'S(TCA)', 'S(TCG)', 'S(AGT)', 'S(AGC)',  # Serine
                      'T(ACT)', 'T(ACC)', 'T(ACA)', 'T(ACG)',  # Threonine
                      'Q(CAA)', 'Q(CAG)',                      # Glutamine
                      'N(AAT)', 'N(AAC)',                      # Asparagine
                      'D(GAT)', 'D(GAC)',                      # Aspartic acid
                      'E(GAA)', 'E(GAG)',                      # Glutamic acid
                      'H(CAT)', 'H(CAC)',                      # Histidine
                      'R(CGT)', 'R(CGC)', 'R(CGA)', 'R(CGG)', 'R(AGA)', 'R(AGG)',  # Arginine
                      'K(AAA)', 'K(AAG)',                      # Lysine
                      '*(TAA)', '*(TAG)', '*(TGA)'             # Stop codons
                      ]
    return pd.DataFrame(index=row_labels, columns=column_values)

def make_dms_matrix(data, score_col, num_aa, wt_seq, min_pos=1, mutant_type='aa'):
    # supplied wt_seq should be dna sequence if mutant_type is 'dna' and aa sequence if mutant_type is 'aa'
    # Set min_pos to first position of the sequence as listed in the data. You can update the labels later in dms_heatmap.

    # drop rows with missing values in the score column
    data = data.dropna(subset=[score_col])

    # scaffold the empty matrix
    matrix = dms_matrix_template(num_aa, min_pos, mutant_type)
    
    if mutant_type == 'aa':
        diff_col = 'aa_seq_diff'
    elif mutant_type == 'dna':
        diff_col = 'codon_diff'

    for index, row in data.iterrows():
        char, col = extract_position(row[diff_col])
        if char in matrix.index and col in matrix.columns:
            matrix.at[char, col] = row[score_col]

    if mutant_type == 'aa':
        for index, amino_acid in enumerate(wt_seq, start=min_pos):
            if amino_acid in matrix.index and index in matrix.columns:
                matrix.at[amino_acid, index] = 'WT'

    if mutant_type == 'dna':
        codon_enumeration = [((i // 3)+1, wt_seq[i:i+3]) for i in range(0, len(wt_seq), 3)]
        for index, codon in codon_enumeration:
            aa = translate_dna(codon)
            codon_index = f'{aa}({codon})'
            if codon_index in matrix.index and index in matrix.columns:
                matrix.at[codon_index, index] = 'WT'
    return matrix

def get_tick_marks(stats_dictionary, suffix=''):
    ticks = [
        stats_dictionary[f'wt_dna_score{suffix}'], 
        stats_dictionary[f'spike_in_score{suffix}'], 
        # stats_dictionary[f'syn_avg{suffix}'],
        # stats_dictionary[f'non_avg{suffix}'], 
        # stats_dictionary[f'mis_aa_avg{suffix}'],
        stats_dictionary[f'wt_dna_score{suffix}'] * 1.25, 
        stats_dictionary[f'wt_dna_score{suffix}'] * 0.5
    ]
    labels = [
        'wt', 
        'E501K',
        # 'syn (avg)', 
        # '* (avg)', 
        # 'miss (avg)',
        '1.25 wt', 
        '0.5 wt'
    ]

    return ticks, labels

def make_col_avg_df(heatmap_df):
    # Calculate the average of each column in the heatmap df and put the results in a new df
    col_avg = heatmap_df.iloc[:-1].mean() # exclude the last row '*'
    col_avg_df = pd.DataFrame(col_avg).transpose()
    return col_avg_df

def fill_wt(dms_matrix, wt_dna_score):
    # For each cell in the DataFrame grab the activity score or identify 'WT'
    heatmap_df = dms_matrix.map(extract_value)

    # Assign WT value
    heatmap_df = heatmap_df.replace('WT', wt_dna_score).astype(float)
    return heatmap_df

def apply_transparency(ax, indices, transparency=0.5):
    # Get the heatmap data from the Axes object
    heatmap = ax.collections[0]
    data = heatmap.get_array().reshape(heatmap.get_array().shape[0], -1)

    # Define the alpha matrix with the same shape as the data matrix
    alpha = np.ones_like(data)

    # Set the transparency for columns not in the specified indices
    for index in range(data.shape[1]):
        if index not in indices:
            alpha[:, index] = transparency  # Set transparency for these columns

    # Apply the alpha channel to the existing heatmap
    heatmap.set_alpha(alpha)

def dms_heatmap(dms_matrix, title, min_pos, max_pos, wt_score, 
                fig_size='small', export=False, output='activity_heatmap.png', 
                tick_values=[], tick_labels=[], motif_indices=None, row_avg=False):
    # calculate number and percent of AA variants without score
    dropout_num, dropout_percent = get_dropout(dms_matrix)

    # Masks
    nan_mask = dms_matrix.isnull()
    wt_mask = dms_matrix == 'WT'

    heatmap_df = fill_wt(dms_matrix, wt_score)

    # get values for column averages to be placed above the heatmap
    col_avg_df = make_col_avg_df(heatmap_df)

    # adjust the size of the figure and subplots
    if fig_size == 'small':
        fig = plt.figure(figsize=(16.5, 12))
        tick_freq = 2
    elif fig_size == 'large':
        fig = plt.figure(figsize=(30, 10))
        tick_freq = 5
    elif fig_size == 'long':
        fig = plt.figure(figsize=(30, 25))
        tick_freq = 5
    
    if row_avg:
        row_avg_df = pd.DataFrame(heatmap_df.mean(axis=1), columns=['Avg'])
        gs = GridSpec(2,3, width_ratios=[1, 35, 1], height_ratios=[1, 45], hspace=0.03, wspace=0.03)
        # Create the subplots using the GridSpec layout
        ax1 = fig.add_subplot(gs[0, 1])  # For the average row
        ax2 = fig.add_subplot(gs[1, 1])  # For the main heatmap
        cax = fig.add_subplot(gs[1, 2])  # For the colorbar
        ax3 = fig.add_subplot(gs[1, 0])
    else:
        gs = GridSpec(2, 2, width_ratios=[35, 1], height_ratios=[1, 20], hspace=0.03, wspace=0.03)
        # # Create the subplots using the GridSpec layout
        ax1 = fig.add_subplot(gs[0, 0])  # For the average row
        ax2 = fig.add_subplot(gs[1, 0])  # For the main heatmap
        cax = fig.add_subplot(gs[1, 1])  # For the colorbar

    # Get the min and max values from the heatmap dataframe
    min = heatmap_df.min().min()
    max = heatmap_df.max().max()
    print(f'Min: {min}, Max: {max}')

    # Use to match scales across subplots
    norm = plt.Normalize(vmin=min, vmax=max)

    # Magma theme - perceptually uniform colormap
    cmap = plt.cm.magma

    # --------------------------------
    # Create the row of positional averages
    sns.heatmap(col_avg_df, annot=False, cmap=cmap, cbar=False, ax=ax1, norm=norm)

    # Create the main heatmap
    ax = sns.heatmap(heatmap_df, cmap=cmap, cbar=False, ax=ax2)
    if nan_mask.any().any():
        sns.heatmap(nan_mask, annot=False, cmap=cmap, cbar=False, mask=~nan_mask, ax=ax2)
    if row_avg:
        sns.heatmap(row_avg_df, cmap=cmap, cbar=False, ax=ax3, norm=norm)
    # --------------------------------
    ### Denote Missing Values
    # Add hatches for values that are missing
    for i in range(len(nan_mask)):
        for j in range(len(nan_mask.columns)):
            if nan_mask.iloc[i, j]:
                ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=True, hatch='//', edgecolor='lightgray', facecolor='white'))
    # --------------------------------
    ### Highlight WT
    # Get the indices where wt_mask is True
    wt_indices = np.where(wt_mask)

    # Add a light dot at each 'WT' location
    ax2.scatter(wt_indices[1] + 0.5, wt_indices[0] + 0.5, color='white', s=30, alpha=0.5)

    ### Apply transparency
    if motif_indices:
        apply_transparency(ax, motif_indices)
    # --------------------------------
    ### Colorbar
    # Define new ticks and labels
    ticks = tick_values
    labels = tick_labels

    # Create the colorbar in the separate column
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, cax=cax, ticks=ticks)
    cbar.ax.set_yticklabels(labels, fontsize=18)
    # ------------------------------
    # Update the x-axis tick labels to match the codon position in GLI2
    x_labels = [str(i) for i in range(min_pos, max_pos+1)]

    ax.set_xticks([i + 0.5 for i in range(0, len(x_labels), tick_freq)])
    ax.set_xticklabels([x_labels[i] for i in range(0, len(x_labels), tick_freq)], rotation=0)

    #-------------------------------
    plot_title = f'{title} - Dropout {dropout_num} variant ({dropout_percent}%)'

    # Set the title with enough space above the subplots
    if fig_size == 'small':
        fig.suptitle(plot_title, fontsize=28, y=0.89)
    elif fig_size == 'large':
        fig.suptitle(plot_title, fontsize=30, y=0.91)

    fig.subplots_adjust(top=0.85, bottom=0.15)

    ax1.set_ylabel('Avg', fontsize=20, labelpad=8, ha='center')
    ax1.set_xticks([])  # Remove x-axis tick marks for the averages subplot
    ax1.set_yticks([])  # Remove y-axis tick marks for the averages subplot
    ax2.set_xlabel('Residue Sequence Number', fontsize=24)
    ax2.set_ylabel('Mutant Amino Acid', fontsize=24)
    ax2.tick_params(axis='both', which='major', labelsize=20)  # Set font size for y-axis tick labels
    if row_avg:
        ax3.set_ylabel('Mutant Amino Acid', fontsize=24, labelpad=8, ha='center')
        ax3.set_yticks(ax2.get_yticks())
        ax3.set_yticklabels(ax2.get_yticklabels())
        ax3.tick_params(axis='y', which='major', labelsize=20)  # Set font size for y-axis tick labels
        ax2.set_ylabel('')
        ax2.set_yticks([])  # Remove y-axis tick marks for the main heatmap
        ax3.set_xticks([])  # Remove x-axis tick marks for the row averages subplot
    #-------------------------------
    if export:
        # Save the heatmap as a PNG file
        plt.savefig(output, dpi=300, format='png', transparent=True)
    else:
        plt.show()


def dms_histogram(data, avgscore_col, export=False, output='activity_histogram.png'):

    # Determine the global range of the data to create consistent bins
    min_score = data[avgscore_col].min()
    max_score = data[avgscore_col].max()
    bins = np.linspace(min_score, max_score, 51)  # Creates 50 bins between the min and max

    unique_values = data['annotate_dna'].unique()
    colors = plt.cm.jet(np.linspace(0, 1, len(unique_values)))  # Using the jet colormap; you can choose another

    # spike-in' and 'wt_dna' annotation
    annotated_values = ['spike-in', 'wt_dna']  # Values to be plotted last
    other_values = [v for v in unique_values if v not in annotated_values]

    # Plot histograms for non-priority values first
    for value in other_values:
        color = colors[unique_values.tolist().index(value)]
        subset = data[data['annotate_dna'] == value]
        plt.hist(subset, bins=bins, color=color, alpha=0.5, label=str(value))

    # Plot spike-in and wt values last, so they appear on top
    for value in annotated_values:
        color = colors[unique_values.tolist().index(value)]
        subset = data[data['annotate_dna'] == value]
        plt.hist(subset, bins=bins, color=color, alpha=0.75, label=str(value), edgecolor='black')

    # Add titles and labels (optional)
    plt.title('Histogram of DNA level Scores')
    plt.xlabel('Activity Score')
    plt.ylabel('Frequency')
    plt.legend(title='annotate_dna')


    if export:
        # Save the heatmap as a PNG file
        plt.savefig(output, dpi=1200, format='png', transparent=True)
    else:
        plt.show()


    # TODO: frequency of each annotation in the 4 bins