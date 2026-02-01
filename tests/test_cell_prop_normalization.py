#!/usr/bin/env python3
"""
Test script for cell proportion normalization feature.

This script demonstrates how to use the new cell proportion normalization
by creating mock data and testing the normalization calculations.
"""

import pandas as pd
import numpy as np
from sortscore.analysis.score import calculate_full_activity_scores

def test_cell_prop_normalization():
    """Test cell proportion normalization with mock data."""
    
    # Create mock count data
    counts = {
        1: {  # Replicate 1
            'A': pd.DataFrame({
                'variant_seq': ['seq1', 'seq2', 'seq3'],
                'count': [100, 200, 150]
            }),
            'B': pd.DataFrame({
                'variant_seq': ['seq1', 'seq2', 'seq3'], 
                'count': [50, 100, 75]
            })
        }
    }
    
    # Mock median GFP values
    median_gfp = {
        1: {'A': 100.0, 'B': 500.0}
    }
    
    # Mock total reads (for sequencing depth normalization)
    total_reads = {
        1: {'A': 1000, 'B': 1000}
    }
    
    # Mock cell proportions - bin A gets 20% of cells, bin B gets 10% of cells
    cell_prop = {
        1: {'A': 0.2, 'B': 0.1}
    }
    
    # Create merged DataFrame
    all_variants = ['seq1', 'seq2', 'seq3']
    merged_df = pd.DataFrame({'variant_seq': all_variants})
    
    # Add count columns
    merged_df['count.r1bA'] = [100, 200, 150]
    merged_df['count.r1bB'] = [50, 100, 75]
    
    print("=== Testing Cell Proportion Normalization ===")
    print("\nInput data:")
    print("Counts:")
    print(f"  Bin A: {[100, 200, 150]} (total reads: 1000, cell prop: 0.2)")
    print(f"  Bin B: {[50, 100, 75]} (total reads: 1000, cell prop: 0.1)")
    print("\nExpected normalization:")
    print("  Without cell prop: reads_per_million = count / total_reads * 1e6")
    print("  With cell prop: normalized = reads_per_million / cell_proportion")
    print(f"  Bin A seq1: (100/1000*1e6)/0.2 = {(100/1000*1e6)/0.2}")
    print(f"  Bin B seq1: (50/1000*1e6)/0.1 = {(50/1000*1e6)/0.1}")
    
    # Test without cell proportion normalization
    print("\n--- Test 1: Without cell proportion normalization ---")
    scores_no_cell_prop = calculate_full_activity_scores(
        counts=counts,
        median_gfp=median_gfp,
        merged_df=merged_df.copy(),
        min_bins=1,
        min_reps=1,
        total_reads=total_reads,
        cell_prop=None  # No cell prop normalization
    )
    
    print("Normalized counts (reads per million):")
    print(f"  norm.count.r1bA seq1: {scores_no_cell_prop.loc[0, 'norm.count.r1bA']}")
    print(f"  norm.count.r1bB seq1: {scores_no_cell_prop.loc[0, 'norm.count.r1bB']}")
    
    # Test with cell proportion normalization  
    print("\n--- Test 2: With cell proportion normalization ---")
    scores_with_cell_prop = calculate_full_activity_scores(
        counts=counts,
        median_gfp=median_gfp,
        merged_df=merged_df.copy(),
        min_bins=1,
        min_reps=1,
        total_reads=total_reads,
        cell_prop=cell_prop  # Include cell prop normalization
    )
    
    print("Normalized counts (with cell proportion correction):")
    print(f"  norm.count.r1bA seq1: {scores_with_cell_prop.loc[0, 'norm.count.r1bA']}")
    print(f"  norm.count.r1bB seq1: {scores_with_cell_prop.loc[0, 'norm.count.r1bB']}")
    
    # Verify the calculation
    expected_a = (100/1000*1e6)/0.2  # 500000
    expected_b = (50/1000*1e6)/0.1   # 500000
    
    actual_a = scores_with_cell_prop.loc[0, 'norm.count.r1bA']
    actual_b = scores_with_cell_prop.loc[0, 'norm.count.r1bB']
    
    print(f"\nVerification:")
    print(f"  Expected bin A: {expected_a}, Actual: {actual_a}")
    print(f"  Expected bin B: {expected_b}, Actual: {actual_b}")
    print(f"  Bin A correct: {abs(actual_a - expected_a) < 1e-6}")
    print(f"  Bin B correct: {abs(actual_b - expected_b) < 1e-6}")
    
    print(f"\nWith cell proportion normalization, seq1 has equal normalized counts")
    print(f"in both bins despite different raw counts and cell proportions.")
    print(f"This corrects for the fact that bin B received fewer cells (10% vs 20%).")
    
    return True

if __name__ == "__main__":
    test_cell_prop_normalization()
    print("\nCell proportion normalization test completed!")