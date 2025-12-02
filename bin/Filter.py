import pandas as pd
import io
import argparse
import sys
import math
from typing import List, Dict, Optional, Union
from tqdm.auto import tqdm

# --- Helper Functions ---

def load_bedgraph(filepath_or_buffer: Union[str, io.StringIO]) -> pd.DataFrame:
    col_names = ['sequence', 'start', 'end', 'coverage', '3prime', '5prime']
    try:
        df = pd.read_csv(
            filepath_or_buffer,
            sep='\t',
            header=None,
            names=col_names,
            comment='#' 
        )
    except FileNotFoundError:
        print(f"Error: Input file not found at {filepath_or_buffer}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error loading bedgraph file: {e}", file=sys.stderr)
        sys.exit(1)

    try:
        df['start'] = df['start'].astype(int)
        df['end'] = df['end'].astype(int)
        df['coverage'] = df['coverage'].astype(int)
        df['3prime'] = df['3prime'].astype(int)
        df['5prime'] = df['5prime'].astype(int)
    except ValueError as e:
        print(f"Error: Could not convert columns to numbers. {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Loaded {len(df)} intervals from bedgraph.")
    return df

def find_high_coverage_regions(df: pd.DataFrame, min_coverage: int = 100) -> pd.DataFrame:
    high_cov_df = df[df['coverage'] >= min_coverage].copy()
    print(f"Found {len(high_cov_df)} intervals with coverage >= {min_coverage}.")
    return high_cov_df

def merge_nearby_regions(df: pd.DataFrame, max_distance: int = 10) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=['sequence', 'start', 'end'])

    df_sorted = df.sort_values(['sequence', 'start'])
    merged_regions = []
    
    for seq, group in df_sorted.groupby('sequence'):
        if group.empty: continue
            
        current_start = group.iloc[0]['start']
        current_end = group.iloc[0]['end']
        
        for i in range(1, len(group)):
            next_region = group.iloc[i]
            gap = next_region['start'] - current_end
            
            if gap <= max_distance:
                current_end = max(current_end, next_region['end'])
            else:
                merged_regions.append({'sequence': seq, 'start': current_start, 'end': current_end})
                current_start = next_region['start']
                current_end = next_region['end']
        
        merged_regions.append({'sequence': seq, 'start': current_start, 'end': current_end})
        
    merged_df = pd.DataFrame(merged_regions)
    print(f"Merged regions into {len(merged_df)} clusters (merge_distance={max_distance}).")
    return merged_df

# --- EXTEND INSTEAD OF FILTER ---

def extend_short_regions_to_min_length(df: pd.DataFrame, min_length: int = 70) -> List[Dict]:
    """
    Checks the length of every merged region.
    If length < min_length, it extends the start and end symmetrically
    to reach the min_length.
    """
    if df.empty:
        return []

    extended_regions = []
    
    # We iterate over rows to calculate extensions
    for _, row in df.iterrows():
        seq = row['sequence']
        start = int(row['start'])
        end = int(row['end'])
        current_len = end - start
        
        if current_len < min_length:
            deficit = min_length - current_len
            # Split the deficit between left and right
            expand_left = math.floor(deficit / 2)
            expand_right = math.ceil(deficit / 2)
            
            new_start = start - expand_left
            new_end = end + expand_right
            
            # Ensure we don't go below 0
            if new_start < 0:
                # If we hit 0 on the left, add the remainder to the right
                diff = 0 - new_start
                new_start = 0
                new_end += diff
            
            extended_regions.append({
                'sequence': seq,
                'region_start': new_start,
                'region_end': new_end
            })
        else:
            # It's already long enough
            extended_regions.append({
                'sequence': seq,
                'region_start': start,
                'region_end': end
            })
            
    print(f"Processed {len(extended_regions)} regions. Short regions were extended to at least {min_length}bp.")
    return extended_regions

# --- EXPAND TO PER-BASE (With Context Padding) ---

def expand_regions_to_perbase(df_bedgraph: pd.DataFrame, selected_regions: List[dict],
                              pad_size: int = 50) -> List[pd.DataFrame]:
    if not selected_regions:
        return []
    
    print(f"\nExpanding regions with context padding (pad={pad_size}bp)...")
    
    all_expanded = []
    regions_df = pd.DataFrame(selected_regions)
    regions_df = regions_df.sort_values(['sequence', 'region_start']).reset_index(drop=True)
    
    for seq_name, seq_regions in tqdm(regions_df.groupby('sequence'), desc="Expanding sequences"):
        seq_regions = seq_regions.reset_index(drop=True)
        seq_intervals = df_bedgraph[df_bedgraph['sequence'] == seq_name].copy()
        
        for idx, region_spec in seq_regions.iterrows():
            region_start = int(region_spec['region_start'])
            region_end   = int(region_spec['region_end'])
            
            prev_end = seq_regions.iloc[idx - 1]['region_end'] if idx > 0 else -1
            next_start = seq_regions.iloc[idx + 1]['region_start'] if idx < len(seq_regions) - 1 else float('inf')
            
            pad_start = int(region_start - pad_size)
            pad_end   = int(region_end + pad_size)

            # FIX: Compare using max/min to avoid casting float('inf') to int
            pad_start = max(pad_start, prev_end)
            pad_end   = min(pad_end, next_start)
            pad_start = max(0, pad_start)

            overlapping = seq_intervals[
                (seq_intervals['start'] < pad_end) &
                (seq_intervals['end'] > pad_start)
            ]
            
            expanded_rows = []
            for pos in range(pad_start, pad_end):
                containing = overlapping[
                    (overlapping['start'] <= pos) &
                    (overlapping['end'] > pos)
                ]
                
                if not containing.empty:
                    row = containing.iloc[0]
                    coverage = row['coverage']
                    prime3 = row['3prime']
                    prime5 = row['5prime']
                else:
                    coverage = 0
                    prime3 = 0
                    prime5 = 0
                
                expanded_rows.append({
                    'sequence': seq_name,
                    'position': pos, 
                    'coverage': coverage,
                    '3prime': prime3,
                    '5prime': prime5
                })
            
            if expanded_rows:
                region_df = pd.DataFrame(expanded_rows)
                # FIX: Changed naming format as requested
                region_df['cluster_id'] = f"{seq_name}_region_{pad_start}_{pad_end}"
                all_expanded.append(region_df)
    
    print(f"Expanded to {len(all_expanded)} regions.")
    return all_expanded

# --- FORMAT ---

def format_expanded_regions(expanded_dfs_list: List[pd.DataFrame]) -> pd.DataFrame:
    if not expanded_dfs_list:
        return pd.DataFrame(columns=[
            'RegionID', 'sequence', 'ori_position', 'position', 
            'coverage', '3prime', '5prime'
        ])
        
    print("Formatting final DataFrame...")
    final_df = pd.concat(expanded_dfs_list, ignore_index=True)
    
    final_df = final_df.rename(columns={
        'cluster_id': 'RegionID',
        'position': 'ori_position' 
    })
    
    final_df = final_df.sort_values(['RegionID', 'ori_position'])
    final_df['position'] = final_df.groupby('RegionID').cumcount() + 1
    
    final_cols = [
        'RegionID', 'sequence', 'ori_position', 'position', 
        'coverage', '3prime', '5prime'
    ]
    return final_df[final_cols]

# --- PIPELINE ---

def process_bedgraph_pipeline(
    filepath_or_buffer: Union[str, io.StringIO],
    min_coverage: int = 100,
    merge_distance: int = 10,
    min_length: int = 70,
    pad_size: int = 50
) -> pd.DataFrame:
    print("--- Starting Bedgraph Processing Pipeline ---")
    
    # 1. Load
    df_bedgraph = load_bedgraph(filepath_or_buffer)
    
    # 2. Find High Coverage
    df_high_cov = find_high_coverage_regions(df_bedgraph, min_coverage)
    
    # 3. Merge Nearby
    df_merged = merge_nearby_regions(df_high_cov, merge_distance)
    
    # 4. Short Regions (Logic Change Here)
    # This ensures every region is at least 'min_length' (70bp)
    selected_regions_list = extend_short_regions_to_min_length(df_merged, min_length)
    
    if not selected_regions_list:
        print("\nNo regions found.")
        return format_expanded_regions([])
        
    # 5. Add Context Padding & Expand
    # This adds the 'pad_size' (50bp) on top of the region boundaries
    expanded_dfs_list = expand_regions_to_perbase(
        df_bedgraph, 
        selected_regions_list, 
        pad_size
    )
    
    # 6. Format
    final_output_df = format_expanded_regions(expanded_dfs_list)
    
    print("--- Pipeline Complete ---")
    return final_output_df

def main():
    parser = argparse.ArgumentParser(description="Find and expand high-coverage regions.")
    parser.add_argument("input_file", type=str, help="Input bedgraph file")
    parser.add_argument("output_file", type=str, help="Output CSV file")
    parser.add_argument("-c", "--min-coverage", type=int, default=20, help="Default: 20")
    parser.add_argument("-d", "--merge-distance", type=int, default=10, help="Default: 10")
    parser.add_argument("-l", "--min-length", type=int, default=70, help="Min length to EXTEND short regions to (default: 70)")
    parser.add_argument("-p", "--pad-size", type=int, default=50, help="Extra context padding (default: 50)")
    
    args = parser.parse_args()
    
    final_df = process_bedgraph_pipeline(
        args.input_file,
        min_coverage=args.min_coverage,
        merge_distance=args.merge_distance,
        min_length=args.min_length,
        pad_size=args.pad_size
    )
    
    if not final_df.empty:
        final_df.to_csv(args.output_file, sep=',', index=False)
        print(f"\nSuccessfully saved {len(final_df)} rows to {args.output_file}")
    else:
        print("\nNo regions found.")

if __name__ == "__main__":
    main()