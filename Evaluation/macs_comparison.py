import pandas as pd
from intervaltree import Interval, IntervalTree
#%%


def read_encode(file_path):

    df = pd.read_csv(file_path)
    header = df.iloc[23]
    new_df = df.iloc[24:]
    new_df.columns = header.map(str)
    new_df.reset_index(drop=True,inplace=True)
    return new_df


def find_overlaps_optimized(df1, df2):
    # Normalize column names and ensure data types are correct
    df1 = df1.rename(columns={'chromosome': 'chr'})
    df2 = df2.rename(columns={'chr': 'chr'})
    df1['start'] = df1['start'].astype(int)
    df1['end'] = df1['end'].astype(int)
    df2['start'] = df2['start'].astype(int)
    df2['end'] = df2['end'].astype(int)

    # Ensure both dataframes are sorted by chromosome first, then start
    df1 = df1.sort_values(by=['chr', 'start'])
    df2 = df2.sort_values(by=['chr', 'start'])

    overlaps = []
    i, n = 0, len(df1)
    # Iterate over each row in df2
    count=0
    for idx2, row2 in df2.iterrows():
        count+=1
        if count%1000==0:
            print(count)
        # Advance the index i in df1 until we find overlapping ranges
        while i < n and (df1.iloc[i]['chr'] < row2['chr'] or (df1.iloc[i]['chr'] == row2['chr'] and df1.iloc[i]['end'] < row2['start'])):
            i += 1
        # Check for overlaps in df1 starting from the current index i
        j = i
        while j < n and df1.iloc[j]['chr'] == row2['chr'] and df1.iloc[j]['start'] <= row2['end']:
            if df1.iloc[j]['end'] >= row2['start']:
                overlaps.append({
                    'chr': row2['chr'],
                    'start_df1': df1.iloc[j]['start'],
                    'end_df1': df1.iloc[j]['end'],
                    'start_df2': row2['start'],
                    'end_df2': row2['end']
                })
            j += 1

    # Convert the list of dictionaries to a DataFrame
    overlaps_df = pd.DataFrame(overlaps)
    return overlaps_df





#%%
macs_control_path = "ENCODE Data/eClip_control_and_target/Control/control_output_peaks.csv"
macs_control_data = read_encode(macs_control_path)
our_data_path = "ENCODE Data/eClip_control_and_target/our_control_peaks.csv"
our_control_data = pd.read_csv(our_data_path)
our_control_data = our_control_data[our_control_data['significant']==True]
#%%

overlap_df = find_overlaps_optimized(our_control_data,macs_control_data)