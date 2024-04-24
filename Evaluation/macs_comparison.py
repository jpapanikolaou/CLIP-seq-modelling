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


def build_confusion_matrix(df1,df2,overlap_df):
   agreed = len(overlap_df['start_df2'])
   in_df1_not_df2 = len(df1) - agreed
   in_df2_not_df1 = len(df2) - agreed

   return{
        "agreed": agreed,
        "in_df1_not_df2": in_df1_not_df2,
        "in_df2_not_df1": in_df2_not_df1
   }