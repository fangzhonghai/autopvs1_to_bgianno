# -*- coding:utf-8 -*-
import pandas as pd
import optparse
import sys


def print_usage(option, opt, value, parser):
    usage_message = r"""
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
    python vep_autopvs1_2_bgi_anno.py -i pvs1.pilot.vep.lof.autopvs1 \
    -o pvs1.pilot.vep.lof.autopvs1.bgi_anno.bed
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
    """
    print(usage_message)
    sys.exit()


def format_trans(input_file):
    df = pd.read_csv(input_file, sep='\t', header=None)
    df.rename(columns={0: 'Name'}, inplace=True)
    df_split = df['Name'].str.split('-', expand=True)
    df_split.columns = ['#CHROM', 'POS', 'REF', 'ALT']
    df_split['Name'] = df['Name']
    df_split['MuType'] = 'delins'
    df_split.loc[(df_split['REF'].str.len() == 1) & (df_split['ALT'].str.len() == 1), 'MuType'] = 'snv'
    df_split.loc[(df_split['REF'].str.len() == 1) & (df_split['ALT'].str.len() > 1), 'MuType'] = 'ins'
    df_split.loc[(df_split['REF'].str.len() > 1) & (df_split['ALT'].str.len() == 1), 'MuType'] = 'del'
    df_split.loc[(df_split['MuType'] == 'ins') & (df_split['REF'].str[0] != df_split['ALT'].str[0]), 'MuType'] = 'delins'
    df_split.loc[(df_split['MuType'] == 'del') & (df_split['REF'].str[0] != df_split['ALT'].str[0]), 'MuType'] = 'delins'
    df_split['#CHROM'] = df_split['#CHROM'].astype('str')
    if len(df_split[df_split['#CHROM'].str.startswith('chr')]):
        df_split['#Chr'] = df_split['#CHROM']
    else:
        df_split['#Chr'] = 'chr' + df_split['#CHROM']
    df_split.loc[df_split['#CHROM'] == 'chrM_NC_012920.1', '#Chr'] = 'chrMT'
    df_split['Start'] = df_split['POS'].astype(int) - 1
    df_split['Stop'] = df_split['POS'].astype(int)
    df_split['Ref'] = df_split['REF']
    df_split['Call'] = df_split['ALT']
    df_split.loc[df_split['MuType'] == 'ins', 'Ref'] = '.'
    df_split.loc[df_split['MuType'] == 'ins', 'Call'] = df_split.loc[df_split['MuType'] == 'ins', 'ALT'].str[1:]
    df_split.loc[df_split['MuType'] == 'del', 'Call'] = '.'
    df_split.loc[df_split['MuType'] == 'del', 'Ref'] = df_split.loc[df_split['MuType'] == 'del', 'REF'].str[1:]
    df_split.loc[df_split['MuType'] == 'ins', 'Start'] = df_split.loc[df_split['MuType'] == 'ins', 'Stop']
    df_split.loc[df_split['MuType'] == 'del', 'Start'] = df_split.loc[df_split['MuType'] == 'del', 'Stop']
    df_split.loc[df_split['MuType'] == 'del', 'Stop'] = df_split.loc[df_split['MuType'] == 'del', 'Stop'] + df_split.loc[df_split['MuType'] == 'del', 'Ref'].str.len()
    df_split.loc[df_split['MuType'] == 'delins', 'Stop'] = df_split.loc[df_split['MuType'] == 'delins', 'Start'] + df_split.loc[df_split['MuType'] == 'delins', 'Ref'].str.len()
    df_merge = pd.merge(df, df_split[['Name', '#Chr', 'Start', 'Stop', 'Ref', 'Call', 'MuType']], on=['Name'])
    return df_merge


if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('-u', '--usage', help='print more info on how to use this script', action="callback", callback=print_usage)
    parser.add_option('-i', '--in', dest='in_file', help='input file', default=None, metavar='file')
    parser.add_option('-o', '--out', dest='out', help='output file', default=None, metavar='string')
    (opts, args) = parser.parse_args()
    in_file = opts.in_file
    out = opts.out
    trans_df = format_trans(in_file)
    trans_df.to_csv(out, sep='\t', index=False, header=None)
