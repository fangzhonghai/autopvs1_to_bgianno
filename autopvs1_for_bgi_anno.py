# -*- coding:utf-8 -*-
import pandas as pd
import optparse
import pyfaidx
import yaml
import sys
import os


def print_usage(option, opt, value, parser):
    usage_message = r"""
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
    python autopvs1_for_bgi_anno.py -c autopvs1.yaml --wkdir /path/to/work -i test.xlsx -o test \
    --file_format excel --out_file_format excel --pvs1 2
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
    """
    print(usage_message)
    sys.exit()


def yaml_read(yaml_file):
    with open(yaml_file, 'r') as y:
        yaml_dic = yaml.load(y, Loader=yaml.FullLoader)
    return yaml_dic


def bgi_anno_format(bed, reference, in_file_format, sheet_no, skip):
    # #CHROM, "Start", "Stop", "Ref", "Call"
    if in_file_format == "excel":
        bed_df = pd.read_excel(bed, sheet_name=sheet_no-1, skiprows=range(skip))
    else:
        bed_df = pd.read_csv(bed, dtype={"#Chr": str}, sep='\t', skiprows=range(skip))
    bed_5_cols_df = bed_df[["#Chr", "Start", "Stop", "Ref", "Call"]].copy()
    if len(bed_5_cols_df[bed_5_cols_df['#Chr'].str.startswith("chr")]):
        bed_5_cols_df['#CHROM'] = bed_5_cols_df['#Chr']
    else:
        bed_5_cols_df['#CHROM'] = "chr" + bed_5_cols_df['#Chr']
    bed_5_cols_df.loc[bed_5_cols_df['#CHROM'] == 'chrMT', '#CHROM'] = 'chrM_NC_012920.1'
    bed_5_cols_df['ID'] = "."
    bed_5_cols_df['QUAL'] = "."
    bed_5_cols_df['FILTER'] = "."
    bed_5_cols_df['INFO'] = "."
    bed_5_cols_df['MuType'] = 'delins'
    bed_5_cols_df.loc[bed_5_cols_df['Ref'] == ".", "MuType"] = 'ins'
    bed_5_cols_df.loc[bed_5_cols_df['Call'] == ".", "MuType"] = 'del'
    bed_5_cols_df.loc[(bed_5_cols_df['Ref'].map(len) == 1) & (bed_5_cols_df['Call'].map(len) == 1) & (bed_5_cols_df['Ref'] != '.')
                      & (bed_5_cols_df['Call'] != '.'), 'MuType'] = 'snp'
    bed_5_cols_df['POS'] = bed_5_cols_df['Stop']
    bed_5_cols_df.loc[bed_5_cols_df['MuType'] == 'del', 'POS'] = bed_5_cols_df.loc[bed_5_cols_df['MuType'] == 'del', 'Start']
    bed_5_cols_df.loc[bed_5_cols_df['MuType'] == 'delins', 'POS'] = bed_5_cols_df.loc[bed_5_cols_df['MuType'] == 'delins', 'Start'] + 1
    bed_5_cols_df['REF'] = bed_5_cols_df['Ref']
    bed_5_cols_df['ALT'] = bed_5_cols_df['Call']
    fa = pyfaidx.Fasta(reference)
    for i in range(bed_5_cols_df.shape[0]):
        if bed_5_cols_df.loc[i, 'MuType'] == 'ins':
            base = str(fa.get_seq(bed_5_cols_df.loc[i, '#CHROM'], bed_5_cols_df.loc[i, 'POS'], bed_5_cols_df.loc[i, 'POS'])).upper()
            bed_5_cols_df.loc[i, 'REF'] = base
            bed_5_cols_df.loc[i, 'ALT'] = base + bed_5_cols_df.loc[i, 'ALT']
        elif bed_5_cols_df.loc[i, 'MuType'] == 'del':
            base = str(fa.get_seq(bed_5_cols_df.loc[i, '#CHROM'], bed_5_cols_df.loc[i, 'POS'], bed_5_cols_df.loc[i, 'POS'])).upper()
            bed_5_cols_df.loc[i, 'ALT'] = base
            bed_5_cols_df.loc[i, 'REF'] = base + bed_5_cols_df.loc[i, 'REF']
    # bed_5_cols_df_rm = bed_5_cols_df[bed_5_cols_df['MuType'] != 'delins'].copy()
    a = bed_5_cols_df[["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]].copy()
    a.sort_values(by=["#CHROM", "POS"], ascending=True, inplace=True)
    b = bed_5_cols_df[["#Chr", "Start", "Stop", "Ref", "Call", "#CHROM", "POS", "REF", "ALT"]].copy()
    return a, b, bed_df


def write_vcf(bed2vcf, prefix, path):
    vcf_file = os.path.join(path, prefix + ".vcf")
    with open(vcf_file, 'w') as f:
        text = r'''##fileformat=VCFv4.2
'''.format(**locals())
        f.write(text)
    bed2vcf.to_csv(vcf_file, mode='a', index=False, sep='\t')


def run_autopvs1(yaml_dic, prefix, path, pvs1_res):
    vcf_file = os.path.join(path, prefix + ".vcf")
    vep_file = vcf_file + ".vep"
    if yaml_dic['lof'] == 'lof':
        lof_file = vep_file + ".lof"
    else:
        lof_file = vep_file + ".all"
    pvs1_file = lof_file + ".autopvs1"
    run_vep = 'sh ' + yaml_dic['vep_sh'] + ' ' + vcf_file + ' > ' + vep_file
    status = os.system(run_vep)
    if status != 0:
        sys.exit(1)
    run_vep_lof = yaml_dic['py3'] + ' ' + yaml_dic['vep_lof'] + ' ' + vep_file
    status = os.system(run_vep_lof)
    if status != 0:
        sys.exit(1)
    run_vep_lof_auto = yaml_dic['py3'] + ' ' + yaml_dic['autopvs1'] + ' ' + lof_file + ' > ' + pvs1_file
    status = os.system(run_vep_lof_auto)
    if status != 0:
        sys.exit(1)
    pvs1_df = pd.read_csv(pvs1_file, sep='\t', header=None)
    if pvs1_res == 2:
        pvs1_df.rename(columns={5: 'AutoPvs1 by decision tree', 6: 'AutoPvs1 by disease mechanism'}, inplace=True)
    else:
        pvs1_df.rename(columns={5: 'AutoPvs1'}, inplace=True)
    pvs1_df_split = pvs1_df[0].str.split('-', expand=True)
    pvs1_df_split.columns = ['#CHROM', 'POS', 'REF', 'ALT']
    pvs1_df_split['#CHROM'] = 'chr' + pvs1_df_split['#CHROM']
    pvs1_df_split['POS'] = pvs1_df_split['POS'].astype('int')
    pvs1_df_split['REF'] = pvs1_df_split['REF'].str.upper()
    pvs1_df_split['ALT'] = pvs1_df_split['ALT'].str.upper()
    if pvs1_res == 2:
        pvs1_df_final = pvs1_df_split.join(pvs1_df[['AutoPvs1 by decision tree', 'AutoPvs1 by disease mechanism']])
    else:
        pvs1_df_final = pvs1_df_split.join(pvs1_df[['AutoPvs1']])
    return pvs1_df_final


if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('-u', '--usage', help='print more info on how to use this script', action="callback", callback=print_usage)
    parser.add_option('-i', '--in', dest='in_file', help='input file', default=None, metavar='file')
    parser.add_option('-c', '--config', dest='config', help='config file', default=None, metavar='file')
    parser.add_option('--skip_row', dest='skip_row', default=0, type=int)
    parser.add_option('--sheet', dest='sheet', default=1, type=int)
    parser.add_option('--file_format', dest='file_format', default='tsv', type='string')
    parser.add_option('-o', '--out', dest='out', help='output file prefix', default=None, metavar='string')
    parser.add_option('--pvs1', dest='pvs1', default=1, type=int)
    parser.add_option('--out_file_format', dest='out_file_format', help='output file format', default='tsv', metavar='string')
    parser.add_option('--wkdir', dest='wkdir', default=None)
    (opts, args) = parser.parse_args()
    in_file = opts.in_file
    file_format = opts.file_format
    sheet = opts.sheet
    skip_row = opts.skip_row
    out = opts.out
    wkdir = opts.wkdir
    config = opts.config
    pvs1_num = opts.pvs1
    out_file_format = opts.out_file_format
    config_dic = yaml_read(config)
    ref = config_dic['ref']
    bed_2_vcf, bed_vcf, bgi_bed = bgi_anno_format(in_file, ref, file_format, sheet, skip_row)
    write_vcf(bed_2_vcf, out, wkdir)
    pvs1 = run_autopvs1(config_dic, out, wkdir, pvs1_num)
    merge_1 = pd.merge(bed_vcf, pvs1, on=['#CHROM', 'POS', 'REF', 'ALT'])
    merge_2 = pd.merge(bgi_bed, merge_1, on=["#Chr", "Start", "Stop", "Ref", "Call"], how='left')
    if pvs1_num == 2:
        merge_2.fillna(value={'AutoPvs1 by decision tree': '.', 'AutoPvs1 by disease mechanism': '.'}, inplace=True)
    else:
        merge_2.fillna(value={'AutoPvs1': '.'}, inplace=True)
    merge_2.drop(columns=['#CHROM', 'POS', 'REF', 'ALT'], inplace=True)
    if out_file_format == 'excel':
        merge_2.to_excel(wkdir + "/" + out + ".vcf.vep." + config_dic['lof'] + ".autopvs1.bgi_anno.xlsx", index=False)
    else:
        merge_2.to_csv(wkdir + "/" + out + ".vcf.vep." + config_dic['lof'] + ".autopvs1.bgi_anno.bed", sep='\t', index=False)
