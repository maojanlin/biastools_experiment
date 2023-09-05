import argparse
import pandas as pd
import numpy as np






if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    #parser.add_argument('-bl', '--list_bed', nargs='+', default=[], help='one or a the list of scanning bed files')
    parser.add_argument('-bed', '--bias_region_bed_file', help='the scanning bed file report')
    parser.add_argument('-tsv', '--bias_tsv_report', help='the simulation bias file')
    args = parser.parse_args()
    
    fn_bed = args.bias_region_bed_file
    fn_tsv = args.bias_tsv_report
    
    pd_bed = pd.read_csv(fn_bed, sep=' ', header=0, index_col=False)
    pd_tsv = pd.read_csv(fn_tsv, sep='\t', header=0)

    assert(pd_bed["#chrom"][0] == pd_tsv["CHR"][0])
    bias_region = list(zip(pd_bed["chromStart"], pd_bed["chromEnd"]))
    list_in = []
    list_out = []
    for site in pd_tsv["HET_SITE"]:
        flag_out = True
        for l_bd, r_bd in bias_region:
            if site >= l_bd and site <= r_bd:
                list_in.append(site)
                flag_out = False
                break
        if flag_out:
            list_out.append(site)

    total_num = len(pd_tsv["HET_SITE"].values)
    print("# Total sites:", total_num)
    print("# sites inside bias region:", len(list_in), len(list_in)/total_num)
    print("# sites outside bias region", len(list_out), len(list_out)/total_num)

