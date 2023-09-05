import argparse
import seaborn as sns
import matplotlib.pyplot as plt

import math
import numpy as np
import pysam



def read_bias_report(fn_bias_report):
    list_bias_SNP = []
    list_bias_gap = []
    f = open(fn_bias_report, 'r')
    header = f.readline()
    for line in f:
        fields = line.split()
        if fields[-1] == '.':
            list_bias_gap.append(fields)
        else:
            list_bias_SNP.append(fields)
    f.close()
    return list_bias_SNP, list_bias_gap


def calculate_SNP_balance(assign_SNP, naive_SNP):
    """
    [[simulate_balance], [map_balance], [assign_balance], [naive_balance]]
    """
    record = [[],[],[],[],[]]
    for idx in range(len(assign_SNP)):
        record[0].append(float(assign_SNP[idx][14]))
        record[1].append(float(assign_SNP[idx][10]))
        record[2].append(float(assign_SNP[idx][5]))
        record[3].append(float(naive_SNP[idx][14]))
        record[4].append(float(naive_SNP[idx][5]))
    return record


def calculate_gap_balance(assign_gap, naive_gap, f_vcf, len_bd):
    list_insert = [ [] for _ in range(len_bd) ]
    list_delete = [ [] for _ in range(len_bd) ]
    for idx in range(len(assign_gap)):
        ref_name  = assign_gap[idx][0]
        var_start = int(assign_gap[idx][1])
        var_segment = f_vcf.fetch(contig=ref_name, start=var_start-1, stop=var_start+1)
        for var in var_segment:
            if var.start+1 != var_start:
                continue
            len_ref = len(var.ref)
            if len(var.alts) == 1:
                len_alt = len(var.alts[0])
            else:
                hap = var.samples[0]['GT']
                if hap[0] != 0:
                    len_alt = len(var.alts[hap[0]-1])
                else:
                    len_alt = len(var.alts[hap[1]-1])
            
            if len_ref > len_alt: # deletion
                diff = min(len_ref - len_alt -1, len_bd-1)
                record = [float(assign_gap[idx][14]), float(assign_gap[idx][10]), float(assign_gap[idx][5]), \
                          float(naive_gap[idx][14]), float(naive_gap[idx][5])]
                list_delete[diff].append(record)
            else:
                diff = min(len_alt - len_ref -1, len_bd-1)
                record = [float(assign_gap[idx][14]), float(assign_gap[idx][10]), float(assign_gap[idx][5]), \
                          float(naive_gap[idx][14]), float(naive_gap[idx][5])]
                list_insert[diff].append(record)
    return list_insert, list_delete


def addlabels(x,y,len_bd):
    for i in range(len(x)):
        plt.text(i-len_bd, y[i], y[i], ha = 'center', va = 'bottom')


def plot_balance(balance_delete, balance_SNP, balance_insert, output_name, len_bd, list_incidents):
    balance_list = [np.zeros(2*len_bd+1), np.zeros(2*len_bd+1), np.zeros(2*len_bd+1), np.zeros(2*len_bd+1), np.zeros(2*len_bd+1)]
    balance_1st  = [np.zeros(2*len_bd+1), np.zeros(2*len_bd+1), np.zeros(2*len_bd+1), np.zeros(2*len_bd+1), np.zeros(2*len_bd+1)]
    balance_3rd  = [np.zeros(2*len_bd+1), np.zeros(2*len_bd+1), np.zeros(2*len_bd+1), np.zeros(2*len_bd+1), np.zeros(2*len_bd+1)]
    for idx in range(len_bd):
        quadra_list = balance_delete[idx]
        np_quadra_list = np.array(quadra_list)
        if len(np_quadra_list) > 0:
            for idy in range(5):
                list_balance = np_quadra_list[:,idy]
                median = np.median(list_balance[~np.isnan(list_balance)])
                balance_list[idy][len_bd-1-idx] = median
                balance_1st [idy][len_bd-1-idx] = median - np.quantile(list_balance[~np.isnan(list_balance)], 0.25)
                balance_3rd [idy][len_bd-1-idx] = np.quantile(list_balance[~np.isnan(list_balance)], 0.75) - median
        else:
            for idy in range(5):
                balance_list[idy][len_bd-1-idx] = np.nan
                balance_1st [idy][len_bd-1-idx] = np.nan
                balance_3rd [idy][len_bd-1-idx] = np.nan
    for idy, list_balance in enumerate(np.array(balance_SNP)):
        median = np.median(list_balance[~np.isnan(list_balance)])
        balance_list[idy][len_bd] = median
        balance_1st [idy][len_bd] = median - np.quantile(list_balance[~np.isnan(list_balance)], 0.25)
        balance_3rd [idy][len_bd] = np.quantile(list_balance[~np.isnan(list_balance)], 0.75) - median
    for idx, quadra_list in enumerate(balance_insert):
        np_quadra_list = np.array(quadra_list)
        if len(np_quadra_list) > 0:
            for idy in range(5):
                list_balance = np_quadra_list[:,idy]
                median = np.median(list_balance[~np.isnan(list_balance)])
                balance_list[idy][idx+len_bd+1] = median
                balance_1st [idy][idx+len_bd+1] = median - np.quantile(list_balance[~np.isnan(list_balance)], 0.25)
                balance_3rd [idy][idx+len_bd+1] = np.quantile(list_balance[~np.isnan(list_balance)], 0.75) - median
        else:
            for idy in range(5):
                balance_list[idy][idx+len_bd+1] = np.nan
                balance_1st[idy][idx+len_bd+1] = np.nan
                balance_3rd[idy][idx+len_bd+1] = np.nan
    
    """
    t = list(range(-len_bd,len_bd+1))
    plt.figure(figsize = (14,5))
    plt.ylim([-0.3, 0.2])
    plt.errorbar(t, (1-balance_list[0]) - (1-balance_list[0]), yerr=(balance_1st[0], balance_3rd[0]), capsize=3, fmt='-o', label="Simulation")
    plt.errorbar(t, (1-balance_list[1]) - (1-balance_list[0]), yerr=(balance_1st[1], balance_3rd[1]), capsize=3, fmt='-o', label="Mapping")
    plt.errorbar(t, (1-balance_list[2]) - (1-balance_list[0]), yerr=(balance_1st[2], balance_3rd[2]), capsize=3, fmt='-o', label="Context-aware assignment")
    plt.errorbar(t, (1-balance_list[4]) - (1-balance_list[3]), yerr=(balance_1st[4], balance_3rd[4]), capsize=3, fmt='-o', label="Naive assignment")
    plt.legend()
    plt.xlabel('Insertion (+) or deletion (-) length')
    plt.ylabel('Fraction of alternate allele')
    plt.axhline(y=0.5, color='gray', linestyle='dashdot', linewidth=0.9)
    if output_name:
        plt.savefig(output_name + '.indel_balance.pdf')
    else:
        plt.show()
    """

    t = list(range(-len_bd,len_bd+1))
    f, (a0, a1, a2) = plt.subplots(3, 1, gridspec_kw={'height_ratios': [3,3,1]})
    f.set_size_inches(15,12)
    #for idx, name in enumerate(list_plot_name):
    #    a0.errorbar(t, (1-balance_list[idx]), yerr=(balance_1st[idx], balance_3rd[idx]), capsize=3, fmt='-o', label=name)
    a0.errorbar(t, (1-balance_list[0]), yerr=(balance_1st[0], balance_3rd[0]), capsize=3, fmt='-o', label="Simulation")
    a0.errorbar(t, (1-balance_list[1]), yerr=(balance_1st[1], balance_3rd[1]), capsize=3, fmt='-o', label="Mapping")
    a0.errorbar(t, (1-balance_list[2]), yerr=(balance_1st[2], balance_3rd[2]), capsize=3, fmt='-o', label="Context-aware assignment")
    a0.errorbar(t, (1-balance_list[4]), yerr=(balance_1st[4], balance_3rd[4]), capsize=3, fmt='-o', label="Naive assignment")
    a0.axhline(y=0.5, color='gray', linestyle='dashdot', linewidth=0.9)
    a0.legend()
    a0.set(ylabel='Fraction of alternate allele')
    #a1.errorbar(t, (1-balance_list[0]) - (1-balance_list[0]), yerr=(balance_1st[0], balance_3rd[0]), capsize=3, fmt='-o', label="Simulation")
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    a1.errorbar(t, (1-balance_list[1]) - (1-balance_list[0]), yerr=(balance_1st[1], balance_3rd[1]), capsize=3, fmt='-o', color=colors[1], label="Mapping")
    a1.errorbar(t, (1-balance_list[2]) - (1-balance_list[0]), yerr=(balance_1st[2], balance_3rd[2]), capsize=3, fmt='-o', color=colors[2], label="Context-aware assignment")
    #a1.errorbar(t, (1-balance_list[4]) - (1-balance_list[3]), yerr=(balance_1st[4], balance_3rd[4]), capsize=3, fmt='-o', label="Naive assignment")
    a1.axhline(y=0, color='gray', linestyle='dashdot', linewidth=0.9)
    a1.legend()
    #a0.set_ylim([0.3, 0.7])
    #a0.axhline(y=0.5, color='gray', linestyle='dashdot', linewidth=0.9)
    a1.set(ylabel='Normalized fraction of alternate allele')
    
    a2.set(xlabel='Insertion (+) or deletion (-) length')
    a2.set(ylabel='# of variants')
    a2.bar(t, list_incidents, align='center', width=0.5, log=True)
    addlabels(t, list_incidents, len_bd)
    #a1.set_yscale('log')
    a2.set_ylim([1, max(list_incidents)*5])
    a2.grid(axis='y', color='gray', linestyle='dashdot', linewidth=0.6)
    if output_name:
        plt.savefig(output_name + '.indel_balance.pdf')
    else:
        plt.show()





if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-ar', '--assignment_report', help='the assignment bias report')
    parser.add_argument('-nr', '--naive_report', help='the naive bias report')
    parser.add_argument('-vcf', '--vcf_report', help='the vcf report for the bias report regions')
    parser.add_argument('-bd', '--boundary', type=int, default=40, help='the boundary indel lengths extend from 0')
    parser.add_argument('-out', '--output_name', help="output file name")
    args = parser.parse_args()

    fn_assign_report = args.assignment_report
    fn_naive_report  = args.naive_report
    fn_vcf = args.vcf_report
    boundary = args.boundary
    output_name = args.output_name
    
    assign_SNP, assign_gap = read_bias_report(fn_assign_report)
    naive_SNP,  naive_gap  = read_bias_report(fn_naive_report)
    
    balance_SNP = calculate_SNP_balance(assign_SNP, naive_SNP)

    f_vcf = pysam.VariantFile(fn_vcf)
    balance_insert, balance_delete = calculate_gap_balance(assign_gap, naive_gap, f_vcf, boundary)
    list_incidents = [len(balance) for balance in balance_delete][::-1] + [len(balance_SNP[0])] + [len(balance) for balance in balance_insert]
    plot_balance(balance_delete, balance_SNP, balance_insert, output_name, boundary, list_incidents)


