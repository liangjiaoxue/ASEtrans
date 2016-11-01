#!/usr/bin/python
# -*- coding: utf-8 -*-

"""This code filter the output from bcftools to get one single sample records"""

import re
from sys import argv
script, input_file, output_file = argv


# dir = os.path.dirname(os.path.realpath(__file__))
# data_file_full = os.path.join(dir,data_file)

OUT = open(output_file, "w")
check_begin = re.compile("^#")
check_data = re.compile("\d")

num_converted_var = 0
num_converted_phase = 0


class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value
# auto


def write_vcf():
    global vcf_hetero
    global vcf_list
    global num_converted_var
    global num_converted_phase
    for key1 in vcf_hetero.keys():
        if key1 == 'S': # scaffold
            for vcf_hetero_line_string in vcf_hetero[key1]:
                vcf_hetero_line = vcf_hetero_line_string.split("\t")
                gt_1_ma_num, gt_1_pa_num, gt_2_ma_num, gt_2_pa_num = get_four_hit(vcf_hetero_line)
                mp_score = gt_1_ma_num - gt_1_pa_num - (gt_2_ma_num - gt_2_pa_num )
                # check score
                check_split = vcf_hetero_line[9].split(':')
                gt_in = check_split[0]
                gt_out = gt_in
                if mp_score < 0:
                    gt_out = gt_in[2] + '/' + gt_in[0]
                elif mp_score == 0:
                    pair_in = [gt_in[0], gt_in[2]]
                    if '0' in pair_in:
                        pos_0 = pair_in.index('0')
                        gt_out = pair_in[1-pos_0] + '/' + pair_in[pos_0]
                    else :
                        # no zero convert to homo
                        gt_out = pair_in[0] + '/' + pair_in[0]
                        num_converted_var += 1
                # score larger than 0 no change
                check_split[0]= gt_out
                vcf_hetero_line[9]= ':'.join(check_split)
                vcf_list.append(vcf_hetero_line)
        else: # P = phased
            # check order of phased SNP
            gt_1_ma_num_phased = 0
            gt_1_pa_num_phased = 0
            gt_2_ma_num_phased = 0
            gt_2_pa_num_phased = 0
            for vcf_phased_line_string in vcf_hetero[key1]:
                vcf_phased_line = vcf_phased_line_string.split("\t")
                gt_1_ma_num, gt_1_pa_num, gt_2_ma_num, gt_2_pa_num = get_four_hit(vcf_phased_line)
                gt_1_ma_num_phased += gt_1_ma_num
                gt_1_pa_num_phased += gt_1_pa_num
                gt_2_ma_num_phased += gt_2_ma_num
                gt_2_pa_num_phased += gt_2_pa_num
            # total score
            mp_score_phased = gt_1_ma_num_phased - gt_1_pa_num_phased - (gt_2_ma_num_phased - gt_2_pa_num_phased)
            change = 0
            if mp_score_phased < 0 :
                change = 1
            elif mp_score_phased == 0 :
                m_zero_num = 0
                p_zero_num = 0
                for vcf_phased_line_string in vcf_hetero[key1]:
                    vcf_phased_line = vcf_phased_line_string.split("\t")
                    if vcf_phased_line[9][0] == '0' : # not good but work
                        m_zero_num += 1
                    if vcf_phased_line[9][2] == '0' :
                        p_zero_num += 1
                if m_zero_num > p_zero_num :
                    change = 1
                elif m_zero_num < p_zero_num :
                    change = 0
                elif  m_zero_num == p_zero_num:
                    # convert to homo
                    change = 2
                else:
                    sys.exit("Phased Error\n")
            if change == 1 :
                for vcf_phased_line_string in vcf_hetero[key1]:
                    vcf_phased_line = vcf_phased_line_string.split("\t")
                    data_split = vcf_phased_line[9].split(':')
                    gt = data_split[0]
                    data_split[0] = gt[2]+'|' + gt[0]
                    vcf_phased_line[9] = ':'.join(data_split)+':'+key1[1:]
                    vcf_phased_line[8] += ':PS'
                    vcf_list.append(vcf_phased_line) # append with change
            elif change == 0: # no change
                for vcf_phased_line_string in vcf_hetero[key1]:
                    vcf_phased_line = vcf_phased_line_string.split("\t")
                    vcf_phased_line[9] = vcf_phased_line[9] + ':' + key1[1:]
                    vcf_phased_line[8] += ':PS'
                    vcf_list.append(vcf_phased_line)
            elif change == 2 :
                num_converted_phase += 1
                for vcf_phased_line_string in vcf_hetero[key1]:
                    vcf_phased_line = vcf_phased_line_string.split("\t")
                    data_split = vcf_phased_line[9].split(':')
                    gt = data_split[0]
                    data_split[0] = gt[0]+'/' + gt[0]
                    vcf_phased_line[9] = ':'.join(data_split)+':'+key1[1:]
                    vcf_phased_line[8] += ':PS'
                    vcf_list.append(vcf_phased_line)                    
            else:
                sys.exit("Error in phase sorting\n")
    ##########################################################
    # data loaded in list
    new_list = sorted(vcf_list, key=lambda x: int(x[1]))
    for line_new in new_list:
        OUT.write("\t".join(line_new)+"\n")
# End od write VCF


def get_four_hit(vcf_line_in):
    gt_in = vcf_line_in[9].split(':')[0]
    len_ref = len(vcf_line_in[3])
    allele_in = [vcf_line_in[3]]
    allele_in.extend(vcf_line_in[4].split(','))
    # print("\t".join(vcf_line_in))
    # print ("\t".join(allele_in) + "allele information" + "\n")
    gt_1 = allele_in[int(gt_in[0])]
    gt_2 = allele_in[int(gt_in[2])]
    gt_1_ma_num = get_snp_hit(gt_1, "m",vcf_line_in[1],len_ref)
    gt_1_pa_num = get_snp_hit(gt_1, "p",vcf_line_in[1],len_ref)
    gt_2_ma_num = get_snp_hit(gt_2, "m",vcf_line_in[1],len_ref)
    gt_2_pa_num = get_snp_hit(gt_2, "p",vcf_line_in[1],len_ref)
    return gt_1_ma_num,gt_1_pa_num,  gt_2_ma_num,  gt_2_pa_num
#######


def get_snp_hit(geno_in, tag_in, site_in, len_ref_in):
    global vcf_ma_pa
    out = 0
    if len(geno_in) == len_ref_in :
        relative_site = 0
        for nt in geno_in:
            site_out = int(site_in)+relative_site
            check_site = str(site_out)+nt
            if vcf_ma_pa[tag_in][check_site] == 1:
                out += 1
            relative_site += 1
    return out
###########################################################


def input_dict(tag, site_in, allele_in, input_gt):
    global vcf_ma_pa
    # allele is list, gt is string
    gt_1 = input_gt[0]
    gt_2 = input_gt[2]
    gt_list = [allele_in[int(gt_1)]]
    if gt_1 != gt_2 :
        gt_list.extend(allele_in[int(gt_2)])
    ref_string = allele_in[0]
    for allele_string in gt_list:
        if len(ref_string)==len(allele_string):
            relative_site = 0
            for nt in allele_string:
                site_new = int(site_in) + relative_site
                site_out = str(site_new)+nt
                vcf_ma_pa[tag][site_out] = 1
                relative_site += 1
    return 1
#################################################

vcf_ma_pa  = AutoVivification()
vcf_list   = []
vcf_hetero = AutoVivification()
tracking   = 0
chr_last   = ''


with open (input_file,"r") as INPUT:
    for line in INPUT:
        if check_begin.match(line):            
            if '#CHROM' == line[0:6]:
                line_list = line.rstrip().split("\t")
                OUT.write("\t".join(line_list[0:10])+"\n")
            else:
                OUT.write(line)
        else:
            line_list = line.rstrip().split("\t")
            check_sample, ma_sample, pa_sample = line_list[9:]
            chr_id = line_list[0]
            site   = line_list[1]
            allele = [line_list[3]]
            allele.extend(line_list[4].split(','))
            format_records = line_list[8].split(':')
            gt_pos = format_records.index('GT')
            ma_gt = ma_sample.split(':')[gt_pos]
            pa_gt = pa_sample.split(':')[gt_pos]
            if chr_id != chr_last:
                print ("Starting:"+chr_id+"\n")
                # new chromosome, check last large record
                if chr_last != '':
                    # not firt line
                    write_vcf() 
                    vcf_ma_pa  = AutoVivification()
                    vcf_list   = []
                    vcf_hetero = AutoVivification()
                ## refresh old ID
                chr_last = chr_id
            # the records need to be load regardless of new or old chromosome ID
            ###########################################
            # put homo in list, hetero in dictionary
            if check_data.search(check_sample):
                check_records = check_sample.split(':')
                check_gt = check_records[gt_pos]
                check_RO = check_records[format_records.index('RO')]
                check_AO = check_records[format_records.index('AO')]
                out1 = line_list[0:10]
                out1[7]='.'
                out1[8] = 'GT:RO:AO'
                out1[9] = check_gt+':'+check_RO+':'+check_AO
                if check_gt[0] == check_gt[2]:
                    vcf_list.append(out1)
                else:
                    if 'PS' in format_records :
                        ps_pos  = format_records.index('PS')
                        ps_data = check_records[ps_pos]
                        vcf_hetero['P' + ps_data]["\t".join(out1)] = 1
                    else :
                        vcf_hetero['S']["\t".join(out1)] = 1
            ############################################
            if check_data.search(ma_gt) :
                input_dict("m",site,allele, ma_gt)
            if check_data.search(pa_gt) :
                input_dict("p",site,allele, pa_gt)
## have last record in memory
write_vcf() 
OUT.close()
print (str(num_converted_var)+" variant are converted \n")
print (str(num_converted_phase)+" phased group are converted\n")
## Author : lxue@uga.edu
