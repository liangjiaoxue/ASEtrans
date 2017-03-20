#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This code sort the phased variants into two haplotypes , 0 first and then 1
Usage:
script vcf_for_sorting.vcf.gz  vcf_p1.vcf.gz, vcf_p2.vcf.gz, vcf_out
Author : lxue@uga.edu
"""

import re
import gzip
import sys
from sys import argv

#################
### FUNCTIONS ###
#################


class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value
# auto


def assign_hetero_zeros(vcf_in):
    sign_list = []
    for line_in in vcf_in:
        if line_in[10] == '1':
            out_string = 'A'
        elif line_in[10] == '-1':
            out_string = 'B'
        elif line_in[10] == '0':
            out_string = '0'
        else :
            sys.exit("Error in sign types\n")
        sign_list.append(out_string)
    sign_string = ''.join(sign_list)
    site_list = []
    p = re.compile("[A|B]+")
    if p.search(sign_string):
        for m in p.finditer(sign_string):
            site_list.append([int(m.start()), m.group()])
        # before the blocks
        sign_out = []
        for idx in range(0,site_list[0][0]):
            sign_out.append(site_list[0][1][0]) # fist string
        # with in the blocks
        for idx in range(len(site_list)-1):
            # letter body
            for short_str in site_list[idx][1]:
                sign_out.append(short_str)
            # zeros
            start_site = site_list[idx][0] + len(site_list[idx][1])
            end_site = site_list[idx + 1][0]
            if site_list[idx][1][-1] == site_list[idx+1][1][0]: # same at two sides
                for idx_zero in range(start_site, end_site):
                    sign_out.append(site_list[idx][1][-1]) # last character of
            else : # different sign at two sides
                zero_buffer = []
                for idx_zero in range(start_site, end_site):
                    zero_buffer.append(vcf_in[idx_zero])
                buff_out = check_buffer(zero_buffer, vcf_in[start_site-1], vcf_in[end_site])
                for line_out in buff_out:
                    sign_out.append(line_out)
        # after the blocks
        for short_str in site_list[-1][1]:
            sign_out.append(short_str)
        start_site = site_list[-1][0] + len(site_list[-1][1])
        for idx in range(start_site, len(vcf_in)):
            sign_out.append(site_list[-1][1][-1])  # last string
        # merge vcf with sign information
        print ("Start_site:"+vcf_in[0][1]+"\t vcf length:"+str(len(vcf_in))+"\t "+"output length:"+str(len(sign_out)))
        for idx in range(len(vcf_in)):
            vcf_in[idx].append(sign_out[idx])
    else :
        for idx in range(len(vcf_in)):
            vcf_in[idx].append('all_zeros')
    return vcf_in


def check_buffer(zero_buffer, start_line, end_line):
    start_site = int(start_line[1])
    end_site = int(end_line[1])
    start_sign = 'A' if start_line[10] == '1' else 'B'
    end_sign = 'A' if end_line[10]=='1' else 'B'
    out_list = []
    gaps = []
    for line_list in zero_buffer:
        gaps.append(int(line_list[1])-start_site)
        start_site = int(line_list[1])
    gaps.append(end_site-start_site)
    max_gap = max(gaps)
    sign_out = start_sign
    for idx in range(len(gaps)-1):
        if gaps[idx] != max_gap:
            out_list.append(sign_out)
        else:
            sign_out = end_sign
            out_list.append(sign_out)
    return out_list


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
    # use global VarTable1
    out = 0
    if len(geno_in) == len_ref_in :
        relative_site = 0
        for nt in geno_in:
            site_out = int(site_in)+relative_site
            check_site = str(site_out)+nt
            if VarTable1.vcf_ma_pa[tag_in][check_site] == 1:
                out += 1
            relative_site += 1
    return out


def check_sliding_windows(record_signed):
    window_size = 10
    record_signed_2 = []
    for idx in range(len(record_signed)):
        # Check the windows for sign
        start_site = idx - window_size
        if start_site < 0:
            start_site = 0
        end_site = idx + window_size + 1
        if end_site > len(record_signed):
            end_site = len(record_signed)
        sign_sum = 0
        for idx_sum in range(start_site, end_site):
            sign_sum += int(record_signed[idx_sum][10])
        record_signed[idx].append(str(sign_sum))
        record_signed_2.append(record_signed[idx])
    return record_signed_2


def check_boundary(record_signed_2):
    ## check boundary
    out_list = []
    for idx in range(len(record_signed_2)):
        sign_check10 = int(record_signed_2[idx][10])
        sign_check11 = int(record_signed_2[idx][11])
        if sign_check11 > 1:
            switch = '1'
        elif sign_check11 < -1:
            switch = '-1'
        else:  # boundary now
            if sign_check10 == sign_check11:
                switch = str(sign_check10)
            else:  # check genome distance
                if (idx - 2) < 0 or (idx + 2) >= len(record_signed_2):
                    switch = str(sign_check10)
                else:
                    dis_before = abs(int(record_signed_2[idx - 2][1]) - int(record_signed_2[idx][1]))
                    dis_after = abs(int(record_signed_2[idx + 2][1]) - int(record_signed_2[idx][1]))
                    if dis_before < dis_after:
                        switch_boundary = int(record_signed_2[idx - 2][11])
                    else:
                        switch_boundary = int(record_signed_2[idx + 2][11])
                    if switch_boundary < 0:
                        switch = '-1'
                    else:
                        switch = '1'
        record_signed_2[idx].append(switch)
        out_list.append(record_signed_2[idx])
    return out_list


class VarTable:
    """Hold the VCF records in memory"""
    def __init__(self):
        self.vcf_hetero = AutoVivification()
        self.vcf_ma_pa = AutoVivification()
        self.vcf_final = []
        self.vcf_hetero_single = []
        self.vcf_hetero_phased = []
        # for output
        self.head = []
        self.summary = 0
        self.total_num = 0
        self.homo_num = 0

    def add_head_line(self,line_in):
        self.head.append(line_in)

    def merge_sort_vcf(self):
        # put hetero_single variants into final
        for line_in in self.vcf_hetero_single :
            out_record = line_in[0:10] # 11 for test
            if line_in[10] == '-1' :
                data_split = out_record[9].split(':')
                gt = data_split[0]
                data_split[0] = gt[2] + '/' + gt[0]
                out_record[9] = ':'.join(data_split)
            self.vcf_final.append(out_record)
        # put hetero phased variants into finsl
        for line_in in self.vcf_hetero_phased :
            out_record = line_in[0:10]  # 12 for test
            if line_in[11] == 'B':
                data_split = out_record[9].split(':')
                gt = data_split[0]
                data_split[0] = gt[2] + '|' + gt[0]
                out_record[9] = ':'.join(data_split)
            self.vcf_final.append(out_record)
        # sort the table
        new_list = sorted(self.vcf_final, key=lambda x: int(x[1]))
        if self.vcf_final[0][0].startswith('Chr'):
            new_list = sorted(new_list, key=lambda x: int(x[0][3:]))
        else :
            new_list = sorted(new_list, key=lambda x: int(x[0][9:]))
        self.vcf_final = new_list

    def vcf_write(self,file_out):
        with open(file_out,"w") as VCFOUT:
            VCFOUT.write("\n".join(self.head)+"\n")
            for line_new in self.vcf_final:
                VCFOUT.write("\t".join(line_new) + "\n")

    def log_write(self,file_out):
        pass
        #log_out = file_out[0:-4] + '.log'
        #with open(log_out,"w") as LOG:
            #LOG.write("Total var No.: {!s}".format(str(self.total_num)) + "\n")
            #LOG.write("Total homo var No.: {!s}".format(str(self.homo_num)) + "\n")
            #description = ['he_single', 'he_phase', 'he_phase_each', 'he_single_convert', 'he_phase_convert',
            #               'he_phase_each_convert']
            #for idx in range(len(self.summary)):
            #    LOG.write(description[idx] + ':' + str(self.summary[idx]) + "\n")


    def read_vcf_file(self,tag_in, file_in, chr_in):
        check_begin = re.compile("^#")
        check_data = re.compile("\d")
        with gzip.open(file_in, "rt") as INVCF:
            for line_vcf in INVCF:
                if not check_begin.match(line_vcf):
                    line_list = line_vcf.rstrip().split("\t")
                    chr_id = line_list[0]
                    if chr_id[0:5] == chr_in:
                        sample_in = line_list[9]
                        site = line_list[1]
                        allele = [line_list[3]]
                        allele.extend(line_list[4].split(','))
                        format_records = line_list[8].split(':')
                        gt_pos = format_records.index('GT')
                        mapa_gt = sample_in.split(':')[gt_pos]
                        if check_data.search(mapa_gt):
                            self.input_dict_function(tag_in, site, allele, mapa_gt)

    def input_dict_function(self, tag, site_in, allele_in, input_gt):
        # allele is list, gt is string
        gt_1 = input_gt[0]
        gt_2 = input_gt[2]
        gt_list = [allele_in[int(gt_1)]]
        if gt_1 != gt_2:
            gt_list.extend(allele_in[int(gt_2)])
        ref_string = allele_in[0]
        for allele_string in gt_list:
            if len(ref_string) == len(allele_string):
                relative_site = 0
                for nt in allele_string:
                    site_new = int(site_in) + relative_site
                    site_out = str(site_new) + nt
                    self.vcf_ma_pa[tag][site_out] = 1
                    relative_site += 1


    def read_main_vcf(self,file_in):
        check_begin = re.compile("^#")
        check_data = re.compile("\d")
        with gzip.open(file_in, "rt") as INPUT:
            add_tag = 0
            for line in INPUT:
                if check_begin.match(line):
                    if line.startswith('##INFO=') and add_tag == 0:
                        self.add_head_line(
                            '##FORMAT=<ID=PS,Number=.,Type=String,Description="Phase set identifier">')
                        self.add_head_line(
                            '##FORMAT=<ID=HN,Number=1,Type=Integer,Description="Haplotype not-assigned">')
                        add_tag = 1
                    self.add_head_line(line.rstrip())
                else:
                    line_list = line.rstrip().split("\t")
                    self.total_num += 1
                    check_sample = line_list[9]
                    # chr_id = line_list[0]
                    # site = line_list[1]
                    allele = [line_list[3]]
                    allele.extend(line_list[4].split(','))
                    format_records = line_list[8].split(':')
                    gt_pos = format_records.index('GT')
                    # put homo in list, hetero in dictionary
                    if check_data.search(check_sample):
                        check_records = check_sample.split(':')
                        check_gt = check_records[gt_pos]
                        check_ro = check_records[format_records.index('RO')]
                        check_ao = check_records[format_records.index('AO')]
                        out1 = line_list
                        out1[7] = '.'
                        out1[8] = 'GT:RO:AO'
                        out1[9] = check_gt + ':' + check_ro + ':' + check_ao
                        if check_gt[0] == check_gt[2]:
                            self.vcf_final.append(out1) # homozygous directly into final
                        else:
                            if 'PS' in format_records:
                                ps_pos = format_records.index('PS')
                                ps_data = check_records[ps_pos]
                                self.vcf_hetero['P' + ps_data]["\t".join(out1)] = 1
                            else:
                                self.vcf_hetero['S']["\t".join(out1)] = 1

    def check_main_vcf(self):
        for key1 in self.vcf_hetero.keys():
            if key1 == 'S':  # scaffold
                for vcf_hetero_line_string in self.vcf_hetero[key1]:
                    vcf_hetero_line = vcf_hetero_line_string.split("\t")
                    gt_1_ma_num, gt_1_pa_num, gt_2_ma_num, gt_2_pa_num = get_four_hit(vcf_hetero_line)
                    mp_score = gt_1_ma_num - gt_1_pa_num - (gt_2_ma_num - gt_2_pa_num)
                    score_out = '1'
                    if mp_score < 0:
                        score_out = '-1'
                    elif mp_score == 0:
                        score_out = '0'
                    vcf_hetero_line.append(score_out)
                    self.vcf_hetero_single.append(vcf_hetero_line)
            else:  # P = phased
                # check order of phased SNP
                record_signed = []
                each_phase = []
                for vcf_phased_line_string in self.vcf_hetero[key1]:
                    vcf_phased_line = vcf_phased_line_string.split("\t")
                    gt_1_ma_num, gt_1_pa_num, gt_2_ma_num, gt_2_pa_num = get_four_hit(vcf_phased_line)
                    mp_score = gt_1_ma_num - gt_1_pa_num - (gt_2_ma_num - gt_2_pa_num)
                    score_out = '1'
                    if mp_score < 0:
                        score_out = '-1'
                    elif mp_score == 0:
                        score_out = '0'
                    vcf_phased_line.append(score_out)
                    vcf_phased_line[9] = vcf_phased_line[9] + ':' + key1[1:]
                    vcf_phased_line[8] += ':PS'
                    if score_out == '0':
                        each_phase.append(vcf_phased_line)
                    else :
                        record_signed.append(vcf_phased_line)
                ##############
                # sort the records
                record_signed = sorted(record_signed, key=lambda x: int(x[1]))
                # check signed with sliding windows
                record_signed_2 = check_sliding_windows(record_signed)
                # check boundary
                record_signed_3 = check_boundary(record_signed_2)
                for line_3_list in record_signed_3:
                    line_3_list_out = line_3_list[0:10]
                    line_3_list_out.append(line_3_list[-1])
                    each_phase.append(line_3_list_out)
                # check zeros in each phase
                each_phase = sorted(each_phase,key=lambda x: int(x[1]))
                phased_out = assign_hetero_zeros(each_phase)
                for line_line in phased_out:
                    self.vcf_hetero_phased.append(line_line)
# End of Class

#################
###   MAIN    ###
#################

if __name__ == "__main__":
    script, input_file, parent1, parent2, chrome_tag, output_file = argv
    if chrome_tag.startswith('scaff'):
        chrome_tag = 'scaff'
    VarTable1 = VarTable()
    # read guides
    VarTable1.read_vcf_file('m', parent1, chrome_tag)
    print("Read m done")
    VarTable1.read_vcf_file('p', parent2, chrome_tag)
    print("Read p done")
    # read main
    VarTable1.read_main_vcf(input_file)
    print("Read main done")
    VarTable1.check_main_vcf()
    VarTable1.merge_sort_vcf()
    VarTable1.vcf_write(output_file)
    #VarTable1.log_write(output_file)
