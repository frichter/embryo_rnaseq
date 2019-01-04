"""Prepare the chromosome matrix.

cd /Users/felixrichter/Dropbox/PhD/embryo_rnaseq/reproseq
python3
"""

repro_f_loc = 'reproseq_data.txt'
repro_chr_loc = 'reproseq_chr_df.txt'
chr_list = ['chr' + str(i) for i in range(1, 23)]
chr_list.extend(['chrX', 'chrY'])
chr_euploid_ct = [2 for i in range(1, 25)]


def clean_mono_list(mono_list, line_dict):
    """Clean the monosomy list."""
    # don't have raw data for 93115_C4_THS_021_BxE8_3_9_17
    if line_dict['Specimen ID'] == '93115_C4_THS_021_BxE8_3_9_17':
        mono_list.remove('chr9')
    if line_dict['Specimen ID'] == '93425_C2_THS_020_BxE2_3_9_17':
        mono_list.remove('chr9')
    if line_dict['Specimen ID'] == '94279_C2_THS_002_BxE1_1/27/17':
        mono_list.remove('chr22')
    return mono_list


def clean_tri_list(tri_list, line_dict):
    """Clean the monosomy list."""
    # don't have raw data for 93115_C4_THS_021_BxE8_3_9_17
    if line_dict['Specimen ID'] == '83738_C7_THS_026_BxE5_3_15_17':
        tri_list.remove('chr11')
    if line_dict['Specimen ID'] == '91122_C1_THS_027_BxE2_3_17_17':
        tri_list.remove('chr3')
    # don't have raw data for 93115_C4_THS_021_BxE8_3_9_17
    if line_dict['Specimen ID'] == '93115_C4_THS_021_BxE8_3_9_17':
        tri_list.remove('chr9')
    if line_dict['Specimen ID'] == '93425_C2_THS_020_BxE2_3_9_17':
        tri_list.remove('chr9')
    if line_dict['Specimen ID'] == '94279_C2_THS_002_BxE1_1/27/17':
        tri_list.remove('chr22')
    return tri_list


count = 0
with open(repro_f_loc, 'r') as repro_f, open(repro_chr_loc, 'w') as out_f:
    header = next(repro_f).strip().split('\t')
    out_line = 'Specimen.ID\t' + '\t'.join(chr_list) + '\n'
    _ = out_f.write(out_line)
    for line in repro_f:
        # reset chr_dict for every patient
        chr_dict = dict(zip(chr_list, chr_euploid_ct))
        line_dict = dict(zip(header, line.strip().split('\t')))
        if line_dict['passed QC?'] == '0':
            continue
        if line_dict['Monosomy'] != '':
            mono_list = ['chr' + i for i in line_dict['Monosomy'].split(',')]
            # print(mono_list)
            mono_list = clean_mono_list(mono_list, line_dict)
            for chr_i in mono_list:
                chr_dict[chr_i] = 1
        if line_dict['Trisomy'] != '':
            tri_list = ['chr' + i for i in line_dict['Trisomy'].split(',')]
            # print(tri_list)
            tri_list = clean_tri_list(tri_list, line_dict)
            for chr_i in tri_list:
                if chr_dict[chr_i] == 1:
                    print(line_dict)
                    print(chr_i)
                    raise ValueError('Cannot have both a monosomy and trisomy')
                chr_dict[chr_i] = 3
        # print(line_dict['Specimen ID'])
        # print(chr_dict)
        chr_ct_str = [str(i) for i in chr_dict.values()]
        out_line = '{}\t{}\n'.format(
            line_dict['Specimen ID'], '\t'.join(chr_ct_str))
        print(out_line)
        _ = out_f.write(out_line)
        count += 1
        # if count > 5:
        #    break

#
#
#
