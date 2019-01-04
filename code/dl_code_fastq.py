"""Download embryo data.

module load python/3.5.0 py_packages/3.5
source ~/venv_ore/bin/activate

cd /sc/orga/projects/chdiTrios/Felix/embryo_rnaseq

python

"""

import re
import os
# from urllib.request import urlretrieve
import requests
from bs4 import BeautifulSoup
import shutil


def download_fastq(f_url, f_local):
    """Download the large fastq files."""
    # response = requests.get(f_url, verify=False, auth=('bobby', 'b0bBy'))
    # urlretrieve(f_url, f_local)
    # with open('filename.txt','w') as fout:
    #    fout.write(response.text)
    # Source: https://stackoverflow.com/a/39217788/10688049
    with requests.get(f_url, stream=True, auth=('bobby', 'b0bBy')) as r:
        with open(f_local, 'wb') as f:
            shutil.copyfileobj(r.raw, f)
        return f_local


home_dir = '/hpc/users/richtf01/embryo_rnaseq/'
# home_dir = '/Users/felixrichter/Dropbox/PhD/embryo_rnaseq/'
html_f = home_dir + 'urls/fastq.html'

dl_dir = '/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/FASTQ/'

# import html document as a BS4 object
with open(html_f, 'r') as html_doc:
    soup = BeautifulSoup(html_doc, 'html.parser')

soup.contents.__len__()  # list
# soup.children # iterator
list(soup.descendants).__len__()  # list
bed_regex = re.compile('fastq.gz$')
f_url_list = []
for linkSoup in soup.find_all('a'):
    f_url = linkSoup['href']
    # print(f_url, f_local)
    if re.search(bed_regex, f_url):
        f_url_list.append(f_url)


len(f_url_list)
count = 0
dl_ct = 0
ct_bin = 0  # 0 75 125
for f_url in f_url_list:
    f_local = dl_dir + re.sub('.*FASTQ/', '', f_url)
    count += 1
    if (count % 10) == 0:
        print(count)
    if (count < ct_bin) or (count > ct_bin + 500):
        continue
    if not os.path.exists(f_local):
        print('Starting', f_url)
        download_fastq(f_url, f_local)
        print(f_local, "downloaded!")
        dl_ct += 1
    else:
        print('Already downloaded', f_local)


# 75888_C4_THS_014_BxE8_2_28_17_S18_L004_R1_001.fastq.gz
# first file should be 3.4 Gb
""" Benchmarking against wget
cd /sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/FASTQ
time wget https://wangy33.u.hpc.mssm.edu/pacbioDataRelease/PBG-X19/\
qGrpzvc9wFQP_pSf-wQePUZ8rRBw33Z7hVC8wkR7xmk/FASTQ/96295_C1_THS_029\
_E17_3_15_17_S41_L008_R2_001.fastq.gz \
/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/FASTQ/\
96295_C1_THS_029_E17_3_15_17_S41_L008_R2_001.fastq.gz


91122_C1_THS_027_BxE3_3_17_17_S37_L008_R2_001.fastq.gz

Comparing downloaded to observed
cd /sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/FASTQ
ls -alh * | tr -s ' ' | cut -d' ' -f5,9

# locally
cd /Users/felixrichter/Dropbox/PhD/embryo_rnaseq/urls
python3

online_dict = {}
with open('fastq_online_data_size.txt', 'r') as f:
    for line in f:
        line_list = line.strip().split('\t')
        online_dict[line_list[0]] = line_list[1]

dl_dict = {}
with open('fastq_dl_data_size.txt', 'r') as f:
    for line in f:
        line_list = line.strip().split(' ')
        dl_dict[line_list[1]] = line_list[0]


count = 0
for k, v in online_dict.items():
    if k in dl_dict:
        # if v != dl_dict[k]:
        onl_size = float(v[:-1])
        dl_size = float(dl_dict[k][:-1])
        # if k == '92182_C5_THS_011_E9_2_23_17_S11_L003_R2_001.fastq.gz':
        #     print(k, v, onl_size, dl_dict[k], dl_size)
        # account for the fact that downloaded files appear to have
        # larger sizes than the online version
        if onl_size < (dl_size - 0.11):
            print(k, v, dl_dict[k])
            # count += 1
            continue
        # check for files where downloaded is smaller than online
        # should only be currently downloading files
        if onl_size > dl_size:
            print(k, v, onl_size, dl_dict[k], dl_size)
    else:
        print(k, 'not downloaded')


len(dl_dict)
len(online_dict)
"""
