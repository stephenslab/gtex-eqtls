import os, sys
import numpy as np
from pandas import DataFrame
from analysis_admin import env
from utils import HPData
def show_merged():
    for item in ['ENSG00000227232.4', 'ENSG00000227232.5', 'ENSG00000227232.6', 'ENSG00000227232.7',  'ENSG00000227232.8', 'ENSG00000227232.9']:
        try:
            data = HPData('merged.h5', item)
            env.log(item)
            # for table in ['beta', 't-stat', 'p-value']:
            for table in ['beta']:
                data.dump(table, True)
        except:
            env.log('%s not found' % item)
    # os.system('cat *.error')

def browse():
    env.error('Press enter to continue')
    raw_input()

env.log("Testing data conversion ....")
os.system('rm -f *.h5; python analysis_admin.py ss_to_h5 test/*.gz --action convert --output .')
for item in ['ENSG00000227232.4', 'ENSG00000227232.5', 'ENSG00000227232.6', 'ENSG00000227232.7',  'ENSG00000227232.8']:
    data = HPData('sumstat1.h5', item)
    data.dump('data', True)
os.system("zcat test/sumstat1.txt.gz | awk '{OFS = \",\"; print $2, $1, $3, $5, $4}'")
browse()
env.log('Testing merge files with unique genes in each file ...')
os.system("rm -f merged.*; python analysis_admin.py ss_to_h5 sumstat1.h5 sumstat2.h5 --action merge --output merged")
show_merged()
browse()
env.log('Testing merge files with partially overlapping SNPs in one gene ...')
os.system("rm -f merged.*; python analysis_admin.py ss_to_h5 sumstat1.h5 sumstat3.h5 --action merge --output merged")
show_merged()
browse()
env.log('Testing merge all files with various entries ...')
os.system("rm -f merged.*; python analysis_admin.py ss_to_h5 sumstat*.h5 --action merge --output merged")
show_merged()
