import os
from utils import HPData
os.system('python analysis_admin.py ss_to_h5 test/Merged.h5 --action max --output Merged.Max.h5')
data = HPData('test/Merged.h5', 'ENSG00000238009.2').dump('t-stat').loc['1_895755_A_AG_b37']
print data
data = HPData('Merged.Max.h5', 'max').dump('t-stat')
print data
