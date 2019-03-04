
import vcf 
import subprocess
import sys 
import os.path
import wave 
import tabix

#REQUIRED#
#pysam


#this will execute the mutationData.py file, which will build a BST from the clinvar.vcf file. Both files must be in the same directory. 
#this will take 3-10 minutes to build. 

import mutationData.py 

#example: code below tests position 86952254 on chromosome 11 for a known diesease caused by mutation. Returns two objects, a boolean and a list of features related to the disease. 
diseaseFound , diseaseInfo = testForDisease(11,86952254)
print diseaseFound
print diseaseInfo