#!/usr/bin/env python
#_*_coding:utf-8_*_
'''determine whether the version of user's python comply with the requirements of  this procedure'''
import sys
if sys.version_info[0] != 3 or sys.version_info[1] != 7:
	print("\nYou are using python" + str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + " , while pHisPred needs python3.7!\n", file =  sys.stderr)
	sys.exit()
import os
sys.path.append(os.getcwd()+'\\codes')
import re
import math
import optparse
import time
import argparse
import re
import numpy as np
from codes import *
import numpy as np
from sklearn import svm
from sklearn import preprocessing

#==============================================================================
# transfer multiple-line fasta format into two-line fasta format
#==============================================================================

def TwoLineFasta (Seq_Array):
    Tmp_sequence_Arr = []
    Tmp_trans_str = ''
    for i in range(len(Seq_Array)):
        if '>' in Seq_Array[i]:
            if i == 0:
                Tmp_sequence_Arr.append(Seq_Array[i])
            else:
                Tmp_sequence_Arr.append(Tmp_trans_str)
                Tmp_sequence_Arr.append(Seq_Array[i])
                Tmp_trans_str = ''
        else:
            if i == len(Seq_Array) - 1:
                Tmp_trans_str = Tmp_trans_str + str(Seq_Array[i])
                Tmp_sequence_Arr.append(Tmp_trans_str)
            else:
                Tmp_trans_str = Tmp_trans_str + str(Seq_Array[i])
    return Tmp_sequence_Arr

#==============================================================================
# construct id array and sequence array
#==============================================================================
   
def Tran_Seq (input_arr):
    label_Arr = []
    FastA_seq_Arr = []
    for n in range(len(input_arr)):
        if n == 0 or n % 2 == 0:
            label = input_arr[n]
            label_Arr.append(label)
        else :
            seq = input_arr[n]
            FastA_seq_Arr.append(seq)
    return (label_Arr,FastA_seq_Arr)

#==============================================================================
# extract peptides surrouding histidine sites
#==============================================================================
   
def Extract_Peptides (input_arr, wz, outfile):
    fout = open(outfile,'w')
    for n in range(len(input_arr)):
        if n == 0 or n % 2 == 0:
            label = input_arr[n]
        else :
            seq = input_arr[n]
            for i in re.finditer('H',seq):
                if i.start() < int(wz/2):
                    segment = seq[0: i.start() + int(wz/2) + 1]
                    m = int(wz/2) - i.start()
                    segment = 'X' * m + segment 
                elif i.start() + int(wz/2) + 1 > len(seq):
                    segment = seq[i.start() - int(wz/2) : len(seq)] 
                    m = int(wz/2) - len(seq) + i.start() + 1
                    segment = segment + 'X' * m
                else:
                    segment = seq[i.start() - int(wz/2) : i.start() + int(wz/2) + 1]
                print((label + "#" + str(i.start()+1)), file = fout)
                print(segment, file = fout)
    fout.close()
    return (outfile)



#==============================================================================
#   Print the final result 
#==============================================================================
def PrintResult(ids,probability,labels,outputfile):
      
    for i in range(len(ids)):
        line = ids[i]
        line += '\t' + str(probability[i])
        if labels[i] == 1:
            line = line + '\t' + 'pHis_site'
        else:
            line = line + '\t' + 'non_pHis_site'
        line += '\n'
        outputfile.write(line)
        
#==============================================================================
# Main program
#==============================================================================

parse=optparse.OptionParser()
parse.add_option('-f','--file',dest='file',action='store',metavar='input files',help="enter proteins in .fasta (.fa) format.")
parse.add_option('-o','--out',dest='outfile',action='store',metavar='output files',help='assign your output file.')
parse.add_option('-t','--type',dest='type',action="store",metavar='class type',default="euka",help="eukaryotes: euka (default); prokaryotes: proka.")
parse.add_option('-w','--window',dest='wz',action="store",metavar='window size', help="specific the window size used for extracting the peptides around histidine sites")
parse.add_option('-m','--model',dest='model',action="store",metavar='trained model', help="specific your own trained model")

(options,args) = parse.parse_args()

#check input and output files
for file in ([options.file,options.outfile,options.type]):
	if not (file):
		parse.print_help()
		sys.exit(0)

inPutFileName = options.file
outPutFileName = options.outfile
ModelType = options.type

Compute_time = time.time()

#==============================================================================
#  data preprocessing
#==============================================================================

inFiles = open(inPutFileName)
inFilesArr = inFiles.read()
inFiles.close()

sequence_Arr = inFilesArr.split('\n')
del inFilesArr
sLen = len(sequence_Arr) - 1
del sequence_Arr[sLen]

flg = 1
mp = int(len(sequence_Arr[1])/2)
if sequence_Arr[1][mp] == 'H':
    flg = 0

if len(sequence_Arr[1]) > 50:
    flg = 1

wz = 31
if ModelType == 'proka':
    wz = 25

if (options.wz):
    wz = int(options.wz)


(filepath, file) = os.path.split(outPutFileName)
(filename, extension) = os.path.splitext(file) 

if flg:
    tmpfile = 'tmp_' + filename + '.fa' 
    tmpfile = filepath + tmpfile
    ARRAY_temp = TwoLineFasta(sequence_Arr)
    Extract_Peptides(ARRAY_temp, wz, tmpfile)
    del ARRAY_temp
else:
    tmpfile = inPutFileName

#==============================================================================
#  feature calcualtion
#==============================================================================

euka_triplets = {'AAA' : ['AER', 'AGV', 'AHF', 'DLH', 'DMP', 'DPV', 'DWP',
                          'EGA', 'EMM', 'FCI', 'FDR', 'FIS', 'FKL', 'FRQ', 
                          'GDW', 'GEM', 'GGA', 'GKD', 'HAW', 'HGD', 'HGE', 
                          'HKI', 'HKY', 'HMN', 'HTF', 'IKN', 'KAT', 'KCP', 
                          'KDI', 'KLP', 'LFI', 'MGQ', 'MNG', 'MPN', 'NAL', 
                          'NVK', 'PTK', 'PYL', 'QGI', 'RHG', 'RYP', 'SGK', 
                          'SMF', 'VCH', 'VFT', 'VNA', 'VYC', 'WGK', 'WPL', 
                          'YGS']}
proka_triplets = {'AAA' : ['FRR', 'GRP', 'HDH', 'HHT', 'HTK', 'RTP', 'SHE',
                           'VVR', 'YHH']}

euka_aaindex = {'POS' : [2, 4, 14, 16, 17, 24, 28]}
proka_aaindex = {'POS' : [5, 10, 11, 12, 13, 14, 17, 19]}

euka_features = ['AAINDEX', 'APAAC', 'BINARY', 'BLOSUM62', 'CKSAAGP', 'CKSAAP', 
                 'CTDC','DPC', 'EAAC', 'EGAAC', 'GAAC', 'GDPC', 'Geary', 'Moran', 
                 'PAAC', 'QSOrder', 'TPC']
euka_features_parameter = {'AAINDEX': '(fastas, **euka_aaindex)', 
                           'APAAC': '(fastas, lambdaValue = 15)', 
                           'BINARY': '(fastas)', 'BLOSUM62': '(fastas)', 
                           'CKSAAGP': '(fastas)', 'CKSAAP': '(fastas)', 
                           'CTDC': '(fastas)','DPC': '(fastas)', 
                           'EAAC': '(fastas)', 'EGAAC': '(fastas)', 
                           'GAAC': '(fastas)', 'GDPC': '(fastas)', 
                           'Geary': '(fastas, nlag = 15)', 
                           'Moran': '(fastas, nlag = 15)', 
                           'PAAC': '(fastas, lambdaValue = 15)', 
                           'QSOrder': '(fastas, nlag = 15)', 
                           'TPC': '(fastas, **euka_triplets)'}

proka_features = ['AAINDEX', 'APAAC', 'BINARY', 'BLOSUM62', 'CKSAAGP', 
                  'CKSAAP', 'CTDC', 'CTDD', 'CTDT', 'CTriad', 'DPC', 'EAAC', 
                  'EGAAC', 'GAAC', 'GDPC', 'GTPC', 'KSCTriad', 'PAAC', 
                  'QSOrder', 'SOCNumber', 'TPC']
proka_features_parameter = {'AAINDEX': '(fastas, **proka_aaindex)', 
                            'APAAC': '(fastas, lambdaValue = 10)', 
                            'BINARY': '(fastas)', 'BLOSUM62': '(fastas)',
                            'CKSAAGP': '(fastas)', 'CKSAAP': '(fastas)', 
                            'CTDC': '(fastas)', 'CTDD': '(fastas)', 
                            'CTDT': '(fastas)', 'CTriad': '(fastas)', 
                            'DPC': '(fastas)', 'EAAC': '(fastas)', 
                            'EGAAC': '(fastas)', 'GAAC': '(fastas)', 
                            'GDPC': '(fastas)', 'GTPC': '(fastas)', 
                            'KSCTriad': '(fastas)', 
                            'PAAC': '(fastas, lambdaValue = 10)', 
                            'QSOrder': '(fastas, nlag = 10)', 
                            'SOCNumber': '(fastas, nlag = 10)', 
                            'TPC': '(fastas, **proka_triplets)'}

euka_feature_index = [9, 21, 37, 92, 120, 175, 180, 237, 251, 273, 327, 
                      429, 518, 633, 1579, 1637, 1667, 1683, 1695, 1715, 
                      1722, 1724, 1737, 1781, 1783, 1801, 1811, 1812, 1896, 
                      1898, 1913, 1914, 1935, 1936, 1944, 1949, 1950, 1952, 
                      1959, 1961, 1979, 1983, 2023, 2028, 2081, 2095, 2105, 
                      2121, 2128, 2228, 2319, 2345, 2747, 2892, 3703, 3743, 
                      4036, 4056, 4656, 5030, 5052, 5055, 5102, 5105, 5127, 
                      5138, 5153, 5177, 5286, 5462, 5550, 5662, 6116, 6533, 
                      6543, 7272, 7278, 7582, 7725, 7901, 7989, 8203, 8223, 
                      8243, 8263, 8302, 8322, 8342, 8604, 8606, 8614, 8626, 
                      8696, 8699, 8721, 8767, 8887, 8967, 9008, 9028, 9067, 
                      9068, 9069, 9070, 9071, 9072, 9073, 9074, 9075, 9076, 
                      9077, 9078, 9079, 9080, 9081, 9082, 9083, 9084, 9085, 
                      9086, 9087, 9088, 9089, 9090, 9091, 9092, 9093, 9094, 
                      9095, 9096, 9097, 9098, 9099, 9100, 9101, 9102, 9103, 
                      9104, 9105, 9106, 9107, 9108, 9109, 9110, 9111, 9112, 
                      9113, 9114, 9115, 9116]
proka_feature_index = [7, 9, 15, 90, 415, 416, 639, 886, 1029, 1065, 1124, 
                       1320, 1909, 1948, 2091, 2148, 2185, 2214, 2244, 2254, 
                       2263, 2347, 2351, 2365, 2387, 2388, 2392, 2414, 2442, 
                       2445, 2447, 2450, 2463, 2567, 2608, 2612, 2617, 2626, 
                       2636, 2645, 2697, 2834, 2925, 2961, 3559, 3825, 4270, 
                       4277, 4297, 4477, 4497, 4517, 4537, 4740, 4890, 4977, 
                       4997, 5037, 5310, 5311, 5322, 5335, 5336, 5346, 5360, 
                       5364, 5371, 5381, 5396, 5406, 5410, 5411, 5419, 5444, 
                       5446, 5585, 5622, 5765, 5985, 6385, 6779, 6785, 7861, 
                       7863, 7864, 7866, 7867, 7870, 7873, 7883, 7887, 7889, 
                       7972, 8052, 8096, 8097, 8099, 8100, 8106, 8119, 8121, 
                       8123, 8457, 8601, 8638, 8781, 8949, 8961, 8981, 9001, 
                       9021, 9041, 9061, 9081, 9223, 9243, 9263, 9312, 9317, 
                       9322, 9327, 9331, 9332, 9337, 9342, 9372, 9382, 9387, 
                       9401, 9402, 9406, 9407, 9418, 9431, 9880, 9899, 9929, 
                       9949, 9956, 9959, 9992, 10008, 10009, 10010, 10011, 
                       10012, 10013, 10014, 10015, 10016]

ka_features = euka_features
ka_features_parameter = euka_features_parameter
ka_feature_index = euka_feature_index
if ModelType == 'proka':
    ka_features = proka_features
    ka_features_parameter = proka_features_parameter
    ka_feature_index = proka_feature_index

	


print('Start calcuating features, ', time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
print('Feature type: AAC')

fastas = readFasta.readFasta(tmpfile)
encodings = eval('AAC.AAC(fastas)')
pids = []
for i in range(len(encodings)):
    pids.append(encodings[i][0])

pids = pids[1:]

for feature in ka_features:
    myFun = feature + '.' + feature + ka_features_parameter[feature]
    print('Feature type: ' + feature)
    tmp_encodings = eval(myFun)
    for i in range(len(encodings)):
        encodings[i].extend(tmp_encodings[i][1:])

for i in range(len(encodings)):
    encodings[i] = [encodings[i][j] for j in ka_feature_index]

encodings = encodings[1:]
outFile  = 'tmp_' + filename + '.tsv' 
outFile  = filepath + outFile 
     
saveCode.savetsv2(encodings, outFile)

print('Finished at ', time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

#==============================================================================
#  Classification 
#==============================================================================
print('Start predicting pHis sites, ', time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

f = 'dat/' + ModelType + '_train.tsv'
if options.model:
	f = options.model

data = np.loadtxt(f, comments='#')
X_train = data[:,1:]
y_train = data[:,0]
scaler = preprocessing.MinMaxScaler().fit(X_train)
X_train = scaler.transform(X_train)
classifier = svm.SVC(kernel='rbf', probability=True, class_weight = 'balanced',gamma='scale').fit(X_train, y_train)


f = open(outFile)
X_test = np.loadtxt(f)  
f.close() 
X_test = scaler.transform(X_test)

prob = classifier.predict_proba(X_test)
labels = classifier.predict(X_test)

Files = open(outPutFileName, 'w')
Tabel = 'Histidine_site' + '\t' + 'Score' + '\t' + 'Class' + '\n'
Files.write(Tabel)
PrintResult(pids,prob[:,1],labels,Files) 
Files.close()

print('Finished at ', time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
print('pHisPred was completely done!')
print("%f seconds for classifying" % (time.time() - Compute_time) + ' ' +
      str(len(pids)) + ' ' + "histidine sites.\n")

if flg:
    os.remove(tmpfile)

os.remove(outFile)
