from Bio import SeqIO
from Bio.Seq import Seq
from CAI import RSCU
from matplotlib import pyplot as plt
import numpy as np

for seq_record in SeqIO.parse("KRas.gb", "genbank"):
	for f in seq_record.features:
		rscu_list = []
		if f.type == 'CDS':
			feature_seq = f.location.extract(seq_record).seq
			coding_result = Seq(str(feature_seq))
			# 新建一个列表,列表才能作为输入数据
			rscu_list.append(coding_result)
			#print(rscu_list)
			codon_pre = RSCU(rscu_list)
			#print(codon_pre)
		# 满足rscu不为空就退出
		if len(rscu_list): 
			break
code = []
value = []

for k in codon_pre:
	code.append(k)
	values = codon_pre[k]
	value.append(values)

plt.switch_backend('Agg')
plt.figure(figsize=(20, 6), dpi=80)
width = 0.8
fig = plt.bar(code, value, width, color ='blue')
plt.ylabel('RSCU')
plt.xlabel('Codons')
plt.title('RSCU of KRas gene')
plt.legend(loc='upper left')
plt.savefig("Codon.png")

#排序输出
#for k in sorted(codon_pre, key = codon_pre.__getitem__, reverse = True):
	#print(k, codon_pre[k])

