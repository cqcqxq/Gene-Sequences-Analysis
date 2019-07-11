from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq


#在genbank文件中提取序列记录
for seq_record in SeqIO.parse("KRas.gb", "genbank"):
	seqs = seq_record.seq
	seq_sum = 0
	'''for i in ('A' ,'T', 'C', 'G'): #该循环可以输出ATCG的数量
		seq_num = seqs.count(i)
		print(i+'碱基的个数是：' + str(seq_num))	
		seq_sum += seq_num
	print('碱基总数是：'+ str(seq_sum))'''
	for f in seq_record.features: #提取特征
		if f.type == 'CDS':	#确定是否为CDS类型
			#确定位置并提取CDS序列
			feature_seq = f.location.extract(seq_record).seq	
			print(feature_seq)
			coding_result = Seq(str(feature_seq), IUPAC.unambiguous_dna)
			#编码
			#print(coding_result.translate())
