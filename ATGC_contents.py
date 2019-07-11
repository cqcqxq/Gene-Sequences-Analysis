from Bio import SeqIO
from Bio.Seq import Seq
from matplotlib import pyplot as plt
plt.switch_backend('Agg')


#在genbank文件中提取序列记录
for seq_record in SeqIO.parse("KRas.gb", "genbank"):
	seqs = seq_record.seq
	A = []
	T = []
	C = []
	G = []
	for i in range(len(seqs)-100):
		window = 100
		moving_area = seqs[0+i:window+i]
		A_ratio = moving_area.count('A')/100
		A.append(A_ratio)
		T_ratio = moving_area.count('T')/100
		T.append(T_ratio)
		C_ratio = moving_area.count('C')/100
		C.append(C_ratio)
		G_ratio = moving_area.count('G')/100
		G.append(G_ratio)
		'''print('A的含量是' + str(A_covers))
		print('T的含量是' + str(T_covers))
		print('C的含量是' + str(C_covers))
		print('G的含量是' + str(G_covers))'''

fig = plt.figure(dpi = 128, figsize = (10, 6))
plt.plot(A, 'black', label="A")
plt.plot(T, 'red', label="T")
plt.plot(C, 'green', label="C")
plt.plot(G, 'blue', label="G")
plt.xlabel('Position', fontsize=16)
plt.ylabel('A/T/G/C %', fontsize=16)
plt.title('Contents of bases', fontsize=20)
plt.legend(loc='upper left')
plt.savefig("ATGC.png")
