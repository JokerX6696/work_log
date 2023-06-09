#!/public/store5/DNA/Test/zhengfuxing/conda/bin/python
import pandas as pd
import re
species = 'tuantoufang'
size_file = '/public/store5/DNA/Project/Reseq/2023/DZOE2023041442-b1/' + species + '/reseq_v3/ref/genome.fa.sizes.chrom'
indel_file = '/public/store5/DNA/Project/Reseq/2023/DZOE2023041442-b1/' + species + '/reseq_v3/result/04.annotation/anno_stat/indel.annotation.xls'
snp_file = '/public/store5/DNA/Project/Reseq/2023/DZOE2023041442-b1/' + species + '/reseq_v3/result/04.annotation/anno_stat/snp.annotation.xls'
out_file_indel = '/public/store5/DNA/Project/Reseq/2023/DZOE2023041442-b1/PA/SNV_Density/' + species + '/' + species + '_InDel_stay.xlsx'
out_file_snp = '/public/store5/DNA/Project/Reseq/2023/DZOE2023041442-b1/PA/SNV_Density/' + species + '/' + species + '_Snp_stay.xlsx'
f_chr = open(size_file, 'r')
chr_file = f_chr.readlines()
f_chr.close

all_dict = {}
all_dict['counts'] = {}
all_dict['length'] = {}
chr_len = {}
pattern = '\./\.|0/0'
for i in chr_file:
    nm = i.split('\t')[0]
    l = i.split('\t')[1].replace('\n','')
    chr_len[nm] = int(l)


df = pd.read_csv('./indel.annotation.xls',sep='\t')
samples_raw = df.columns.tolist()[5:]
samples = []
for i in samples_raw:
    if i == 'Func.refGene' or i == 'Gene.refGene' or i == 'GeneDetail.refGene' or i == 'ExonicFunc.refGene' or i == 'AAChange.refGene':
        break
    else:
        samples.append(i)
del(samples_raw)
# 只保留指定染色体
df = df[df['chromosome'].isin(chr_len)]

for sample in samples:  
    df_temp = df[~df[sample].str.contains(pattern)] # 返回 bool 正则匹配后逻辑取反 使用 ~
    chr_counts = df_temp['chromosome'].value_counts().to_dict() # 统计不同样本 染色体数目 出现次数 存入字典
    all_dict['counts'][sample] = chr_counts
    all_dict['length'][sample] = chr_len

counts_df = pd.DataFrame(all_dict['counts'])
len_df = pd.DataFrame(all_dict['length'])
new_df = len_df.div(counts_df)

writer = pd.ExcelWriter(out_file_indel, engine='xlsxwriter')
new_df.to_excel(writer, sheet_name='Density', index=True)
counts_df.to_excel(writer, sheet_name='counts', index=True)
len_df.to_excel(writer, sheet_name='chr_length', index=True)
writer.save()



all_dict = {}
all_dict['counts'] = {}
all_dict['length'] = {}
chr_len = {}
pattern = '\./\.|0/0'
for i in chr_file:
    nm = i.split('\t')[0]
    l = i.split('\t')[1].replace('\n','')
    chr_len[nm] = int(l)


df = pd.read_csv('./snp.annotation.xls',sep='\t')
samples_raw = df.columns.tolist()[5:]
samples = []
for i in samples_raw:
    if i == 'Func.refGene' or i == 'Gene.refGene' or i == 'GeneDetail.refGene' or i == 'ExonicFunc.refGene' or i == 'AAChange.refGene':
        break
    else:
        samples.append(i)
del(samples_raw)
# 只保留指定染色体
df = df[df['chromosome'].isin(chr_len)]

for sample in samples:  
    df_temp = df[~df[sample].str.contains(pattern)] # 返回 bool 正则匹配后逻辑取反 使用 ~
    chr_counts = df_temp['chromosome'].value_counts().to_dict() # 统计不同样本 染色体数目 出现次数 存入字典
    all_dict['counts'][sample] = chr_counts
    all_dict['length'][sample] = chr_len

counts_df = pd.DataFrame(all_dict['counts'])
len_df = pd.DataFrame(all_dict['length'])
new_df = len_df.div(counts_df)

writer = pd.ExcelWriter(out_file_snp, engine='xlsxwriter')
new_df.to_excel(writer, sheet_name='Density', index=True)
counts_df.to_excel(writer, sheet_name='counts', index=True)
len_df.to_excel(writer, sheet_name='chr_length', index=True)
writer.save()




