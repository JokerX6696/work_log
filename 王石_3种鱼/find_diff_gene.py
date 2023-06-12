#!/public/store5/DNA/Test/zhengfuxing/conda/bin/python
# 挑选不同的排列组合
import re
import pandas as pd

def stat(t,h,l):        
    h = 'F2-1';l = 'F1-1'
    samples = [h,l]

    temp_df = df[df[h] != df[l]]
    all_list = last_name + samples + next_name
    temp_df = temp_df[all_list]
    outname = h + '_' + l + '_' + t + '_' + '_diff.xls'
    outlist = h + '_' + l + '_' + t + '_' +'_diff_genelist.xls'
    out = open(outlist, 'w')
    temp_df.to_csv(outname, index=False)
    gene_list_temp = temp_df['Gene.refGene'].to_list()
    # 修改基因名称
    pattern = 'gene-'
    gene_list = [i.replace(pattern,'') for i in list(set(gene_list_temp)) ]
    for i in gene_list:
        print(i,file=out)
    out.close

group_high = ['F2-1', 'F2-2', 'F3-3']
group_low = ['F1-1', 'F1-2', 'F1-3']


indel = 'indel.annotation.xls'
snp = 'snp.annotation.xls'
for type in [snp, indel]:
    df = pd.read_csv(type, sep='\t')
    df = df.replace("\s\|\s.+", "", regex=True)
    df['Gene.refGene'] = df['Gene.refGene'].replace("Gene-", "")
    colnames = df.columns.to_list()
    last_name = colnames[0:5]
    next_name = colnames[-5:]
    for h in group_high: #  h = 'F2-1';l = 'F1-1'
        for l in group_low:
            stat(type,h,l)
