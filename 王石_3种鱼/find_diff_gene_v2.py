#!/public/store5/DNA/Test/zhengfuxing/conda/bin/python
# 挑选不同的排列组合
import re
import pandas as pd

def stat(t,h,l):        
    samples = h + l
    h0, h1, h2 = h
    l0, l1, l2 = l
    condition1 = df[h0] == df[h1]
    condition2 = df[h0] == df[h2] 
    condition3 = df[l0] == df[l1] 
    condition4 = df[l0] == df[l2]

    condition5 = df[h0] != df[l0]
    condition6 = df[l0] != './.' 
    condition7 = df[h0] != './.' 


    temp_df = df[condition1 & condition2 & condition3 & condition4 & condition5 & condition6 & condition7]
    all_list = last_name + samples + next_name
    temp_df = temp_df[all_list]
    outname = 'F1_F2' + '_' + t.replace('.annotation.xls','') + '_diff.xls'
    outlist = 'F1_F2' + '_' + t.replace('.annotation.xls','') +'_diff_genelist.xls'
    out = open(outlist, 'w')
    temp_df.to_csv(outname, index=False, sep='\t')
    gene_list_temp = temp_df['Gene.refGene'].to_list()
    # 修改基因名称
    pattern = 'gene-'
    gene_list = [i.replace(pattern,'') for i in list(set(gene_list_temp)) ]
    g_list = []
    for gene in gene_list:
        split_gene = re.split(',|;', gene)
        for b in split_gene:
            if b not in g_list:
                g_list.append(b)
    for i in g_list:
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
            stat(type,group_high,group_low)
