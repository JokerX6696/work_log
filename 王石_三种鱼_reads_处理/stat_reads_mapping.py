#!/public/store5/DNA/Test/zhengfuxing/conda/bin/python


samples_prefix = ['CB-3', 'F1-2', 'F1-3', 'F2-1', 'F2-2', 'F3-3','F1-1', 'NCRC-F1-1']
samples = []
for prefix in samples_prefix:
    samples.append('/public/store5/DNA/Project/Reseq/2023/DZOE2023041442-b1/PA/reads_name/' + prefix)


for sample in samples:
    fj = open((sample + '_jiyu_mapped.txt'),'r')
    jiyu = set(fj.readlines())
    fj.close

    fl = open((sample + '_liyu_mapped.txt'),'r')
    liyu = set(fl.readlines())
    fl.close

    ft = open((sample + '_tuantoufang_mapped.txt'),'r')
    tuantoufang = set(ft.readlines())
    ft.close

    # 三品种共有
    sample2 = sample.replace('/public/store5/DNA/Project/Reseq/2023/DZOE2023041442-b1/PA/reads_name/','')
    f_name = sample2 + '_common_readsname.xls'
    f = open(f_name, 'w')
    cm = jiyu & liyu & tuantoufang
    for read in cm:
        print(read,end="",file=f)
    f.close
    del(cm)

    # 交集

    f_name = sample2 + '_j_l_readsname.xls'
    f = open(f_name, 'w')
    j_l = (jiyu & liyu) - tuantoufang
    for read in j_l:
        print(read,end="",file=f)
    f.close
    del(j_l)


    f_name = sample2 + '_j_t_readsname.xls'
    f = open(f_name, 'w')
    j_t = (jiyu & tuantoufang) - liyu
    for read in j_t:
        print(read,end="",file=f)
    f.close
    del(j_t)

    f_name = sample2 + '_l_t_readsname.xls'
    f = open(f_name, 'w')
    l_t = (liyu & tuantoufang) - jiyu
    for read in l_t:
        print(read,end="",file=f)
    f.close
    del(l_t)



