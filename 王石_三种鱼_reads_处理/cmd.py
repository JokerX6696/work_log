#!/bin/env python3
samples = ['CB-3', 'F1-2', 'F1-3', 'F2-1', 'F2-2', 'F3-3','F1-1', 'NCRC-F1-1']
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
    fo = open((sample + '_stat.txt'),'w')
    # 三品种共有
    cm = len(jiyu & liyu & tuantoufang)
    print('三样本共有reads数目:',cm,sep=' ',file=fo)
    # 特有
    uniq_J = len(jiyu - (liyu | tuantoufang))
    print('jiyu 特有 reads数目:',uniq_J,sep=' ',file=fo)

    uniq_L = len(liyu - (jiyu | tuantoufang))
    print('liyu 特有 reads数目:',uniq_L,sep=' ',file=fo)

    uniq_T = len(tuantoufang - (jiyu | liyu))
    print('tuantoufang 特有 reads数目:',uniq_T,sep=' ',file=fo)

    # 交集

    j_l = len((jiyu & liyu) - tuantoufang)
    print('maped to liyu jiyu  reads数目:',j_l,sep=' ',file=fo)
    j_t = len((jiyu & tuantoufang) - liyu)
    print('maped to jiyu tuantoufang  reads数目:',j_t,sep=' ',file=fo)
    l_t = len((liyu & tuantoufang) - jiyu)
    print('maped to liyu tuantoufang  reads数目:',l_t,sep=' ',file=fo)
    fo.close
