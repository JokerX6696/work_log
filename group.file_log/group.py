##
#!/data/USER/zhengfuxing/conda/after-sale/bin/python
import re
import pandas as pd

file = 'D:/desk/make_Script/DZOE2023041063_2023052905.xlsx'

grouped_file = 'D:/desk/make_Script/config.ini'

#################排除已分析过的组  grouped 已经分析过得到常规分组
fd = open(grouped_file, 'r',encoding='UTF-8')
df_fd = fd.read()
fd.close
grouped = re.search('group_vs\s=\s?(.+)\n', df_fd).group(1).replace(' ','').replace('\t','').split(';')
##########################################

all_samples = pd.read_excel(file,sheet_name='样本信息')
all_samples = [i for i in all_samples['样本分析名称']]

########################################
df = pd.read_excel(file,sheet_name='样本分组信息')
# 将分组信息储存到字典 key 为组 value 样本的列表
dict_group = {}

for j in range(0, df.shape[0]):
    temp_samples = df.iloc[j]['samples'].split(',')
    group_name = df.iloc[j]['group']
    for l in temp_samples:
        if group_name in dict_group:
            dict_group[group_name].append(l)
        else:
            dict_group[group_name] = [l]

# 确定补充拓展分组的数目
group_vs = pd.read_excel(file,sheet_name='差异比较信息')
group_vs_new = []
for i in group_vs['groups']:
    k = i.replace(' ','').replace('\t','')
    if k in grouped:
        continue
    else:
        group_vs_new.append(k)

group_num = len(group_vs_new) # 待分组的数目
group_al = []#已完成的分组
data = ""
split_char = '################################################\n'
while( group_num > 0 ):
    ll = set()
    jj = False
    group_temp = []
    for j in group_vs_new:
        if j not in group_al:
            groups_raw = j.replace(' ','').replace('\t','')
            groups = groups_raw.split(',')
            for group in groups:
                if ll & set(dict_group[group]) == set():  # 无交集
                    ll = ll | set(dict_group[group])
                    for q in dict_group[group]:
                        data += q + '\t' + group + '\n'
                        #print(q,group,sep='\t')
                else:  # 有交集
                    jj = True
                    break
        else: # 如果已经被分组了 则跳过 
            continue
        if jj:
            None
        else:
            group_al.append(j)
            group_num -= 1
            group_temp.append(groups_raw)
        #print('------------------------------------------------')
    data += "#  "  
    for g in group_temp:
        data += g + ';'
    data += '\n' +split_char
final_data = data.split(split_char)


