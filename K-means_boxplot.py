import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

# 生成示例数据
np.random.seed(0)
X = np.random.rand(100, 2)

# 设置 K 值范围
k_values = range(2, 11)

# 存储每个 K 值下的聚类误差（SSW）
ssw_values = {}

# 计算每个 K 值下的聚类误差（SSW）
for k in k_values:
    kmeans = KMeans(n_clusters=k, random_state=0)
    kmeans.fit(X)
    
    # 获取每个样本点到聚类中心的距离平方
    keyname = 'k = ' + str(k) 
    distances = kmeans.transform(X)
    cluster_errors = list(np.min(distances, axis=1) ** 2)
    
    # 计算聚类误差（SSW）
    ssw = np.sum(cluster_errors)
    ssw_values[keyname] = cluster_errors
# 提取 x 轴刻度和对应的数据
x_labels = list(ssw_values.keys())
data_values = list(ssw_values.values())
# 绘制箱线图
plt.boxplot(data_values, labels=x_labels)
# 设置标题和轴标签
plt.title('Boxplot by Categories')
plt.xlabel('Categories')
plt.ylabel('Values')
# 显示图形
plt.show()

# 打印每个 K 值下的聚类误差（SSW）
# for k, ssw in zip(k_values, ssw_values):
#     print(f"K={k}: SSW={ssw}")