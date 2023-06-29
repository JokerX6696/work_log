#!D:/application/anaconda/python.exe
import numpy as np
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
# delta K
# 生成示例数据
np.random.seed(0)
X = np.random.rand(100, 2)

# 初始化变量
max_clusters = 10
distances = []

# 计算不同聚类数下的聚类结果和平均距离
for k in range(1, max_clusters+1):
    kmeans = KMeans(n_clusters=k, random_state=0)
    kmeans.fit(X)
    avg_distance = np.mean(np.min(kmeans.transform(X), axis=1))
    distances.append(avg_distance)

# 对平均距离进行对数变换
log_distances = np.log(distances)

# 计算相邻聚类数对应的对数变换值之间的差异
differences = np.diff(log_distances)  # 

# 找到具有最大差异值的聚类数
optimal_clusters = np.argmax(differences) + 2  # +2是因为从2开始聚类数递增

# 绘制Delta K曲线
plt.plot(range(2, max_clusters+1), differences)
plt.xlabel('Number of clusters')
plt.ylabel('Delta K')
plt.title('Delta K Method')
plt.show()

print("Optimal number of clusters:", optimal_clusters)


