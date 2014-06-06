topology-layout
===============

calculate 3d coordinates for  topology graph.

已知无向连通拓图（G,E)，根据E给G的每个元素计算一个3维坐标，达到以下目标：

    1、有连接的节点在3维空间中尽量靠近
    2、所有节点在空间中分布尽量均匀

算法整体描述：

    1、在G中选出若干节点，称“核心节点”，其集合为“核心节点集合”，记为C，设C的大小为N
    2、对每个核心节点，以之为根，广度优先遍历G，得到G中每个节点到各个核心节点的跳数（共N个），组成该节点的N维空间坐标向量
    3、运用主成分分析法（PCA）将每个节点的坐标向量从N维降为3维，得到每个节点的3维坐标
    
坐标节点的选取：

    1、选取连接数最多的节点u1作为第1个核心节点，加入C
    2、对G-C中的节点u，定义u到C的距离为u到C中已有节点的距离的最小值
    3、将当前G-C中与C距离最远的节点作为下一个核心节点加入C，更新G-C所有节点到C的距离
    4、将3重复N-1次，一共得到N个核心节点. N大小取50足够.
    
    
详见 《Graph Drawing by High-Dimensional Embedding》
http://www.emis.de/journals/JGAA/accepted/2004/HarelKoren2004.8.2.pdf

    
