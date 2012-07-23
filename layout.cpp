//============================================================================
// Name        : layout.cpp
// Author      : Chen Kai
// Version     :
// Copyright   : Your copyright notice
// Description : layout for bgp。
//============================================================================

#include <iostream>
#include <stdio.h>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define MD 			 	50
#define ADJ  			300.0				//坐标扩张系数
#define EOR  			0.0001				//迭代法计算特征向量时迭代结束的界限
#define MAX_RUN  		20			//迭代法计算特诊向量时迭代次数限制

using namespace std;

struct LinkStruct
{
	int peer;
	LinkStruct* next;
};

struct Node
{
	int l;
	LinkStruct* links;
};

struct NodeSort
{
	int l;
	int index;
};

typedef double Array50[50];
typedef double Array3[3];

int compare (const void * a, const void * b)
{  //compare的返回值应表示a>b 或a==b 或 a<b 你可用正数、0、负数表示，只要返回值包含了这三种取值就行了，一般情况下，常返回两数相减的结果
  return ( ((NodeSort*)a)->l < ((NodeSort*)b)->l );
}

int main(int argc, char* argv[]) {

	clock_t  clockBegin, clockEnd,t1, t2;
	clockBegin = clock();

	int total_nodes,peer;
	double delta = 10.0;
	double adj  = 1.0 ;
	string inputFile = "input.txt";
	string outputFile = "output.txt";

	if (argc > 1)
		inputFile = argv[1];
	if (argc > 2)
		outputFile = argv[2];
	if (argc > 3)
		adj = (double)atoi(argv[3]);
	if (argc > 4)
		delta = (double)atoi(argv[4]);

	LinkStruct* pLink = NULL;
	LinkStruct* pLast = NULL;

	freopen(inputFile.c_str(),"r",stdin);

	cin >> total_nodes;
	cout << total_nodes << endl;

	Node* nodes = new Node[total_nodes];
	memset(nodes, 0, sizeof(nodes));

	for (int i = 0; i < total_nodes; i++)
	{
		cin >> nodes[i].l;
		if (!nodes[i].l)
			continue;
		cin >> peer;
		pLink = new LinkStruct;
		pLink->peer = peer;
		pLink->next = NULL;
		nodes[i].links = pLink;
		pLast = pLink;

		for (int j = 1, l = nodes[i].l; j < l; j++)
		{
			cin >> peer;
			pLink = new LinkStruct;
			pLink->peer = peer;
			pLink->next = NULL;
			pLast->next = pLink;
			pLast = pLink;
		}
	}

	NodeSort* nodeSort = new NodeSort[total_nodes];
	for (int i = 0; i < total_nodes; i++)
	{
		nodeSort[i].l = nodes[i].l;
		nodeSort[i].index = i;
	}

	qsort (nodeSort, total_nodes, sizeof(NodeSort), compare);

	int* nodeHash = new int[total_nodes];
	for (int i = 0; i < total_nodes; i++)
		nodeHash[nodeSort[i].index] = i;

	for (int i = 0; i < total_nodes; i++)
	{
		pLink = nodes[i].links;
		for (int j = 0, l = nodes[i].l; j < l; j++)
		{
			pLink->peer = nodeHash[pLink->peer];
			pLink = pLink->next;
		}
	}

	delete []nodeHash;
	delete []nodeSort;


	int h,t,p;  				 						// 队列头与对列尾
	double ss[MD][MD];			 						//协方差矩阵ss

	Array50* xx = new Array50[total_nodes];		   		// 高维坐标数据
	int* deg = new int[total_nodes];					// 到当前BFS中心节点的距离
	int* d = new int[total_nodes];					 	// 到所有中心节点的最短距离
	int* q = new int[total_nodes];						// BFS搜索 队列


	//初始化
	for (int i = 0; i< total_nodes; i++)
		d[i] = 100000;

	//memset(xx,0,sizeof(xx));      ===========ERROR!!!!!=========
	for (int i = 0; i < total_nodes; i++)
		for (int j = 0; j < MD; j++)
			xx[i][j] = 0;

	memset(ss,0,sizeof(ss));
	memset(deg,0,sizeof(deg));

	cout << "BFS.......";

	t1 = clock();

	//p = r; //50个中心节点的第一个
	p = rand() % 10;
	for (int k = 0; k < MD; k++)
	{
		//广度优先搜索
		memset(q,0,sizeof(q));
		for (int i = 0; i< total_nodes; i++)
			deg[i] = -1;

		h = 0, t = 0;
		q[t++] = p;
		deg[p] = 0;
		while ( h < t){

			pLink = nodes[ q[h] ].links;
			for (int j = 0,l = nodes[ q[h] ].l; j < l; j++){
				int tmp = pLink->peer;
				pLink = pLink->next;
				if (deg[tmp] > -1) continue;
				q[t++] = tmp;
				deg[tmp] = deg[ q[h] ] + 1;
			}
			h++;
		}

		//构造原始高维坐标数据，确定下一个坐标
		int  far = 0;
		for (int i = 0; i < total_nodes; i ++)
		{
			xx[i][k] = deg[i];

			if (deg[i] < d[i])
				d[i]  = deg[i];

			if (d[i] > far){
				far  = d[i];
				p = i;
			}
		}
	}

	t2 = clock();
	cout << "BFS time : " << (t2 - t1)/1000.0 << "s" << endl;

	//释放内存
	for (int i = 0; i < total_nodes; i++)
	{
		pLink = nodes[i].links;
		while (pLink != NULL){
			pLast = pLink->next;
			delete(pLink);
			pLink = pLast;
		}
	}

	delete []nodes;
	delete []deg;
	delete []d;
	delete []q;

	cout << "average......." << endl;

	double average;
	//预处理，去平均值
	for ( int i = 0; i < MD; i++){

		average = 0.0;
		for (int j = 0; j < total_nodes; j++)
			average += xx[j][i] ;
		average = average / total_nodes;

		for (int j = 0; j < total_nodes; j++)
			xx[j][i] -= average;
	}

	cout << "cov.......";

	t1 = clock();

	//计算协方差矩阵ss
	for ( int i = 0; i < MD; i ++)
		for (  int j = 0; j < MD; j++ ){
			if (i > j ) {
				ss[i][j] = ss[j][i];
				continue;
			}
			ss[i][j] = 0;
			for (int v = 0; v < total_nodes; v++)
				ss[i][j] += xx[v][i] * xx[v][j];

		}

	t2 = clock();
	cout << "cov time : " << (t2 - t1)/1000.0 << "s" << endl;

	//转换为相关系数矩阵
	for ( int i = 0; i < MD; i++)
		for ( int j = 0; j < MD; j++)
			ss[i][j] = ss[i][j] / sqrt( (double)ss[i][i] * ss[j][j]);

	//幂迭代法求协方差矩阵ss的特征值最大的3个特征向量
	double uuu[3][50];
	double uu[50];
	double ml,ee;
	time_t ts;
	srand((unsigned) time(&ts));

	memset(uuu,0,sizeof(uuu));
	memset(uu,0,sizeof(uu));

	cout << "vec......." << endl;

	for (int i = 0; i < 3; i++){

		//生成迭代的初始向量
		ml = 0.0;

		for (int j = 0; j < MD; j++){
			uu[j] = rand() / (double)RAND_MAX;
			ml += uu[j] * uu[j];
		}

		//初始向量标准化
		ml = sqrt(ml);
		for (int j = 0; j < MD; j++)
			uu[j] = uu[j] / ml;

		int run = 0;
		for (;;)
		{
			run ++;
			for (int j = 0; j< MD; j++)
				uuu[i][j] = uu[j];

			//与已生成的主成分向量正交化
			for (int j = 0; j < i; j ++){
				double sss = 0.0;
				for (int k = 0; k < MD; k ++)
					sss += uuu[i][k] * uuu[j][k];
				for (int k = 0; k < MD; k ++)
					uuu[i][k] -= sss * uuu[j][k];
			}

			//迭代计算下一个
			ml = 0.0;
			for (int j = 0; j < MD; j ++ ){
				uu[j] = 0.0;
				for (int k = 0; k < MD; k++)
					uu[j] += ss[j][k] * uuu[i][k];
				ml += uu[j] * uu[j];
			}

			//标准化
			ml = sqrt(ml);
			for (int j = 0; j < MD; j++)
				uu[j] = uu[j] / ml;

			//计算前后改变量，决定是否终止。
			ml = 0.0;
			for (int j = 0; j < MD; j++)
				ml += uu[j] * uuu[i][j];
			ee = 1- ml;

			if (ee < EOR || run > MAX_RUN)
				break;
		}

		for (int j = 0; j< MD; j++)
				uuu[i][j] = uu[j];
	}

	//投影
	double aX = 0.0, aY = 0.0, aZ = 0.0;
	double bX = 0.0, bY = 0.0, bZ = 0.0;

	Array3* pos = new Array3[total_nodes];
	//memset(pos,0,sizeof(pos));                 ==============!!!!ERROR!!!!!=============
	for (int i = 0; i < total_nodes; i++)
		for (int j = 0; j < 3; j++)
			pos[i][j] = 0;

	cout << "pos......." << endl;

	for (int i = 0; i < total_nodes; i ++){
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < MD; k++ )
				pos[i][j] += xx[i][k] * uuu[j][k];
		}

		double tmp = pos[i][0];
		pos[i][0] = pos[i][1];
		pos[i][1] = -tmp;
	}

	for (int i = 0; i < total_nodes; i ++){

		aX += abs(pos[i][0]);
		aY += abs(pos[i][1]);
		aZ += abs(pos[i][2]);
		bX += pos[i][0];
		bY += pos[i][1];
		bZ += pos[i][2];
	}

	cout << "adj......." << endl;

	//整体调整
	aX = aX/aZ;
	aY = aY/aZ;
	bX = bX/(double)total_nodes;
	bY = bY/(double)total_nodes;
	bZ = bZ/(double)total_nodes;

	//delta 随机偏移量，解决有些节点的连接情况完全相同时位置重合的问题。

	for (int i = 0; i < total_nodes; i ++){
		pos[i][0] = (pos[i][0] - bX) * ADJ * adj + delta* rand()/(double)RAND_MAX;
		pos[i][1] = (pos[i][1] - bY) * ADJ / aY * aX / 1.5 * adj+ delta*rand()/(double)RAND_MAX;
		pos[i][2] = (pos[i][2] - bZ) * ADJ * aX * adj + delta*rand()/(double)RAND_MAX;
	}

	clockEnd = clock();
	cout << "total time : " << (clockEnd - clockBegin)/1000.0 << "s" << endl;

	freopen(outputFile.c_str(),"w",stdout);
	for (int i = 0; i < total_nodes; i ++)
		cout << pos[i][0] << " " << pos[i][1] << " " << pos[i][2] << endl;

	delete []pos;
	delete []xx;

	clockEnd = clock();
	cout << clockEnd - clockBegin << endl;

	return 0;
}
