//
//  maxFlow.hpp
//  genic
//
//  Created by sanshanxiashi on 2017/3/23.
//  Copyright © 2017年 sanshanxiashi. All rights reserved.
//

#ifndef maxFlow_hpp
#define maxFlow_hpp

#include<iostream>
#include<algorithm>
#include<string>
#include<sstream>
#include<set>
#include<vector>
#include<stack>
#include<map>
#include<queue>
#include<deque>
#include<cstdlib>
#include<stdio.h>
#include<cstring>
#include<cmath>
#include<ctime>
#include<functional>
#define maxn 10000
#define INF 2147483647
using namespace std;


int maxFlow(struct MCMF& hua);

struct Edge{
    int from, to ,cap ,flow, cost;
};

struct Node{
    int f,indx;
    double averCost;
};


struct MCMF{
    
    int n,m;
    int s , t , k;
    int serverCost, netNode , consumeNode,numOfserver;
    
    vector<Edge> edges;
    vector<int> G[maxn];
    struct Node degree[maxn]; //某个网络点的出度
    int consumeEdge[505];  //与消费点相连的边的 编号
    int consumeFlow[505];  //消费点的容量
    int isNetToCons[1005]; //判断网络点i是否与某个消费点相连，-1为不相连，>=0为消费点的编号
    int inq[maxn];
    int d[maxn];
    int p[maxn];
    int a[maxn];
    
    int costNow;
    int consumeMaxFlow; //最大的 消费点需求容量
    
    //for search path
    int fa[1505];
    int find_flow;
    int p_edge[1505]; //点对应的边的编号
    int ok;
    int find_vis[1505];

    
    //ceshi
    int conflow[505], conflowNeed[505];
    
    //输出
    vector<int> res[maxn];  // p条输出路径
    vector<int> outFlow;  //p条路径的 各自流量
    int outCnt; //记录有多少条输出路径
    
    
    void init(int n,int netNode,int consumeNode,int serverCost);
    void AddEdge(int from , int to, int cap , int cost) ;
    bool BellmanFord(int s, int t , int& flow, int& cost, int& restF);
    int Mincost(int s, int t, int k);
    void search_path();
    void update_flow(int minflow , int fa[]);
    void dfs_findPath(int cur, int minflow);
    void OutputPath(string& str);
};

int getServersAndGetMincost(vector<int>& serverList, int& serverNum ,struct MCMF hua, struct MCMF& huaReal, string& str, int printPath);
//用来排序的判断函数
bool cmp(Edge a, Edge b);

int DirectLink(vector<int>& serverList, int& serverNum ,struct MCMF hua ,struct MCMF& huaReal , string& str);

void  minDisofXY(struct MCMF hua ,int i ,int j, int cmin[1005][505]);

int newGetServersAndGetMincost(vector<int>& serverList,vector<int> consumeList,int& serverNum,MCMF hua,MCMF& huaReal);
int newGetServersAndGetMincost2(vector<int>& serverList,vector<int> consumeList,int& serverNum,MCMF& hua, string& str);
#endif /* maxFlow_hpp */
