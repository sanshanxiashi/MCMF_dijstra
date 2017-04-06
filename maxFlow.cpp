//
//  maxFlow.cpp
//  genic
//
//  Created by sanshanxiashi on 2017/3/23.
//  Copyright © 2017年 sanshanxiashi. All rights reserved.
//

#include "maxFlow.hpp"
#include <iostream>
#include <algorithm>
#include <string>
#include <sstream>
#include <set>
#include <vector>
#include <stack>
#include <map>
#include <queue>
#include <deque>
#include <stdio.h>
#include <cstring>
#include <cmath>
#include <ctime>
#include <stdlib.h>
#include <functional>
#define maxn 10000

using namespace std;


    
    
    void MCMF:: init(int n,int netNode,int consumeNode,int serverCost){
        this->n = n;
        this->netNode= netNode;
        this->consumeNode= consumeNode;
        this->serverCost= serverCost;
        this->numOfserver = 0;
        
        this->s= netNode + consumeNode;
        this->t= this->s + 1;
        
        
        for(int i=0; i<n; i++) G[i].clear(), res[i].clear(), outFlow.clear() ;
        memset(isNetToCons, -1 ,sizeof(isNetToCons));
        edges.clear();
        memset(a,0,sizeof(a));
    }
    
    void MCMF:: AddEdge(int from , int to, int cap , int cost) {
        edges.push_back((Edge){from , to ,cap , 0, cost});
        edges.push_back((Edge){to , from ,0 , 0, -cost});
        m = (int)edges.size();
        G[from].push_back(m-2);
        G[to].push_back(m-1);
    }
    
    //沿最短路增广的过程
    bool MCMF:: BellmanFord(int s, int t , int& flow, int& cost, int& restF){
        
        for(int i=0; i<n; i++) d[i] = INF;
        memset(inq , 0, sizeof(inq));
        d[s] = 0; inq[s] = 1; p[s] = 0; a[s] = INF;
        
        queue<int> Q;
        Q.push(s);
        while(!Q.empty()){
            int u = Q.front(); Q.pop();
            
            //printf("in u=%d\n",u);
            
            inq[u] = 0;
            for(int i=0; i<G[u].size(); i++){
                Edge& e= edges[G[u][i]];
                if(e.cap >e.flow && d[e.to]>d[u]+e.cost){
                    d[e.to] = d[u] + e.cost;
                    p[e.to] = G[u][i];
                    a[e.to] = min(a[u], e.cap - e.flow);
                    if(!inq[e.to])  { Q.push(e.to); inq[e.to] = 1; }
                }
            }
            //   printf("out\n");
        }
        
        //printf("d[t]=%d\n",d[t]);
        if(d[t] == INF) return false; //s-t不连通，失败退出
        //找到增广路径
        a[t]= min(a[t], restF);
        restF = restF - a[t];
        flow += a[t];
        cost += d[t] * a[t];
        
        //printf("a[t]= %d , restF=%d , flow = %d , cost = %d \n",a[t], restF, flow, cost);
        
        int u=t;
        while(u!=s){
            if(u==t)
            {
                conflow[edges[p[u]].from - netNode] +=a[t];
              //  printf("%d : %d\n",edges[p[u]].from - netNode , conflow[edges[p[u]].from - netNode]);
            }
            //记录当前的增广路径
            res[outCnt].push_back(u);
          //  printf("%d ",u);
            edges[p[u]].flow += a[t];
            edges[p[u]^1].flow -= a[t];
            u = edges[p[u]].from;
        }
        
        outCnt++;
        outFlow.push_back(a[t]);
      //  printf(": %d\n",a[t]);
        
        
        
        if(!restF) return false;
        else return true;
    }
    
    
    //最小费用流算法
    int MCMF:: Mincost(int s, int t, int k){
        int flow = 0, cost = 0, restF = k;
        outCnt=0;
        memset(conflow, 0, sizeof(conflow));
        while(BellmanFord(s, t, flow, cost, restF));
      //  printf("\nPathCost= %d",cost);
        cost += serverCost*numOfserver;
      //  printf("        totCost= %d\n",cost);
        
       /*
        for(int i=0; i<consumeNode; i++){
            printf("%d: real =%d  ,now =%d ------ ",i, conflowReal[i], conflow[i]);
            if(conflowReal[i]!=conflow[i]) puts("caooooo");
            else printf("\n");
        }*/
        
      //  printf("flow = %d , k =%d\n",flow, k);
        if(flow < k) cost = -1;   //流量不能满足要求
        
        this->costNow = cost;
        return cost ;
    }
    
    
    void MCMF:: OutputPath(string& str){
        char c[50];
        
        printf("outCnt = %d\n\n", outCnt);
        sprintf(c, "%d\n\n", outCnt);
        str += c;
        
        for(int i=0; i< outCnt ;i++){
            for(int j= (int)res[i].size()-1 ; j>0; j--)
            {
                printf("%d ",res[i][j]>=netNode?(res[i][j]-netNode):res[i][j]);
                sprintf(c,"%d ",res[i][j]>=netNode?(res[i][j]-netNode):res[i][j]);
                str +=c;
            }
            printf("%d\n",outFlow[i]);
            sprintf(c, "%d\n",outFlow[i]);
            str +=c;
        }
    }



    void MCMF:: update_flow(int minflow , int fa[]){
    //从s开始
    int cnt = outCnt;
   // printf("cnt=%d: minflow=%d\n",cnt, minflow);
    int u=t;
    while(u!=s){
        
        //流量更新，减掉minflow
      //  if(p[u]==s) printf("p_edeg=%d ,s.flow=%d ",p_edge[u], edges[p_edge[u]].flow);
        edges[p_edge[u]].flow -=minflow;
        if(u!=s) res[cnt].push_back(u);
      //   if(p[u]==s) printf("^^^^ , s.flow=%d \n",edges[p_edge[u]].flow);
        u = fa[u];
    }
    outFlow.push_back(minflow);

        
    find_flow += minflow;
    outCnt++;
}


//cur表示当前点
void MCMF:: dfs_findPath(int cur, int minflow){
    
    if(ok) return;
    
    //从s正方向找
    if(cur == t){
        
        //更新整条路上的 流量
        update_flow(minflow, fa);
        ok=1;
        return ;
    }
    
    for(int i=0; i<G[cur].size(); i++){
        Edge e = edges[G[cur][i]];
      //  printf("e.to=%d e.flow =%d ,vis=%d\n",e.to, e.flow , find_vis[e.to]);
        if(e.flow<=0) continue; //反向边不搜索
        if(find_vis[e.to]) continue;
        
        fa[e.to] = cur;
        p_edge[e.to] = G[cur][i];
        find_vis[e.to] = 1;
        dfs_findPath(e.to , min(minflow, e.flow));
        find_vis[e.to] = 0;
        
        if(ok) return;
    }
    
}




    void MCMF::search_path(){

        printf("search start\n");
        outCnt = 0;
        for(int i=0; i<1500;i++)
        {
            res[i].clear();
            outFlow.clear();
        }
        find_flow = 0;
        fa[s] = -1;
        while(find_flow< k ){
            
            for(int i=0; i<(int)G[s].size(); i++){
                //puts("ss");
                Edge e = edges[G[s][i]];
               // printf("i=%d, e.flow =%d e.from =%d, e.to =%d\n",i, e.flow , e.flow , e.to);
                fa[e.to] = s;
                p_edge[e.to] = G[s][i];
                
                if(e.flow>0) {  //搜路径 ,初始为服务点
                    ok=0;
                    memset(find_vis, 0, sizeof(find_vis));
                    find_vis[s] = 1;
                    dfs_findPath(e.to , INF);
                   // printf("find_flow =%d , k=%d\n",find_flow, k);
                }
                
                
            }
        }
        


    }



// -------------------------------------------------------------



    bool cmp(Edge a, Edge b){
        return a.cap > b.cap;
    }



//当前方案的最大流 ,hua为实参
    int maxFlow(struct MCMF& hua){
        
        hua.outCnt = 0;
        int res = hua.Mincost(hua.s, hua.t, hua.k);
        return res;
    }



//-------求每个网络点到消费点的最小花费

void  minDisofXY(struct MCMF hua ,int i ,int j, int cmin[1005][505]){ //不修改 实参数的值
    
    
    //源点与服务点连
    hua.AddEdge(hua.s, hua.edges[hua.consumeEdge[i]].to , INF, 0);
    //消费点与汇点相连
    int Ccap = hua.edges[hua.consumeEdge[j]].cap;
    hua.AddEdge(j + hua.netNode, hua.t, Ccap, 0);
    
    //int tmp = hua.Mincost( i, j+hua.netNode, e.cap);
    int tmp = hua.Mincost( hua.s, hua.t, Ccap);
    //tmp = !tmp?INF:tmp; tmp为-1表示不可到达或不能满足流量
    printf("#%d %d : %d\n",i , j ,tmp);
    
    cmin[i][j] = tmp;
    
}



//获得服务点 并返回最小费用 ,hua为形参 ,huaReal为实参
int getServersAndGetMincost(vector<int>& serverList, int& serverNum ,struct MCMF hua ,struct MCMF& huaReal , string& str,  int printPath){
    int len,ret;
    
    //int servers[]={2};
    //int servers[]={7,13,15,22,37,38,43}; //case0
    //int servers[]={6,7,13,17,35,41,48};  //case1
    //int servers[]={10,22,26,29,35};  //case3
    //int servers[]={0,1,24};  //demo
    //memcpy(serverList, servers, sizeof(servers));
    
    // len=3;
    //服务点设置为 与消费点最近的网络点
    /*    len = hua.consumeNode;
     for(int i=0; i< len; i++){
     int consume =i + hua.netNode;
     Edge e = hua.edges[hua.G[consume][0]];
     serverList[i] = e.to ;
     }
     */
    len = (int)serverList.size();
    
    hua.numOfserver=len;
    //建立 源点 与 服务点的连线
    for(int i=0; i< len ; i++)
    {
        //源点与服务点连
        hua.AddEdge(hua.s, serverList[i], INF, 0);
        //消费点与汇点相连
     //   int Ccap = hua.edges[hua.consumeEdge[serverList[i]]].cap;
     //   hua.AddEdge(serverList[i] + hua.netNode, hua.t, Ccap, 0);
    }
    
    
    
    //清空之前数据
    huaReal.outFlow.clear();
    for(int i=0; i<hua.outCnt; i++)
        hua.res[i].clear();
    
    
    //对当前结果 判优，并进行保存
    ret = maxFlow(hua);

    if(printPath)
    {
        hua.search_path();
        hua.OutputPath(str);
    }
    
    return ret;
}



//添加了消费点列表
//获得服务点 并返回最小费用 ,hua为形参 ,huaReal为实参
int newGetServersAndGetMincost(vector<int>& serverList,vector<int> consumeList,int& serverNum,MCMF hua,MCMF& huaReal)
{
    int len,ret;
    len = serverNum;
    
    hua.numOfserver= (int)serverList.size();
    //建立 源点 与 服务点的连线
    for(int i=0; i< serverList.size() ; i++)
    {
        //源点与服务点连
        hua.AddEdge(hua.s, serverList[i], INF, 0);
        
    }
    
    //记录t点的总流量
    int tot_k=0;
    //消费点与汇点相连
    for(int i=0; i<consumeList.size(); i++){
        int Ccap = hua.edges[hua.consumeEdge[consumeList[i]]].cap;
        tot_k +=Ccap;
      //  printf("cL: %d , %d\n",consumeList[i], Ccap);
        hua.AddEdge(consumeList[i] + hua.netNode, hua.t, Ccap, 0);
    }
    hua.k=tot_k;
    
    //清空之前数据
    huaReal.outFlow.clear();
    for(int i=0; i<hua.outCnt; i++)
        hua.res[i].clear();
    
    
    //对当前结果 判优，并进行保存
    ret = maxFlow(hua);
    
    
    if(ret>0){ //正确解
        if(huaReal.costNow<0 || huaReal.costNow>ret){
            // printf("^^^^^^ huaReal.costNow =%d\n",huaReal.costNow);
            huaReal.outFlow.assign(hua.outFlow.begin(), hua.outFlow.end());
            for(int i=0; i<hua.outCnt; i++)
                huaReal.res[i].assign(hua.res[i].begin(), hua.res[i].end());
            huaReal.outCnt = hua.outCnt;
            huaReal.costNow = hua.costNow;
        }
    }
    
    //  huaReal.OutputPath(str);
    
    return ret;
}



int newGetServersAndGetMincost2(vector<int>& serverList,vector<int> consumeList,int& serverNum,MCMF& hua , string& str)
{
    int len,ret;
    len = serverList.size();
    
    hua.numOfserver= (int)serverList.size();
    //建立 源点 与 服务点的连线
    for(int i=0; i< serverList.size() ; i++)
    {
        //源点与服务点连
        hua.AddEdge(hua.s, serverList[i], INF, 0);
        
    }
    
    //记录t点的总流量
    int tot_k=0;
    //消费点与汇点相连
    for(int i=0; i<consumeList.size(); i++){
        int Ccap = hua.edges[hua.consumeEdge[consumeList[i]]].cap;
        tot_k +=Ccap;
        //  printf("cL: %d , %d\n",consumeList[i], Ccap);
        hua.AddEdge(consumeList[i] + hua.netNode, hua.t, Ccap, 0);
    }
    hua.k=tot_k;
    
    //清空之前数据
    for(int i=0; i<hua.outCnt; i++)
        hua.res[i].clear();
    
    
    //对当前结果 判优，并进行保存
    ret = maxFlow(hua);
 /*   if(ret>0)
    {
        hua.search_path();
        hua.OutputPath(str);
    }else{
        for(int i=0; i<hua.consumeNode; i++)
        printf("%d : needflow= %d  realf= %d\n",i ,hua.conflowNeed[i], hua.conflow[i]);
    }
  
  */
    return ret;
}








//--------------- 纯粹直连方式
int DirectLink(vector<int>& serverList, int& serverNum ,struct MCMF hua ,struct MCMF& huaReal , string& str){
    int len,ret;
    
    //int servers[]={2};
    //int servers[]={7,13,15,22,37,38,43}; //case0
    //int servers[]={6,7,13,17,35,41,48};  //case1
    //int servers[]={10,22,26,29,35};  //case3
    //int servers[]={0,1,24};  //demo
    //memcpy(serverList, servers, sizeof(servers));
    
    // len=3;
    //服务点设置为 与消费点最近的网络点
     len = hua.consumeNode;
     for(int i=0; i< len; i++){
     int consume =i + hua.netNode;
     Edge e = hua.edges[hua.G[consume][0]];
     serverList.push_back(e.to) ;
     }
    
    
    hua.numOfserver=len;
    //建立 源点 与 服务点的连线
    for(int i=0; i< len ; i++)
    {
        hua.AddEdge(hua.s, serverList[i], INF, 0);
    }
    
    int tot_k =0;
    //消费点与汇点相连
    for(int i=0; i<hua.consumeNode; i++){
        int Ccap = hua.edges[hua.consumeEdge[i]].cap;
        tot_k +=Ccap;
        //  printf("cL: %d , %d\n",consumeList[i], Ccap);
        hua.AddEdge(i + hua.netNode, hua.t, Ccap, 0);
    }
    hua.k=tot_k;
    
    //清空之前数据
    huaReal.outFlow.clear();
    for(int i=0; i<hua.outCnt; i++)
        hua.res[i].clear();
    
    //对当前结果 判优，并进行保存
    ret = maxFlow(hua);
    /*
    if(ret>0){ //正确解
        if(huaReal.costNow<0 || huaReal.costNow>ret){
             // printf("^^^^^^ huaReal.costNow =%d\n",huaReal.costNow);
            huaReal.outFlow.assign(hua.outFlow.begin(), hua.outFlow.end());
            for(int i=0; i<hua.outCnt; i++)
                huaReal.res[i].assign(hua.res[i].begin(), hua.res[i].end());
            huaReal.outCnt = hua.outCnt;
            huaReal.costNow = hua.costNow;
        }
    }
    */
    //huaReal.OutputPath(str);
    
    
    printf("纯直连方式，minCost =  %d\n\n",ret);
    return ret;
}






