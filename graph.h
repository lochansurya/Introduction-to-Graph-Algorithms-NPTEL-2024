#ifndef GRAPH_H
#define GRAPH_H
#include<vector>
#include<utility>
#include<limits>

using namespace std;

const int INF = numeric_limits<int>::max();

class graph{
    private:
        int V;
        vector<tuple<int, int, int>> E;
    public:
        graph(int v):V(v){}
        void addEdge(int u, int v, int weight);
        void makeUndirected();
        void printAdjList(int node);
        vector<vector<int>> getAdjMatrix();
        vector<vector<int>> getWeightMatrix();
        void dijkstraSSSP(int src, vector<int> &dist);
        void bellmanFordSSSP(int src, vector<int> &dist);
        void bellmanFordAPSP(vector<vector<int>> &dist);
        void floydWarshallAPSP(vector<vector<int>> &dist);
        void johnsonAPSP(vector<vector<int>> &dist);
        int find(vector<int> &parent, int i);
        vector<int> primMST(int source);
        void kruskalMST();
        void printMST();
        friend void printGraph(const graph &G);
};

#endif