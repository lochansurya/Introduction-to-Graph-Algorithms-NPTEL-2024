#include"graph.h"
#include<iostream>
#include<vector>
#include<limits>
#include<tuple>
#include<algorithm>
#include<utility>
#include<set>
#include<numeric> //for iota

using namespace std;
void printGraph(const graph &G)
{
    cout << endl;
    for (auto edge : G.E)
    {
        cout << get<0>(edge) << " " << get<1>(edge) << " " << get<2>(edge) << "\n";
    }
    cout << endl;
}

void printMatrix(vector<vector<int>> matrix){
    cout << endl;
    for (auto row : matrix)
    {
        for(auto value: row){
            if(value==INF){
                cout<<"INF ";
            }
            else{
                cout<<value<<" ";
            }
        }
        cout<<endl;
    }
    cout << endl;
}

void graph::addEdge(int u, int v, int weight)
{
    E.push_back(make_tuple(u, v, weight));
}

void graph::makeUndirected()
{
    vector<tuple<int, int, int>> undirectedEdges;
    set<pair<int, int>> addedEdges;

    for (auto edge : E)
    {
        int u = get<0>(edge);
        int v = get<1>(edge);
        int weight = get<2>(edge);

        if (addedEdges.find({u, v}) == addedEdges.end() && addedEdges.find({v, u}) == addedEdges.end())
        {
            undirectedEdges.push_back(edge);
            addedEdges.insert({u, v});
            undirectedEdges.push_back(make_tuple(v, u, weight));
        }
    }

    E = undirectedEdges;
}

void graph::printAdjList(int node){
    cout << endl;
    for (auto edge : E)
    {
        if(get<0>(edge)==node){
            cout<<get<1>(edge)<<" "<<get<2>(edge)<<endl;
        }
    }
    cout << endl;
}

int graph::find(vector<int> &parent, int i)
{
    if (parent[i] != i)
    {
        parent[i] = find(parent, parent[i]); // Path compression
    }
    return parent[i];
}

vector<vector<int>> graph::getWeightMatrix()
{
    vector<vector<int>> weightMatrix(V, vector<int>(V, INF));
    for (auto edge : E)
    {
        int u = get<0>(edge);
        int v = get<1>(edge);
        int weight = get<2>(edge);
        weightMatrix[u][v] = weight;
    }

    // Set diagonal entries to 0
    for (int i = 0; i < V; ++i)
    {
        weightMatrix[i][i] = 0;
    }

    return weightMatrix;
}

void graph::dijkstraSSSP(int src, vector<int> &dist){
    dist[src] = 0;
    vector<bool> visited(V, false);
    for(int i=0; i<V-1; i++){
        int u = -1;
        for(int j=0; j<V; j++){
            if(!visited[j] && (u==-1 || dist[j]<dist[u])){
                u = j;
            }
        }
        visited[u] = true;
        for(auto edge: E){
            if(get<0>(edge)==u){
                int v = get<1>(edge);
                int weight = get<2>(edge);
                if(dist[u]!=INF && dist[u]+weight<dist[v]){
                    dist[v] = dist[u]+weight;
                }
            }
        }
    }
}

void graph::bellmanFordSSSP(int src, vector<int> &dist){
    dist[src] = 0;
    for(int i=1; i<V; i++){
        for(auto edge: E){
            int u = get<0>(edge);
            int v = get<1>(edge);
            int weight = get<2>(edge);
            if(dist[u]!=INF && dist[u]+weight<dist[v]){
                dist[v] = dist[u]+weight;
            }
        }
    }
    for(auto edge: E){
        int u = get<0>(edge);
        int v = get<1>(edge);
        int weight = get<2>(edge);
        if(dist[u]!=INF && dist[u]+weight<dist[v]){
            cout<<"Graph contains negative weight cycle"<<endl;
            return;
        }
    }
}

void graph::bellmanFordAPSP(vector<vector<int>>& distMatrix){
    for(int i=0; i<V; i++){
        distMatrix[i][i] = 0;
    }
    for(int i=1; i<V; i++){
        for(auto edge: E){
            int u = get<0>(edge);
            int v = get<1>(edge);
            int weight = get<2>(edge);
            if(distMatrix[u][i-1]!=INF && distMatrix[u][i-1]+weight<distMatrix[v][i]){
                distMatrix[v][i] = distMatrix[u][i-1]+weight;
            }
        }
    }
    for(auto edge: E){
        int u = get<0>(edge);
        int v = get<1>(edge);
        int weight = get<2>(edge);
        if(distMatrix[u][V-1]!=INF && distMatrix[u][V-1]+weight<distMatrix[v][V]){
            cout<<"Graph contains negative weight cycle"<<endl;
            return;
        }
    }
}

void graph::floydWarshallAPSP(vector<vector<int>> &dist){
    for(int i=0; i<V; i++){
        for(int j=0; j<V; j++){
            dist[i][j] = INF;
        }
    }
    for(auto edge: E){
        int u = get<0>(edge);
        int v = get<1>(edge);
        int weight = get<2>(edge);
        dist[u][v] = weight;
    }
    for(int i=0; i<V; i++){
        dist[i][i] = 0;
    }
    for(int k=0; k<V; k++){
        for(int i=0; i<V; i++){
            for(int j=0; j<V; j++){
                if(dist[i][k]!=INF && dist[k][j]!=INF && dist[i][k]+dist[k][j]<dist[i][j]){
                    dist[i][j] = dist[i][k]+dist[k][j];
                }
            }
        }
    }
}

void graph::johnsonAPSP(vector<vector<int>> &dist)
{
    vector<int> h(V + 1, INF);
    h[V] = 0; // The additional vertex will have a distance of 0 from itself.

    // Step 1: Add an extra vertex with index V and edges of weight 0 from this vertex to all other vertices.
    for (int i = 0; i < V; ++i)
    {
        E.push_back(make_tuple(V, i, 0));
    }

    // Step 2: Run Bellman-Ford from the new vertex V
    bellmanFordSSSP(V, h);

    // Step 3: Reweight all the edges in the original graph
    for (auto &edge : E)
    {
        int u = get<0>(edge);
        int v = get<1>(edge);
        int weight = get<2>(edge);
        if (u != V)
        { // Skip edges connected to the added vertex
            edge = make_tuple(u, v, weight + h[u] - h[v]);
        }
    }

    // Remove the last added vertex
    E.erase(E.end() - V, E.end());

    // Step 4: Construct and print the reweighted matrix
    cout << "Reweighted matrix:" << endl;
    vector<vector<int>> reweightedMatrix = getWeightMatrix();
    printMatrix(reweightedMatrix);

    // Step 5: Use Dijkstra's algorithm for each vertex on the reweighted graph
    for (int i = 0; i < V; ++i)
    {
        vector<int> dijkstraDist(V, INF);
        dijkstraSSSP(i, dijkstraDist);

        for (int j = 0; j < V; ++j)
        {
            if (dijkstraDist[j] != INF)
            {
                dist[i][j] = dijkstraDist[j] + h[j] - h[i];
            }
        }
    }

    // Step 6: Set diagonal entries to 0
    for (int i = 0; i < V; ++i)
    {
        dist[i][i] = 0;
    }

    // Step 7: Restore the original graph weights
    for (auto &edge : E)
    {
        int u = get<0>(edge);
        int v = get<1>(edge);
        int weight = get<2>(edge);
        edge = make_tuple(u, v, weight - h[u] + h[v]);
    }

    // print the weight matrix for confirmation
    cout << "Weight matrix after Johnson's algorithm:" << endl;
    vector<vector<int>> weightMatrix = getWeightMatrix();
    printMatrix(weightMatrix);

}

vector<int> graph::primMST(int source){
    vector<int> parent(V, -1);
    vector<int> key(V, INF);
    vector<bool> inMST(V, false);
    key[source] = 0;
    for(int i=0; i<V-1; i++){
        int u = -1;
        for(int j=0; j<V; j++){
            if(!inMST[j] && (u==-1 || key[j]<key[u])){
                u = j;
            }
        }
        inMST[u] = true;
        for(auto edge: E){
            if(get<0>(edge)==u){
                int v = get<1>(edge);
                int weight = get<2>(edge);
                if(!inMST[v] && weight<key[v]){
                    key[v] = weight;
                    parent[v] = u;
                }
            }
        }
    }
    cout << "Parent Array: ";
    for (int i = 0; i < V; i++)
    {
        if(parent[i]!=-1){
            cout<<parent[i]<<" "<<i<<endl;
        }
    }
    return parent;
}

void graph::kruskalMST()
{
    vector<int> parent(V);
    iota(parent.begin(), parent.end(), 0); // Initialize parent

    vector<pair<int, pair<int, int>>> edges;
    for (auto edge : E)
    {
        edges.push_back({get<2>(edge), {get<0>(edge), get<1>(edge)}});
    }
    sort(edges.begin(), edges.end());

    vector<int> mstEdges;
    for (auto edge : edges)
    {
        int u = edge.second.first;
        int v = edge.second.second;
        int weight = edge.first;

        int setU = find(parent, u);
        int setV = find(parent, v);
        if (setU != setV)
        {
            parent[setU] = setV;
            mstEdges.push_back(u);
            mstEdges.push_back(v);
            cout << u << " " << v << endl;
        }
    }
}

int main()
{
    graph G(5);
    G.addEdge(0, 1, 3);
    G.addEdge(0, 2, 8);
    G.addEdge(0, 4, -4);
    G.addEdge(1, 3, 1);
    G.addEdge(1, 4, 2);
    G.addEdge(2, 1, 4);
    G.addEdge(3, 0, 2);
    G.addEdge(3, 2, -5);
    G.addEdge(4, 3, 6);

    // Print the original graph
    cout << "Original Graph:" << endl;
    printGraph(G);

    // Print the weight matrix
    vector<vector<int>> weightMatrix = G.getWeightMatrix();
    cout << "Weight matrix:" << endl;
    printMatrix(weightMatrix);

    // Run and print results of Dijkstra’s algorithm from node 0
    vector<int> dijkstraDist(5, INF);
    G.dijkstraSSSP(0, dijkstraDist);
    cout << "Shortest paths from node 0 using Dijkstra's algorithm:" << endl;
    for (int i = 0; i < dijkstraDist.size(); ++i)
    {
        cout << "Distance to node " << i << ": " << (dijkstraDist[i] == INF ? "INF" : to_string(dijkstraDist[i])) << endl;
    }

    // Run and print results of Bellman-Ford’s algorithm from node 0
    vector<int> bellmanFordDist(5, INF);
    G.bellmanFordSSSP(0, bellmanFordDist);
    cout << "Shortest paths from node 0 using Bellman-Ford algorithm:" << endl;
    for (int i = 0; i < bellmanFordDist.size(); ++i)
    {
        cout << "Distance to node " << i << ": " << (bellmanFordDist[i] == INF ? "INF" : to_string(bellmanFordDist[i])) << endl;
    }

    // Run and print results of Floyd-Warshall’s algorithm
    vector<vector<int>> floydWarshallDist(5, vector<int>(5, INF));
    G.floydWarshallAPSP(floydWarshallDist);
    cout << "All-pairs shortest paths using Floyd-Warshall algorithm:" << endl;
    printMatrix(floydWarshallDist);

    // Run and print results of Johnson’s algorithm
    vector<vector<int>> johnsonDist(5, vector<int>(5, INF));
    G.johnsonAPSP(johnsonDist);
    cout << "All-pairs shortest paths using Johnson's algorithm:" << endl;
    printMatrix(johnsonDist);

    // Make the graph undirected before calling MST algorithms
    G.makeUndirected();

    // Run and print results of Prim’s MST algorithm from node 0
    cout << "Minimum Spanning Tree (MST) edges using Prim’s algorithm:" << endl;
    vector<int> parent = G.primMST(0);

    // print the parent array
    cout << "Parent array after calling primMST: ";
    for (int i = 0; i < parent.size(); i++)
    {
        cout << parent[i] << " ";
    }
    cout << endl;

    // Run and print results of Kruskal’s MST algorithm
    cout << "Minimum Spanning Tree (MST) edges using Kruskal’s algorithm:" << endl;
    G.kruskalMST();

    
    return 0;
}
