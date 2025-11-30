#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <stdexcept>
#include <functional>
#include <utility>
#include <queue>
#include <algorithm>

using namespace std;

// Tasks
// 2. Run Dijkstra to generate SPT and write SPT to file
// 3. Run Kruskal to generate MST and write MST to file

struct Edge {
  int destination;
  double distance;
  double cost;
};

struct MSTEdge {
  int source;
  int destination;
  double cost;
};

struct SPT {
  vector<int> parent;
  vector<double> distance;
};

// Type alias for readability
using Graph = vector<vector<Edge>>;
using EdgeList = vector<MSTEdge>;

/**
 * Reads power network graph from text file and builds an adjacency list.
 *
 * The input file format is:
 *   N plantID
 *   u v distance cost
 *   u v distance cost
 *   ...
 *
 * The graph is treated as undirected: for each line (u, v, ...),
 * edges are added in both directions in the adjacency list.
 *
 * @param fileName  Path to the input graph file.
 * @param plantID   Output parameter that will receive the plant node ID.
 * @return          Adjacency list representation of the graph.
 *
 * @throws std::runtime_error if the file cannot be opened or the format is invalid.
 */
Graph readGraphFromFile(const string& fileName, int& plantID){
  ifstream file(fileName);
  if (!file) {
    throw runtime_error("Could not open input file: " + fileName);
  }

  int N;
  if (!(file >> N >> plantID)) {
    throw runtime_error("Failed to read N and plantID from file: " + fileName);
  }

  Graph G(N);

  int s1, s2;
  double dist, cost;

  while (file >> s1 >> s2 >> dist >> cost){
    if (s1 < 0 || s1 >= N || s2 < 0 || s2 >= N) {
      throw runtime_error("Node index out of range in file: " + fileName);
    }

    Edge e1;
    e1.destination = s2;
    e1.distance = dist;
    e1.cost = cost;
    G[s1].push_back(e1);

    Edge e2;
    e2.destination = s1;
    e2.distance = dist;
    e2.cost = cost;
    G[s2].push_back(e2);
  }

  file.close();

  return G;
}

/**
 * Computes the shortest-path tree (SPT) from the plant using Dijkstra's algorithm.
 *
 * The SPT is computed using Euclidean distance stored in each Edge. This produces
 * the minimum-distance routing tree from the generation plant to all substations.
 * A min-priority queue (distance, node) pair structure is used, and stale queue
 * entries are ignored using a visited[] guard.
 *
 * @param G       Undirected weighted adjacency list representing the network.
 * @param plant   Node ID of the generation plant (SPT root).
 *
 * @return A pair (parent, dist):
 *         - parent[i] gives the parent of node i in the shortest-path tree
 *           (parent[plant] = -1).
 *         - dist[i] gives the shortest-path distance from plant to node i.
 *
 * The caller is responsible for writing the resulting tree to a file.
 */
SPT generateShortestPathTree(const Graph& G, int plant){
  int N = G.size();
  vector<int> parent(N, -1);
  vector<double> dist(N, numeric_limits<double>::infinity());
  vector<bool> visited(N, false);
  // Each entry in the priority queue is of form (distance, destination)
  priority_queue<
    pair<double, int>,
    vector<pair<double, int> >,
    greater<pair<double, int> >
  > pq;

  dist[plant] = 0;
  pq.push({0.0, plant});

  while(!pq.empty()){
    pair<double, int> top = pq.top();
    int u = top.second;
    double d = dist[u];
    pq.pop();
    if (visited[u]){
	continue;
      }
    
    for (Edge e : G[u]){
      int v = e.destination;
      double newDist = d + e.distance;
      if (dist[v] > newDist){
	parent[v] = u;
	dist[v] = newDist;
	pq.push({newDist, v});
      }
    }
    visited[u] = true;
  }

  SPT tree;
  tree.parent = parent;
  tree.distance = dist;
  return tree;
}

void writeSPTtoFile(const SPT& tree, const string& fileName){
  ofstream file(fileName);

  if (!file) {
    throw runtime_error("Could not open file for writing SPT: " + fileName);
  }

  int N = tree.parent.size();

  for (int i = 0; i < N; i++){
    file << i << " " << tree.parent[i] << " " << tree.distance[i] << "\n";
  }

  file.close();
}


/**
 * @brief Generates a deduplicated and cost-sorted edge list from an adjacency-list graph.
 *
 * This function converts an undirected adjacency-list representation into a flat list
 * of MSTEdge records. Because the graph is undirected and each edge appears twice
 * (u→v and v→u), only edges with u < v are retained to avoid duplicates.
 *
 * The returned list is sorted in non-decreasing order of edge cost, making it suitable
 * for use in Kruskal's minimum spanning tree algorithm.
 *
 * @param G   The adjacency-list representation of the graph. G[u] contains all edges
 *            outgoing from u, each with a destination and cost.
 *
 * @return A vector of MSTEdge objects, each representing an undirected edge (u, v, cost),
 *         sorted by increasing cost.
 */
EdgeList generateSortedEdgeList(const Graph& G){
  int N = G.size();
  EdgeList edgeList;
  for (int u = 0; u < N; u++){
     for (const Edge& e : G[u]){
       int v = e.destination;
       if (u < v){
	 MSTEdge newEdge;
	 newEdge.source = u;
	 newEdge.destination = e.destination;
	 newEdge.cost = e.cost;
	 edgeList.push_back(newEdge);
       }
    }
  }

  // Sort edges in increasing order by cost
  sort(edgeList.begin(), edgeList.end(),
         [](const MSTEdge& a, const MSTEdge& b) {
             return a.cost < b.cost;
         });

  return edgeList;
}

/**
 * @brief Merges two disjoint-set components using union by rank.
 *
 * This function assumes that @p u and @p v are the root representatives of their
 * respective sets as returned by Find(). The parent array is updated so that one
 * root becomes the parent of the other, preserving the union-by-rank heuristic to 
 * maintain near-constant amortized time complexity.
 *
 * @param u      Root of the first set.
 * @param v      Root of the second set.
 * @param rank   Rank array tracking approximate tree heights for each root.
 * @param parent Parent array representing the disjoint-set forest structure.
 */
void Union(int u, int v, vector<int>& rank, vector<int>& parent){
  if (rank[u] == rank[v]){
    parent[v] = u;
    rank[u]++;
  }
  else if (rank[u] > rank[v]){
    parent[v] = u;
  }
  else {
    parent[u] = v;
  }
}


/**
 * @brief Finds the root of the set containing a given element, with path compression.
 *
 * This function recursively finds the root of the element @p u in the disjoint-set
 * forest. Path compression is applied: each node visited on the find path is directly
 * linked to the root, dramatically improving performance for future queries.
 *
 * @param u       The element whose representative (root) is to be found.
 * @param parent  Parent array representing the disjoint-set forest. This array is
 *                mutated to perform path compression.
 *
 * @return The root representative of the set containing @p u.
 */
int Find(int u, vector<int>& parent){
  if (parent[u] != u){
    parent[u] = Find(parent[u], parent);
  }
  return parent[u];
}

/**
 * @brief Computes the minimum spanning tree (MST) using Kruskal's algorithm.
 *
 * This function assumes the edge list @p E is sorted in non-decreasing order of cost.
 * It processes edges in that order, using union–find with union-by-rank and path 
 * compression to efficiently test whether an edge connects two previously disjoint 
 * components. Edges that do not create cycles are added to the MST.
 *
 * @param E  A sorted list of undirected edges (source, destination, cost).
 * @param N  The total number of nodes in the graph (labeled 0 to N-1).
 *
 * @return A vector of MSTEdge objects representing the minimum spanning tree.
 *         The MST will contain exactly N-1 edges if the graph is connected.
 */
EdgeList generateMinimumCostTree(const EdgeList& E, int N){
  EdgeList MST;
  vector<int> parent(N);
  vector<int> rank(N, 0);
  for (int i = 0; i < N; i++){
    parent[i] = i;
  }

  for (const MSTEdge& e : E){
    int rootU = Find(e.source, parent);
    int rootV = Find(e.destination, parent);
    if (rootU == rootV){
      // Cycle detection heuristic
      continue;
    }
    MST.push_back(e);
    Union(rootU, rootV, rank, parent);
  }
  return MST;
}


/**
 * @brief Writes the minimum spanning tree edge list to a text file.
 *
 * The output format is one edge per line:
 *     source destination cost
 *
 * @param MST       The list of MST edges generated by Kruskal's algorithm.
 * @param fileName  Path to the output text file.
 *
 * @throws std::runtime_error if the file cannot be opened.
 */
void writeMSTtoFile(const EdgeList& MST, const string& fileName){
  ofstream file(fileName);

  if (!file) {
    throw runtime_error("Could not open file for writing MST: " + fileName);
  }


  for (const MSTEdge& e : MST){
    file << e.source << " " << e.destination << " " << e.cost << "\n";
  }

  file.close();
}


int main(int argc, char** argv){
  try {
    int plantID;
    Graph G = readGraphFromFile(argv[1], plantID);

    SPT tree = generateShortestPathTree(G, plantID);
    writeSPTtoFile(tree, "shortestPathTree.txt");

    EdgeList sortedEdgeList = generateSortedEdgeList(G);
    EdgeList MST = generateMinimumCostTree(sortedEdgeList, G.size());
    writeMSTtoFile(MST, "minimumCostTree.txt");
    
  }
  catch (const exception& e){
    cerr << "Error: " << e.what() << endl;
    return 1;
  }


  return 0;
}
