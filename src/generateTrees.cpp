#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <stdexcept>
#include <functional>
#include <utility>
#include <queue>

using namespace std;

// Tasks
// 2. Run Dijkstra to generate SPT and write SPT to file
// 3. Run Kruskal to generate MST and write MST to file

struct Edge {
  int destination;
  double distance;
  double cost;
};

struct SPT {
  vector<int> parent;
  vector<double> distance;
};

// Type alias for readability
using Graph = vector<vector<Edge>>;

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







// Is this needed?
bool compareCosts(const Edge& e1, const Edge& e2){
  return e1.cost < e2.cost;
}





int main(int argc, char** argv){
  try {
    int plantID;
    Graph G = readGraphFromFile(argv[1], plantID);

    SPT tree = generateShortestPathTree(G, plantID);
    writeSPTtoFile(tree, "shortestPathTree.txt");
  }
  catch (const exception& e){
    cerr << "Error: " << e.what() << endl;
    return 1;
  }


  return 0;
}
