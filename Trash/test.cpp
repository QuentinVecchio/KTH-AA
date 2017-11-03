#include <iostream>
#include <vector>
#include <cmath>
#include <time.h>
#include <deque>
#include <climits>
#include <algorithm>
#include <ctime>

#define TIME_LIMIT 1.8

using namespace std;

typedef struct {
    int id;
    double x;
    double y;
} Node;

typedef struct {
    int u;
    int v;
    int d;
} Edge;

// Generals Methods
void printSolution(const vector<Node> &solutionInit);
int euclidianDist(double x1, double y1, double x2, double y2);
int euclidianDistBetweenNodes(Node n1, Node n2);
int fitness(const vector<Node> &solution);
double random0and1();

// Kruskal Method
int find(vector<int> parent, int i);
vector<int> union(vector<int> parent, int x, int y);
vector<Edge> getEdgesSorted(vector<vector<int> > matrix);
vector<Node> getMinimumSpanningTree(vector<Node> G, vector<vector<int> > matrix);


int main() {
    auto begin = clock();

    // Variables initialisation
    int N;
    srand (time(NULL));
    vector<Node> nodes;
    vector<Node> optimaleSolution;
    vector<vector<int> > distanceMatrix;

    // Get number of nodes
    cin >> N;

    // Get the N nodes
    for(int i=0;i<N;i++) {
        double x, y;
        cin >> x >> y;
        Node n = {i,x,y};
        nodes.push_back(n);
    }

    // Create matrix distance
    distanceMatrix = createDistanceMatrix(nodes);
    


    // Print solution
    printSolution(optimaleSolution);

    return 0;
}

vector<Node> metricTSP(vector<Node> G, vector<vector<int> > matrix) {
    vector<Node> minimumSpanningTree = getMinimumSpanningTree(G, matrix);
    
}

vector<Node> getMinimumSpanningTree(vector<Node> G, vector<vector<int> > matrix) {
    vector<Edge> edgesSorted = getEdgesSorted(matrix);
    vector<Edge> subset;
    vector<int> parents = vector<int>(matrix.size(), -1);
    for (int v=0;v<matrix.size();++v) {
        parents[v] = v;
    }

    for(int i=0;i<edgesSorted.size();i++) {
        int x = find(parents);
        int y = find(parents);

        if(x != y) {
            subset.push_back(edgesSorted[i]);
            parents = union(parents, x, y)
        }

        if(subset.size() == G.size() - 1) {
            break;
        }
    }
    return subset;
}

vector<Edge> getEdgesSorted(vector<vector<int> > matrix) {
    vector<Edge> edges;

    for(int i=0;i<matrix.send();i++) {
        for(int j=0;j<matrix[i].send();j++) {
            edges.push_back({i, j, matrix[i][j]});
        }
    }

    sort(edges.begin( ), edges.end( ), [ ]( const Edge& lhs, const Edge& rhs ) {
       return lhs.d < rhs.d;
    });

    return edges;
}

vector<int> union(vector<int> parent, int x, int y) {
    int xset = find(parent, x);
    int yset = find(parent, y);
    parent[xset] = yset;
    return parent;
}

int find(vector<int> parent, int i) {
    if (parent[i] == -1)
        return i;
    return find(parent, parent[i]);
}

int fitness(const vector<Node> &solution) {
    int f = 0;
    for(int i=0;i<solution.size()-1;++i) {
        f += euclidianDistBetweenNodes(solution[i], solution[i+1]);
    }
    if(solution.size()-1 > 1) {
        f += euclidianDistBetweenNodes(solution[solution.size()-1], solution[0]);
    }
    return f;
}

void printSolution(const vector<Node> &solutionInit) {
    for(int i=0;i<solutionInit.size();++i) {
        cout << solutionInit[i].id << endl;
    }
}

vector<vector<int> > createDistanceMatrix(const vector<Node> & solution) {
    vector<vector<int> > matrix = vector<vector<int> >(solution.size());
    for(int i=0;i<solution.size();++i) {
        for(int j=0;j<solution.size();++j) {
            if(i != j) {
                matrix[i].push_back(euclidianDistBetweenNodes(solution[i], solution[j]));
            } else {
                matrix[i].push_back(INT_MAX);
            }
        }
    }
    return matrix;
}

double random0and1() {
    return static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
}

int euclidianDistBetweenNodes(Node n1, Node n2) {
    return euclidianDist(n1.x, n1.y, n2.x, n2.y);
}

int euclidianDist(double x1, double y1, double x2, double y2) {
    return round(sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2)));
}