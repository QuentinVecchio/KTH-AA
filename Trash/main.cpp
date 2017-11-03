#include <iostream>
#include <vector>
#include <cmath>
#include <time.h>
#include <deque>
#include <climits>
#include <algorithm>

using namespace std;

typedef struct {
    int id;
    double x;
    double y;
} Node;

void printSolution(const vector<Node> &solutionInit);
int euclidianDist(double x1, double y1, double x2, double y2);
int euclidianDistBetweenNodes(Node n1, Node n2);
vector<Node> greedyTour(const vector<Node> &solutionInit, int start, vector<vector<int> > matrix);
vector<Node> optimalGreedyTour(const vector<Node> &solutionInit, vector<vector<int> > matrix);
int fitness(const vector<Node> &solution);
vector<Node> twoOpt(const vector<Node> & solution, vector<vector<int> > matrix);
vector<vector<int> > createDistanceMatrix(const vector<Node> & solution);
vector<Node> reverseGraph(const vector<Node> & solution, int n1, int n2);

int main() {
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
    
    // Call greedy tour
    optimaleSolution = greedyTour(nodes, 0, distanceMatrix);
    //printSolution(optimaleSolution);
    optimaleSolution = twoOpt(optimaleSolution, distanceMatrix);
    //cout << fitness(optimaleSolution) << endl;
    //optimaleSolution = optimalGreedyTour(nodes);

    // Print solution
    printSolution(optimaleSolution);

    return 0;
}

vector<Node> optimalGreedyTour(const vector<Node> &solutionInit, vector<vector<int> > matrix) {
    vector<Node> best;
    int bestFitness = INT_MAX;
    int borne = 50;
    vector<bool> alreadyVisited = vector<bool>(solutionInit.size(), false);
    if(borne >= solutionInit.size()) {
        borne = solutionInit.size();
    }
    for(int i=0;i<borne;i++) {
        int index = i;
        if(borne != solutionInit.size()) {
            do {
                index = rand() % solutionInit.size();
            } while(alreadyVisited[index]);
            alreadyVisited[index] = true;
        }
        vector<Node> solution = greedyTour(solutionInit, index, matrix);
        int fit = fitness(solution);
        if(fit < bestFitness) {
            bestFitness = fit;
            best = solution;
        }
    }
    return best;
}

vector<Node> greedyTour(const vector<Node> &solutionInit, int start, vector<vector<int> > matrix) {
    vector<Node> solution = vector<Node>(solutionInit.size());
    vector<bool> used = vector<bool>(solutionInit.size(), false);
    solution[start] = solutionInit[start];
    used[start] = true;
    for(int i=0;i<solutionInit.size();i++) {
        if(i != start) {
            int best = -1;
            for(int j=0;j<solutionInit.size();j++) {
                if(!used[j] && (best == -1 || matrix[solution[i-1].id][solutionInit[j].id] < matrix[solution[i-1].id][solutionInit[best].id])) {
                    best = j;
                }
            }
            solution[i] = solutionInit[best];
            used[best] = true;
        }
    }
   return solution;
}

/*vector<Node> twoOpt(const vector<Node> & solution, vector<vector<int> > matrix) {
    vector<Node> solution = vector<Node>(solution.size());
    vector<bool> used = vector<bool>(solution.size(), false);

    bool localOptimal = false;
    while(!localOptimal) {
        localOptimal = true;
        for (u_i = 0, v_i = 1; u_i < solution.size(); ++u_i, ++v_i) {
            u = solution[u_i].id;
            v = tour[v_i % solution.size()].id;

            for (size_t k = 0; k < solution.size(); ++k) {
                w_i = position[neighbor[u][k]];
            }
        }
    }


    for (u_i = 0, v_i = 1; u_i < N; ++u_i, ++v_i) {
            u = tour[u_i];
            v = tour[v_i % N];

            // For each edge wz (w k:th closest neighbor of u).
            for (size_t k = 0; k < neighbor.cols(); ++k) {
                w_i = position[neighbor[u][k]];
                z_i = w_i + 1;
                w = tour[w_i];
                z = tour[z_i % N];

                if (v == w || z == u) {
                    continue; // Skip adjacent edges.
                }

                // d[u][w] + min is a lower bound on new length.
                // d[u][v] + max is an upper bound on old length.
                if (d[u][w] + min > d[u][v] + max) {
                    break; // Go to next edge uv.
                }

                if (d[u][w] + d[v][z] < d[u][v] + d[w][z]) {
                    //   --u w--        --u-w->
                    //      X     ===>
                    //   <-z v->        <-z-v--
                    reverse(tour, v_i % N, w_i, position);
                    max = maximum(max, d[u][w], d[v][z]);
                    locallyOptimal = false;
                    break;
                }
            }
        }

}*/

vector<Node> twoOpt(const vector<Node> & solution, vector<vector<int> > matrix) {
    vector<Node> newSolution(solution);
    vector<Node> bestSolution(solution);
    int bestFitness = fitness(newSolution);
    //cout << "Fitness at begining : " << bestFitness << endl;
    bool improvement = true;
    while(improvement) {
        improvement = false;
        for(int i=0;i<newSolution.size();i++) {
            int nodeI = newSolution[i].id;
            int indexNextNodeToI = 0;
            if(i < (newSolution.size()-1)) {
                indexNextNodeToI = i+1;
            }
            int nextNodeToI = newSolution[indexNextNodeToI].id;
            for(int j=0;j<newSolution.size();j++) {
                if(i != j) {
                    int nodeJ = newSolution[j].id;
                    int indexNextNodeToJ = 0;
                    if(j < (newSolution.size()-1)) {
                        indexNextNodeToJ = j+1;
                    }
                    int nextNodeToJ = newSolution[indexNextNodeToJ].id;
                    if((matrix[nodeI][nextNodeToI] + matrix[nodeJ][nextNodeToJ]) > (matrix[nodeI][nextNodeToJ] + matrix[nodeJ][nextNodeToI])) {
                        newSolution = reverseGraph(newSolution, i, j);
                        int fit = fitness(newSolution);
                        if(fit < bestFitness) {
                            bestFitness = fit;
                            bestSolution = newSolution;
                            //cout << "Fitness : " << bestFitness << endl;
                            improvement = true;
                            //return bestSolution;
                        } else {
                            newSolution = bestSolution;
                        }
                    }
                }
            }
        }
        //improvement = false;
    }
    //cout << "Fitness : " << bestFitness << endl;
    return bestSolution;
}

vector<Node> reverseGraph(const vector<Node> & solution, int n1, int n2) {
    vector<Node> newSolution;
    vector<Node> reversePath;

    int i=0;
    while(i != n1 && i != n2) {
        newSolution.push_back(solution[i]);
        i++;
    }

    reversePath.push_back(solution[i]);i++;
    // From start to first point
    while(i != n1 && i != n2) {
        reversePath.push_back(solution[i]);
        i++;
    }
    reversePath.push_back(solution[i]);i++;

    for(int x=reversePath.size()-1;x>=0;x--) {
        newSolution.push_back(reversePath[x]);
    }

    while(i<solution.size()) {
        newSolution.push_back(solution[i]);
        i++;
    }

    return newSolution;
}

/*vector<Node> twoOpt(const vector<Node> & solution, vector<vector<int> > matrix) {
    vector<Node> newSolution;
    vector<Node> reversePath;
    int bestFitness = fitness(solution);
    /*int p1 = 0;
    int p2 = 0;
    int p3 = 0;
    int p4 = 0;
    int p1 = rand() % solution.size();
    int p2 = p1;
    int best = INT_MAX;
    for(int i=0;i<solution.size();++i) {
        if(i != p1 && i != p1+1 && i != p1-1) {
            if(!((i == solution.size() && p1 == 0) || (i == 0 && p1 == solution.size() ))) {
                int d = euclidianDistBetweenNodes(solution[p1], solution[i]);
                if(d < best) {
                    best = d;
                    p2 = i;
                }
            }
        }
    }

    //cout << "P1 : " << p1 << endl;
    //cout << "P2 : " << p2 << endl;
    /*do {
        p1 = rand() % solution.size();
        p2 = rand() % solution.size();
        if(p1 < solution.size()-2) {
            p3 = 0;
        } else {
            p3 = p1+1;
        }
        if(p2 < solution.size()-2) {
            p4 = 0;
        } else {
            p4 = p2+1;
        }
        ite++;
    } while((p1 == p2 || p1 == p4 || p3 == p2) && ite < 20);
    if(ite == 20) {
        return solution;
    }

    // From start to first point
    int i=0;
    while(i != p1 && i != p2) {
        newSolution.push_back(solution[i]);
        i++;
    }
    newSolution.push_back(solution[i]);

    // From start to first point
    i++;
    while(i != p1 && i != p2) {
        reversePath.push_back(solution[i]);
        i++;
    }

    std::reverse(reversePath.begin(),reversePath.end());

    for(int x=0;x<reversePath.size();x++) {
        newSolution.push_back(reversePath[x]);
    }

    while(i<solution.size()) {
        newSolution.push_back(solution[i]);
        i++;
    }

    if(fitness(newSolution) > bestFitness) {
        return newSolution;
    } else {
        return solution;
    }
}
*/

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

int euclidianDistBetweenNodes(Node n1, Node n2) {
    return euclidianDist(n1.x, n1.y, n2.x, n2.y);
}

int euclidianDist(double x1, double y1, double x2, double y2) {
    return round(sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2)));
}