#include <iostream>
#include <vector>
#include <cmath>
#include <time.h>
#include <deque>
#include <climits>
#include <algorithm>
#include <ctime>

#define TIME_LIMIT 0.9
#define K 12 // K-Nearest neighbor

using namespace std;

// Structure of Node
typedef struct {
    int id;
    double x;
    double y;
} Node;

// Print a solution in Kattis Format (One node/line)
void printSolution(const vector<Node> &solutionInit);

// Calculating the euclidian distance between 2 pair of coordinate
int euclidianDist(double x1, double y1, double x2, double y2);

// Calculating the euclidian distance between 2 nodes
int euclidianDistBetweenNodes(Node n1, Node n2);

// Return a random number between 0 and 1
double random0and1();

// Calculating a fitness of a solution (sum of distances between nodes + distance from last node to the first node)
int fitness(const vector<Node> &solution);

// Create a matrix with distances
vector<vector<int> > createDistanceMatrix(const vector<Node> & solution);

// Function that reverse a graph (reverse a path between two nodes)
vector<Node> reverseGraph(const vector<Node> & solution, int n1, int n2);

// Simple algorithm to search an optimal path when we begin from node 0
vector<Node> greedyTour(const vector<Node> &solutionInit, int start, vector<vector<int> > matrix);

// Algorithm twoOpt optimal (hard to compute ~0.4s)
vector<Node> twoOpt(const vector<Node> & solution, vector<vector<int> > matrix);

// Algorithm twoOpt fast (~0.01s)
vector<Node> twoOptFast(const vector<Node> & solution, vector<vector<int> > matrix);

// Algorithm twoOpt with K nearest neighbord # TODO 
vector<Node> KtwoOpt(const vector<Node> & solution, vector<vector<int> > matrix, int k);

// Algorithm twoHalfOpt with K nearest neighbord # TODO {La fonction est bonne comme ça mais peut etre améliore en prenant réellement les K plus proches voisins}
vector<Node> KtwoHalfOpt(const vector<Node> & solution, vector<vector<int> > matrix, int k);

// META-Heuristique SA
vector<Node> simulatedAnnealing(const vector<Node>& solutionInit, vector<vector<int> > matrix, double temperatureInit, double temperatureMin, int energyMin, double delta = 0.5);

// Function which returns a neighbor solution (inverse two nodes randomly)
vector<Node> neighborhood(const vector<Node>& solution);

// Function which switch nodes in a graph
vector<Node> moveNodeInGraph(const vector<Node> & solution, int n1, int n2, int n3);

// Algorithm KthreeOpt with K nearest neighbord # TODO 
vector<Node> KthreeOpt(const vector<Node> & solution, vector<vector<int> > matrix, int k);


int main() {
    auto begin = clock();

    // Variables initialisation
    int N;
    srand (time(NULL));
    vector<Node> nodes;
    vector<Node> optimaleSolution;
    vector<vector<int> > distanceMatrix;

    double temperatureMin = 0;
    double temperatureInit = 30; // 60 iterations
    int energyMin = 0;

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
    //cout << "Fitness Greddy : " << fitness(optimaleSolution) << endl;

    optimaleSolution = twoOpt(optimaleSolution, distanceMatrix);
    //cout << "Fitness 2-Opt optimal : " << fitness(optimaleSolution) << endl;

    while (double(clock() - begin) / CLOCKS_PER_SEC < TIME_LIMIT) {
        optimaleSolution = simulatedAnnealing(optimaleSolution, distanceMatrix, temperatureInit , temperatureMin, energyMin);
        //cout << "Fitness SA : " << fitness(optimaleSolution) << endl;
        
        optimaleSolution = twoOptFast(optimaleSolution, distanceMatrix);
        //cout << "Fitness fast 2-Opt : " << fitness(optimaleSolution) << endl;
    }
    optimaleSolution = twoOpt(optimaleSolution, distanceMatrix);
    //cout << "Fitness 2-Opt optimal: " << fitness(optimaleSolution) << endl;
    optimaleSolution = KtwoHalfOpt(optimaleSolution, distanceMatrix, K);
    //cout << "Fitness K-2.5-Opt : " << fitness(optimaleSolution) << endl;

    // Print solution
    printSolution(optimaleSolution);

    return 0;
}

vector<Node> simulatedAnnealing(const vector<Node>& solutionInit, vector<vector<int> > matrix, double temperatureInit, double temperatureMin, int energyMin, double delta) {
    double t = temperatureInit;
    vector<Node> s(solutionInit);
    vector<Node> sBest(solutionInit);
    int e = fitness(s);
    //cout << "Temperature Initial : " << e << endl;
    double eBest = fitness(s);
    bool change = false;
    while(t > temperatureMin && e > energyMin) {
        vector<Node> sN(neighborhood(s));
        int eN = fitness(sN);
        //cout << "Temperature Neighbor : " << eN << endl;
        if((eN < e) || random0and1() < exp(-abs(double(eN) - double(e))/t) ) {
            s = sN;
            e = eN;
            if(eN < eBest) {
                eBest = eN;
                sBest = sN;//twoOpt(sN, matrix);
            }
            change = true;
        }
        t -= delta;
    }
    //cout << "Temperature Best : " << eBest << endl;
    return sBest;
}

vector<Node> neighborhood(const vector<Node>& solution) {
    vector<Node> newSolution(solution);
    int i = rand() % (newSolution.size()-2) + 1;
    int p1 = i-1;
    int p2 = i+1;

    Node mem = newSolution[p1];
    newSolution[p1] = newSolution[p2];
    newSolution[p2] = mem;
    return newSolution;
}

vector<Node> greedyTour(const vector<Node> &solutionInit, int start, vector<vector<int> > matrix) {
    vector<Node> solution = vector<Node>(solutionInit.size());
    vector<bool> used = vector<bool>(solutionInit.size(), false);
    solution[start] = solutionInit[start];
    used[start] = true;
    for(int i=0;i<solutionInit.size();i++) {
        int nodeI = solution[solutionInit.size()-1].id;
        if(i > 0) {
            nodeI = solution[i-1].id;
        }
        if(i != start) {
            int best = -1;
            for(int j=0;j<solutionInit.size();j++) {
                if(!used[j] && (best == -1 || matrix[nodeI][solutionInit[j].id] < matrix[nodeI][solutionInit[best].id])) {
                    best = j;
                }
            }
            solution[i] = solutionInit[best];
            used[best] = true;
        }
    }
   return solution;
}

vector<Node> twoOptFast(const vector<Node> & solution, vector<vector<int> > matrix) {
    vector<Node> newSolution(solution);
    vector<Node> bestSolution(solution);
    int bestFitness = fitness(newSolution);
    //cout << "Fitness at begining : " << bestFitness << endl;
    bool improvement = true;
    bool start_over = false;
    while(improvement) {
        improvement = false;
        for(int i=1;i<bestSolution.size() && !start_over;i++) {
            for(int j=i+1;j<bestSolution.size()-1 && !start_over;j++) { 
            //solution[i-1]->dist(solution[j]) + solution[i]->dist(solution[j+1]) < solution[i-1]->dist(solution[i]) + solution[j]->dist(solution[j+1]       
                if((matrix[solution[i-1].id][solution[j].id] + matrix[solution[i].id][solution[j+1].id]) < (matrix[solution[i-1].id][solution[i].id] + matrix[solution[j].id][solution[j+1].id])) {
                    newSolution = reverseGraph(bestSolution, i, j);
                    int fit = fitness(newSolution);
                    if(fit < bestFitness) {
                        bestFitness = fit;
                        bestSolution = newSolution;
                        //cout << "Fitness : " << bestFitness << endl;
                        improvement = true;
                        start_over = true;
                    } else {
                        //newSolution = bestSolution;
                        start_over = false;
                    }
                }
            }
        }
        if (!start_over)
            break;
    }
    //cout << "Fitness : " << bestFitness << endl;
    return bestSolution;
}

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
                    if((matrix[nodeI][nodeJ] < matrix[nodeI][nextNodeToI] || matrix[nodeI][nextNodeToI] < matrix[nodeJ][nextNodeToJ]) &&
                        (matrix[nodeI][nextNodeToI] + matrix[nodeJ][nextNodeToJ]) > (matrix[nodeI][nextNodeToJ] + matrix[nodeJ][nextNodeToI])) {
                        newSolution = reverseGraph(newSolution, i, j);
                        int fit = fitness(newSolution);
                        if(fit < bestFitness) {
                            bestFitness = fit;
                            bestSolution = newSolution;
                            //cout << "Fitness : " << bestFitness << endl;
                            improvement = true;
                        } else {
                            newSolution = bestSolution;
                        }
                    }
                }
            }
        }
    }
    //cout << "Fitness : " << bestFitness << endl;
    return bestSolution;
}

vector<Node> KtwoOpt(const vector<Node> & solution, vector<vector<int> > matrix, int k) {
    vector<Node> newSolution(solution);
    vector<Node> bestSolution(solution);
    int bestFitness = fitness(newSolution);
    if(k > newSolution.size()) {
        k = newSolution.size();
    }
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
            for(int j=0;j<k;j++) {
                if(i != j) {
                    int nodeJ = newSolution[j].id;
                    int indexNextNodeToJ = 0;
                    if(j < (newSolution.size()-1)) {
                        indexNextNodeToJ = j+1;
                    }
                    int nextNodeToJ = newSolution[indexNextNodeToJ].id;
                    if((matrix[nodeI][nodeJ] < matrix[nodeI][nextNodeToI] || matrix[nodeI][nextNodeToI] < matrix[nodeJ][nextNodeToJ]) &&
                        (matrix[nodeI][nextNodeToI] + matrix[nodeJ][nextNodeToJ]) > (matrix[nodeI][nextNodeToJ] + matrix[nodeJ][nextNodeToI])) {
                        newSolution = reverseGraph(newSolution, i, j);
                        int fit = fitness(newSolution);
                        if(fit < bestFitness) {
                            bestFitness = fit;
                            bestSolution = newSolution;
                            //cout << "Fitness : " << bestFitness << endl;
                            improvement = true;
                        } else {
                            newSolution = bestSolution;
                        }
                    }
                }
            }
        }
     }
    //cout << "Fitness : " << bestFitness << endl;
    return bestSolution;
}

vector<Node> KthreeOpt(const vector<Node> & solution, vector<vector<int> > matrix, int k) {
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
            int borneMax = k;
            if(k > matrix.size()-1) {
                borneMax = matrix.size()-1;
            }
            for(int j=0;j<borneMax;j++) {
                if(i != j) {
                    int nodeJ = newSolution[j].id;
                    int indexNextNodeToJ = 0;
                    if(j < (newSolution.size()-1)) {
                        indexNextNodeToJ = j+1;
                    }
                    int nextNodeToJ = newSolution[indexNextNodeToJ].id;
                    int borneMax2 = k;
                    if(k > matrix.size()-1) {
                        borneMax2 = matrix.size()-1;
                    }
                    for(int l=0;l<borneMax2;l++) {
                        if(i != l && l != j) {
                            int nodeL = newSolution[l].id;
                            int indexNextNodeToL = 0;
                            if(l < (newSolution.size()-1)) {
                                indexNextNodeToL = l+1;
                            }
                            int nextNodeToL = newSolution[indexNextNodeToL].id;
                            if((matrix[nodeI][nodeJ] < matrix[nodeI][nextNodeToI] || matrix[nodeI][nextNodeToI] < matrix[nodeJ][nextNodeToJ]) &&
                                (matrix[nodeI][nextNodeToI] + matrix[nodeJ][nextNodeToJ]) > (matrix[nodeI][nextNodeToJ] + matrix[nodeJ][nextNodeToI])) {
                                newSolution = reverseGraph(newSolution, i, j);
                                newSolution = reverseGraph(newSolution, j, l);
                                int fit = fitness(newSolution);
                                if(fit < bestFitness) {
                                    bestFitness = fit;
                                    bestSolution = newSolution;
                                    //cout << "Fitness : " << bestFitness << endl;
                                    improvement = true;
                                } else {
                                    newSolution = bestSolution;
                                }
                            }
                            /*if((matrix[nodeI][nextNodeToI] + matrix[nodeJ][nextNodeToJ]) > (matrix[nodeI][nextNodeToJ] + matrix[nodeJ][nextNodeToI]) &&
                                matrix[nodeJ][nextNodeToJ] + matrix[nodeL][nextNodeToL] > (matrix[nodeJ][nextNodeToL] + matrix[nodeL][nextNodeToJ])) {*/
                                
                                //
                                
                            //}
                        }
                    }
                }
            }
        }
    }
    //cout << "Fitness : " << bestFitness << endl;
    return bestSolution;
}

vector<Node> KtwoHalfOpt(const vector<Node> & solution, vector<vector<int> > matrix, int k) {
    vector<Node> newSolution(solution);
    vector<Node> bestSolution(solution);
    int bestFitness = fitness(newSolution);
    //cout << "Fitness at begining : " << bestFitness << endl;
    bool improvement = true;
    while(improvement) {
        improvement = false;
        for(int i=0;i<newSolution.size();i++) {
            int borneMax = i+3+k;
            if(borneMax > newSolution.size()) {
                borneMax = newSolution.size();
            }
            for(int j=i+3;j<borneMax;j++) {
                int opt1_cur = matrix[newSolution[i].id][newSolution[i+1].id] +
                                matrix[newSolution[i+1].id][newSolution[i+2].id] +
                                matrix[newSolution[j-1].id][newSolution[j].id];

                int option1 = matrix[newSolution[i].id][newSolution[i+2].id] +
                                matrix[newSolution[j-1].id][newSolution[i+1].id] +
                                matrix[newSolution[i+1].id][newSolution[j].id];

                if(option1 < opt1_cur){
                    Node temp = newSolution[i+1];
                    for(int m=i+2;m<j;m++) {
                        newSolution[m-1] = newSolution[m];
                    }
                    newSolution[j-1] = temp;
                    improvement = true;
                }

                int opt2_cur = matrix[newSolution[i].id][newSolution[i+1].id] +
                                matrix[newSolution[j-2].id][newSolution[j-1].id] +
                                matrix[newSolution[j-1].id][newSolution[j].id];
                int option2 = matrix[newSolution[i].id][newSolution[j-1].id] +
                                matrix[newSolution[j-1].id][newSolution[i+1].id] +
                                matrix[newSolution[j-2].id][newSolution[j].id];
                
                if (option2 < opt2_cur) {
                    Node temp = newSolution[j-1];
                    for(int m=j-2;m>i;m--) {
                        newSolution[m+1] = newSolution[m];
                    }
                    newSolution[i+1] = temp;
                    improvement = true;
                }

                int fit = fitness(newSolution);

                if(bestFitness > fit) {
                    bestFitness = fit;
                    bestSolution = newSolution;
                }

            }
        }
    }
    //cout << "Fitness : " << bestFitness << endl;
    return bestSolution;
}

vector<Node> moveNodeInGraph(const vector<Node> & solution, int n1, int n2, int n3) {
    if(n1 == n2 || n1 ==n3 || n2 == n3) {
        return solution;
    }
    vector<Node> newSolution(solution);

    if (n1 < n2) {
        for (int i = n1; i < n2; i++) {
            newSolution[i] = newSolution[i+1];
        }
        newSolution[n2] = newSolution[n1];
    } else {
        for (int i = n1; i > n3; i--) {
            newSolution[i] = newSolution[i-1];
        }
        newSolution[n3] = newSolution[n1];
    }
    return newSolution;
}

vector<Node> reverseGraph(const vector<Node> & solution, int n1, int n2) {
    if(n1 == n2) {
        return solution;
    }
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