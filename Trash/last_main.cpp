#include <iostream>
#include <vector>
#include <cmath>
#include <time.h>
#include <deque>
#include <climits>
#include <algorithm>

using namespace std;

int euclidianDist(double x1, double y1, double x2, double y2);
vector<pair<int, pair<double,double> > > annealingSimulated(double temperatureInit, const vector<pair<int, pair<double,double> > >& solutionInit, double temperatureMin, int energyMin, double delta = 0.9);
vector<pair<int, pair<double,double> > > tabuSearch(const vector<pair<int, pair<double,double> > >& solutionInit, int maxIteration, int sizeTabuMemory);
vector<pair<int, pair<double,double> > > neighborhood(const vector<pair<int, pair<double,double> > >& solution);
vector<vector<pair<int, pair<double,double> > > > neighborhood2(const vector<pair<int, pair<double,double> > >& solution);
vector<pair<int, pair<double,double> > > neighborhood2opt(const vector<pair<int, pair<double,double> > >& solution);
vector<vector<pair<int, pair<double,double> > > > neighborhood2opt2(const vector<pair<int, pair<double,double> > >& solution);
int fitness(const vector<pair<int, pair<double,double> > > &solution);
double random0and1();
void printSolution( const vector<pair<int, pair<double,double> > >& solutionInit);
bool solutionInMemory(const deque<vector<pair<int, pair<double,double> > > >& memory, const vector<pair<int, pair<double,double> > > solution);
vector<pair<int, pair<double,double> > > initialSolution(vector<pair<int, pair<double,double> > > points);

int main() {
    // Variables init
    srand (time(NULL));
    int N;
    vector<pair<int, pair<double,double> > > points;
    vector<pair<int, pair<double,double> > > initialSol;
    vector<pair<int, pair<double,double> > > optimaleSolution;

    // Get number of points
    cin >> N;

    // Get the N points
    for(int i=0;i<N;i++) {
        double x, y;
        cin >> x >> y;
        pair<double,double> value(x,y);
        pair<int,pair<double,double> > p(i,value);
        points.push_back(p);
    }

    initialSol = initialSolution(points);
    //printSolution(initialSol);
    //cout << "---------" << endl;
    // Call the metha heuristic
    
    double temperatureMin = 0;
    double temperatureInit = 100000;
    int energyMin = 0;

    /*printSolution(points);
    cout << endl;
    printSolution(initialSol);
    cout << endl;*/
    optimaleSolution = annealingSimulated(temperatureInit, initialSol, temperatureMin, energyMin);
    
    
    int maxIteration = 5000;
    int sizeTabuMemory = 5;
    //optimaleSolution = tabuSearch(initialSol, maxIteration, sizeTabuMemory);
    
    // Print solution
    printSolution(optimaleSolution);
    //printSolution(points);
}

vector<pair<int, pair<double,double> > > initialSolution(vector<pair<int, pair<double,double> > > points) {
    int start = rand() % points.size(); // Starting point
    vector<bool> nodeVisited = vector<bool>(points.size(), false);
    vector<pair<int, pair<double,double> > > initSolution;
    nodeVisited[start] = true;
    initSolution.push_back(points[start]);
    bool continued = true;

    while(continued) {
        int bestPoint = -1;
        int bestPath = INT_MAX;
        for(int i=0;i<points.size();i++) {
            if(points[i].first != initSolution[initSolution.size()-1].first && !nodeVisited[points[i].first]) {
                int dist = euclidianDist(initSolution[initSolution.size()-1].second.first, initSolution[initSolution.size()-1].second.second, points[i].second.first, points[i].second.second);
                if(dist<bestPath) {
                    bestPath = dist;
                    bestPoint = i;
                }
            }
        }
        if(bestPoint != -1) {
            nodeVisited[bestPoint] = true;
            initSolution.push_back(points[bestPoint]);
        } else {
            continued = false;
        }
    }
    return initSolution;
}

int fitness(const vector<pair<int, pair<double,double> > > &solution) {
    int f = 0;
    for(int i=0;i<solution.size()-1;++i) {
        f += euclidianDist(solution[i].second.first, solution[i].second.second, solution[i+1].second.first, solution[i+1].second.second);
    }
    if(solution.size()-1 > 1) {
        f += euclidianDist(solution[solution.size()-1].second.first, solution[solution.size()-1].second.second, solution[0].second.first, solution[0].second.second);
    }
    return f;
}

double random0and1() {
    return static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
}

vector<vector<pair<int, pair<double,double> > > > neighborhood2(const vector<pair<int, pair<double,double> > >& solution) {
    vector<vector<pair<int, pair<double,double> > > > solutions;
    for(int t=0;t<10;t++) {
        solutions.push_back(neighborhood(solution));
    }
    return solutions;
}

vector<pair<int, pair<double,double> > > neighborhood(const vector<pair<int, pair<double,double> > >& solution) {
    vector<pair<int, pair<double,double> > > newSolution(solution);
    int p1 = 0;
    int p2 = 0;
    int ite = 0;
    do {
        p1 = rand() % newSolution.size();
        p2 = rand() % newSolution.size();
        ite++;
    } while(p1 == p2 && ite < 30);
    if(ite == 30) {
        if((p1+1)<newSolution.size()) {
            p2 = p1+1;
        } else {
            p2 = p1-1;
        }  
    }

    pair<int, pair<double,double> > mem = newSolution[p1];
    newSolution[p1] = newSolution[p2];
    newSolution[p2] = mem;
    return newSolution;
}

vector<pair<int, pair<double,double> > > neighborhood2opt(const vector<pair<int, pair<double,double> > >& solution) {
    vector<pair<int, pair<double,double> > > newSolution;
    vector<pair<int, pair<double,double> > > reversePath;
    int p1 = 0;
    int p2 = 0;
    int p3 = 0;
    int p4 = 0;
    do {
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
    } while(p1 == p2 || p1 == p4 || p3 == p2);

    //cout << "We choose : " << p1 << " " << p2 << endl;

    // From start to first point
    int i=0;
    while(i != p1 && i != p2) {
        newSolution.push_back(solution[i]);
        i++;
    }
    newSolution.push_back(solution[i]);
    //cout << newSolution.size() << endl;
    // From start to first point
    i++;
    while(i != p1 && i != p2) {
        reversePath.push_back(solution[i]);
        i++;
    }
    //cout << reversePath.size() << endl;
    std::reverse(reversePath.begin(),reversePath.end());
    //cout << reversePath.size() << endl;
    for(int x=0;x<reversePath.size();x++) {
        newSolution.push_back(reversePath[x]);
    }
    //cout << newSolution.size() << endl;

    while(i<solution.size()) {
        newSolution.push_back(solution[i]);
        i++;
    }
    //cout << newSolution.size() << endl;
    if(newSolution.size() != solution.size()) {
        cout << "COUILLLLLLLLE" << endl;
    }
    return newSolution;
}

vector<vector<pair<int, pair<double,double> > > > neighborhood2opt2(const vector<pair<int, pair<double,double> > >& solution) {
    vector<vector<pair<int, pair<double,double> > > > solutions;
    for(int t=0;t<10;t++) {
        solutions.push_back(neighborhood2opt(solution));
    }
    return solutions;
}

vector<pair<int, pair<double,double> > > annealingSimulated(double temperatureInit, const vector<pair<int, pair<double,double> > >& solutionInit, double temperatureMin, int energyMin, double delta) {
    double t = temperatureInit;
    vector<pair<int, pair<double,double> > > s(solutionInit);
    vector<pair<int, pair<double,double> > > sBest(solutionInit);
    int e = fitness(s);
    double eBest = fitness(s);
    //cout << "Temperature Initial : " << e << endl;
    while(t > temperatureMin && e > energyMin) {
        vector<pair<int, pair<double,double> > > sN(neighborhood(s));
        int eN = fitness(sN);
        //cout << "Temperature Neighbor : " << eN << endl;
        if((eN < e) || random0and1() < exp(-abs(double(eN) - double(e))/t) ) {
            s = sN;
            e = eN;
            if(eN < eBest) {
                eBest = eN;
                sBest = sN;
            }
        }
        t -= delta;
    }
    //cout << "Temperature Best : " << eBest << endl;
    return sBest;
}

void printSolution( const vector<pair<int, pair<double,double> > >& solutionInit) {
    for(int i=0;i<solutionInit.size();++i) {
        cout << solutionInit[i].first << endl;
    }
    //cout << endl;
}

vector<pair<int, pair<double,double> > > tabuSearch(const vector<pair<int, pair<double,double> > >& solutionInit, int maxIteration, int sizeTabuMemory) {
    vector<pair<int, pair<double,double> > > s(solutionInit);
    vector<pair<int, pair<double,double> > > sBest(solutionInit);
    int currentIteration = 0;
    deque<vector<pair<int, pair<double,double> > > > memory;
    memory.push_back(s);
    //cout << "Fitness Initial : " << fitness(sBest) << endl;
    while(maxIteration > currentIteration) {
        vector<vector<pair<int, pair<double,double> > > > solutions(neighborhood2(s));
        //vector<pair<int, pair<double,double> > > sN(neighborhood(s));
        vector<pair<int, pair<double,double> > > sN(solutions[0]);

        for(int i=0;i<solutions.size();++i) {
            if(!solutionInMemory(memory, solutions[i]) && fitness(solutions[i]) < fitness(sN)) {
                sN = solutions[i];
            }
        }

        if(fitness(sN) < fitness(sBest)) {
            sBest = sN;
        }
        memory.push_back(sN);

        if(memory.size() > sizeTabuMemory) {
            memory.pop_front();
        }

        currentIteration++;
    }
    //cout << "Fitness Best : " << fitness(sBest) << endl;
    return sBest;
}

bool solutionInMemory(const deque<vector<pair<int, pair<double,double> > > >& memory, const vector<pair<int, pair<double,double> > > solution) {
    bool found = false;
    for(int i=0;i<memory.size();i++) {
        found = true;
        for(int j=0;j<memory[i].size();j++) {
            if(solution[j] != memory[i][j]) {
                found = false;
                break;
            }
        }
        if(found)
            break;
    }
    return found;
}

int euclidianDist(double x1, double y1, double x2, double y2) {
    return round(sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2)));
}

