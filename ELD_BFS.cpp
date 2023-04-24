#include <iostream>
#include <queue>
#include <vector>
#include <cmath>

using namespace std;

const int MAXN = 100;
const int INF = 1e9;

int N; // number of generating units
int C[MAXN]; 
 double P_MIN[2] = {200.0, 150.0};  
 double P_MAX[2] = {450.0, 350.0};


int calc_cost(vector<int>& x) {
    int cost = 0;
    cost = (0.004*pow(x[0],2)+5.3*x[0]+500)+(0.006*pow((x[1]),2)+5.5*(x[1])+400)+(0.009*pow((975-x[0]-x[1]),2)+5.8*(975-x[0]-x[1])+200);   

    return cost;
}

int economic_load_dispatch() {
    queue<vector<int>> q;
    // initialize the queue with the state where all generating units are operating at their minimum output levels
    vector<int> s(P_MIN,P_MIN+2);
    q.push(s);
    
    int min_cost = INF;
    
    while (!q.empty()) {
        vector<int> u = q.front(); 
        q.pop();
        int cost = calc_cost(u);
        if (cost < min_cost) {
            min_cost = cost;
            cout << u[0]<<" "<<u[1]<<" operating cost: " << calc_cost(u) << endl;
   
        }
        // explore the neighbors of the current state by increasing the output of each generator by 1
        
        for (int i = 0; i < N; i++) {
          vector<int> v = u; 
            if (u[i] < P_MAX[i]) {
                v[i]++;

                if (calc_cost(v) < calc_cost(u)) {
                    q.push(v);
                }
            }
        }
    }
    return min_cost;
}

int main() {
     N=2;
    
    // solve the economic load dispatch problem using the BFS algorithm
    int min_cost = economic_load_dispatch();
    // output the minimum operating cost
    cout << "Minimum operating cost: " << min_cost << endl;
    return 0;
}
