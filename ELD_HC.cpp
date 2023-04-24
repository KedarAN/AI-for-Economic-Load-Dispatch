#include <iostream>
#include <cmath>
#include <random>
 
using namespace std;
 
const int MAX_GEN = 10000;  // Maximum number of iterations
const int N = 2;  // Number of generators
const double PL = 975.0;  // Load demand
//const double P_MIN[N] = {200.0, 150.0, 100.0};  // Minimum power output of generators
//const double P_MAX[N] = {450.0, 350.0, 225.0};  // Maximum power output of generators
const double P_MIN[N] = {200.0, 150.0};  
const double P_MAX[N] = {450.0, 350.0};
double calc_cost2(double P[N]) {
    double cost = 0.0;
    cost = (0.004*pow(P[0],2)+5.3*P[0]+500)+(0.006*pow((P[1]),2)+5.5*(P[1])+400)+(0.009*pow((975-P[0]-P[1]),2)+5.8*(975-P[0]-P[1])+200);   
    return cost;
}
 
// Perform Hill Climbing optimization
void hill_climbing(double P[N]) {
    cout<<"in hc algo"<<endl;
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dist(0, N-1);
    uniform_int_distribution<> dist2(-100, 100);
 
 
    double P_best[N];
            for (int i = 0; i < N; i++) {
                    P_best[i]=P[i]; 
                }
    
    double cost = calc_cost2(P);
    double best_cost = cost;
 
    for (int k = 0; k < MAX_GEN; k++) {
        
        for (int i = 0; i < N; i++) {
            P[i]=P_best[i]; 
        }
 
        int l = dist(gen);  // Randomly select a generator to adjust
        int j = dist2(gen);
        double P_step = (j* P[l])/100;  // Step size for power output adjustment
        P[l] = P[l] + P_step;
                
 
        // Check if power output exceeds maximum limit
        if (P[l] > P_MAX[l]) {
            P[l] = P_MAX[l];
        }
        if (P[l] < P_MIN[l]) {
            P[l] = P_MIN[l];
        }
        cost = calc_cost2(P);
         cout<<P[0]<<" "<<P[1]<<" "<< 975-P[0]-P[1]<<" "<< cost << endl;
 
        // Check if power balance constraint is satisfied
            if (cost < best_cost) {
                best_cost = cost;
                for (int i = 0; i < N; i++) {
                    P_best[i]=P[i]; 
                }
            }
 
    }
    cout << "Optimal power outputs: ";
    for (int i = 0; i < N; i++) {
        cout << P_best[i] << " ";
    }
    cout<< 975-P_best[0]-P_best[1];
    cout << endl;
    cout << "Total cost of power generation: " << calc_cost2(P_best) << endl;
}
 
int main() {
    // Initial power outputs (randomly generated)
    double P[N] = {200, 150};
 
    hill_climbing(P);
    return 0;
}
