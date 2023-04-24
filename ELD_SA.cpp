#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
using namespace std;
 
double objective_function(std::vector<double> x) {
    double sum = 0.0;
    sum = (0.004*pow(x[0],2)+5.3*x[0]+500)+(0.006*pow((x[1]),2)+5.5*(x[1])+400)+(0.009*pow((975-x[0]-x[1]),2)+5.8*(975-x[0]-x[1])+200);   
    //cout<<" x1: "<<x[0]<<" x2: "<<x[1]<<" cost: "<<sum<<endl;
        
    
    return sum;
 
}
 
int main() {
    // Define the initial solution
    srand(time(NULL));
    std::vector<double> x_initial(2);
    for (int i = 0; i < x_initial.size(); i++) {
 
         if(i==0)
            x_initial[i] = abs(((double) rand() / RAND_MAX)) * 250 + 200;
            else
            x_initial[i] = abs(((double) rand() / RAND_MAX)) * 200 + 150;
          
    }
 
 
    // Define the temperature schedule
    double T_initial = 150000000000;
    double T_final = 150000000;
    double cooling_rate = 0.9999999;
 
    // Define the maximum number of iterations
    auto max_iterations = 10000000000;
 
    // Initialize the current solution and cost
    std::vector<double> x_current = x_initial;
    double cost_current = objective_function(x_current);
 
    // Initialize the best solution and cost
    std::vector<double> x_best = x_initial;
    double cost_best = objective_function(x_best);
    cout <<"running..."<<endl;
    // Simulated annealing algorithm
    for (int i = 0; i < max_iterations; i++) {
        // Generate a new solution
        std::vector<double> x_new(x_current.size());
        
        for (int j = 0; j < x_new.size(); j++) {
            x_new[j] = x_current[j] + abs(((double) rand() / RAND_MAX))* 2000 - 1000;
            if(j==0){
            if(x_new[j] >450) x_new[j] =450;
            if(x_new[j] <200) x_new[j] =200;
            }
            else{
            if(x_new[j] >350) x_new[j] =350;
            if(x_new[j] <150) x_new[j] =150;
            }
        }
 
        double cost_new = objective_function(x_new);
    
        // Calculate the acceptance probability
        double delta_cost = cost_new - cost_current;
        double acceptance_probability = exp(-delta_cost / T_initial);
 
        // Decide whether to accept the new solution
        if (acceptance_probability > (double) rand() / RAND_MAX) {
            x_current = x_new;
            cost_current = cost_new;
        }
 
        // Update the best solution and cost
        if (cost_current < cost_best) {
            x_best = x_current;
            cost_best = cost_current;
            cout<<" x1: "<<x_best[0]<<" x2: "<<x_best[1]<<" cost: "<<cost_best<<endl;
        }
 
        //Decrease the temperature
        T_initial *= cooling_rate;
        if (T_initial < T_final) {
            break;
        }
    }
 
    // Print the results
    std::cout << "Final solution: ";
    for (int i = 0; i < x_best.size(); i++) {
        std::cout << x_best[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "Minimum cost: " << cost_best << std::endl;
 
    return 0;
}
