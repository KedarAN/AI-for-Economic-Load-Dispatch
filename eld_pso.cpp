#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <math.h>
 
using namespace std;
 
 
double evaluate(vector<double> x) {
    // Evaluate the fitness of the particle's position using the Rosenbrock function
    double rank,ans;
    
     ans = (0.004*pow(x[0],2)+5.3*x[0]+500)+(0.006*pow((x[1]),2)+5.5*(x[1])+400)+(0.009*pow((975-x[0]-x[1]),2)+5.8*(975-x[0]-x[1])+200);
    
    return ans;
}
class Particle {
public:
    vector<double> position; // Particle's position in the search space
    vector<double> velocity; // Particle's velocity in the search space
    vector<double> best_position; // Particle's personal best position
    double best_fitness; // Particle's personal best fitness
 
    // Constructor to initialize the particle's position and velocity
    Particle(int dim) {
        position.resize(dim);
        velocity.resize(dim);
        best_position.resize(dim);
        for (int i = 0; i < dim; i++) {
            if(i==0)
            position[i] = abs(((double) rand() / RAND_MAX)) * 250 + 200;
            else
            position[i] = abs(((double) rand() / RAND_MAX)) * 200 + 150;
          
             // Initialize the position to a random value between -5 and 5
            velocity[i] = 0; // Initialize the velocity to zero
            best_position[i] = position[i]; // Initialize the personal best position to the initial position
        }
        best_fitness = INFINITY; // Initialize the personal best fitness to infinity
    }
 
    // Update the particle's velocity based on the global best position and the particle's personal best position
    void update_velocity(vector<double> global_best_position, double omega, double phi_p, double phi_g) {
        for (int i = 0; i < position.size(); i++) {
            double r_p = ((double) rand() / RAND_MAX); // Generate a random number between 0 and 1
            double r_g = ((double) rand() / RAND_MAX); // Generate a random number between 0 and 1
            velocity[i] = omega * velocity[i] + phi_p * r_p * (best_position[i] - position[i]) + phi_g * r_g * (global_best_position[i] - position[i]); // Update the velocity based on the PSO equation
        }
    }
 
    // Update the particle's position based on its velocity
    void update_position() {
        for (int i = 0; i < position.size(); i++) {
            position[i] += velocity[i]; // Update the position based on the velocity
            if(i==0){
            if(position[i]>450) position[i]=450;
            if(position[i]<200) position[i]=200;
            }
            else{
            if(position[i]>350) position[i]=350;
            if(position[i]<150) position[i]=150;
            }
        }
 
    }
};
 
vector<double> pso(int dim, int num_particles, int num_iterations) {
    
    // Initialize the particles and global best position
    vector<Particle*> particles(num_particles);
 
    for(int i=0;i<num_particles;i++){
        Particle *p = new Particle(dim);
        particles[i]= p;
            
    }
    
  
    vector<double> global_best_position(dim);
    double global_best_fitness = INFINITY;
    
 
    // Run the PSO algorithm for the specified number of iterations
    for (int iter = 0; iter < num_iterations; iter++) {
       
        // Update the personal best position and fitness of each particle
        for (int i = 0; i < num_particles; i++) {
            double fitness = evaluate(particles[i]->position);
            // Evaluate the fitness of the particle's
            if (fitness < particles[i]->best_fitness) { // Update the personal best position and fitness if the fitness is better than the current personal best fitness
                particles[i]->best_position = particles[i]->position;
                particles[i]->best_fitness = fitness;
            }
            if (fitness < global_best_fitness) { // Update the global best position and fitness if the fitness is better than the current global best fitness
                global_best_position = particles[i]->position;
                global_best_fitness = fitness;
            }
    }
 
    // Update the velocity and position of each particle
        for (int i = 0; i < num_particles; i++) {
            particles[i]->update_velocity(global_best_position, 0.5, 0.5, 0.5); // Update the velocity using the PSO equation with omega = 0.5, phi_p = 0.5, and phi_g = 0.5
            particles[i]->update_position(); // Update the position based on the velocity
        }
    }
 
return global_best_position; // Return the global best position
}
 
 
int main() {
    // Set the random seed
    //srand(time(NULL));
    // Set the dimension of the search space, the number of particles, and the number of iterations
    int dim = 2;
    int num_particles = 10;
    int num_iterations = 100;
 
    // Run the PSO algorithm and print the global best position and fitness
    vector<double> global_best_position = pso(dim, num_particles, num_iterations);
    double global_best_fitness = evaluate(global_best_position);
    cout << "Global best position: ";
    for (int i = 0; i < dim; i++) {
        cout << global_best_position[i] << " ";
    }
    cout << endl;
    cout << "Global best fitness: " << global_best_fitness << endl;
 
    return 0;
}
