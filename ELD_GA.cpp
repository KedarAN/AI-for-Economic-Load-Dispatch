#include <iostream>
#include <math.h>
#include <random>
#include <vector> 
#include <algorithm>
 
 
const int SAMPLE_SIZE = 1000; // How many top solutions to select
const int NUM = 60000;        // Number of solutions per generation
 
struct Solution
{
    double rank, x,y,ans;
    void fitness()
    {
         ans = (0.004*pow(x,2)+5.3*x+500)+(0.006*pow((y),2)+5.5*(y)+400)+(0.009*pow((975-x-y),2)+5.8*(975-x-y)+200);
        
		rank = (ans == 0) ? 9999 : std::abs(10000000/ans);
        
 
    }
};
 
int main()
{
 
    // Create initial random solutions
    //
    std::random_device device;
    std::uniform_real_distribution<double> unif1(200,450);
    std::uniform_real_distribution<double> unif2(150,350);
    std::uniform_int_distribution<int> cross(0,SAMPLE_SIZE-1);
    std::uniform_real_distribution<double> m(0.99,1.01);
 
    std::vector<Solution> solutions;
    std::vector<Solution> sample;
    
    for(int i = 0; i < NUM; i++)
        solutions.push_back(Solution{
            0, 
            unif1(device),
            unif2(device)
        });
    int a = 200;
    while(a > 0)
    {
        // Run our fitness function
        //
        for(auto& s : solutions) { s.fitness(); }
 
        // Sort our solutions by rank (optimize by using partial sort)
        //
        // std::partial_sort(
        //     solutions.begin(), 
        //     solutions.begin() + SAMPLE_SIZE, 
        //     solutions.end(), 
        // [](const auto& lhs,const auto& rhs){
        //     return lhs.rank > rhs.rank;
        // });
        
        // Sort our solutions by rank (slower compared to partial sort)
        //
        std::sort(
            solutions.begin(), 
            solutions.end(), 
        [](const auto& lhs,const auto& rhs){
            return lhs.rank > rhs.rank;
        });
        
        // Print top solutions
        //
        std::for_each(
            solutions.begin(), 
            solutions.begin()+1 , [=](const auto& s){
            std::cout << std::fixed 
                << "Rank " << static_cast<int>(s.rank)
                << "\nx:" << s.x <<" " << " y:" << s.y<<" z:"<<975-s.x-s.y<<" ans:"<<s.ans<<" \n";
        });
 
        // Take top solutions
        //
        
        
        std::copy(
            solutions.begin(), 
            solutions.begin() + SAMPLE_SIZE,
            std::back_inserter(sample)
        );
        solutions.clear();
 
        // Mutate the top solutions by %
 
        std::for_each(sample.begin(), sample.end(), [&](auto& s){
            s.x *= m(device);
            s.y *= m(device);
        });
 
        // Cross over
 
        for(int i = 0; i < NUM; i++)
        {
            solutions.push_back(Solution{
                0,
                (sample[cross(device)].x),
                (sample[cross(device)].y)
            });
        }
        sample.clear();
        
         std::for_each(
            solutions.begin(), 
            solutions.end() , [&](auto& s){
            if(s.x>450) s.x=450;
            if(s.x<200) s.x=200;
            if(s.y>350) s.y=350;
            if(s.y<150) s.y=150;
        });
        
	a-=1;
    }
 
}
