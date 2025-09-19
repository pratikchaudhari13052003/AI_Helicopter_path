#include "solver.h"
#include <iostream>
#include <chrono>
#include <random>
#include <algorithm>


using namespace std;
// You can add any helper functions or classes you need here.

/**
 * @brief The main function to implement your search/optimization algorithm.
 * * This is a placeholder implementation. It creates a simple, likely invalid,
 * plan to demonstrate how to build the Solution object. 
 * * TODO: REPLACE THIS ENTIRE FUNCTION WITH YOUR ALGORITHM.
 */

 // global or static engine
mt19937 rng(random_device{}()); 


int randomInt(int low, int high) {
    uniform_int_distribution<int> dist(low, high);
    return dist(rng);
}

double randomDouble(double low, double high) {
    uniform_real_distribution<double> dist(low, high);
    return dist(rng);
}

void possible_villages(const ProblemData& problem,vector<int> possible[]){

    double d_max=problem.d_max;
    
    for(auto helicopter:problem.helicopters){
        int heli_id=helicopter.id;
        int heli_home_city_id=problem.helicopters[heli_id-1].home_city_id;
        
        double distance_capacity=helicopter.distance_capacity;

    for(auto village:problem.villages){

        if(2*distance(problem.cities[helicopter.home_city_id-1],village.coords) 
                <=distance_capacity){
                    possible[heli_id-1].push_back(village.id);
                }


    }
}



}

//Here i have given one trip to each helicopter 
// and randomly assign village to each helicopter 
// and also food is generated randomly
//Improve this by making multiple trips and make it valid

Solution randomGenerator(const ProblemData& problem,vector<int> possible[]){

    Solution sol;

    vector<int> food_delivered(problem.villages.size(),0);
    vector<int> food_delivered_other(problem.villages.size(),0);

    //I am writing it for one trip only have to make multiple trips but good for debuging

    for (const auto & helicopter:problem.helicopters){

        
        HelicopterPlan plan;
        // sol[0].push_back(helicopter.id);
        
        plan.helicopter_id=helicopter.id;
        
        int heli_id=helicopter.id;
        int home_city_id=helicopter.home_city_id;
        
        if (possible[heli_id-1].empty()) {
           sol.push_back(plan);
           continue;
       }

        double weight_capacity=helicopter.weight_capacity;
        
    

       int maxVisits = min((int)possible[heli_id-1].size(), 3); // Limit to 3 for feasibility
        int numVisits = randomInt(1, maxVisits); 

        // Shuffle possible villages so visits are random
        vector<int> shuffled = possible[heli_id-1];
        random_shuffle(shuffled.begin(), shuffled.end());
        Trip trip;

        // int last_coords=helicopter.home_city_id;
        double total_dry_food=0;
        double total_peri_food=0;
        double total_other_food=0;

        
        double distance_capacity=helicopter.distance_capacity;

        double total_distance=0;
        Point last_coords = problem.cities[home_city_id-1];
        
        for(int i=1;i<=numVisits;i++){

            int village_id=shuffled[i-1];

            total_distance += distance(problem.villages[village_id-1].coords , last_coords);
            last_coords=problem.villages[village_id-1].coords;

            if(total_distance + distance(problem.cities[home_city_id-1] , last_coords) >=  distance_capacity){
                break;
            }
            
            Drop d;
            
            int total_population=problem.villages[village_id-1].population;
            
            double dry_food=0;
            double peri_food=0;
            double other_food=0;
            
            peri_food = randomInt(0.0,(9*total_population - food_delivered[village_id-1]));
            dry_food = randomInt(0.0,(9*total_population - food_delivered[village_id-1]-dry_food));
            other_food = randomInt(0.0,(total_population -food_delivered_other[village_id-1]));
            
            food_delivered[village_id-1]+=dry_food;
            food_delivered[village_id-1]+=peri_food;
            food_delivered_other[village_id-1]+=other_food;
            
            d.village_id=i;
            d.dry_food=dry_food;
            d.perishable_food=peri_food;
            d.other_supplies=other_food;
            
            trip.drops.push_back(d);
            
        }
 // Set pickup amounts based on total drops
        trip.dry_food_pickup = 0;
        trip.perishable_food_pickup = 0;
        trip.other_supplies_pickup = 0;
        
        for (const auto& drop : trip.drops) {
            trip.dry_food_pickup += drop.dry_food;
            trip.perishable_food_pickup += drop.perishable_food;
            trip.other_supplies_pickup += drop.other_supplies;
        }
        
        if (!trip.drops.empty()) {
            plan.trips.push_back(trip);
        }
        
        sol.push_back(plan);
    }

   


    return sol;

}


double evaluateSolution(const ProblemData& problem, Solution& current) {
    const double INVALID_SCORE = -1e9;
    
    // Global tracking arrays - persist across entire solution
    vector<double> food_delivered(problem.villages.size() + 1, 0.0);
    vector<double> other_delivered(problem.villages.size() + 1, 0.0);
    vector<double> helicopter_total_distances(problem.helicopters.size() + 1, 0.0);
    
    double total_trip_cost = 0.0;
    double total_value = 0.0;

    for (const auto &plan : current) {
        int heli_id = plan.helicopter_id; 
                
        if (heli_id <= 0 || heli_id > (int)problem.helicopters.size()) return INVALID_SCORE;
        const Helicopter &heli = problem.helicopters[heli_id - 1];

        for (const auto &trip : plan.trips) {
            // Check weight constraint
            double trip_pickup_weight = trip.dry_food_pickup * problem.packages[0].weight
                                      + trip.perishable_food_pickup * problem.packages[1].weight
                                      + trip.other_supplies_pickup * problem.packages[2].weight;
            if (trip_pickup_weight > heli.weight_capacity + 1e-9) return INVALID_SCORE;

            // Track trip distance and drops
            int home_city_idx = heli.home_city_id;
            if (home_city_idx <= 0 || home_city_idx > (int)problem.cities.size()) return INVALID_SCORE;
            Point current_location = problem.cities[home_city_idx - 1];

            double trip_distance = 0.0;
            int total_drop_dry = 0;
            int total_drop_peri = 0;
            int total_drop_other = 0;

            for (const auto &drop : trip.drops) {
                int village_id = drop.village_id;
                if (village_id <= 0 || village_id > (int)problem.villages.size()) return INVALID_SCORE;
                const Village &village = problem.villages[village_id - 1];
                
                // Calculate trip distance
                trip_distance += distance(current_location, village.coords);
                current_location = village.coords;

                // Value Capping Logic - EXACTLY like format_checker
                double max_food_needed = village.population * 9.0;
                double food_room_left = max(0.0, max_food_needed - food_delivered[village_id]);
                double food_in_this_drop = drop.dry_food + drop.perishable_food;
                double effective_food_this_drop = min(food_in_this_drop, food_room_left);
                
                // Prioritize perishable food (like format_checker)
                double effective_vp = min((double)drop.perishable_food, effective_food_this_drop);
                double value_from_p = effective_vp * problem.packages[1].value;
                double remaining_effective_food = effective_food_this_drop - effective_vp;
                double effective_vd = min((double)drop.dry_food, remaining_effective_food);
                double value_from_d = effective_vd * problem.packages[0].value;
                
                total_value += value_from_p + value_from_d;

                // Other supplies logic
                double max_other_needed = village.population * 1.0;
                double other_room_left = max(0.0, max_other_needed - other_delivered[village_id]);
                double effective_vo = min((double)drop.other_supplies, other_room_left);
                total_value += effective_vo * problem.packages[2].value;

                // Update global delivery tracking
                food_delivered[village_id] += food_in_this_drop;
                other_delivered[village_id] += drop.other_supplies;

                // Track total drops for validation
                total_drop_dry += drop.dry_food;
                total_drop_peri += drop.perishable_food;
                total_drop_other += drop.other_supplies;
            }

            // Return to home city
            trip_distance += distance(current_location, problem.cities[home_city_idx - 1]);

            // Check constraints
            if (trip_distance > heli.distance_capacity + 1e-9) return INVALID_SCORE;
            if (total_drop_dry > trip.dry_food_pickup) return INVALID_SCORE;
            if (total_drop_peri > trip.perishable_food_pickup) return INVALID_SCORE;
            if (total_drop_other > trip.other_supplies_pickup) return INVALID_SCORE;

            helicopter_total_distances[heli_id] += trip_distance;

            // Trip cost - only charge if villages visited (like format_checker)
            double trip_cost = (trip.drops.size() > 0) ? 
                              (heli.fixed_cost + (heli.alpha * trip_distance)) : 0;
            total_trip_cost += trip_cost;
        }

        // Check DMax constraint
        if (helicopter_total_distances[heli_id] > problem.d_max + 1e-9) return INVALID_SCORE;
    }

    return total_value - total_trip_cost;
}
void optimizePackagePickup(Trip& trip, double weight_capacity, const ProblemData& problem) {
    // Priority: perishable > dry > other (based on values)
    double remaining_weight = weight_capacity;
    
    // Reset pickup amounts
    trip.dry_food_pickup = 0;
    trip.perishable_food_pickup = 0;
    trip.other_supplies_pickup = 0;
    
    // Calculate total demand from drops
    int total_dry_needed = 0, total_perishable_needed = 0, total_other_needed = 0;
    for (const auto& drop : trip.drops) {
        total_dry_needed += drop.dry_food;
        total_perishable_needed += drop.perishable_food;
        total_other_needed += drop.other_supplies;
    }
    
    // Pick up needed amounts, prioritizing by value/weight ratio
    double dry_ratio = problem.packages[0].value / problem.packages[0].weight;
    double perishable_ratio = problem.packages[1].value / problem.packages[1].weight;
    double other_ratio = problem.packages[2].value / problem.packages[2].weight;
    
    // Sort by value/weight ratio
    vector<pair<double, int>> ratios = {{dry_ratio, 0}, {perishable_ratio, 1}, {other_ratio, 2}};
    sort(ratios.rbegin(), ratios.rend());
    
    vector<int> needed = {total_dry_needed, total_perishable_needed, total_other_needed};
    vector<int*> pickup_ptrs = {&trip.dry_food_pickup, &trip.perishable_food_pickup, &trip.other_supplies_pickup};
    
    for (auto& ratio_pair : ratios) {
        int type = ratio_pair.second;
        double weight = problem.packages[type].weight;
        int max_possible = int(remaining_weight / weight);
        int to_pickup = min(max_possible, needed[type]);
        
        *pickup_ptrs[type] = to_pickup;
        remaining_weight -= to_pickup * weight;
    }
}

Solution randomNeighbor(const ProblemData& problem,Solution& current,vector<int> possible[]){
 Solution neighbor = current;
    
    if (neighbor.empty()) {
        return randomGenerator(problem,possible);
    }
    
    int action = randomInt(1, 8); // Expanded action space
    
    // Select random helicopter and trip
    int heli_idx = randomInt(0, neighbor.size() - 1);
    HelicopterPlan& plan = neighbor[heli_idx];
    
    if (plan.trips.empty()) {
        // Add a new trip
        Trip new_trip;
        int village_id = randomInt(1, problem.villages.size());
        
        Drop drop;
        drop.village_id = village_id;
        drop.dry_food = randomInt(1, 10);
        drop.perishable_food = randomInt(1, 10);
        drop.other_supplies = randomInt(1, 5);
        
        new_trip.drops.push_back(drop);
        optimizePackagePickup(new_trip, problem.helicopters[plan.helicopter_id - 1].weight_capacity, problem);
        
        plan.trips.push_back(new_trip);
        return neighbor;
    }
    
    int trip_idx = randomInt(0, plan.trips.size() - 1);
    Trip& trip = plan.trips[trip_idx];
    
    if (action == 1) {
        // Swap villages between two drops in the same trip
        if (trip.drops.size() > 1) {
            int i = randomInt(0, trip.drops.size() - 1);
            int j = randomInt(0, trip.drops.size() - 1);
            if (i != j) {
                swap(trip.drops[i].village_id, trip.drops[j].village_id);
            }
        }
    }
    else if (action == 2) {
        // Add a new village to existing trip
        if (trip.drops.size() < 4) { // Limit villages per trip
            Drop new_drop;
            new_drop.village_id = randomInt(1, problem.villages.size());
            new_drop.dry_food = randomInt(0, 20);
            new_drop.perishable_food = randomInt(0, 20);
            new_drop.other_supplies = randomInt(0, 10);
            trip.drops.push_back(new_drop);
            optimizePackagePickup(trip, problem.helicopters[plan.helicopter_id - 1].weight_capacity, problem);
        }
    }
    else if (action == 3) {
        // Remove a village from trip
        if (!trip.drops.empty()) {
            int drop_idx = randomInt(0, trip.drops.size() - 1);
            trip.drops.erase(trip.drops.begin() + drop_idx);
            if (!trip.drops.empty()) {
                optimizePackagePickup(trip, problem.helicopters[plan.helicopter_id - 1].weight_capacity, problem);
            }
        }
    }
    else if (action == 4) {
        // Modify drop amounts for a village
        if (!trip.drops.empty()) {
            int drop_idx = randomInt(0, trip.drops.size() - 1);
            Drop& drop = trip.drops[drop_idx];
            
            int change = randomInt(-5, 5);
            int food_type = randomInt(0, 2);
            
            if (food_type == 0) {
                drop.dry_food = max(0, drop.dry_food + change);
            } else if (food_type == 1) {
                drop.perishable_food = max(0, drop.perishable_food + change);
            } else {
                drop.other_supplies = max(0, drop.other_supplies + change);
            }
            
            optimizePackagePickup(trip, problem.helicopters[plan.helicopter_id - 1].weight_capacity, problem);
        }
    }
    else if (action == 5) {
        // Add a completely new trip
        Trip new_trip;
        
        // Select 1-2 random villages
        int num_villages = randomInt(1, 2);
        for (int i = 0; i < num_villages; i++) {
            Drop drop;
            drop.village_id = randomInt(1, problem.villages.size());
            
            // Calculate reasonable drop amounts based on village population
            const Village& village = problem.villages[drop.village_id - 1];
            int food_needed = village.population * 9;
            int other_needed = village.population;
            
            drop.dry_food = randomInt(0, food_needed);
            drop.perishable_food = randomInt(0, food_needed - drop.dry_food);
            drop.other_supplies = randomInt(0, other_needed);
            
            new_trip.drops.push_back(drop);
        }
        
        optimizePackagePickup(new_trip, problem.helicopters[plan.helicopter_id - 1].weight_capacity, problem);
        plan.trips.push_back(new_trip);
    }
    else if (action == 6) {
        // Remove a trip
        if (plan.trips.size() > 1) {
            plan.trips.erase(plan.trips.begin() + trip_idx);
        }
    }
    else if (action == 7) {
        // Optimize drop amounts for entire trip based on village needs
        if (!trip.drops.empty()) {
            for (auto& drop : trip.drops) {
                const Village& village = problem.villages[drop.village_id - 1];
                int food_needed = village.population * 9;
                int other_needed = village.population;
                
                // Redistribute current amounts more intelligently
                int total_food = drop.dry_food + drop.perishable_food;
                if (total_food > 0) {
                    // Prefer perishable food
                    double perishable_ratio = 0.8;
                    drop.perishable_food = min(int(total_food * perishable_ratio), food_needed);
                    drop.dry_food = min(total_food - drop.perishable_food, food_needed - drop.perishable_food);
                }
                drop.other_supplies = min(drop.other_supplies, other_needed);
            }
            optimizePackagePickup(trip, problem.helicopters[plan.helicopter_id - 1].weight_capacity, problem);
        }
    }
    else if (action == 8) {
        // Reorder villages in trip for better routing
        if (trip.drops.size() > 1) {
            random_shuffle(trip.drops.begin(), trip.drops.end());
        }
    }
    
    return neighbor;
    

}


Solution solve(const ProblemData& problem) {
    
    vector<int> possible[problem.helicopters.size()+1];

   possible_villages(problem,possible);
cout << "Starting improved solver with simulated annealing..." << endl;
    
    Solution best_solution;
    double best_score = -1e9;
    
    using namespace std::chrono;
    auto start = steady_clock::now();
    auto time_limit = duration<double>(problem.time_limit_minutes * 60); // Convert to seconds
    
    int restart_count = 0;
    const int max_restarts = 20;
    double time_per_restart = problem.time_limit_minutes * 60.0 / max_restarts;
    
    while (steady_clock::now() - start < time_limit && restart_count < max_restarts) {
        cout << "Restart " << restart_count + 1 << "/" << max_restarts << endl;
        
        // Generate initial solution
        Solution current = randomGenerator(problem,possible);
        double current_score = evaluateSolution(problem, current);
        
        Solution local_best = current;
        double local_best_score = current_score;
        
        // Simulated annealing parameters
        double initial_temp = 1000.0;
        double final_temp = 0.1;
        double cooling_rate = 0.99;
        double temperature = initial_temp;
        
        auto restart_start = steady_clock::now();
        auto restart_limit = duration<double>(time_per_restart);
        
        int iterations = 0;
        int improvements = 0;
        
        while (steady_clock::now() - restart_start < restart_limit) {
            Solution neighbor = randomNeighbor(problem, current,possible);
            double neighbor_score = evaluateSolution(problem, neighbor);
            
            // Accept if better or with probability based on temperature
            bool accept = false;
            if (neighbor_score > current_score) {
                accept = true;
                improvements++;
            } else if (temperature > final_temp) {
                double delta = neighbor_score - current_score;
                double probability = exp(delta / temperature);
                accept = (randomDouble(0.0, 1.0) < probability);
            }
            
            if (accept) {
                current = neighbor;
                current_score = neighbor_score;
                
                if (current_score > local_best_score) {
                    local_best = current;
                    local_best_score = current_score;
                }
            }
            
            // Cool down
            temperature *= cooling_rate;
            iterations++;
            
            // If no valid solution found for a while, restart early
            if (iterations > 1000 && current_score < -1e8) {
                break;
            }
        }
        
        cout << "Restart " << restart_count + 1 << " completed. Iterations: " << iterations 
             << ", Improvements: " << improvements << ", Score: " << local_best_score << endl;
        
        if (local_best_score > best_score) {
            best_solution = local_best;
            best_score = local_best_score;
            cout << "New best score: " << best_score << endl;
        }
        
        restart_count++;
    }
    
    // Final optimization pass
    if (best_score > -1e8) {
        cout << "Performing final optimization..." << endl;
        Solution current = best_solution;
        double current_score = best_score;
        
        auto final_start = steady_clock::now();
        auto final_limit = duration<double>(min(30.0, problem.time_limit_minutes * 60.0 * 0.1));
        
        while (steady_clock::now() - final_start < final_limit) {
            Solution neighbor = randomNeighbor(problem, current,possible);
            double neighbor_score = evaluateSolution(problem, neighbor);
            
            if (neighbor_score > current_score) {
                current = neighbor;
                current_score = neighbor_score;
                
                if (current_score > best_score) {
                    best_solution = current;
                    best_score = current_score;
                }
            }
        }
    }
    
    cout << "Solver finished. Best score = " << best_score << endl;
    return best_solution;
}








//random restart have not been applied properly

// Solution solve(const ProblemData& problem) {
//     cout << "Starting solver..." << endl;


//    Solution current=randomGenerator(problem);

//     double currentScore=evaluateSolution(problem,current);

//     Solution bestSolution=current;

//     double bestScore=currentScore;

//     // cout<<bestScore<<endl;

    
//     using namespace std::chrono;

//     auto start = steady_clock::now();
//     auto time_limit = duration<double>(problem.time_limit_minutes );

//     while (steady_clock::now() - start < time_limit){
       
//         Solution neighbor = randomNeighbor(problem,current);
//         double neighborScore=evaluateSolution(problem,neighbor);

        
//         if(neighborScore >= bestScore){
//             current=neighbor;
            
//             bestScore = neighborScore;
            
//             bestSolution=neighbor;
            

//         }
       

//         if(bestScore < -1e8){
//             current = randomGenerator(problem);
//             bestScore = evaluateSolution(problem,current);
//         }



//     }




//     cout << "Solver finished. Best score = " << bestScore << endl;
//     return bestSolution;
// }
