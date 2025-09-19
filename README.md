# Disaster Relief Helicopter Routing

This repository formulates and solves a disaster-relief helicopter routing problem as a search/optimization task. The objective is to maximize the delivered aid value (prioritizing perishable food) while minimizing logistical trip costs, subject to helicopter weight, per-trip distance, and total-distance (DMax) constraints.

**What I implemented:** core optimization routine (simulated annealing with random-restart), neighborhood operators, and a solution evaluator that enforces all constraints and computes objective = delivered value − trip costs.

# Disaster Relief Helicopter Routing (AI Assignment)

## Overview
During floods, timely delivery of food and supplies is critical. This project models the delivery problem as a combinatorial search/optimization problem:
- Villages and cities are points in 2D; villages have populations (demand).
- Packages: dry food, perishable food, other supplies (weights & values).
- Helicopters: home city, weight capacity, per-trip distance capacity, and total distance cap (DMax).
- Each trip starts and ends at the helicopter's home city. Each trip costs `F + α * distance`.
- Objective: maximize `Total Value = delivered_value − total_trip_cost`, while respecting capacities and distance constraints.

This repo contains a working solver that searches for helicopter trip plans and package allocations using **simulated annealing with random restarts** and a set of local neighborhood moves.

## Key files
- `solver.h` / `solver.cpp` (or `solver.cc`) — main solver implementation (core functions described below).
- `README.md` — this file.
- (Runner / grader files provided by course) — call `solve(problem)` to obtain a `Solution`.

## My contribution (what I implemented)
I implemented the core optimization module (C++):
- `solve(const ProblemData& problem)` — Simulated Annealing + random restarts driver, time-limited.
- `randomGenerator(...)` — Creates an initial candidate solution (feasible-first, randomized).
- `randomNeighbor(...)` — Extensive neighborhood operator set (add/remove trips, add/remove villages in a trip, reorder drops, modify drop amounts).
- `evaluateSolution(...)` — Validates constraints and computes objective (delivered value capped by village needs; perishable prioritized; trip costs).
- `optimizePackagePickup(...)` — Heuristic to set package pickup counts per trip based on drop needs and value/weight ratio.
- `possible_villages(...)` — Quickly prunes villages unreachable by per-trip helicopter range (2× distance to/from home city).

## Algorithm / Heuristics
- Simulated annealing with temperature decay and restart strategy (multiple restarts to escape local optima).
- Neighborhood moves include: swap villages, add/remove village drops, change drop amounts, add/remove entire trips, reorder drop sequence, and a per-trip reoptimization of pickup amounts.
- Objective respects real constraints: trip weight, per-trip distance, total helicopter distance (DMax), and value capping (max food = 9 × population; other = 1 × population).

## How to run
This code is designed to plug into the assignment runner that constructs `ProblemData` and expects `solve(problem)` to return a `Solution`. Typical flow:
1. Build with your course's build system (or `g++` / `cmake` if provided).
2. Run the provided grader or input driver which will call `solve`.
3. The solver respects `problem.time_limit_minutes` — it stops when the time budget elapses.

> Note: If you want a CLI harness, add a small `main()` that reads JSON/problem file, constructs `ProblemData`, calls `solve()` and prints the `Solution`.

## Results & notes
- The solver produces feasible solutions and improves them by repeated restarts and local search.
- The implementation prioritizes perishable food when assigning pickup amounts and uses value/weight ratio to allocate limited helicopter capacity.
- Current implementation uses one trip per helicopter by default in the initial generator; simulated annealing evolves the solution to add/remove trips.

## Suggested improvements (future work)
- Implement a proper initial solution that uses a greedy cluster-by-proximity (not purely random).
- Integrate a TSP solver or greedy nearest-neighbor to order villages better within each trip.
- Consider multi-trip planning per helicopter and explicit load-splitting across trips.
- Track and log best objective per iteration so you can plot convergence.
- Add unit tests for `evaluateSolution` and small fixed-instance sanity checks.

## License
Add the license you prefer (MIT recommended for student projects).
