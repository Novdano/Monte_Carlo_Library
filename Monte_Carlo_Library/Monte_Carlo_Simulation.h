#pragma once

#include <math.h>;
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"

using namespace std;
using namespace tbb;

#define PARALLEL_RUN = 1;
#define SEQUENTIAL_RUN = 0;

class Monte_Carlo_Simulation
{
private:
	double S_0;
	double K;
	double r;
	double sigma;
	double T;
	int M;
	int N;


public:
	Monte_Carlo_Simulation(double init, double rate, double vol, double mat, int num_sim, int num_steps) {};
	virtual double underlying_evolution(double price) {};
	virtual double volatility_evolution(double vol) {};
	virtual double payoff_function(double S) {};
	double run_simulations(int run_type) {}
};


