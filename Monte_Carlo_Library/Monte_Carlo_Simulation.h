#pragma once

#include <math.h>
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"

using namespace std;
using namespace tbb;


class applySimulations
{
private:
	double S_0;
	double K;
	double r;
	double sigma;
	double T;
	int M;
	int N;
	double(*evolve_under)(double, double, double, double, int);
	double(*evolve_vol)(double, double, double, double, int);
	double(*payoff_fn)(double, double);

public:
	applySimulations(double init, double strike, double rate, double vol, double mat, int num_sim, int num_steps,
		double(*under_change)(double, double, double, double, int), double(*vol_change)(double, double, double, double, int),
		double(*calc_payoff)(double, double));
	void reset(double init, double strike, double rate, double vol, double mat, int num_sim, int num_steps);
	void operator() (const tbb::blocked_range<int>& q) const;
};


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
	double(*evolve_under)(double, double, double, double, int);
	double(*evolve_vol)(double, double, double, double, int);
	double(*payoff_fn)(double, double);
	applySimulations par_app;


public:
	Monte_Carlo_Simulation(double init, double strike, double rate, double vol, double mat, applySimulations simul, int num_sim, int num_steps);
	double run_simulations(int run_type);
	void reset_input(double init, double strike, double rate, double vol, double mat);
};


