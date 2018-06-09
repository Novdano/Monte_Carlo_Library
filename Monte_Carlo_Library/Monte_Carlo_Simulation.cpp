#include "Monte_Carlo_Simulation.h"

#include <math.h>
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/cache_aligned_allocator.h"


class applySimulations {
private:
	double S_0;
	double r;
	double sigma;
	double T;
	double *result;
	int M;
	int N;
	double (Monte_Carlo_Simulation::*evolve_stock)(double);
	double (Monte_Carlo_Simulation::*evolve_vol)(double);
	double (Monte_Carlo_Simulation::*payoff_fn)(double);

public:
	applySimulations(double *res_cont,  double init, double rate, double vol, double mat, int num_sim, int num_steps, 
		double (Monte_Carlo_Simulation::)(double), double (Monte_Carlo_Simulation::*vol_change)(double), double (Monte_Carlo_Simulation::*calc_payoff)(double)) :
		result(res_cont), S_0(init), r(rate), sigma(vol), T(mat),M(num_sim), N(num_steps), 
		evolve_stock(stock_change), evolve_vol(vol_change), payoff_fn(calc_payoff) {}

	void operator() (const tbb::blocked_range<int>&q) const {
		tbb::atomic<double> total_payoff = 0.0;
		for (int i = q.begin(); i < q.end(); i++) {
			double S = S_0;
			double vol = sigma;
			for (int n = 0; i < N; i++) {
				S = (*Monte_Carlo_Simulation::underlying_evolution)(S);
				vol = (Monte_Carlo_Simulation::*evolve_stock)(vol);
			}
			double payoff = (Monte_Carlo_Simulation::*payoff_fn)(S);
			while (total_payoff.compare_and_swap(total_payoff + payoff, total_payoff) != total_payoff) {}
		}
		*result = (total_payoff / M) * std::pow((1 +  r * (T / N)), -N);
	}



};


Monte_Carlo_Simulation::Monte_Carlo_Simulation(double init, double rate, double vol,double mat, int num_sim, int num_steps)
{
	S_0 = init;
	r = rate;
	sigma = vol;
	r = rate;
	T = mat;
	M = num_sim;
	N = num_steps;
}

double Monte_Carlo_Simulation::run_simulations(int run_type) {
	tbb::atomic<double> total_payoff = 0;
	if (run_type == SEQUENTIAL_RUN) {
		for (int i = 0, i < M; i++) {
			double S = S_0;
			double vol = sigma;
			for (int n = 0; j < N; j++) {
				S = underlying_evolution(S, r, vol, T);
				vol = volatility_evolution(S, r, vol, T);
			}
			total_payoff = total_payoff+ payoff_function(S);
		}
		return (total_payoff / M) * pow((1 + r * (T / N)), N);
	}
	else if (run_type == PARALLEL_RUN) {
		double *result = (double*) malloc(sizeof(double));
		parallel_for(blocked_range<int>(0, M), applySimulations(result, S_0, r, sigma, T, M, N, 
			&underlying_evolution, &volatility_evolution, &payoff_function));
	}
}

