#include "Monte_Carlo_Simulation.h"

#include <math.h>
#include <iostream>
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/cache_aligned_allocator.h"
#include "tbb/mutex.h"

tbb::mutex payoff_mutex;
tbb::atomic<double> total_payoff = 0;

double average(double *list, int size) {
	double res = 0;
	for (int i=0; i < size; i++) {
		res = res + list[i];
	}
	return res / size;
}


applySimulations::applySimulations(double init, double strike, double rate, double vol, double mat, int num_sim, int num_steps,
	double(*under_change)(double, double, double, double, int), double(*vol_change)(double, double, double, double, int),
	double(*calc_payoff)(double, double)) {
	S_0 = init; 
	K = strike;
	r = rate;
	sigma = vol; 
	T = mat;
	M = num_sim;
	N = num_steps;
	evolve_under = under_change;
	evolve_vol = vol_change;
	payoff_fn = calc_payoff;
}

void applySimulations::reset(double init, double strike, double rate, double vol, double mat, int num_sim, int num_steps) {
	S_0 = init;
	K = strike;
	r = rate;
	sigma = vol;
	T = mat;
	M = num_sim;
	N = num_steps;
}

void applySimulations::operator() (const tbb::blocked_range<int>& q) const {
	for (int i = q.begin(); i != q.end(); ++i) {
		double S = S_0;
		double vol = sigma;
		for (int n = 0; n < N; n++) {
			S = (*evolve_under)(S, vol, r, T, N);
			vol = (*evolve_vol)(S, vol, r, T, N);
		}
		double payoff = (*payoff_fn)(S, K);
		total_payoff = total_payoff + payoff;
	}
}

Monte_Carlo_Simulation::Monte_Carlo_Simulation(double init, double strike, double rate, double vol,double mat, applySimulations simul, int num_sim, int num_steps):
	par_app(simul)
{
	S_0 = init;
	K = strike;
	r = rate;
	sigma = vol;
	T = mat;
	M = num_sim;
	N = num_steps;
}




double Monte_Carlo_Simulation::run_simulations(int run_type) {
	total_payoff = 0;
	if (run_type == 0) {
		//cout << S_0 <<" "<< sigma<<" " << K <<" "<< r <<" "<< T <<" "<< endl;
		for (int i = 0; i < M; i++) {
			double S = S_0;
			double vol = sigma;
			for (int n = 0; n < N; n++) {
				S = (*evolve_under)(S, vol, r, T, N);
				vol = (*evolve_vol)(S, vol, r, T, N);
			}
			total_payoff = total_payoff+ (*payoff_fn)(S,K);
			/*cout << S << K<<endl;
			cout << (*payoff_fn)(S, K) << endl;
			Sleep(100);*/
		}
	}
	else if (run_type == 1) {
		par_app.reset(S_0, K, r, sigma, T, M, N);
		parallel_for(blocked_range<int>(0, M), par_app);
	}
	double delta_t = T / N;
	/*cout << total_payoff / M << endl;
	cout << (total_payoff / M) * pow((1 + r * delta_t), -N) << endl;
	system("pause");*/
	return (total_payoff / M) * pow((1 + r * delta_t), -N);
}

void Monte_Carlo_Simulation::reset_input(double init, double strike, double rate, double vol, double mat)
{
	S_0 = init;
	K = strike;
	r = rate;
	sigma = vol;
	T = mat;
}
