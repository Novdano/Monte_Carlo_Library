#include <math.h>
#include <random>
#include <iostream>
#include <fstream>
#include "Monte_Carlo_Simulation.h"
#include "boost/math/distributions/normal.hpp"

double evolve_stock_BS(double S, double sigma, double r, double T, int N) {
	std::random_device rd{};
	std::mt19937 gen{ rd() };
	std::normal_distribution<double> d{ 0, 1 };
	double delta = T / N;
	double res = S + S*r*delta + S* sigma * pow(delta, 0.5) * d(gen);
	return res;
}

double evolve_vol_BS(double S, double sigma, double r, double T, int N) {
	return sigma;
}

double evolve_stock_SV(double S, double sigma, double r, double T, int N) {
	std::random_device rd{};
	std::mt19937 gen{ rd() };
	std::normal_distribution<double> d{ 0, 1 };
	double delta = T / N;
	double res = S + S * r*delta + S * pow(sigma,0.5) * pow(delta, 0.5) * d(gen);
	return res;
}

double evolve_vol_SV(double S, double sigma, double r, double T, int N) {
	std::random_device rd{};
	std::mt19937 gen{ rd() };
	std::normal_distribution<double> d{ 0, 1 };
	sigma = sigma + 1 * (0.004 - sigma) * (T / N) + 0.001*d(gen) * pow(sigma, 0.5) * pow(T / N, 0.5);
	return sigma;
}

double payoff_call(double S, double K) {
	if (S >= K) {
		return S - K;
	}
	return 0;
}
double BS_Call(double S_0, double K, double r, double sigma, double T) {
	boost::math::normal N(0, 1);
	double d1 = (std::log(S_0 / K) + ((r + (0.5 * std::pow(sigma, 2))) * T)) / (sigma * std::pow(T, 0.5));
	double d2 = (std::log(S_0 / K) + ((r - (0.5 * std::pow(sigma, 2))) * T)) / (sigma * std::pow(T, 0.5));
	return S_0 * cdf(N, d1) - cdf(N, d2) * K * std::exp(-r * T);
}

int main() {
	ofstream file;
	ofstream runtime_file;
	file.open("BSvsSV.dat");
	runtime_file.open("runtime.txt");
	applySimulations sim_BS = applySimulations(50, 30, 0.175, 0.02, 1, 10000, 10000, &evolve_stock_BS, &evolve_vol_BS, &payoff_call);
	applySimulations sim_SV = applySimulations(50, 30, 0.175, 0.02, 1, 10000, 10000, &evolve_stock_SV, &evolve_vol_SV, &payoff_call);
	Monte_Carlo_Simulation MCBS = Monte_Carlo_Simulation(50, 30, 0.0175, 0.02, 1, sim_BS, 10000, 10000);
	Monte_Carlo_Simulation MCSV = Monte_Carlo_Simulation(50, 30, 0.0175, 0.02, 1, sim_SV, 10000, 10000);
	file << "Strike BS_Price SV_Price" << endl;
	const clock_t start = clock();
	for (int K = 30; K <= 50; K++) {
		cout << K << endl;
		double BSCall = BS_Call(50, K, 0.0175, 0.02, 1);
		MCBS.reset_input(50, K, 0.0175, 0.02, 1);
		MCSV.reset_input(50, K, 0.0175, 0.02, 1);
		double BS = MCBS.run_simulations(1);
		double SV = MCSV.run_simulations(1);
		file << K << " " << BS << " " << SV << " "<<BSCall<<endl;
	}
	const clock_t end = clock();
	runtime_file << float(end - start)/CLOCKS_PER_SEC << endl;
	runtime_file.close();
	file.close();
	return 0;
}