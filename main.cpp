#include<fstream>
#include "SpinSystem.h"
//#include <armadillo>


int main() {
	System1D single_spin_system(1);
	single_spin_system.add_spin(vector<double>{0.95, sqrt(0.0975), 0});
	double dt = 1e-15;
	int counts = 100;
	vector<vector<double>> spin_vector(counts, vector<double>(3));
	spin_vector[0] = single_spin_system.get_spins()[0];
	for (int i = 1; i < counts; i++) {
		single_spin_system.evolve_Heun(dt);
		spin_vector[i] = single_spin_system.get_spins()[i];
	}

	ofstream fileStream("output.txt");

	for (int i = 0; i < spin_vector.size(); i++) {
		fileStream << spin_vector[i][0] << '\t' << spin_vector[i][1] << '\t' << spin_vector[i][2] << '\n';
	}

	return 0;
}

