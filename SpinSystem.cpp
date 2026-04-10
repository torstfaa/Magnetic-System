#include "SpinSystem.h"

vector<double> operator+(const vector<double>& v1, const vector<double>& v2) {
	vector<double> v3(v1.size());
	for (int i = 0; i < v1.size(); i++) {
		v3[i] = v1[i] + v2[i];
	}

	return v3;
};

vector<double> operator-(const vector<double>& v1, const vector<double>& v2) {
	vector<double> v3(v1.size());
	for (int i = 0; i < v1.size(); i++) {
		v3[i] = v1[i] - v2[i];
	}

	return v3;
};

vector<double> operator*(const double& c, const vector<double>& v) {
	vector<double> u(v.size());
	for (int i = 0; i < v.size(); i++) {
		u[i] = c*v[i];
	}
	return u;
}

vector<double> operator/(const vector<double>& v, const double& c) {
	vector<double> u(v.size());
	for (int i = 0; i < v.size(); i++) {
		u[i] = v[i]/c;
	}
	return u;
}

double operator*(const vector<double>& v1, const vector<double>& v2) {
    double return_value = 0;
    for (int i = 0; i < v1.size(); i++) {
        return_value += v1[i]*v2[i];
    }
    return return_value;
}

vector<double> cross(const vector<double>& v1, const vector<double>& v2) {
    vector<double> return_vec(3);
    return_vec[0] = v1[1]*v2[2] - v1[2]*v2[1];
    return_vec[1] = v1[2]*v2[0] - v1[0]*v2[2];
    return_vec[2] = v1[0]*v2[1] - v1[1]*v2[0];
    return return_vec;
}

void System1D::add_spin(vector<double> spin) {
    spins.push_back(spin);
}

inline vector<vector<double>> System1D::get_spins() {
    return spins;
}

vector<double> System1D::coupling_force(int i) {
    vector<double> force = {0.0, 0.0, 0.0};
    if ((i != 0) && (i != spins.size()-1)) {
        vector<double> force = coupling*(spins[i-1] + spins[i+1])/mu_B;
    } else if (periodic) {
        if (i == 0) {
            vector<double> force = coupling*(spins[1] + spins[spins.size()-1])/mu_B;
        } else {
            vector<double> force = coupling*(spins[0] + spins[spins.size()-2])/mu_B;
        }
    } else if (i == 0) {
        vector<double> force = coupling*spins[1]/mu_B;
    } else {
        vector<double> force = coupling*spins[spins.size()-2]/mu_B;
    }
    return force;
}

inline vector<double> System1D::anisotropy_force(int i) {
    return vector<double>{0, 0, 2*d_z*spins[i][2]};
}

inline vector<double> System1D::mag_force() {
    return vector<double>{0, 0, b_0};
}

vector<double> System1D::thermal_force(double dt, vector<double> random_vector) {
    double thermal_constant = sqrt((2*alpha*k_B*temperature)/(gamma*mu_B*dt));
    default_random_engine generator; //Move elsewere?
    normal_distribution<double> distribution_x(0, 1.0);
    normal_distribution<double> distribution_y(0, 1.0);
    normal_distribution<double> distribution_z(0, 1.0);
    vector<double> thermal_force = {distribution_x(generator), distribution_y(generator), 
        distribution_z(generator)};
    return thermal_constant*thermal_force;
}

vector<double> System1D::total_force(int i, double dt) {
    vector<double> return_vector = coupling_force(i) + mag_force();
}

vector<double> System1D::time_evolution(vector<double> spin) {
    vector<double> return_vector = (-gamma/(1 + pow(alpha, 2)))*(cross(spin, mag_force()) + 
    alpha*cross(spin, cross(spin, mag_force()))); //temporary
}

void System1D::evolve_Heun(double dt) {
    for (int i = 0; i < spins.size(); i++) {
        vector<double> intermediate_spin = spins[i] + dt*time_evolution(spins[i]);
        spins[i] = spins[i] + dt*(time_evolution(spins[i]) + time_evolution(intermediate_spin))/2;
    }
}