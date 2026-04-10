#include<vector>
#include<cmath>
#include<random>
#include<iostream>

using namespace std;

constexpr double pi = 3.14159265;
constexpr double mu_B = 9.274e-24;
constexpr double gamma = 1.76e11;
constexpr double mu_0 = pi*4e-7;
constexpr double k_B = 1.3807e-23;

vector<double> operator+(const vector<double>& v1, const vector<double>& v2);
vector<double> operator-(const vector<double>& v1, const vector<double>& v2);
vector<double> operator*(const double& c, const vector<double>& v);
vector<double> operator/(const vector<double>& v, const double& c);
double operator*(const vector<double>& v1, const vector<double>& v2);

vector<double> cross(const vector<double>& v1, const vector<double>& v2);

class System1D {
    private:
        double alpha;
        double b_0;
        double d_z;
        double coupling;
        double temperature;
        bool periodic;
        vector<vector<double>> spins;
    public:
        void add_spin(vector<double> spin);
        vector<vector<double>> get_spins();
        vector<double> coupling_force(int i);
        vector<double> anisotropy_force(int i);
        vector<double> mag_force();
        vector<double> thermal_force(double dt, vector<double> random_vector);
        vector<double> total_force(int i, double dt);
        vector<double> time_evolution(vector<double> spin);
        void evolve_Heun(double dt);
        System1D(double alpha = 0, double b_0 = 0, double d_z = 0, double coupling = 0, 
        double temperature = 0, bool periodic = false);
};

class System2D {
    private:
        double alpha;
        double b_0;
        double d_z;
};