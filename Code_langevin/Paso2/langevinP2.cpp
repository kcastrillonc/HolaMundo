// langevin.cpp :
#include <functional>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fstream>
#include <string>
#include <sstream>

#ifndef RAND_LF_H_
#define RAND_LF_H_

class particle {
private:
public:
    double R;
    double x, y;
    double c_x, c_y; //For oscillation forces
    double vx, vy;
    double Dvx, Dvy;
    double Fx, Fy;
    double W;

    particle(double Ri, double xi, double c_xi, double yi, double c_yi, double vxi, double vyi,
        double Fxi, double Fyi)
    {
        R = Ri; x = xi; c_x = c_xi; y = yi; c_y = c_yi; vx = vxi; vy = vyi; Fx = Fxi; Fy = Fyi; Dvx = 0; Dvy = 0; W = 0;
    }
    void reset_F() { Fx = Fy = 0; }
    void add_F(double Fx0, double Fy0);
};

void particle::add_F(double Fx0, double Fy0) {
    Fx += Fx0;
    Fy += Fy0;
}

class lf_integrator {
public:
    double dt;
    double f, T;
    bool bonded;
    lf_integrator(double dti, double fi, double Ti)
    {
        dt = dti; f = fi; T = Ti; bonded = false;
    }
    void initialize_v(std::vector<particle>& particles);
    void update_F(std::vector<particle>& particles, double k);
    void update_v1(std::vector<particle>& particles);
    void impulse_Dv(std::vector<particle>& particles, gsl_rng* r);
    void update_x(std::vector<particle>& particles, double k);
    void update_v2(std::vector<particle>& particles);
    void check_bonded(std::vector<particle>& particles, double k, double thres_bond, double thres_unbond);
    double x_avg(std::vector<particle>& particles);
    double x_std(std::vector<particle>& particles);
    double vx_avg(std::vector<particle>& particles);
    double vx_std(std::vector<particle>& particles);
};

void lf_integrator::initialize_v(std::vector<particle>& particles) {
    for (auto& p : particles) {
        p.vx += 0.5 * dt * p.Fx;
        p.vy += 0.5 * dt * p.Fy;
    }
}

void lf_integrator::update_F(std::vector<particle>& particles, double k) {
    double F;
    for (auto& p : particles) {
        p.reset_F();
        F = -k * (p.x - p.c_x);
        p.add_F(F, 0);
    }
}

void lf_integrator::update_v1(std::vector<particle>& particles) {
    for (auto& p : particles) {
        p.vx += dt * p.Fx;
        p.vy += dt * p.Fy;
    }
}

void lf_integrator::impulse_Dv(std::vector<particle>& particles, gsl_rng* r) {
    double eta;
    for (auto& p : particles) {
        eta = gsl_ran_gaussian(r, 1.);
        p.Dvx = -f * p.vx + sqrt(f * (2 - f) * T) * eta;
    }
}

void lf_integrator::update_x(std::vector<particle>& particles, double k) {
    double dx, dy;
    for (auto& p : particles) {
        dx = dt * (p.vx + 0.5 * p.Dvx);
        dy = dt * (p.vy + 0.5 * p.Dvy);
        p.x += dx;
        p.y += dy;
        p.W += -k * (p.x - p.c_x) * dx;
    }
}

void lf_integrator::update_v2(std::vector<particle>& particles) {
    for (auto& p : particles) {
        p.vx += p.Dvx;
        p.vy += p.Dvy;
    }
}

void lf_integrator::check_bonded(std::vector<particle>& particles, double k, double thres_bond, double thres_unbond) {
    double dist = particles[1].x - particles[0].x;
    double F;
    if (fabs(dist) <= thres_bond) {
        bonded = true;
        F = -k * (dist - thres_bond);
        particles[1].add_F(F, 0);
        particles[0].add_F(-F, 0);
    }
    else if (fabs(dist) <= thres_unbond && bonded == true) {
        F = -k * (dist - thres_bond);
        particles[1].add_F(F, 0);
        particles[0].add_F(-F, 0);
    }
    else if (fabs(dist) > thres_unbond && bonded == true) {
        bonded = false;
    }
}

double lf_integrator::x_avg(std::vector<particle>& particles) {
    double avg = 0;
    for (auto& p : particles) {
        avg += p.x;
    }
    avg /= particles.size();
    return avg;
}

double lf_integrator::x_std(std::vector<particle>& particles) {
    double avg = this->x_avg(particles);
    double std = 0;
    for (auto& p : particles) {
        std += pow((p.x - avg), 2);
    }
    std /= particles.size() - 1;
    std = sqrt(std);
    return std;
}

double lf_integrator::vx_avg(std::vector<particle>& particles) {
    double avg = 0;
    for (auto& p : particles) {
        avg += p.vx;
    }
    avg /= particles.size();
    return avg;
}

double lf_integrator::vx_std(std::vector<particle>& particles) {
    double avg = this->vx_avg(particles);
    double std = 0;
    for (auto& p : particles) {
        std += pow((p.vx - avg), 2);
    }
    std /= particles.size() - 1;
    std = sqrt(std);
    return std;
}

#endif // RAND_LF_H_

double k1_ini = 0.05, k1_fin = 0.2, k2_fin = 1.0;
int num = 20, t_end = 2000;

void ciclo(std::vector<particle>& particles, lf_integrator osc, gsl_rng* r, int teq);

template <typename T>
std::string ToString(T val)
{
    std::stringstream stream;
    stream << val;
    return stream.str();
}

std::ofstream outfile;

int main() {
    int SEED = 0;
    gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, SEED);

    particle pini(0.0, /*x*/ 0.0, 0.0, 0.0, 0.0, /*vx*/ 0.0, 0.0, 0.0, 0.0);

    std::vector<particle> bs(10000,pini); /*1000 partículas con
                                           condiciones iniciales pini*/
    //prueba(bs);

    double gamma = 2.1, dt = 1/*/(10*gamma)*/, t_end = 1000/**dt*/;
    double k =0.05, f=1-exp(-dt*gamma), T=0.5; /*kboltzmann=1, masa=1*/

    /*outfile << "k=" << " " << k << "\n"
            << "gamma=" << " " << gamma << "\n"
            << "T=" << " " << T << "\n";*/

    lf_integrator osc(dt, f, T);
    
    osc.update_F(bs, k);
    osc.initialize_v(bs);

    /*std::string ss = "datosP1_k"+ToString(k)+
                        "_gamma"+ToString(gamma)+
                            "_T"+ToString(T)+".dat";

    outfile.open(ss.c_str());
    outfile.close();*/

    ciclo(bs, osc, r, t_end/40);

    gsl_rng_free(r);

    std::cout << "Esta cosa compila." << "\n";

    return 0;
}

void ciclo(std::vector<particle>& particles, lf_integrator osc, gsl_rng* r, int teq){
    outfile.open("parametros.dat");

    /*Compresión isotérmica*/
    //outfile << "Etapa 1" << "\n";

    double k1 = k1_ini;
    double k1_add = (k1_fin - k1_ini)/num;

    for (int t = 1; t <= t_end/4; t ++) {
        if(t%teq == 0){
            k1 = k1 + k1_add;

            outfile << k1 << " " << osc.T << "\n";
        }
        osc.update_F(particles, k1);

        osc.update_v1(particles);
        osc.impulse_Dv(particles, r);
        osc.update_x(particles, k1);
        osc.update_v2(particles);
    }

    double alpha = (osc.T * osc.T)/k1;  // Cantidad proporcional a la entropía

    /*Compresión adiabática*/
    //outfile << "Etapa 2" << "\n";

    double k2 = k1;
    double k2_add = (k2_fin - k1_fin)/num;

    for (int t = 1; t <= t_end/4; t ++) {
        if(t%teq == 0){
            k2 = k2 + k2_add;
            osc.T = sqrt(alpha*k2);

            outfile << k2 << " " << osc.T << "\n";
        }
        osc.update_F(particles, k2);

        osc.update_v1(particles);
        osc.impulse_Dv(particles, r);
        osc.update_x(particles, k2);
        osc.update_v2(particles);
    }

    /*Expansión isotérmica*/
    //outfile << "Etapa 3" << "\n";
    double k3_fin = (osc.T * osc.T)/5;

    double k3 = k2;
    double k3_add = (k3_fin - k2_fin)/num;

    for (int t = 1; t <= t_end/4; t ++) {
        if(t%teq == 0){
            k3 = k3 + k3_add;

            outfile << k3 << " " << osc.T << "\n";
        }
        osc.update_F(particles, k3);

        osc.update_v1(particles);
        osc.impulse_Dv(particles, r);
        osc.update_x(particles, k3);
        osc.update_v2(particles);
    }

    alpha = (osc.T * osc.T)/k3;  // Cantidad proporcional a la entropía

    /*Expansión adiabática*/
    //outfile << "Etapa 4" << "\n";

    double k4 = k3;
    double k4_add = (k1_ini - k3_fin)/num;

    for (int t = 1; t <= t_end/4; t ++) {
        if(t%teq == 0){
            k4 = k4 + k4_add;
            osc.T = sqrt(alpha*k4);

            outfile << k4 << " " << osc.T << "\n";
        }
        osc.update_F(particles, k4);

        osc.update_v1(particles);
        osc.impulse_Dv(particles, r);
        osc.update_x(particles, k4);
        osc.update_v2(particles);
    }

    outfile.close();
}
