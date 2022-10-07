#include <iostream>
#include <cmath>
#include <string>
#include <tgmath.h>
#include <fstream>
#include <vector>

using namespace std;

double pi = 3.14159;
double hbar = 6.626 * pow(10, -34) / (2 * pi);
double e = 1.602 * pow(10.0, -19);
double m = 9.11 * pow(10, -31);
double m_perp = .198 * m;
double m_z = .92 * m;
double omega0 = .008 * e / hbar; //8 meV
double E_VS = .33 * e / 1000; // valley splitting in meV, for 21 use .75 meV, for 31 use .33 meV
double g = 1.998;
double r = 1.1 * pow(10.0, -9); //dipole size, 1.1 nm
double rho = 2200; //silicon mass density 2200 kg/m^3 in SiO2 (Peihao's paper)
double v_t = 3750; //3750, 5420
double v_l = 5900; //5900, 9330
double T = .15;
double R = 2000.0;
double R_k = 26000.0;
double l0 = 100.0 * pow(10.0, -9);

double Rashba = 45;
double Dressel = 0;

double xi_d = 5.0 * e; //dilation deformation 5 eV
double xi_u = 8.77 * e; //uniaxial sheer deformation 8.77 eV



double Zeeman_Energy(double B) {
    double ub = e * hbar / (2 * m);
    //cout << g * ub * B / hbar << endl;
    return g * ub * B;
}

double sin_gamma_sq(double B) {
    double e3 = (E_VS - Zeeman_Energy(B));
    double delta23 = r * m_perp * (E_VS) * (Dressel + Rashba) / (sqrt(2) * hbar);
    double e3_tilde = sqrt(e3 * e3 + delta23 * delta23);
    double gamma = atan(delta23 / e3_tilde);
    return (sqrt(e3 * e3 + delta23 * delta23) + e3) / (2 * sqrt(e3 * e3 + delta23 * delta23));

}

double cos_gamma_sq(double B) {
    double e3 = (E_VS - Zeeman_Energy(B));
    double delta23 = r * m_perp * (E_VS) * (Dressel + Rashba) / (sqrt(2) * hbar);
    double e3_tilde = sqrt(e3 * e3 + delta23 * delta23);
    double gamma = atan(delta23 / e3_tilde);
    return (sqrt(e3 * e3 + delta23 * delta23) - e3) / (2 * sqrt(e3 * e3 + delta23 * delta23));

}

double Single_Electron_Energy(double B, int spin, int valley) {
    //double omega_c = e * B / m_perp;
    //double Omega = sqrt(.25 * omega_c * omega_c + omega0 * omega0);
    //double E = hbar * Omega * 1000 / e;
    double E = valley * E_VS + spin * Zeeman_Energy(B);
    return E;
}

double F_SV(double B) { //maybe SV should involve valleys??
    //double delta = 2 * m_perp * E_VS * pow(Dressel + Rashba, 2.0) * r / hbar;
    double delta = r* m_perp* (E_VS)* (Dressel + Rashba) / (sqrt(2) * hbar);
    return 1 - 1.0 / sqrt(1 + delta * delta / pow(E_VS - Zeeman_Energy(B), 2.0));

}

// order: spin1, valley1 are for lower energy state
          //spin2, valley2 are for upper energy state
double spin_valley(double B, int spin1, int valley1, int spin2, int valley2) {
    double Energy_diff = Single_Electron_Energy(B, spin2, valley2) - Single_Electron_Energy(B, spin1, valley1);
    double omega_z = Zeeman_Energy(B) / hbar;
    double pre = 4.0 * pi * (R / R_k) * omega_z * F_SV(B) * r * r *pow(l0, -2.0);
    return pre/ tanh(hbar * omega_z / (2 * 1.38 * pow(10, -23) * T));

}


double intravalley_SO(double B) {
    double omega_z = Zeeman_Energy(B) / hbar;
    double pre = 4.0 *  pi * (R/R_k) * pow(Dressel + Rashba, 2.0) * pow(omega_z, 3.0) * pow(omega0, -4.0) * pow(l0, -2.0);
    return pre / tanh(hbar * omega_z / (2 * 1.38 * pow(10, -23) * T));
}


int main() {
    double data_points = 120;
    double beginning = 1.0;
    double ending = 10.0;
    int spin1 = -1.0;
    int valley1 = -1.0;
    int spin2 = 1.0;
    int valley2 = -1.0;
    int spin3 = -1.0;
    int valley3 = 1.0;

    ofstream myfile;
    myfile.open("relaxation.txt");
    for (double i = beginning * data_points; i <= ending * data_points; i++) {
        double B = i / data_points;
        myfile << B << " ";
        if (E_VS >= Zeeman_Energy(B)) {
            myfile << spin_valley(B, spin1, valley1, spin2, valley2) + intravalley_SO(B) << " ";
            myfile << spin_valley(B, spin1, valley1, spin2, valley2) << " ";
            myfile << intravalley_SO(B) << endl;

        }
        else {
            myfile << spin_valley(B, spin1, valley1, spin3, valley3) + intravalley_SO(B) << " ";
            myfile << spin_valley(B, spin1, valley1, spin3, valley3) << " ";
            myfile << intravalley_SO(B) << endl;

        }
    }
    myfile.close();
    return 0;
}