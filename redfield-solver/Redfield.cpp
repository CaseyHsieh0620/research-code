#include <iostream>
#include <Eigen/Dense>
using namespace Eigen;
#include <complex>
using namespace std::complex_literals;
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <boost/numeric/odeint.hpp>
#include <fstream>
#include <string>
#include "parameters.cpp"
#include<cmath>

MatrixXcd Energybasis(MatrixXcd matrix, MatrixXcd energy) {
    Eigen::MatrixXcd m = energy.adjoint() * matrix * energy;
    return m;
}

MatrixXcd regularbasis(MatrixXcd matrix, MatrixXcd energy) {
    Eigen::MatrixXcd m  = energy * matrix * energy.adjoint();
    return m;
}

MatrixXcd timepropogation(MatrixXcd ERDM, MatrixXcd ESm, MatrixXcd Eham, BathParam bath) {
    // Here we are solving for our frequencies so that we can calculate the spectral density
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(Eham);
    Eigen::VectorXd eigen = solver.eigenvalues();
    Eigen::MatrixXcd w = Eigen::MatrixXcd::Zero(bath.states, bath.states);
    for (int i = 0; i < bath.states; i++) {
        for (int j = 0; j < bath.states; j++) {
            double wij = eigen(j) - eigen(i);
            w(i,j) = wij;
        }
    }

    // Here we calculate Spectral Density
    Eigen::MatrixXcd jw = Eigen::MatrixXcd::Zero(bath.states, bath.states);
    for(int i = 0; i < bath.states; i++) {
        for (int j = 0; j < bath.states; j++) {
            jw(i,j) = bath.eta * w(i,j) * std::exp(-std::abs(w(i,j))/(bath.w));
        }
    }

    //calculate the rate constant divided by 2 
    Eigen::MatrixXcd rate = Eigen::MatrixXcd::Zero(bath.states, bath.states);
    for (int i = 0; i < bath.states; i++) {
        for (int j = 0; j < bath.states; j++) {
            if (jw(i,j).real() == 0) {
                rate(i,j) = 0;
            } else {
                rate(i, j) = M_PI * jw(i,j) * (1.0 + (1.0 / tanh(w(i,j) / (2 * bath.temp))));
            }
        }
    }

    // Redfield Tensor calculation time
    Eigen::MatrixXcd redfield = Eigen::MatrixXcd::Zero(bath.states, bath.states);
    for (int i = 0; i < bath.states; i++) {
        for(int j = 0; j < bath.states; j++) {
            redfield(i,j) = ESm(i,j) * rate(i,j);
        }
    }

    //Calculate disappation in arbitrary units we dont' take h bar 
    Eigen::MatrixXcd x = Eigen::MatrixXcd::Zero(bath.states, bath.states);
    Eigen::MatrixXcd disspation = Eigen::MatrixXcd::Zero(bath.states, bath.states);
    x = redfield * ERDM - ERDM * redfield.adjoint();
    disspation = ESm * x - x * ESm;

    //These are here from doing euler method or whatever for it. 
    //Calculate the first term in the redfield equation
    //Eigen::MatrixXcd m = Eigen::MatrixXcd::Zero(bath.states, bath.states);
    //m = -1.0i * (Eham * ERDM - ERDM * Eham);

    // Putting it all together and adding dp/dt = 
    //Eigen::MatrixXcd dp = Eigen::MatrixXcd::Zero(bath.states, bath.states);
    //dp = m - disspation;

    
    return -disspation;
}

int main() {
    BathParam bath = setparam();
    Eigen::MatrixXcd RDM = GetRDM(bath);
    Eigen::MatrixXcd Sm = GetSm(bath);
    Eigen::MatrixXcd ham = GetH(bath);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(ham);
    Eigen::MatrixXcd energy = solver.eigenvectors();
    // write it all in energy basis
    Eigen::MatrixXcd ERDM = Energybasis(RDM, energy);
    Eigen::MatrixXcd ESm = Energybasis(Sm,energy);
    Eigen::MatrixXcd Eham = Energybasis(ham,energy);
    //Write the header of the output file
    std::ofstream outFile("Output.txt");
    outFile << std::scientific << std::setprecision(6);
    if (outFile.is_open()) {
        outFile << "Time" << "\t";
        for (int i = 0; i < bath.states; i++) {
            for (int j = 0; j <= i; j++) {
                if (i == j) {
                    outFile << "rho" << i + 1 << "_" << j + 1 <<"_real" << "\t";
                } else {
                    // find both real and im 
                    outFile << "rho" << i + 1 << "_" << j + 1 <<"_real" << "\t";
                    outFile << "rho" << i + 1 << "_" << j + 1 <<"imag" << "\t";
                }
            }
        }
        outFile << "norm" << "\n";
    }
    Eigen::MatrixXcd dp = Eigen::MatrixXcd::Zero(bath.states,bath.states);
    Eigen::MatrixXcd ORDM = regularbasis(ERDM,energy);
    Eigen::VectorXd eigen = solver.eigenvalues();
    for (double k = 0; k < (bath.nprop * bath.dt); k = k + bath.dt) {
        outFile << k << "\t";
        for(int i = 0; i < bath.states; i++) {
            for (int j = 0; j <= i; j++) {
                if (i == j) {
                    outFile << ORDM(i,j).real() << "\t";
                } else {
                    outFile <<  ORDM(i,j).real() << "\t";
                    outFile <<  ORDM(i,j).imag() << "\t";
                }
            }
        }
        outFile<<ORDM.trace().real() << "\t" ;
        for (int i = 0; i <bath.states; i++) {
            for (int j = 0; j <bath.states; j++) {
                double dE = eigen(i)-eigen(j);
                ERDM(i,j) *= std::exp(-1.0i * dE * bath.dt);
            }
        }

        dp = timepropogation(ERDM, ESm, Eham, bath);
        ERDM += dp * bath.dt;
        ORDM = regularbasis(ERDM,energy);
        outFile << "\n";
    }
    
    // here we have p(t+dt) = p(t) + dp/dt * dt 



    return 0;
}
