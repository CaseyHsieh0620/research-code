#include <iostream>
#include <complex>
#include <fstream>
#include <string>
#include <Eigen/Dense>
 using namespace Eigen;
#include <complex>
using namespace std::complex_literals;


struct BathParam {
    double eta;
    double w;
    double temp;
    double states;
    double nprop;
    double dt;
};

BathParam setparam() {
    std::ifstream file("parameters.txt");
    if (!file.is_open()) {
        throw std::runtime_error("Where is parameters.txt??");
    }
    std::string currentline;
    BathParam bath;
    int count = 0;
    while (std::getline(file, currentline)) {
        size_t eqpos = currentline.find('=');
        double numbervalue = 0.0;
        if (eqpos != std::string::npos) {
            std::string stringvalue = currentline.substr(eqpos + 1);
            try {
                numbervalue = std::stod(stringvalue);
            } catch (std::exception& e) {
                std::cout << "parameters.txt not formatted correctly" << std::endl;
            }
            if (currentline.find("Coupling Strength") != std::string::npos) {
                bath.eta = numbervalue; 
                count = count + 1;
            }
            if (currentline.find("Frequency") != std::string::npos) {
                bath.w = numbervalue;
                count = count + 1;
            }
            if (currentline.find("Temperature") != std::string::npos) {
                bath.temp = numbervalue;
                count = count + 1;
            }
            if (currentline.find("Number of States") != std::string::npos) {
                bath.states = numbervalue;
                count = count + 1;
            }
            if (currentline.find("Number of Steps") != std::string::npos)  {
                bath.nprop = numbervalue;
                count = count + 1;
            }
            if (currentline.find("dt") != std::string::npos)  {
                bath.dt = numbervalue;
                count = count + 1;
            }
            if (count == 6) {
                break;
            }
        }
    }
    return bath; 
};

MatrixXcd GetSm(BathParam baths) {
    int N = baths.states;
    Eigen::MatrixXcd Sm = Eigen::MatrixXcd::Zero(N, N);    
    std::ifstream file("Sm.txt");
    if (!file.is_open()) {
        throw std::runtime_error("You're missing the System Operator Matrix");
    }
    std::string line;
    int linecount  = 0;
    while (std::getline(file, line)) {
        if(line.find('#') == std::string::npos && line.empty() == false) {
            std::stringstream ss(line);
            std::complex<double> value;
            for (int i = 0; i < N; i++) {
                if (ss >> value) {
                    Sm(linecount, i) = value; 
                } else {
                    throw std::runtime_error("Line  " + std::to_string(linecount) + " has the incorrect # of columns in Sm");
                }
            }
            linecount++;
            if (linecount == N) {
                break;
            }
        }
    }
    if (linecount != N) {
        throw std::runtime_error("Your System Operator Matrix is formated incorrectly, You have " + std::to_string(linecount) + " amount of lines but it should be " + std::to_string(N));
    }
    return Sm;
}

MatrixXcd GetH(BathParam baths) {
    int N = baths.states;
    Eigen::MatrixXcd Sm = Eigen::MatrixXcd::Zero(N, N);    
    std::ifstream file("Hamiltonian.txt");
    if (!file.is_open()) {
        throw std::runtime_error("You're missing the Hamiltonian Matrix");
    }
    std::string line;
    int linecount  = 0;
    while (std::getline(file, line)) {
        if(line.find('#') == std::string::npos && line.empty() == false) {
            std::stringstream ss(line);
            std::complex<double> value;
            for (int i = 0; i < N; i++) {
                if (ss >> value) {
                    Sm(linecount, i) = value; 
                } else {
                    throw std::runtime_error("Line  " + std::to_string(linecount) + " has the incorrect # of columns in Hamiltonian");
                }
            }
            linecount++;
            if (linecount == N) {
                break;
            }
        }
    }
    if (linecount != N) {
        throw std::runtime_error("Your Hamiltonian Matrix is formated incorrectly, You have " + std::to_string(linecount) + " amount of lines but it should be " + std::to_string(N));
    }
    return Sm;
}

MatrixXcd GetRDM(BathParam baths) {
    int N = baths.states;
    Eigen::MatrixXcd Sm = Eigen::MatrixXcd::Zero(N, N);    
    std::ifstream file("IntRDM.txt");
    if (!file.is_open()) {
        throw std::runtime_error("You're missing the IntRDM file");
    }
    std::string line;
    int linecount  = 0;
    while (std::getline(file, line)) {
        if(line.find('#') == std::string::npos && line.empty() == false) {
            std::stringstream ss(line);
            std::complex<double> value;
            for (int i = 0; i < N; i++) {
                if (ss >> value) {
                    Sm(linecount, i) = value; 
                    if (linecount == i) {
                        if (Sm(linecount, i) != value.real()) {
                            std::cout << "Made Diaganol Real" << std::endl;
                        }
                        Sm(linecount, i) = value.real();
                    }
                } else {
                    throw std::runtime_error("Line  " + std::to_string(linecount) + " has the incorrect # of columns in IntRDM should be " + std::to_string(N));
                }
            }
            linecount++;
            if (linecount == N) {
                break;
            }
        }
    }
    if (linecount != N) {
        throw std::runtime_error("Your IntRDM Matrix is formated incorrectly, You have " + std::to_string(linecount) + " amount of lines but it should be " + std::to_string(N));
    }
    return Sm;
}





