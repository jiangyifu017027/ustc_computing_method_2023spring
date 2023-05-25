#include<iostream>
#include<vector>
#include<complex>
#include<cmath>
#include<random>
#include<algorithm>
#include<fstream>
#include<sstream>

using namespace std;

const double PI = 3.1415926535;

vector<complex<double>> FFT(const vector<complex<double>>& f) {
    int n = f.size();
    if (n == 1) return f;

    complex<double> wn = exp(-2 * PI / n * complex<double>(0.0, 1.0));
    complex<double> w = 1.0;
    vector<complex<double>> f0;
    vector<complex<double>> f1;
    for (unsigned int i = 0; i < f.size(); i += 2) {
        f0.push_back(f[i]);
        f1.push_back(f[i + 1]);
    }
    vector<complex<double>> g0 = FFT(f0);
    vector<complex<double>> g1 = FFT(f1);
    vector<complex<double>> g(n, complex<double>(0.0, 0.0));
    for (int k = 0; k < n / 2; ++k) {
        g[k] = (g0[k] + w * g1[k]) / 2.0;
        g[k + n / 2] = (g0[k] - w * g1[k]) / 2.0;
        w = w * wn;
    }
    return g;
}

vector<complex<double>> IFFT(const vector<complex<double>>& f) {
    int n = f.size();
    if (n == 1) return f;

    complex<double> wn = exp(2 * PI / n * complex<double>(0.0, 1.0));
    complex<double> w = 1.0;
    vector<complex<double>> f0;
    vector<complex<double>> f1;
    for (unsigned int i = 0; i < f.size(); i += 2) {
        f0.push_back(f[i]);
        f1.push_back(f[i + 1]);
    }
    vector<complex<double>> g0 = IFFT(f0);
    vector<complex<double>> g1 = IFFT(f1);
    vector<complex<double>> g(n, complex<double>(0.0, 0.0));
    for (int k = 0; k < n / 2; ++k) {
        g[k] = g0[k] + w * g1[k];
        g[k + n / 2] = g0[k] - w * g1[k];
        w = w * wn;
    }
    return g;
}

complex<double> f1(double t) {
    return complex<double>(0.7 * sin(4 * PI * t) + sin(10 * PI * t), 0.0);
}

complex<double> f2(double t) {
    uniform_real_distribution<double> uid{0, 1};
    random_device rd;
    default_random_engine dre{rd()}; //uid(dre)
    return complex<double>(0.7 * sin(4 * PI * t) + sin(10 * PI * t) + 0.3 * uid(dre), 0.0);
}

int main() {
    //考虑f1,n=2^4的情况
    double n = pow(2.0, 4);
    vector<complex<double>> f;
    for (double i = 0.0; fabs(n - i) >=  1e-5; i += 1.0) {
        f.push_back(f1(i / n));
    }
    vector<complex<double>> g = FFT(f);
    ofstream dataG1;
    dataG1.open("C:\\ustc_computing_method\\dataG1.txt");
    for (int i = 0; i < int(n); ++i) {
        double temp = fabs(g[i]);
        dataG1 << i << " " << temp << endl;
    }
    dataG1.close();
    /*
    std::cout << "取f1,n=2^4:" << endl;
    std::cout << "g:" << endl;
    std::cout << "实部"  << "  " << "虚部" << endl;
    for (int i = 0; i < g.size(); ++i) {
        std::cout << real(g[i]) << "  " << imag(g[i]) << endl;
    }
    */
    //考虑f1,n=2^7的情况
    n = pow(2.0, 7);
    f.empty();
    for (double i = 0.0; fabs(n - i) >=  1e-5; i += 1.0) {
        f.push_back(f1(i / n));
    }
    g.empty();
    g = FFT(f);
    ofstream dataG2;
    dataG2.open("C:\\ustc_computing_method\\dataG2.txt");
    for (int i = 0; i < int(n); ++i) {
        double temp = fabs(g[i]);
        dataG2 << i << " " << temp << endl;
    }
    dataG2.close();
    /*
    std::cout << "取f1,n=2^7:" << endl;
    std::cout << "g:" << endl;
    std::cout << "实部"  << "  " << "虚部" << endl;
    for (int i = 0; i < g.size(); ++i) {
        std::cout << real(g[i]) << "  " << imag(g[i]) << endl;
    }
    */
    //考虑f2,n=2^7的情况
    f.empty();
    for (double i = 0.0; fabs(n - i) >=  1e-5; i += 1.0) {
        f.push_back(f2(i / n));
    }
    g.empty();
    g = FFT(f);
    ofstream dataG3;
    dataG3.open("C:\\ustc_computing_method\\dataG3.txt");
    for (int i = 0; i < int(n); ++i) {
        double temp = fabs(g[i]);
        dataG3 << i << " " << temp << endl;
    }
    dataG3.close();
    /*
    std::cout << "取f2,n=2^7:" << endl;
    std::cout << "g:" << endl;
    std::cout << "实部"  << "  " << "虚部" << endl;
    for (int i = 0; i < g.size(); ++i) {
        std::cout << real(g[i]) << "  " << imag(g[i]) << endl;
    }
    */
    return 0;
}