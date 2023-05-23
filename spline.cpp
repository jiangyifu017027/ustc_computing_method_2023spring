#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>
#include<sstream>

using namespace std;

int main() {
    vector<vector<double>> coordi;
    ifstream pointdata;
    pointdata.open("C:\\ustc_computing_method\\point.txt");
    string line;
    while (getline(pointdata, line)) {
        stringstream ss; 
        vector<double> tempdata;
        ss << line; 
        if (!ss.eof()) {
            double temp;
            while (ss >> temp) {
                tempdata.push_back(temp);
            }
        }
        coordi.push_back(tempdata);
    }
    pointdata.close();
    //测试读入是否有误
    /*
    for (int i = 0; i < coordi.size(); ++i) {
        for (int j = 0; j < coordi[0].size(); ++j) {
            std::cout << coordi[i][j] << "  ";
        }
        std::cout << endl;
    }
    */
    //考虑将压铁移动后的结果
    //coordi[10][1] = 10.0;
    vector<double> h_i;
    for (unsigned int i = 0; i < coordi.size() - 1; ++i) {
        h_i.push_back(coordi[i + 1][0] - coordi[i][0]);
    }
    //测试h_i是否正确
    
    std::cout << h_i.size() << endl;
    for (unsigned int i = 0; i < h_i.size(); ++i) {
        std::cout << h_i[i] << " ";
    }
    
    vector<double> lamda;
    for (unsigned int i = 1; i < h_i.size(); ++i) {
        lamda.push_back(h_i[i] / (h_i[i] + h_i[i - 1]));
    }
    //测试lambda是否正确
    std::cout << lamda.size() << endl;
    for (unsigned int i = 0; i < lamda.size(); ++i) {
        std::cout << lamda[i] << "  ";
    }
    std::cout << endl;
    vector<double> mu;
    for (unsigned int i = 0; i < lamda.size(); ++i) {
        mu.push_back(1.0 - lamda[i]);
    }
    vector<double> d_i;
    for (unsigned int i = 1; i < coordi.size() - 1; ++i) {
        d_i.push_back(6.0 / (h_i[i] + h_i[i - 1]) * ((coordi[i + 1][1] - coordi[i][1]) / h_i[i] - (coordi[i][1] - coordi[i - 1][1]) / h_i[i - 1]));
    }
    std::cout << d_i.size() << endl;
    for (unsigned int i = 0; i < d_i.size(); ++i) {
        std::cout << d_i[i] << "  ";
    }
    std::cout << endl;
    int n = d_i.size();
    std::cout << "n:" << n << endl;
    //考虑采用追赶法求解三对角方程
    vector<double> a(n, 2.0);
    vector<double> b(n);
    for (unsigned int i = 0; i < b.size() - 1; ++i) {
        b[i] = lamda[i];
    }
    b[n - 1] = 0.0;
    vector<double> c(n, 0.0);
    for (unsigned int i = 1; i < c.size(); ++i) {
        c[i] = mu[i];
    }
    vector<double> f(d_i.begin(), d_i.end());
    vector<double> u(n, 0.0);
    vector<double> v(n, 0.0);
    vector<double> x(n, 0.0);
    vector<double> y(n, 0.0);
    for (int k = 0; k < n; ++k) {
        if (k == 0) u[k] = a[k];
        else u[k] = a[k] - c[k] * v[k - 1];
        v[k] = b[k] / u[k];
        if (k == 0) y[k] = f[k] / u[k];
        else y[k] = (f[k] - c[k] * y[k - 1]) / u[k];
    }
    for (int i = n - 1; i >= 0; --i) {
        if (i == n - 1) x[i] = y[i];
        else x[i] = y[i] - v[i] * x[i + 1];
    }
    std::cout << x.size() << endl;
    std::cout << "M:" << endl;
    for (unsigned int i = 0; i < x.size(); ++i) {
        std::cout << x[i] << "  ";
    }
    std::cout << endl;
    vector<double> M(n + 2, 0.0);
    for (int i = 1; i <= n; ++i) {
        M[i] = x[i - 1];
    }
    //输出结果
    for (int i = 1; i <= 20; ++i) {
        std::cout << "S_" << i << "(x):";
        double C = coordi[i - 1][1] / h_i[i - 1] - h_i[i - 1] * M[i - 1] / 6.0;
        double D = coordi[i][1] / h_i[i - 1] - h_i[i - 1] * M[i] / 6.0;
        std::cout << M[i - 1] / (6.0 * h_i[i - 1]) << "(" << coordi[i][0] << "-x)^3+" << M[i] / (6.0 * h_i[i - 1]) << "(x-(" << coordi[i][0] << "))^3+" << C << "(" << coordi[i][0] << "-x)+" << D << "(x-(" << coordi[i - 1][0] << "))" << endl;
    }
    return 0;
}