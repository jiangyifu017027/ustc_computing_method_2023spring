#include<iostream>
#include<cmath>
#include<fstream>
#include<vector>
#include<sstream>
#include<string>
#include<iomanip>

using namespace std;

vector<double> x_P;
vector<double> y_P;
vector<double> x_Q;
vector<double> y_Q;

double distance(double x1, double y1, double x2, double y2) {
    return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
}

double dotproduct(double x1, double y1, double x2, double y2) {
    return x1 * x2 + y1 * y2;
}

double f1(double x0, double y0, double disOP, double disOQ, double dotPQ) {
    return (disOP * disOP * x0 * x0 + 2 * dotPQ * x0 * y0 + disOQ * disOQ * y0 * y0 - 1);
}

double f2(double x0, double y0, double disOQ, double dotPQ) {
    return (2 * dotPQ * x0 * y0 + 2 * disOQ * disOQ * y0 * y0 - y0 - 1 + x0);
}

double J11(double x0, double y0, double disOP, double dotPQ) {
    return (2 * disOP * disOP * x0 + 2 * dotPQ * y0);
}

double J12(double x0, double y0, double disOQ, double dotPQ) {
    return (2 * disOQ * disOQ * y0 + 2 * dotPQ * x0);
}

double J21(double x0, double y0, double dotPQ) {
    return (2 * dotPQ * y0 + 1);
}

double J22(double x0, double y0, double disOQ, double dotPQ) {
    return (2 * dotPQ * x0 + 4 * disOQ * disOQ * y0 - 1);
}



int main() {
    x_P.clear();
    x_Q.clear();
    y_P.clear();
    y_Q.clear();
    ofstream mirrorcup_answer;
    mirrorcup_answer.open("C:\\ustc_computing_method\\mirrorcup_answer.txt2");
    ifstream read_mirrorcup;
    read_mirrorcup.open("C:\\ustc_computing_method\\mirrorcup.txt");
    string line;
    while (getline(read_mirrorcup, line)) {
        stringstream ss;
        ss << line;
        if (!ss.eof()) {
            string temp;
            int order = 1;
            while (ss >> temp) {
                double tmp = stod(temp);
                if (order == 1) {
                    x_P.push_back(tmp);
                    order++;
                }
                else if (order == 2) {
                    y_P.push_back(tmp);
                    order++;
                }
                else if (order == 3) {
                    x_Q.push_back(tmp);
                    order++;
                }
                else if (order == 4) {
                    y_Q.push_back(tmp);
                }
            }
        }
    } 
    read_mirrorcup.close();
    const double epsilon = 1e-6;
    int len = x_P.size();
    for (int i = 0; i < len; ++i) {
        double x0 = -1 / x_P[i], y0 = 0;
        double disOP = -x_P[i];
        double disOQ = distance(x_Q[i], y_Q[i], 0, 0);
        double dotPQ = dotproduct(x_P[i], y_P[i], x_Q[i], y_Q[i]);
        double delx = 1, dely = 1;
        while (max(delx, dely) >= epsilon) {
            double f_1 = f1(x0, y0, disOP, disOQ, dotPQ);
            double f_2 = f2(x0, y0, disOQ, dotPQ);
            double J_11 = J11(x0, y0, disOP, dotPQ);
            double J_12 = J12(x0, y0, disOQ, dotPQ);
            double J_21 = J21(x0, y0, dotPQ);
            double J_22 = J22(x0, y0, disOQ, dotPQ);
            delx = (-f_1 * J_22 + f_2 * J_12) / (J_11 * J_22 - J_21 * J_12);
            dely = (-f_1 * J_21 + f_2 * J_11) / (J_21 * J_12 - J_22 * J_11);
            x0 += delx;
            y0 += dely;
        }
        double x_T = x0 * x_P[i] + y0 * x_Q[i];
        double y_T = x0 * y_P[i] + y0 * y_Q[i];
        double disPT = distance(x_P[i], y_P[i], x_T, y_T);
        double disQT = distance(x_Q[i], y_Q[i], x_T, y_T);
        double x_R = (x_T - x_P[i]) * disQT / disPT + x_T;
        double y_R = y_T * ((disPT + disQT) / disPT);
        mirrorcup_answer<< setprecision(8) << x_T<< " "<< setprecision(8) << y_T <<" "<< setprecision(8) << x_R<< " "<< setprecision(8) << y_R<< endl;
    }
    mirrorcup_answer.close();
    return 0;
}