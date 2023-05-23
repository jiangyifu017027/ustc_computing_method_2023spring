#include<iostream>
#include<cmath>
#include<fstream>
#include<vector>
#include<sstream>
#include<string>

using namespace std;

vector<double> x_P;
vector<double> y_P;
vector<double> x_Q;
vector<double> y_Q;

double ThetaMax(double x_Q, double y_Q) {
    return 3.1415926 - atan(y_Q / x_Q);
}

double distance(double x1, double y1, double x2, double y2) {
    return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
}

double f(double x, double x_P, double y_P, double x_Q, double y_Q) {
    double disOQ = distance(0, 0, x_Q, y_Q);
    double disOP = distance(0, 0, x_P, y_P);
    double x_T = cos(x);
    double y_T = sin(x);
    double disTQ = distance(x_T, y_T, x_Q, y_Q);
    double disTP = distance(x_T, y_T, x_P, y_P);
    return(((pow(disTQ, 2) + 1.0 - pow(disOQ, 2)) * disTP - (pow(disTP, 2) + 1.0 + pow(disOP, 2)) * disTQ));
} // 靠近max则为负数

int main() {
    x_P.clear();
    x_Q.clear();
    y_P.clear();
    y_Q.clear();
    ofstream mirrorcup_answer;
    mirrorcup_answer.open("C:\\ustc_computing_method\\mirrorcup_answer.txt1");
    ifstream read_mirrorcup;
    read_mirrorcup.open("C:\\custc_computing_method\\mirrorcup.txt");
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
    const double epsilon = 1e-8;
    int len = x_P.size();
    for (int i = 0; i < len; ++i) {
        double thetaMax = 3.1415926;
        double thetaMin = ThetaMax(x_Q[i], y_Q[i]);
        while ((thetaMax - thetaMin) >= epsilon) {
            double midpoint = (thetaMax + thetaMin) / 2;
            double temp = f(midpoint, x_P[i], y_P[i], x_Q[i], y_Q[i]);
            if (temp >= 0) {
                thetaMin = midpoint;
            }
            else {
                thetaMax = midpoint;
            }
        }
        double theta = (thetaMax + thetaMin) / 2;
        double x_T = cos(theta);
        double y_T = sin(theta);
        mirrorcup_answer<<x_T<<" "<< y_T<<endl;
    }
    mirrorcup_answer.close();
    return 0;
}