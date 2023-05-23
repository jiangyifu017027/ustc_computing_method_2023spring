#include<iostream>
#include<vector>
#include<cmath>

using namespace std;

//精确解函数
double f(double epsilon, double a, double x) {
    return ((1 - a) / (1 - exp(-1 / epsilon))) * (1 - exp(-x / epsilon)) + a * x;
}

//计算精确解并将其输出
void exactSolution(double epsilon, double a, double n) {
    int len = int(n);
    vector<double> result(len - 1, 0);
    double interval = 1 / n;
    cout << "精确解为:" << endl; 
    for (double i = 1; i < n; ++i) {
        double x = i * interval;
        result[i - 1] = f(epsilon, a, x);
        cout << result[i - 1] << "  ";
        if (int(i) % 10 == 0) cout << endl;
    }
    cout << endl;
    return;
}

//列主元Gauss消元法 计算并输出结果
void Guass(double n, double epsilon, double a) {
    int len = static_cast<int>(n);
    double coeMatrix[len - 1][len] = {0};
    double h = 1 / n;
    cout << "Gauss列主元法解得:" << endl;
    for (int i = 0; i < len - 1; ++i) {
        if (i == 0) {
            coeMatrix[i][0] = -(2 * epsilon + h);
            coeMatrix[i][1] = epsilon + h;
            coeMatrix[i][len - 1] = a * h * h;
        }
        else if (i == len - 2) {
            coeMatrix[i][len - 3] = epsilon;
            coeMatrix[i][len - 2] = - (2 * epsilon + h);
            coeMatrix[i][len - 1] = a * h * h - epsilon - h;
        }
        else {
            coeMatrix[i][i - 1] = epsilon;
            coeMatrix[i][i] = - (2 * epsilon + h);
            coeMatrix[i][i + 1] = epsilon + h;
            coeMatrix[i][len - 1] = a * h * h;
        }
    }
    for (int i = 0; i < len - 1; ++i) {
        int k = i;
        for (int j = i + 1; j < len - 1; ++j) {
            if (fabs(coeMatrix[k][i]) < fabs(coeMatrix[j][i])) {
                k = j;
            }
        }
        for (int j = i; j < len; ++j) {
            double temp = coeMatrix[i][j];
            coeMatrix[i][j] = coeMatrix[k][j];
            coeMatrix[k][j] = temp;
        }
        for (int j = i + 1; j < len - 1; ++j) {
            coeMatrix[j][i] = coeMatrix[j][i] / coeMatrix[i][i];
            for (int k = i + 1; k < len; ++k) {
                coeMatrix[j][k] = coeMatrix[j][k] - coeMatrix[j][i] * coeMatrix[i][k];
            }
            coeMatrix[j][i] = 0;
        }
    }
    for (int i = len - 2; i >= 0; --i) {
        for (int j = i + 1; j < len - 1; ++j) {
            coeMatrix[i][len - 1] = coeMatrix[i][len - 1] - coeMatrix[i][j] * coeMatrix[j][len - 1];
        }
        coeMatrix[i][len - 1] = coeMatrix[i][len - 1] / coeMatrix[i][i];
    }
    for (int i = 0; i < len - 1; ++i) {
        cout << coeMatrix[i][len - 1] << " ";
        if ((i + 1) % 10 == 0) cout << endl;
    }
    cout << endl;
}

double norm(double x1[], double x2[], int len) {
    double diffvalue = 0;
    for (int i = 0; i < len - 1; ++i) {
        if (fabs(x1[i] - x2[i]) > diffvalue) {
            diffvalue = fabs(x1[i] - x2[i]);
        }
    }
    return diffvalue;
}

//Gauss-Seidel迭代法 求解并输出答案
void Gauss_Seidel(double n, double epsilon, double a) {
    int len = static_cast<int>(n);
    double aMatr[len - 1][len - 1] = {0};
    double b[len - 1] = {0};
    double x1[len - 1] = {1};
    double x2[len - 1] = {0};
    double h = 1 / n;
    cout << "Gauss-Seidel迭代算法计算得到:" << endl;
    for (int i = 0; i < len - 1; ++i) {
        if (i == 0) {
            aMatr[i][0] = -(2 * epsilon + h);
            aMatr[i][1] = epsilon + h;
        }
        else if (i == len - 2) {
            aMatr[i][len - 3] = epsilon;
            aMatr[i][len - 2] = -(2 * epsilon + h);
        }
        else {
            aMatr[i][i - 1] = epsilon;
            aMatr[i][i] = - (2 * epsilon + h);
            aMatr[i][i + 1] = epsilon + h;
        }
    }
    for (int i = 0; i < len - 1; ++i) {
        if (i == len - 2) {
            b[i] = a * h * h - epsilon - h;
        }
        else {
            b[i] = a * h * h;
        }
    }
    double error = 1e-6;
    while (norm(x1, x2, len) > error) {
        for (int i = 0; i < len - 1; ++i) {
            x1[i] = x2[i];
        }
        for (int i = 0; i < len - 1; ++i) {
            double s = 0;
            for (int j = 0; j < len - 1; ++j) {
                s = s + aMatr[i][j] * x2[j];
            }
            x2[i] = (b[i] - s + aMatr[i][i] * x2[i]) / aMatr[i][i];
        }
    }
    for (int i = 0; i < len - 1; ++i) {
        cout << x1[i] << " ";
        if ((i + 1) % 10 == 0) {
            cout << endl;
        }
    }
    cout << endl;
} 

int main() {
    double n = 100;
    double epsilon = 1;
    double a = 0.5;
    exactSolution(epsilon, a, n);
    //Guass(n, epsilon, a);
    Gauss_Seidel(n, epsilon, a);
    epsilon = 0.1;
    exactSolution(epsilon, a, n);
    Guass(n, epsilon, a);
    Gauss_Seidel(n, epsilon, a);
    epsilon = 0.01;
    exactSolution(epsilon, a, n);
    Guass(n, epsilon, a);
    Gauss_Seidel(n, epsilon, a);
    epsilon = 0.001;
    exactSolution(epsilon, a, n);
    Guass(n, epsilon, a);
    Gauss_Seidel(n, epsilon, a);
    return 0;
}
