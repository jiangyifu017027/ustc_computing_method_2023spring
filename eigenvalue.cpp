#include<iostream>
#include<vector>
#include<cmath>

using namespace std;

void Doolittle(vector<vector<double>> L, vector<vector<double>> U, vector<double>& Y, vector<double>& X) {
    int n = Y.size();
    vector<double> temp(n, 0);
    for (int i = 0; i < n; ++i) {
        double sum = 0;
        for (int j = 0; j < i; ++j) {
            sum += L[i][j] * temp[j];
        }
        temp[i] = Y[i] - sum;
    }
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0;
        for (int j = i + 1; j < n; ++j) {
            sum += U[i][j] * X[j];
        }
        X[i] = (temp[i] - sum) / U[i][i];
    }
}

void eigenvalue(vector<vector<double>> coeffiMatrix, vector<double> vecquan) {
    int n = coeffiMatrix.size();
    vector<vector<double>> L(n, vector<double>(n, 0));
    vector<vector<double>> U(n, vector<double>(n, 0));
    // Doolittle分解 求解L,U矩阵 经检验L,U的求解代码没有问题
    for (int i = 0; i < n; ++i) {
        L[i][i] = 1;
    }
    for (int k = 0; k < n; ++k) {
        for (int j = k; j < n; ++j) {
            double temp = 0;
            for (int r = 0; r < k; ++r) {
                temp += L[k][r] * U[r][j];
            }
            U[k][j] = coeffiMatrix[k][j] - temp;
        }
        for (int i = k + 1; i < n; ++i) {
            double temp = 0;
            for (int r = 0; r < k; ++r) {
                temp += L[i][r] * U[r][k];
            }
            L[i][k] = (coeffiMatrix[i][k] - temp) / U[k][k];
        }
    }
    /*
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << L[i][j] << " ";
        }
        cout << endl;
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << U[i][j] << " ";
        }
        cout << endl;
    }
    */
    vector<vector<double>> resultY;
    vector<vector<double>> resultX;
    vector<double> resultLamda;
    vector<double> Y(vecquan.begin(), vecquan.end());
    vector<double> X(vecquan.begin(), vecquan.end());
    resultY.push_back(Y);
    Doolittle(L, U, Y, X);
    resultX.push_back(X);
    double LamdaCur = Y[0] / X[0];
    resultLamda.push_back(LamdaCur);
    double LamdaPre = LamdaCur;
    do {
        LamdaPre = LamdaCur;
        double Max = X[0];
        for (int i = 0; i < X.size(); ++i) {
            if (X[i] >= Max) Max = X[i];
        }
        for (int i = 0; i < X.size(); ++i) {
            Y[i] = X[i] / Max;
        }
        resultY.push_back(Y);
        Doolittle(L, U, Y, X);
        resultX.push_back(X);
        LamdaCur = Y[0] / X[0];
        resultLamda.push_back(LamdaCur);
    } while (fabs(LamdaPre - LamdaCur) >= 1e-5);
    cout << "迭代过程中的Y:" << endl;
    for (int i = 0; i < resultY.size(); ++i) {
        for (int j = 0; j < resultY[0].size(); ++j) {
            cout << resultY[i][j] << "  ";
        }
        cout << endl;
    }
    cout << "迭代过程中的X:" << endl;
    for (int i = 0; i < resultX.size(); ++i) {
        for (int j = 0; j < resultX.size(); ++j) {
            cout << resultX[i][j] << "  ";
        }
        cout << endl;
    }
    cout << "迭代过程中的Lamda:" << endl;
    for (int i = 0; i < resultLamda.size(); ++i) {
        cout << resultLamda[i] << "  ";
    }
    cout << endl;
    int len = resultY.size() - 1;
    cout << "特征向量为:" << endl;
    for (int i = 0; i < resultY[0].size(); ++i) {
        cout << resultY[len][i] << "  ";
    }
    cout << endl;
}

int main() {
    vector<vector<double>> coeffiMatrix1(5, vector<double>(5, 0));
    vector<vector<double>> coeffiMatrix2{{4, -1, 1, 3}, {16, -2, -2, 5}, {16, -3, -1, 7}, {6, -4, 2, 9}};
    for (double i = 0; i < 5; ++i) {
        for (double j = 0; j < 5; ++j) {
            coeffiMatrix1[i][j] = 1 / (9 - i - j);
        }
    }
    vector<double> vecquan1(5, 1);
    vector<double> vecquan2(4, 1);
    //eigenvalue(coeffiMatrix1, vecquan1);
    //eigenvalue(coeffiMatrix2, vecquan2);
    eigenvalue(coeffiMatrix1, vecquan1);
    eigenvalue(coeffiMatrix2, vecquan2);
    return 0;
}