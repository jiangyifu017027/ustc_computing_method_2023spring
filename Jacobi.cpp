//计算方法第四次实验
//Jacobi方法求解对称阵特征值

#include<iostream>
#include<vector>
#include<cmath>
#include<random>
#include<algorithm>
#include<fstream>
#include<sstream>

using namespace std;

//选取非对角线按模最大元素
double matrixMax(vector<vector<double>> a, int & poix, int & poiy) {
    int len = a.size();
    double result = 0.0;
    for (int i = 0; i < len; ++i) {
        for (int j = 0; j < len; ++j) {
            if (i != j && fabs(a[i][j]) > result) {
                result = fabs(a[i][j]);
                poix = i;
                poiy = j;
            }
        }
    }
    return result;
}

//计算非对角线元素平方和
double matrixSum(vector<vector<double>> a) {
    int len = a.size();
    double result = 0.0;
    for (int i = 0; i < len; ++i) {
        for (int j = 0; j < len; ++j) {
            if (i != j) {
                result += a[i][j] * a[i][j];
            }
        }
    }
    return result;
}

int main() {
    uniform_real_distribution<double> uid{0, 1};
    random_device rd;
    default_random_engine dre{rd()}; //uid(dre)
    vector<vector<double>> ran(4, vector<double>(3, 0.0));
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            ran[i][j] = uid(dre);
        }
    }
    //输出随机矩阵A
    std::cout << "random matrix A is:" << endl;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            std::cout << ran[i][j] << "  ";
        }
        std::cout << endl;
    }
    //计算A的转置
    vector<vector<double>> ranT(3, vector<double>(4, 0.0));
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            ranT[i][j] = ran[j][i];
        }
    }
    //计算乘积
    vector<vector<double>> ranA(4, vector<double>(4, 0.0));
    vector<vector<double>> ranAT(3, vector<double>(3, 0.0));
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            for (int k = 0; k < 3; ++k) {
                ranA[i][j] += ran[i][k] * ranT[k][j];
            }
        }
    }
    std::cout << "random matrix A*A^T is:" << endl;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            std::cout << ranA[i][j] << "  ";
        }
        std::cout << endl;
    }
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 4; ++k) {
                ranAT[i][j] += ranT[i][k] * ran[k][j];
            }
        }
    }
    std::cout << "random matrix A^T*A is:" << endl;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            std::cout << ranAT[i][j] << "  ";
        }
        std::cout << endl;
    }
    double epsilon = 1e-6;
    vector<vector<double>> Q(ranA.size(), vector<double>(ranA.size(), 0.0));
    for (int i = 0; i < Q.size(); ++i) {
        Q[i][i] = 1.0;
    }
    while (matrixSum(ranA) >= epsilon) {
        std::cout << "非对角元素平方和为：" << matrixSum(ranA) << endl;
        int poix = 0, poiy = 0;
        double maxA = matrixMax(ranA, poix, poiy);
        double s = (ranA[poiy][poiy] - ranA[poix][poix]) / (2 * ranA[poix][poiy]);
        double t = 0.0;
        if (s == 0.0) {
            t = 1.0;
        }
        else {
            double t1 = -s - sqrt(s * s + 1.0);
            double t2 = -s + sqrt(s * s + 1.0);
            t = (fabs(t1) > fabs(t2)) ? t2 : t1;
        }
        double c =  1.0 / sqrt(t * t + 1.0);
        double d = t / sqrt(t * t + 1.0);
        vector<vector<double>> ranA1(ranA);
        for (int i = 0; i < ranA.size(); ++i) {
            if (i == poix || i == poiy) continue;
            ranA[i][poix] = c * ranA1[poix][i] - d * ranA1[poiy][i];
            ranA[poix][i] = ranA[i][poix];
            ranA[i][poiy] = c * ranA1[poiy][i] + d * ranA1[poix][i];
            ranA[poiy][i] = ranA[i][poiy];
        }
        ranA[poix][poix] = ranA1[poix][poix] - t * ranA1[poix][poiy];
        ranA[poiy][poiy] = ranA1[poiy][poiy] + t * ranA1[poix][poiy];
        ranA[poix][poiy] = 0.0;
        ranA[poiy][poix] = 0.0;
        vector<vector<double>> Q1(Q.size(), vector<double>(Q.size(), 0.0));
        for (int i = 0; i < Q.size(); ++i) Q1[i][i] = 1.0;
        Q1[poix][poix] = c;
        Q1[poix][poiy] = d;
        Q1[poiy][poix] = -d;
        Q1[poiy][poiy] = c;
        decltype(Q1) Q2 = Q;
        for (int i = 0; i < Q.size(); ++i) {
            for (int j = 0; j < Q.size(); ++j) {
                double temp = 0.0;
                for (int k = 0; k < Q.size(); ++k) {
                    temp += Q2[i][k] * Q1[k][j];
                }
                Q[i][j] = temp;
            }
        }
    }
    vector<double> eigen(ranA.size(), 0.0);
    for (int i = 0; i < ranA.size(); ++i) {
        eigen[i] = ranA[i][i];
    }
    //sort(eigen.begin(), eigen.end(), greater<double>());
    vector<vector<double>> QT(ranAT.size(), vector<double>(ranAT.size(), 0.0));
    for (int i = 0; i < QT.size(); ++i) {
        QT[i][i] = 1.0;
    }
    while (matrixSum(ranAT) >= epsilon) {
        int poix = 0, poiy = 0;
        double maxA = matrixMax(ranAT, poix, poiy);
        double s = (ranAT[poiy][poiy] - ranAT[poix][poix]) / (2 * ranAT[poix][poiy]);
        double t = 0.0;
        if (s == 0.0) {
            t = 1.0;
        }
        else {
            double t1 = -s - sqrt(s * s + 1.0);
            double t2 = -s + sqrt(s * s + 1.0);
            t = (fabs(t1) > fabs(t2)) ? t2 : t1;
        }
        double c =  1.0 / sqrt(t * t + 1.0);
        double d = t / sqrt(t * t + 1.0);
        vector<vector<double>> ranA1(ranAT);
        for (int i = 0; i < ranAT.size(); ++i) {
            if (i == poix || i == poiy) continue;
            ranAT[i][poix] = c * ranA1[poix][i] - d * ranA1[poiy][i];
            ranAT[poix][i] = ranAT[i][poix];
            ranAT[i][poiy] = c * ranA1[poiy][i] + d * ranA1[poix][i];
            ranAT[poiy][i] = ranAT[i][poiy];
        }
        ranAT[poix][poix] = ranA1[poix][poix] - t * ranA1[poix][poiy];
        ranAT[poiy][poiy] = ranA1[poiy][poiy] + t * ranA1[poix][poiy];
        ranAT[poix][poiy] = 0.0;
        ranAT[poiy][poix] = 0.0;
        vector<vector<double>> Q1(QT.size(), vector<double>(QT.size(), 0.0));
        for (int i = 0; i < QT.size(); ++i) Q1[i][i] = 1.0;
        Q1[poix][poix] = c;
        Q1[poix][poiy] = d;
        Q1[poiy][poix] = -d;
        Q1[poiy][poiy] = c;
        decltype(Q1) Q2 = QT;
        for (int i = 0; i < QT.size(); ++i) {
            for (int j = 0; j < QT.size(); ++j) {
                double temp = 0.0;
                for (int k = 0; k < QT.size(); ++k) {
                    temp += Q2[i][k] * Q1[k][j];
                }
                QT[i][j] = temp;
            }
        }
    }
    vector<double> eigenT(ranAT.size(), 0.0);
    for (int i = 0; i < ranAT.size(); ++i) {
        eigenT[i] = ranAT[i][i];
    }
    std::cout << "A*A^T的特征值是:" << endl;
    for (auto lamda : eigen) {
        std::cout << lamda << "  "; 
    }
    std::cout << endl;
    std::cout << "A^T*A的特征值是:" << endl;
    for (auto lamda : eigenT) {
        std::cout << lamda << "  "; 
    }
    std::cout << endl;
    //输出特征向量
    std::cout << "A*A^T对应的特征向量为:" << endl;
    for (int i = 0; i < Q.size(); ++i) {
        for (int j = 0; j < Q.size(); ++j) {
            std::cout << Q[i][j] << "  ";
        }
        std::cout << endl;
    }
    std::cout << "A^T*A对应的特征向量为:" << endl;
    for (int i = 0; i < QT.size(); ++i) {
        for (int j = 0; j < QT.size(); ++j) {
            std::cout << QT[i][j] << "  ";
        }
        std::cout << endl;
    }
    //对特征值和特征向量进行排序
    for (int i = 0; i < eigen.size() - 1; ++i) {
        for (int j = eigen.size() - 1; j > i; --j) {
            if (eigen[j] > eigen[j - 1]) {
                double temp = eigen[j];
                eigen[j] = eigen[j - 1];
                eigen[j - 1] = temp;
                for (int k = 0; k < eigen.size(); ++k) {
                    double temp = Q[k][j];
                    Q[k][j] = Q[k][j - 1];
                    Q[k][j - 1] = temp;
                }
            }
        }
    }
    for (int i = 0; i < eigenT.size() - 1; ++i) {
        for (int j = eigenT.size() - 1; j > i; --j) {
            if (eigenT[j] > eigenT[j - 1]) {
                double temp = eigenT[j];
                eigenT[j] = eigenT[j - 1];
                eigenT[j - 1] = temp;
                for (int k = 0; k < eigenT.size(); ++k) {
                    double temp = QT[k][j];
                    QT[k][j] = QT[k][j - 1];
                    QT[k][j - 1] = temp;
                }
            }
        }
    }
    std::cout << "A*A^T的特征值是:" << endl;
    for (auto lamda : eigen) {
        std::cout << lamda << "  "; 
    }
    std::cout << endl;
    std::cout << "A^T*A的特征值是:" << endl;
    for (auto lamda : eigenT) {
        std::cout << lamda << "  "; 
    }
    std::cout << endl;
    std::cout << "A*A^T对应的特征向量为:" << endl;
    for (int i = 0; i < Q.size(); ++i) {
        for (int j = 0; j < Q.size(); ++j) {
            std::cout << Q[i][j] << "  ";
        }
        std::cout << endl;
    }
    std::cout << "A^T*A对应的特征向量为:" << endl;
    for (int i = 0; i < QT.size(); ++i) {
        for (int j = 0; j < QT.size(); ++j) {
            std::cout << QT[i][j] << "  ";
        }
        std::cout << endl;
    }
    //求sigma矩阵
    vector<vector<double>> sigma(eigen.size(), vector<double>(eigenT.size(), 0.0));
    for (int i = 0; i < eigenT.size(); ++i) {
        sigma[i][i] = eigen[i];
    }
    std::cout << "sigma矩阵为:" << endl;
    for (int i = 0; i < eigen.size(); ++i) {
        for (int j = 0; j < eigenT.size(); ++j) {
            std::cout << sigma[i][j] << "  ";
        }
        std::cout << endl;
    }
    //通过U矩阵确定V矩阵的符号
    for (int i = 0; i < eigenT.size(); ++i) {
        double temp1 = 0.0, temp2 = 0.0;
        for (int j = 0; j < eigenT.size(); ++j) {
            temp1 += ran[i][j] * QT[j][i];
            temp2 += Q[i][j] * sigma[j][i];
        }
        if (fabs(temp1 - temp2) >= 1e-2) {
            for (int k = 0; k < eigenT.size(); ++k) {
                QT[k][i] = -QT[k][i];
            }
        }
    }
    std::cout << "V^T矩阵为:" << endl;
    for (int i = 0; i < eigenT.size(); ++i) {
        for (int j = 0; j < eigenT.size(); ++j) {
            std::cout << QT[j][i] << "  ";
        }
        std::cout << endl;
    }
    std::cout << "此时A*A^T对应的特征向量即为U矩阵" << endl;
    //考虑另外一个问题的解决
    //读取数据
    vector<vector<double>> iris;
    ifstream irisdata;
    irisdata.open("C:\\c\\c_study\\iris.txt");
    string line;
    while (getline(irisdata, line)) {
        stringstream ss; 
        vector<double> tempdata;
        ss << line; 
        if (!ss.eof()) {
            double temp;
            while (ss >> temp) {
                tempdata.push_back(temp);
                char ch;
                if (!(ss >> ch)) break;
            }
        }
        iris.push_back(tempdata);
    }
    irisdata.close();
    std::cout << "读入数据为:" << endl;
    for (int i = 0; i < iris.size(); ++i) {
        for (int j = 0; j < iris[0].size(); ++j) {
            std::cout << iris[i][j] << "  ";
        }
        std::cout << endl;
    } //确定数据读入成功
    //去中心化
    vector<double> middle(iris.size(), 0.0);
    for (int i = 0; i < iris.size(); ++i) {
        double temp = 0.0;
        for (int j = 0; j < iris[0].size() - 1; ++j) {
            temp += iris[i][j];
        }
        middle[i] = temp / 4;
    }
    for (int i = 0; i < iris.size(); ++i) {
        for (int j = 0; j < iris[0].size() - 1; ++j) {
            iris[i][j] -= middle[i];
        }
    }
    std::cout << "将数据去中心化后为:" << endl;
    for (int i = 0; i < iris.size(); ++i) {
        for (int j = 0; j < iris[0].size(); ++j) {
            std::cout << iris[i][j] << "  ";
        }
        std::cout << endl;
    }
    //计算协方差矩阵
    int m = iris.size(), n = iris[0].size() - 1;
    vector<vector<double>> covar(n, vector<double>(n, 0.0));
    for (int i = 0; i < covar.size(); ++i) {
        for (int j = 0; j < covar[0].size(); ++j) {
            for (int k = 0; k < m; ++k) {
                covar[i][j] += iris[k][i] * iris[k][j];
            }
            covar[i][j] /= m;
        }
    }
    std::cout << "协方差矩阵为:" << endl;
    for (int i = 0; i < covar.size(); ++i) {
        for (int j = 0; j < covar[0].size(); ++j) {
            std::cout << covar[i][j] << "  ";
        }
        std::cout << endl;
    }
    //Jacobi方法求其特征值特征向量
    vector<vector<double>> QX(covar.size(), vector<double>(covar.size(), 0.0));
    for (int i = 0; i < QX.size(); ++i) {
        QX[i][i] = 1.0;
    }
    while (matrixSum(covar) >= epsilon) {
        int poix = 0, poiy = 0;
        double maxA = matrixMax(covar, poix, poiy);
        double s = (covar[poiy][poiy] - covar[poix][poix]) / (2 * covar[poix][poiy]);
        double t = 0.0;
        if (s == 0.0) {
            t = 1.0;
        }
        else {
            double t1 = -s - sqrt(s * s + 1.0);
            double t2 = -s + sqrt(s * s + 1.0);
            t = (fabs(t1) > fabs(t2)) ? t2 : t1;
        }
        double c =  1.0 / sqrt(t * t + 1.0);
        double d = t / sqrt(t * t + 1.0);
        vector<vector<double>> ranA1(covar);
        for (int i = 0; i < covar.size(); ++i) {
            if (i == poix || i == poiy) continue;
            covar[i][poix] = c * ranA1[poix][i] - d * ranA1[poiy][i];
            covar[poix][i] = covar[i][poix];
            covar[i][poiy] = c * ranA1[poiy][i] + d * ranA1[poix][i];
            covar[poiy][i] = covar[i][poiy];
        }
        covar[poix][poix] = ranA1[poix][poix] - t * ranA1[poix][poiy];
        covar[poiy][poiy] = ranA1[poiy][poiy] + t * ranA1[poix][poiy];
        covar[poix][poiy] = 0.0;
        covar[poiy][poix] = 0.0;
        vector<vector<double>> Q1(QX.size(), vector<double>(QX.size(), 0.0));
        for (int i = 0; i < QX.size(); ++i) Q1[i][i] = 1.0;
        Q1[poix][poix] = c;
        Q1[poix][poiy] = d;
        Q1[poiy][poix] = -d;
        Q1[poiy][poiy] = c;
        decltype(Q1) Q2 = QX;
        for (int i = 0; i < QX.size(); ++i) {
            for (int j = 0; j < QX.size(); ++j) {
                double temp = 0.0;
                for (int k = 0; k < QX.size(); ++k) {
                    temp += Q2[i][k] * Q1[k][j];
                }
                QX[i][j] = temp;
            }
        }
    }
    vector<double> eigenX(covar.size(), 0.0);
    for (int i = 0; i < covar.size(); ++i) {
        eigenX[i] = covar[i][i];
    }
    std::cout << "协方差矩阵的特征值是:" << endl;
    for (auto lamda : eigenX) {
        std::cout << lamda << "  "; 
    }
    std::cout << endl;
    //对特征值特征向量进行排序
    for (int i = 0; i < eigenX.size() - 1; ++i) {
        for (int j = eigenX.size() - 1; j > i; --j) {
            if (eigenX[j] > eigenX[j - 1]) {
                double temp = eigenX[j];
                eigenX[j] = eigenX[j - 1];
                eigenX[j - 1] = temp;
                for (int k = 0; k < eigenX.size(); ++k) {
                    double temp = QX[k][j];
                    QX[k][j] = QX[k][j - 1];
                    QX[k][j - 1] = temp;
                }
            }
        }
    }
    std::cout << "协方差矩阵的特征值是:" << endl;
    for (auto lamda : eigenX) {
        std::cout << lamda << "  "; 
    }
    std::cout << endl;
    std::cout << "对应的特征向量为:" << endl;
    for (int i = 0; i < QX.size(); ++i) {
        for (int j = 0; j < QX.size(); ++j) {
            std::cout << QX[i][j] << "  ";
        }
        std::cout << endl;
    }
    //抽取主方向特征向量
    vector<double> x1, x2;
    for (int i = 0; i < QX.size(); ++i) {
        x1.push_back(QX[i][0]);
        x2.push_back(QX[i][1]);
    }
    //计算坐标
    vector<vector<double>>  project(iris.size(), vector<double>(3, 0.0));
    for (int i = 0; i < iris.size(); ++i) {
        project[i][2] = iris[i][4];
    }
    for (int i = 0; i < iris.size(); ++i) {
        double temp1 = 0.0;
        for (int j = 0; j < 4; ++j) {
            temp1 += x1[j] * iris[i][j];
        }
        project[i][0] = temp1;
        double temp2 = 0.0;
        for (int j = 0; j < 4; ++j) {
            temp2 += x2[j] * iris[i][j];
        }
        project[i][1] = temp2;
    }
    std::cout << "PCA后的坐标是:" << endl;
    for (int i = 0; i < project.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            std::cout << project[i][j] << "  ";
        }
        std::cout << endl;
    }
    //将结果写入txt文件，尝试用python进行可视化分析
    ofstream projectdata;
    projectdata.open("C:\\c\\c_study\\irisdata.txt");
    for (int i = 0; i < project.size(); ++i) {
        projectdata<< project[i][0] << " " << project[i][1] << " " << project[i][2] << endl;
    }
    projectdata.close();
    return 0;
}