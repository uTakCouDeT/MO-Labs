#include <iostream>
#include <vector>
#include <iomanip>


/// Преобразование в матрицу
std::vector<std::vector<double>>
ConvertToMatrix(std::vector<double> c, std::vector<std::vector<double>> A, std::vector<double> b, double F = 0) {
    std::vector<std::vector<double>> SimplexTable = A;
    int k = 0;
    for (std::vector<double> &i: SimplexTable) {
        i.push_back(b[k]);
        ++k;
    }
    SimplexTable.push_back(c);
    SimplexTable[SimplexTable.size() - 1].push_back(F);
    return SimplexTable;
}

/// Вывод условия в каноническом виде
void CanonicalPrint(std::vector<double> c, std::vector<std::vector<double>> A, std::vector<double> b,
                    std::string extremum = "max") {

    std::vector<int> RowXmas;
    std::vector<int> ColXmas;
    bool Max;

    // Заполнение массива переменных
    for (int i = 0; i < c.size(); ++i) {
        RowXmas.push_back(i + 1);
    }
    for (int i = c.size(); i < b.size() + c.size(); ++i) {
        ColXmas.push_back(i + 1);
    }

    // Инвертирование коэфициентов функции
    for (double &i: c) {
        i = -i;
    }

    // Минимум или максимум
    if (extremum == "max") {
        Max = true;
    } else if (extremum == "min") {
        Max = false;
    } else {
        std::cerr << "Enter: max/min ";
        exit(1);
    }

    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    std::cout << "F   =   ";
    std::cout << std::fixed << std::setprecision(2) << -c[0] << " * x" << RowXmas[0];
    for (int i = 1; i < c.size(); ++i) {
        if (c[i] <= 0) { std::cout << "   +   "; } else { std::cout << "   -   "; }
        std::cout << std::fixed << std::setprecision(2) << -c[i] << " * x" << RowXmas[i];
    }
    std::cout << "  -->  ";
    if (Max) { std::cout << "max"; } else { std::cout << "min"; }
    std::cout << "\n------------------------------------------------------------------\n";
    auto Matrix = ConvertToMatrix(c, A, b);
    for (int i = 0; i < Matrix.size() - 1; ++i) {
        std::cout << std::fixed << std::setprecision(2) << Matrix[i][0] << " * x" << RowXmas[0];
        for (int j = 1; j < Matrix[i].size() - 1; ++j) {
            if (Matrix[i][j] >= 0) { std::cout << "   +   "; } else { std::cout << "   -   "; }
            std::cout << std::fixed << std::setprecision(2) << std::abs(Matrix[i][j])
                      << " * x" << RowXmas[j];
        }
        std::cout << "   <=  ";
        if (Matrix[i][Matrix[i].size() - 1] >= 0) { std::cout << ' '; }
        std::cout << Matrix[i][Matrix[i].size() - 1];
        std::cout << std::endl;
    }
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
}

double FunctionValue(std::vector<double> c, std::vector<double> X) {
    double F = 0;
    if (c.size() == X.size()) {
        for (int i = 0; i < c.size(); ++i) {
            F += c[i] * X[i];
        }
    }
    return F;
}

bool Condition(std::vector<std::vector<double>> A, std::vector<double> b, std::vector<double> X) {
    bool res = true;
    for (int i = 0; i < A.size(); ++i) {
        if (FunctionValue(A[i], X) > b[i]) {
            res = false;
        }
    }
    return res;
}


void BruteForce(std::vector<double> c, std::vector<std::vector<double>> A, std::vector<double> b,
                std::string extremum = "max") {

    // Минимум или максимум
    bool Max = true;
    if (extremum == "max") {
        Max = true;
    } else if (extremum == "min") {
        Max = false;
    } else {
        std::cerr << "Enter: max/min ";
        exit(1);
    }
    double X1, X2, X3 = 0;
    double F = FunctionValue(c, {X1, X2, X3});

    if (c.size() == 3) {
        std::cout << std::endl;
        for (double x1 = 0; x1 < 100; ++x1) {
            for (double x2 = 0; x2 < 100; ++x2) {
                for (double x3 = 0; x3 < 100; ++x3) {
                    if (Condition(A, b, {x1, x2, x3})) {
                        double f = FunctionValue(c, {x1, x2, x3});
                        std::cout << std::setprecision(0) << "(" << x1 << ", " << x2 << ", " << x3
                                  << ") F = " << f << std::endl;
                        if (Max && f > F || !Max && f < F) {
                            F = f;
                            X1 = x1;
                            X2 = x2;
                            X3 = x3;
                        }
                    }
                }
            }
        }
    }
    std::cout << "\nOptimal integer solution:\n" << "(" << X1 << ", " << X2 << ", " << X3 << ") F = " << F << std::endl;
}

int main() {
    std::vector<double> c = {1, 5, 5};
    std::vector<std::vector<double>> A = {{4, 1,   1},
                                          {1, 4,   0},
                                          {0, 0.5, 4}};
    std::vector<double> b = {5, 7, 6};

    CanonicalPrint(c, A, b, "max");
    BruteForce(c, A, b, "max");

}
