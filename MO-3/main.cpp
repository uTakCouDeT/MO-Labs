#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <set>
#include <cmath>

//Variant 23

class Simplex {
    std::vector<double> c;
    std::vector<std::vector<double>> A;
    std::vector<double> b;
    std::vector<int> RowXmas;
    std::vector<int> ColXmas;
    double F;
    int r;
    int k;
    bool Max;
    bool SolutionExists;

public:

    Simplex() : F(0), r(0), k(0), Max(true), SolutionExists(true) {}

    /// Ввод данных
    Simplex(std::vector<double> _c, std::vector<std::vector<double>> a,
            std::vector<double> _b, const std::string &extremum = "max") : c(std::move(_c)), A(std::move(a)),
                                                                           b(std::move(_b)), F(0), r(0), k(0),
                                                                           SolutionExists(true) {
        // Проверка входных данных
        if (A.size() == b.size()) {
            for (auto &i: A) {
                if (i.size() != c.size()) {
                    std::cerr << "Invalid arguments!";
                    exit(1);
                }
            }
        } else {
            std::cerr << "Invalid arguments!";
            exit(1);
        }

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
    }

    bool GetSolutionExists() const {
        return SolutionExists;
    }

    double GetFunction() const {
        return F;
    }

    std::vector<double> GetVecb() const {
        return b;
    }

    std::vector<int> GetXmas() const {
        return ColXmas;
    }

    /// Преобразование в матрицу
    std::vector<std::vector<double>> ConvertToMatrix() {
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

    /// Поиск опорного решения
    bool CheckSolution() {
        bool isSolution = true;

        // Проверка свободных членов
        for (int i = 0; i < b.size(); ++i) {
            if (b[i] < 0) {
                isSolution = false;
                r = i;
                for (int j = r + 1; j < b.size(); ++j) {
                    if (b[j] < 0 && std::abs(b[j]) > std::abs(b[r])) {
                        r = j;
                    }
                }
                break;
            }
        }

        if (!isSolution) {

            // Поиск столбца с отрицательным элементом
            SolutionExists = false;
            for (int i = 0; i < A[r].size(); ++i) {
                if (A[r][i] < 0) {
                    k = i;
                    SolutionExists = true;
                    for (int j = k + 1; j < A[r].size(); ++j) {
                        if (A[r][j] < 0 && std::abs(A[r][j]) >= std::abs(A[r][k])) {
                            k = j;
                        }
                    }
                    break;
                }
            }

            // Поиск минимального положительного отношения
            double Ratio;
            for (int i = 0; i < A.size(); ++i) {
                if (A[i][k] != 0 && b[i] / A[i][k] > 0) {
                    Ratio = b[i] / A[i][k];
                    r = i;
                    for (int j = r + 1; j < A.size(); ++j) {
                        if (A[j][k] != 0 && b[j] / A[j][k] > 0 && b[j] / A[j][k] < Ratio) {
                            Ratio = b[j] / A[j][k];
                            r = j;
                        }
                    }
                    break;
                }
            }
        }
        return isSolution;
    }

    /// Поиск оптимального решения
    bool CheckOptimality() {
        bool isOptimal = true;

        // Поиск максимального по модулю элемента неподходящего знака
        for (int i = 0; i < c.size(); ++i) {
            if ((Max && c[i] <= 0) || (!Max && c[i] >= 0)) {
                isOptimal = false;
                k = i;
                for (int j = i + 1; j < c.size(); ++j) {
                    if (((Max && c[j] <= 0) || (!Max && c[j] >= 0)) && std::abs(c[j]) >= std::abs(c[k])) {
                        k = j;
                    }
                }
                break;
            }
        }

        if (!isOptimal) {

            // Поиск минимального положительного отношения
            for (int i = 0; i < A.size(); ++i) {
                if (A[i][k] != 0 && b[i] / A[i][k] > 0) {
                    double Ratio = b[i] / A[i][k];
                    r = i;
                    for (int j = i + 1; j < A.size(); ++j) {
                        if (A[j][k] != 0 && b[j] / A[j][k] > 0 && b[j] / A[j][k] < Ratio) {
                            r = j;
                            Ratio = b[j] / A[j][k];
                        }
                    }
                    break;
                }
            }
        }
        return isOptimal;
    }

    /// Преобразование симплекс таблицы
    void SimplexTableFormation(bool Silence = false) {
        if (!Silence) {
            std::cout << "The resolving row: x" << ColXmas[r]
                      << "\nThe resolving col: x" << RowXmas[k]
                      << "\n=====================\n" << std::endl;
        }

        // Преобразование в матрицу
        auto Matrix = ConvertToMatrix();
        auto NewMatrix = Matrix;

        // Преобразование матрицы
        for (int i = 0; i < Matrix.size(); ++i) {
            for (int j = 0; j < Matrix[i].size(); ++j) {
                if (i == r && j == k) {
                    NewMatrix[i][j] = 1 / Matrix[r][k];
                } else if (i == r) {
                    NewMatrix[i][j] = Matrix[i][j] / Matrix[r][k];
                } else if (j == k) {
                    NewMatrix[i][j] = -Matrix[i][j] / Matrix[r][k];
                } else {
                    NewMatrix[i][j] = Matrix[i][j] - Matrix[i][k] * Matrix[r][j] / Matrix[r][k];
                }
            }
        }

        // Преобразование из матрицы
        for (int i = 0; i < c.size(); ++i) {
            c[i] = NewMatrix[NewMatrix.size() - 1][i];
        }
        F = NewMatrix[NewMatrix.size() - 1][NewMatrix[NewMatrix.size() - 1].size() - 1];
        for (int i = 0; i < b.size(); ++i) {
            b[i] = NewMatrix[i][NewMatrix[i].size() - 1];
        }
        for (int i = 0; i < A.size(); ++i) {
            for (int j = 0; j < A[i].size(); ++j) {
                A[i][j] = NewMatrix[i][j];
            }
        }
        std::swap(RowXmas[k], ColXmas[r]);
        r = k = 0;
    }

    /// Симплекс метод
    void SimplexMethod(bool Silence = false) {
        std::cout << std::endl;
        CanonicalPrint();
        if (!Silence) {
            std::cout << std::endl;
            Print();
        }

        // Нахождение опорного решения
        while (SolutionExists) {
            if (CheckSolution()) {
                if (!Silence) {
                    std::cout << "A reference solution is found\n======================================"
                              << std::endl;
                }
                break;
            } else {
                SimplexTableFormation(Silence);
                if (!Silence) { Print(); }
            }
        }

        // Вывод при отсутствии опорного решения
        if (!SolutionExists) {
            std::cout << "No solution exists!";
        }

        // Нахождение оптимального решения
        while (SolutionExists) {
            if (CheckOptimality()) {
                std::cout << "\nF = " << F << std::endl;
                for (int i = 0; i < c.size(); ++i) {
                    bool zero = true;
                    for (int j = 0; j < ColXmas.size(); ++j) {
                        if (ColXmas[j] == (i + 1)) {
                            std::cout << "x" << i + 1 << " = " << b[j] << std::endl;
                            zero = false;
                            break;
                        }
                    }
                    if (zero) {
                        std::cout << "x" << i + 1 << " = " << 0 << std::endl;
                    }
                }
                break;
            } else {
                SimplexTableFormation(Silence);
                if (!Silence) { Print(); }
            }
        }
    }

    /// Вывод симплекс таблицы
    void Print() {
        auto Matrix = ConvertToMatrix();
        std::cout << "======================================\n\t";
        for (int i = 0; i < Matrix[0].size() - 1; ++i) {
            std::cout << " x" << RowXmas[i] << "\t";
        }
        std::cout << " Si" << std::endl;
        for (int i = 0; i < Matrix.size() - 1; ++i) {
            std::cout << "x" << ColXmas[i];
            for (double j: Matrix[i]) {
                if (j == 0) { j = std::abs(j); }
                std::cout << '\t';
                if (j >= 0) { std::cout << ' '; }
                std::cout << std::fixed << std::setprecision(2) << j;
            }
            std::cout << std::endl;
        }
        std::cout << "F";
        for (double i: Matrix[Matrix.size() - 1]) {
            if (i == 0) { i = std::abs(i); }
            std::cout << '\t';
            if (i >= 0) { std::cout << ' '; }
            std::cout << std::fixed << std::setprecision(2) << i;
        }
        std::cout << "\n======================================" << std::endl;
    }

    /// Вывод условия в каноническом виде
    void CanonicalPrint() {
        std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
        std::cout << "F   =   ";
        for (int i = 0; i < c.size() - 1; ++i) {
            std::cout << std::fixed << std::setprecision(2)
                      << -c[i] << " * x" << RowXmas[i] << "   +   ";
        }
        std::cout << std::fixed << std::setprecision(2)
                  << -c[c.size() - 1] << " * x" << RowXmas[c.size() - 1] << "  -->  ";
        if (Max) { std::cout << "max"; } else { std::cout << "min"; }
        std::cout << "\n------------------------------------------------------------------\n";
        auto Matrix = ConvertToMatrix();
        for (int i = 0; i < Matrix.size() - 1; ++i) {
            for (int j = 0; j < Matrix[i].size() - 2; ++j) {
                std::cout << std::fixed << std::setprecision(2)
                          << Matrix[i][j] << " * x" << RowXmas[j] << "   +   ";
            }
            std::cout << std::fixed << std::setprecision(2) << Matrix[i][Matrix[i].size() - 2]
                      << " * x" << RowXmas[Matrix[i].size() - 2]
                      << "   <=   " << Matrix[i][Matrix[i].size() - 1] << std::endl;
        }
        std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    }

    bool operator==(const Simplex &simplex) const {
        return (c == simplex.c && A == simplex.A && b == simplex.b && F == simplex.F &&
                RowXmas == simplex.RowXmas && ColXmas == simplex.ColXmas && Max == simplex.Max);
    }

    bool operator<(const Simplex &simplex) const {
        return F < simplex.F;
    }

};

class MVG {
    std::vector<double> c;
    std::vector<std::vector<double>> A;
    std::vector<double> b;
    double F;
    bool Max;
    bool SolutionExists;
    std::set<Simplex> Solutions;
public:

    MVG() : F(0), Max(true), SolutionExists(true) {}

    /// Ввод данных
    MVG(std::vector<double> _c, std::vector<std::vector<double>> a,
        std::vector<double> _b, const std::string &extremum = "max") : c(std::move(_c)), A(std::move(a)),
                                                                       b(std::move(_b)), F(0),
                                                                       SolutionExists(true) {
        // Проверка входных данных
        if (A.size() == b.size()) {
            for (auto &i: A) {
                if (i.size() != c.size()) {
                    std::cerr << "Invalid arguments!";
                    exit(1);
                }
            }
        } else {
            std::cerr << "Invalid arguments!";
            exit(1);
        }

        // Минимум или максимум
        if (extremum == "max") { Max = true; }
        else if (extremum == "min") { Max = false; }
        else {
            std::cerr << "Enter: max/min ";
            exit(1);
        }
    }

    void BranchAndBoundaryMethod(bool Silence = true) {
        std::cout << "\nFirst Simplex Method:";
        std::string extremum;
        if (Max) { extremum = "max"; } else { extremum = "min"; }
        Simplex T = Simplex(c, A, b, extremum);
        T.SimplexMethod(Silence);
        if (T.GetSolutionExists()) {
            bool IsIntegerSolution = true;
            int BranchingVariable = 0;
            double BranchingVariableValue = 0;
            std::vector<double> Variables = T.GetVecb();
            for (int i = 0; i < Variables.size(); ++i) {
                if (Variables[i] != floor(Variables[i])) {
                    IsIntegerSolution = false;
                    BranchingVariable = T.GetXmas()[i];
                    BranchingVariableValue = Variables[i];
                    break;
                }
            }
            if (IsIntegerSolution) { Solutions.insert(T); }
            else {
                std::cout << "\n==================================================================\n";
                std::cout << "First branching variable: x" << BranchingVariable << " = "
                          << BranchingVariableValue << std::endl;
                std::cout << "------------------------------------------------------------------\n";
            }
        }
    }

    void PrintSet() {
        for (const Simplex &i: Solutions) {
            std::cout << i.GetFunction();
        }
    }

};

int main() {
    std::vector<double> c = {1, 5, 5};
    std::vector<std::vector<double>> A = {{4, 1,   1},
                                          {1, 4,   0},
                                          {0, 0.5, 4}};
    std::vector<double> b = {5, 7, 6};

    auto mvg = MVG(c, A, b, "max");
    //mvg.BranchAndBoundaryMethod();


    //std::cout << "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nExample for 2x3:\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";

    std::vector<double> c1 = {-12, 1};
    std::vector<std::vector<double>> A1 = {{6, -1},
                                           {2, 5},};
    std::vector<double> b1 = {12, 20};

    auto mvg1 = MVG(c1, A1, b1, "min");
    mvg1.BranchAndBoundaryMethod();

    return 0;
}