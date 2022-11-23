#include <iostream>
#include <vector>
#include <string>
#include <iomanip>

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
    bool solutionExists;

public:

    Simplex() : F(0), r(0), k(0), Max(true), solutionExists(true) {}

    /// Ввод данных
    Simplex(std::vector<double> _c, std::vector<std::vector<double>> a,
            std::vector<double> _b, const std::string &extremum = "max") : c(std::move(_c)), A(std::move(a)),
                                                                           b(std::move(_b)), F(0), r(0), k(0),
                                                                           solutionExists(true) {
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
            solutionExists = false;
            for (int i = 0; i < A[r].size(); ++i) {
                if (A[r][i] < 0) {
                    k = i;
                    solutionExists = true;
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
                if (std::abs(c[i]) >= std::abs(c[k])) {
                    k = i;
                }
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
    void SimplexTableFormation(const std::string &variable = "x", bool Silence = false) {
        if (!Silence) {
            std::cout << "The resolving row: " << variable << ColXmas[r]
                      << "\nThe resolving col: " << variable << RowXmas[k]
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
    void SimplexMethod(const std::string &function = "F", const std::string &variable = "x", bool Silence = false) {
        std::cout << std::endl;
        CanonicalTransformation(function, variable);
        std::cout << std::endl;
        if (!Silence) { Print(function, variable); }

        // Нахождение опорного решения
        while (solutionExists) {
            if (CheckSolution()) {
                if (!Silence) {
                    std::cout << "A reference solution is found\n======================================" << std::endl;
                }
                break;
            } else {
                SimplexTableFormation(variable, Silence);
                if (!Silence) { Print(function, variable); }
            }
        }

        // Вывод при отсутствии опорного решения
        if (!solutionExists){
            std::cout << "No solution exists!";
        }

        // Нахождение оптимального решения
        while (solutionExists) {
            if (CheckOptimality()) {
                if (Silence) { Print(function, variable); }
                std::cout << "\n" << function << " = " << F << std::endl;
                for (int i = 0; i < c.size(); ++i) {
                    bool zero = true;
                    for (int j = 0; j < ColXmas.size(); ++j) {
                        if (ColXmas[j] == (i + 1)) {
                            std::cout << variable << i + 1 << " = " << b[j] << std::endl;
                            zero = false;
                            break;
                        }
                    }
                    if (zero) {
                        std::cout << variable << i + 1 << " = " << 0 << std::endl;
                    }
                }
                break;
            } else {
                SimplexTableFormation(variable, Silence);
                if (!Silence) { Print(function, variable); }
            }
        }
    }

    /// Вывод симплекс таблицы
    void Print(const std::string &function = "F", const std::string &variable = "x") {
        auto Matrix = ConvertToMatrix();
        std::cout << "======================================\n\t";
        for (int i = 0; i < Matrix[0].size() - 1; ++i) {
            std::cout << ' ' << variable << RowXmas[i] << "\t";
        }
        std::cout << " Si" << std::endl;
        for (int i = 0; i < Matrix.size() - 1; ++i) {
            std::cout << variable << ColXmas[i];
            for (double j: Matrix[i]) {
                if (j == 0) { j = std::abs(j); }
                std::cout << '\t';
                if (j >= 0) { std::cout << ' '; }
                std::cout << std::fixed << std::setprecision(2) << j;
            }
            std::cout << std::endl;
        }
        std::cout << function;
        for (double i: Matrix[Matrix.size() - 1]) {
            if (i == 0) { i = std::abs(i); }
            std::cout << '\t';
            if (i >= 0) { std::cout << ' '; }
            std::cout << std::fixed << std::setprecision(2) << i;
        }
        std::cout << "\n======================================" << std::endl;
    }

    /// Вывод условия в каноническом виде
    void CanonicalPrint(const std::string &function = "F", const std::string &variable = "x", bool Dual = false) {
        std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
        std::cout << function << "   =   ";
        for (int i = 0; i < c.size() - 1; ++i) {
            std::cout << std::fixed << std::setprecision(2) << -c[i] << " * " << variable << RowXmas[i] << "   +   ";
        }
        std::cout << std::fixed << std::setprecision(2) << -c[c.size() - 1]
                  << " * " << variable << RowXmas[c.size() - 1] << "  -->  ";
        if (Max) { std::cout << "max"; } else { std::cout << "min"; }
        std::cout << "\n------------------------------------------------------------------\n";

        auto Matrix = ConvertToMatrix();
        for (int i = 0; i < Matrix.size() - 1; ++i) {
            std::cout << std::fixed << std::setprecision(2) << Matrix[i][0] << " * " << variable << RowXmas[0];
            for (int j = 1; j < Matrix[i].size() - 1; ++j) {
                if (Matrix[i][j] >= 0) { std::cout << "   +   "; } else { std::cout << "   -   "; }
                std::cout << std::fixed << std::setprecision(2) << std::abs(Matrix[i][j])
                          << " * " << variable << RowXmas[j];
            }
            std::cout << "   <=   " << Matrix[i][Matrix[i].size() - 1] << std::endl;
        }
        std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    }

    /// Вывод преобразований в каноническом виде
    void CanonicalTransformation(const std::string &function = "F", const std::string &variable = "x") {
        CanonicalPrint(function, variable);
        auto Matrix = ConvertToMatrix();
        for (int i = 0; i < Matrix.size() - 1; ++i) {
            std::cout << std::fixed << std::setprecision(2) << Matrix[i][0] << " * " << variable << RowXmas[0];
            for (int j = 1; j < Matrix[i].size() - 1; ++j) {
                if (Matrix[i][j] >= 0) { std::cout << "   +   "; } else { std::cout << "   -   "; }
                std::cout << std::fixed << std::setprecision(2) << std::abs(Matrix[i][j])
                          << " * " << variable << RowXmas[j];
            }
            std::cout << "   +   " << variable << ColXmas[i] << "   =   "
                      << Matrix[i][Matrix[i].size() - 1] << std::endl;
        }
        std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

        for (int i = 0; i < Matrix.size() - 1; ++i) {
            std::cout << variable << ColXmas[i] << "   =   " << Matrix[i][Matrix[i].size() - 1] << "   -   (";
            std::cout << std::fixed << std::setprecision(2)
                      << Matrix[i][0] << " * " << variable << RowXmas[0];
            for (int j = 1; j < Matrix[i].size() - 1; ++j) {
                if (Matrix[i][j] >= 0) { std::cout << "   +   "; } else { std::cout << "   -   "; }
                std::cout << std::fixed << std::setprecision(2) << std::abs(Matrix[i][j]) << " * " << variable
                          << RowXmas[j];
            }
            std::cout << ")" << std::endl;
        }
        std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    }

    /// Транспонирование марицы
    static std::vector<std::vector<double>> Transpone(std::vector<std::vector<double>> a) {
        std::vector<std::vector<double> > B(a[0].size(), std::vector<double>(a.size()));
        for (int i = 0; i < a.size(); ++i) {
            for (int j = 0; j < a[i].size(); ++j) {
                B[j][i] = -a[i][j];
            }
        }
        return B;
    }

    /// Формулировка и решение двойственной задачи
    void DualTask(bool Silence = false) {
        Max = !Max;
        std::swap(c, b);
        for (double &i: c) {
            i = -i;
        }
        ColXmas = {};
        RowXmas = {};
        for (int i = 0; i < c.size(); ++i) {
            RowXmas.push_back(i + 1);
        }
        for (int i = c.size(); i < b.size() + c.size(); ++i) {
            ColXmas.push_back(i + 1);
        }
        A = Transpone(A);
        SimplexMethod("P", "y", Silence);
    }
};


int main() {
    std::vector<double> c = {1, 5, 5};
    std::vector<std::vector<double>> A = {{4, 1,   1},
                                          {1, 4,   0},
                                          {0, 0.5, 4}};
    std::vector<double> b = {5, 7, 6};

    Simplex SimplexTable = Simplex(c, A, b, "max");
    SimplexTable.SimplexMethod("F", "x");
    Simplex DualSimplexTable = Simplex(c, A, b, "max");
    DualSimplexTable.DualTask();

    //std::cout << "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nExample for 4x2:\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";

    std::vector<double> c1 = {-4, -18, -30, -5};
    std::vector<std::vector<double>> A1 = {{3,  1,  -4, -1},
                                           {-2, -4, -1, 1}};
    std::vector<double> b1 = {-3, -3};

    Simplex SimplexTable1 = Simplex(c1, A1, b1, "max");
    //SimplexTable1.SimplexMethod();
    Simplex SimplexTable2 = Simplex(c1, A1, b1, "max");
    //SimplexTable2.DualTask();

    return 0;
}
