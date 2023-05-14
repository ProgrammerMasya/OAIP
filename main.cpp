#include <iostream>
#include <vector>
#include <cmath>

// Функция для разложения Холецкого
void choleskyDecomposition(const std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& L) {
    int n = A.size();

    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0.0;
            if (j == i) {
                for (int k = 0; k < j; k++) {
                    sum += std::pow(L[j][k], 2);
                }
                L[j][j] = std::sqrt(A[j][j] - sum);
            } else {
                for (int k = 0; k < j; k++) {
                    sum += L[i][k] * L[j][k];
                }
                L[i][j] = (A[i][j] - sum) / L[j][j];
            }
        }
    }
}

// Функция для решения системы уравнений методом квадратного корня
std::vector<double> solveEquations(const std::vector<std::vector<double>>& L, const std::vector<double>& b) {
    int n = L.size();
    std::vector<double> y(n, 0.0);
    std::vector<double> x(n, 0.0);

    // Решение Ly = b
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < i; j++) {
            sum += L[i][j] * y[j];
        }
        y[i] = (b[i] - sum) / L[i][i];
    }

    // Решение L^T * x = y
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += L[j][i] * x[j];
        }
        x[i] = (y[i] - sum) / L[i][i];
    }

    return x;
}

int main() {
    // Ввод размерности системы уравнений
    int n;
    std::cout << "Введите размерность системы уравнений: ";
    std::cin >> n;

    // Ввод матрицы коэффициентов A
    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
    std::cout << "Введите матрицу коэффициентов A:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cin >> A[i][j];
        }
    }

    // Ввод вектора свободных членов b
    std::vector<double> b(n, 0.0);
    std::cout << "Введите вектор свободных членов b:\n";
    for (int i = 0; i < n; i++) {
        std::cin >> b[i];
    }

    // Создание и инициализация матрицы L
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));

    // Разложение Холецкого
    choleskyDecomposition(A, L);

    // Решение системы уравнений
    std::vector<double> solution = solveEquations(L, b);

    // Вывод результата
    std::cout << "Решение системы уравнений:\n";
    for (int i = 0; i < n; i++) {
        std::cout << "x[" << i << "] = " << solution[i] << "\n";
    }

    return 0;
}