#include <iostream>
#include <vector>
#include <cmath>
#include <limits>

using namespace std;

const int N = 9; // Количество точек
const int M = 5; // Степень многочлена + 1

// Функция для решения СЛАУ методом Гаусса
vector<double> gaussianElimination(vector<vector<double> >& A) {
    int n = A.size();
    vector<double> x(n, 0.0);

    for (int i = 0; i < n - 1; i++) {
        // Поиск максимального элемента в текущем столбце
        for (int k = i + 1; k < n; k++) {
            if (fabs(A[i][i]) < numeric_limits<double>::epsilon()) {
                cerr << "Matrix is singular or nearly singular!" << endl;
                return vector<double>(M, numeric_limits<double>::quiet_NaN());
            }

            double c = -A[k][i] / A[i][i];
            for (int j = i; j < n + 1; j++) {
                A[k][j] += c * A[i][j];
            }
        }
    }

    for (int i = n - 1; i >= 0; i--) {
        x[i] = A[i][n] / A[i][i];
        for (int k = i - 1; k >= 0; k--) {
            A[k][n] -= A[k][i] * x[i];
        }
    }
    return x;
}

int main() {
    vector<double> x_values(N);
    double x_vals[N] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    for (int i = 0; i < N; ++i) {
        x_values[i] = x_vals[i];
    }

    vector<double> y_values(N);
    double y_vals[N] = {1.0, 0.5, 0.333333, 0.25, 0.2, 0.166667, 0.142857, 0.125, 0.111111};
    for (int i = 0; i < N; ++i) {
        y_values[i] = y_vals[i];
    }

    vector<vector<double> > A(M, vector<double>(M + 1, 0.0));
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            for (int k = 0; k < N; k++) {
                A[i][j] += pow(x_values[k], i + j);
            }
        }
        for (int k = 0; k < N; k++) {
            A[i][M] += pow(x_values[k], i) * y_values[k];
        }
    }

    vector<double> coefficients = gaussianElimination(A);

    cout << "Coefficients of the polynomial:" << endl;
    for (int i = 0; i < M; i++) {
        cout << "a" << i << " = " << coefficients[i] << endl;
    }

    return 0;
}
