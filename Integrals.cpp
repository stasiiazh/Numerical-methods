#include <iostream>
#include <cmath>

using namespace std;

// Тип для указателя на функцию одной переменной
typedef double (*FuncPtr)(double);

// Тип для указателя на функцию двух переменных
typedef double (*Func2Ptr)(double, double);

double trapezoidal(FuncPtr f, double a, double b, int n, double eps) {
    double integral_n, integral_2n, h_n, h_2n, sum_n, sum_2n;

    do {
        h_n = (b - a) / n;
        sum_n = 0.5 * (f(a) + f(b));
        for (int i = 1; i < n; i++) {
            sum_n += f(a + i * h_n);
        }
        integral_n = h_n * sum_n;

        h_2n = (b - a) / (2 * n);
        sum_2n = 0.5 * (f(a) + f(b));
        for (int i = 1; i < 2 * n; i++) {
            sum_2n += f(a + i * h_2n);
        }
        integral_2n = h_2n * sum_2n;

        if (fabs(integral_n - integral_2n) >= 3 * eps) {
            n *= 2;
        }
    } while (fabs(integral_n - integral_2n) >= 3 * eps);

    return integral_2n; // Возвращаем значение
}

double simpsons(FuncPtr f, double a, double b, int n, double eps) {
    double integral, integral_prev;
    integral = integral_prev = 0.0;

    do {
        integral_prev = integral;
        n = (n % 2 == 0) ? n : n + 1; // Убедимся, что n чётное
        double h = (b - a) / n;
        double sum = f(a) + f(b);

        for (int i = 1; i < n; i += 2) {
            sum += 4 * f(a + i * h);
        }
        for (int i = 2; i < n; i += 2) {
            sum += 2 * f(a + i * h);
        }

        integral = h / 3 * sum;
        n *= 2;
    } while (fabs(integral - integral_prev) > 15 * eps);

    return integral_prev;
}

double simpsons2(Func2Ptr f, double ax, double bx, int nx, double ay, double by, int ny) {
    double hx = (bx - ax) / nx;
    double hy = (by - ay) / ny;
    double sum = 0.0;

    for (int i = 0; i <= nx; i++) {
        for (int j = 0; j <= ny; j++) {
            double x = ax + i * hx;
            double y = ay + j * hy;
            double weight = (i == 0 or i == nx) ? 1 : (i % 2 == 1) ? 4 : 2;
            weight *= (j == 0 or j == ny) ? 1 : (j % 2 == 1) ? 4 : 2;
            sum += weight * f(x, y);
        }
    }
    return hx * hy / 9 * sum;
}

//одной переменной
double func(double x) {
    return 1 / sqrt(x*x*x - 1);
}

//двух переменных
double func2(double x, double y) {
    return 1 / sqrt(x*x*x - 1);
}

int main() {
    double a = 1.3, b = 2.621;
    int n = 100;
    double eps = 0.0001; 

    double integralTrapezoidal = trapezoidal(func, a, b, n, eps);
    cout << "Integral (Trapezoidal): " << integralTrapezoidal << endl;

    double integralSimpson = simpsons(func, a, b, n, eps);
    cout << "Integral (Simpson): " << integralSimpson << endl;

    double ax = a, bx = b, ay = a, by = b;
    int nx = 10, ny = 10; 

    double integral2 = simpsons2(func2, ax, bx, nx, ay, by, ny);
    cout << "Double integral (Simpson's rule): " << integral2 << endl;

    return 0;
}

