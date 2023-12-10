#include <iostream>
#include <cmath>
using namespace std;

// Тип для указателя на функцию одной переменной
typedef double (*FuncPtr)(double);

// Тип для указателя на функцию двух переменных
typedef double (*Func2Ptr)(double, double);

double trapezoidal(FuncPtr f, double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.5 * (f(a) + f(b));
    for (int i = 1; i < n; i++) {
        sum += f(a + i * h);
    }
    return h * sum;
}

double simpsons(FuncPtr f, double a, double b, int n) {
    double h = (b - a) / n;
    double sum = f(a) + f(b);

    for (int i = 1; i < n; i += 2) {
        sum += 4 * f(a + i * h);
    }

    for (int i = 2; i < n - 1; i += 2) {
        sum += 2 * f(a + i * h);
    }

    return h / 3 * sum;
}

double simpsons2(Func2Ptr f, double ax, double bx, int nx, double ay, double by, int ny) {
    double hx = (bx - ax) / nx;
    double hy = (by - ay) / ny;
    double sum = 0.0;

    for (int i = 0; i <= nx; i++) {
        for (int j = 0; j <= ny; j++) {
            double x = ax + i * hx;
            double y = ay + j * hy;
            double weight = (i == 0 || i == nx) ? 1 : (i % 2 == 1) ? 4 : 2;
            weight *= (j == 0 || j == ny) ? 1 : (j % 2 == 1) ? 4 : 2;
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

    double integralTrapezoidal = trapezoidal(func, a, b, n);
    double integralSimpson = simpsons(func, a, b, n);

    cout << "Integral (Trapezoidal): " << integralTrapezoidal << endl;
    cout << "Integral (Simpson): " << integralSimpson << endl;

    double ax = a, bx = b, ay = a, by = b;
    int nx = 10, ny = 10; 

    double integral2 = simpsons2(func2, ax, bx, nx, ay, by, ny);

    cout << "Double integral (Simpson's rule): " << integral2 << endl;

    return 0;
}
