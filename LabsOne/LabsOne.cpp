#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
using namespace std;
// Функция f(x, v) для системы уравнений
double f(double x, double v) {
    return pow(x, 4) + 5 * pow(x, 2) - 4; //любая другая функция
}

// Функция для метода Рунге-Кутта 4-го порядка
void runge_kutta_4th_order(double x0, double v0, double h, int n_steps) {
    double xn = x0;
    double vn = v0;

    // Векторы для хранения результатов
    vector<double> x_vals;
    vector<double> v_vals;

    x_vals.push_back(xn);
    v_vals.push_back(vn);

    for (int i = 0; i < n_steps; i++) {
        // Вычисление коэффициентов k1, k2, k3, k4
        double k1 = f(xn, vn);
        double k2 = f(xn + h / 2.0, vn + h / 2.0 * k1);
        double k3 = f(xn + h / 2.0, vn + h / 2.0 * k2);
        double k4 = f(xn + h, vn + h * k3);

        // Обновление значений переменных x и v
        xn = xn + h;  // Обновляем x на один шаг
        vn = vn + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);  // Обновляем v

        // Сохраняем значения в векторах
        x_vals.push_back(xn);
        v_vals.push_back(vn);
    }

    // Вывод результатов
    cout << "t\tx(t)\tv(t)\n";
    for (int i = 0; i <= n_steps; i++) {
        cout << i * h << "\t" << x_vals[i] << "\t" << v_vals[i] << "\n";
    }

}

int main() {
    setlocale(LC_ALL, "RUS");
    // Начальные условия
    double x0;  // x(0)
    double v0;// v(0)
    cout << "Введите начальные условия:" << endl;
    cout << "x0 = ";
    cin >> x0;
    cout << "v0 = ";
    cin >> v0;
    cout << endl;

    // Шаг и количество шагов
    double h;  // Шаг по времени
    int n_steps;  // Количество шагов
    cout << "Введите шаг и количество шагов : " << endl;
    cout << "h = ";
    cin >> h;
    cout << "n = ";
    cin >> n_steps;
    cout << endl;
    // Вызываем метод Рунге-Кутта
    runge_kutta_4th_order(x0, v0, h, n_steps);

    return 0;
}
