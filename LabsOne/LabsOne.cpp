//наша первая лабораторная работа по численным методам
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
using namespace std;
// Функция f(x, v) для системы уравнений

double u_test(double x, double u0) {
    return u0 * exp(x);
}

double test_f(double x, double v) {
    return v;
}

double f_first_task(double x, double v) {
    return (x/(pow(x,2)+1))*pow(v,2)+v-pow(v,3)*sin(10*x); 
}

pair<double, double> step_of_the_method(double(*f)(double, double),double xn, double vn, double h) {
    // Вычисление коэффициентов k1, k2, k3, k4
    double k1 = f(xn, vn);
    double k2 = f(xn + h / 2.0, vn + h / 2.0 * k1);
    double k3 = f(xn + h / 2.0, vn + h / 2.0 * k2);
    double k4 = f(xn + h, vn + h * k3);

    // Обновление значений переменных x и v
    xn = xn + h;  // Обновляем x на один шаг
    vn = vn + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);  // Обновляем v
    return {xn, vn};
}

double calculate_S(double vn_plus_1, double v2n, double epsilon) {
    double S = abs(v2n - vn_plus_1) / 15.0; //2^4 - 1 = 15;
    return S;
}

int check_the_point(double vn_plus_1, double v2n, double *h, double epsilon) {
    double S = calculate_S(vn_plus_1, v2n, epsilon);
    if ((epsilon / 15.0) <= S && S < epsilon) { //хорошая точка
        return 1; //значит, что точка хорошая и мы продолжаем счет
    }
    if (S < (epsilon / 15.0)) {
        *h = *h * 2;
        return 2; //аналогично с предыдущим случаем, только еще поменяли шаг
    }
    if (S >= epsilon) {//точка плохая, меняем шаг и начинаем той же точки (xn,vn)
        *h = *h / 2.0;
        return 0;
    }
    return 0;
}

// Функция для метода Рунге-Кутта 4-го порядка
vector<vector<double>> runge_kutta_4th_order(double(*f)(double, double),double x0, double v0,//передаю функцию, которая отвечает за правую часть
                                                        double h, int n_steps, double epsilon, int need_epsilon,
                                                        vector<pair<double,double>> *changes_step, vector<pair<double, double>> *test_task) {
    double xn = x0;
    double vn = v0;
    double hn = h;
    double v2n;
    int c1=0, c2=0, is_task_test= (f == test_f), counter_for_c=0;
    //создаю пары вне цикла, чтобы сто раз не выделялась память, а просто перезаписывалось значение
    pair<double, double> half_step_point;
    pair<double, double> new_point_with_half_step;
    pair<double, double>new_point;

    // Вектор для хранения результатов
    vector<vector<double>> numerical_solution;

    if (need_epsilon) { // с контролем локальной погрешности
        vector<double> new_raw = { xn, vn, NAN, NAN, NAN,hn };//положили первую строку сразу
        numerical_solution.push_back(new_raw);
        for (int i = 1; i < n_steps+1; i++) {
            //добавить выход за границу
            while (1) {//пока не найдем хорошую точку
                c1 = 0, c2 = 0;
                new_point = step_of_the_method(f, xn, vn, hn); //поcчитали новую точку с обычным шагом

                half_step_point = step_of_the_method(f,xn, vn, hn / 2.0);//считаем точку с половинным шагом
                new_point_with_half_step = step_of_the_method(f,half_step_point.first, half_step_point.second, hn / 2.0);
                v2n = new_point_with_half_step.second;

                int olp = check_the_point(new_point.second, new_point_with_half_step.second, &hn, epsilon); // 1 - точка хорошая, шаг тот же,
                                                                                                            //2 - точка хорошая, шаг в два раза больше,
                                                                                                            //3 - точка плохая, шаг в два раза меньше
                if (olp == 2) {//счетчик удвоения шага
                    c2 += 1;
                }
                else if (olp == 0) {//счетчик деления шага на два
                    c1 += 1;
                }

                if (olp) {//если точка хорошая
                    xn = new_point.first;
                    vn = new_point.second;//обновляем точку
                    vector<double> new_raw = {xn, vn, v2n, (vn-v2n), calculate_S(vn,v2n,epsilon),hn};//наш результат за этот шаг
                    numerical_solution.push_back(new_raw);
                    if (is_task_test) {//для тестовой задачи еще аналитическое решение
                        test_task->push_back({ u_test(xn,v0),abs(u_test(xn,v0) - vn) });
                    }
                    break;
                }
            }
            changes_step->push_back({ (*changes_step)[i-1].first + c1,(*changes_step)[i-1].second + c2});//добавили счетчик изменения шага для этого шага
        }
    }
    else { //без контроля локальной погрешности
        vector<double> new_raw = { xn, vn, hn };//положили первую строку сразу
        numerical_solution.push_back(new_raw);
        for (int i = 0; i < n_steps; i++) {
            //добавить выход за границу
            new_point = step_of_the_method(f, xn, vn, hn);
            xn = new_point.first;
            vn = new_point.second;
            if (is_task_test) {
            vector<double> new_raw = { xn, vn, hn};//наш результат за этот шаг
            numerical_solution.push_back(new_raw);
            test_task->push_back({ u_test(xn,v0),abs(u_test(xn,v0) - vn) });
            }
            else {
                vector<double> new_raw = { xn, vn, hn };//наш результат за этот шаг
                numerical_solution.push_back(new_raw);
            }
        }
    }
    return numerical_solution;
}

int main() {
    setlocale(LC_ALL, "RUS");
    // Начальные условия
    double x0, b;  // [x(0), b]
    double v0;// v(0)
    int C1=0, C2=0;
    vector<pair<double, double>> changes_of_the_step = { {C1,C2} };
    vector<pair<double, double>>for_test_task;
    cout << "Введите левую и правую границы:" << endl;
    cout << "a = ";
    cin >> x0;
    cout << "b = ";
    cin >> b;
    cout << endl;

    cout << "Введите начальные условия:" << endl;
    cout << "v0 = ";
    cin >> v0;
    cout << endl;


    // Шаг и количество шагов
    double h;  // Шаг по времени
    int n_steps;  // Количество шагов
    cout << "Введите шаг и количество шагов: " << endl;
    cout << "h = ";
    cin >> h;
    h = abs(h);
    cout << "n = ";
    cin >> n_steps;
    n_steps = abs(n_steps);
    cout << endl;

    double epsilon=0.1;
    cout << "Нужно ли контролировать локальную погрешность? Если да - введите 1, если нет - введите 0: " << endl;
    int need_epsilon = 1;
    cin >> need_epsilon;
    cout << endl;
    if (need_epsilon) {
        cout << "Введите параметр контроля погрешности: "<<endl;
        cout << "e = ";
        cin >> epsilon;
        cout << endl;
        epsilon = abs(epsilon);
    }

    double(*f)(double, double)  = f_first_task;
    // Вызываем метод Рунге-Кутта
    vector<vector<double>> result = runge_kutta_4th_order(f,x0, v0, h, n_steps, epsilon, need_epsilon, &changes_of_the_step, &for_test_task);

    //вывод

    if (f == test_f && need_epsilon==0) {
        cout << "i   xn   vn   hn   un   un-vn" << endl;
        for (int i = 0; i < n_steps; i++) {
            cout << i << "   ";
            for (int j=0; j < result[i].size(); ++j) {
                cout << result[i][j] << "   ";
            }
            cout << for_test_task[i].first << "   " << for_test_task[i].second;
            cout << endl;
        }
    }
    else if (f == test_f && need_epsilon != 0) {
        cout << "i   xn   vn   v2n   vn-v2n   hn   S   c1   c2  un   un-vn" << endl;
        for (int i = 0; i < n_steps; i++) {
            cout << i << "   ";
            for (int j = 0; j < result[i].size(); ++j) {
                cout << result[i][j] << "   ";
            }
            cout << changes_of_the_step[i].first << "   " << changes_of_the_step[i].second<<"   ";
            cout << for_test_task[i].first << "   " << for_test_task[i].second;
            cout << endl;
        }
    }
    else if (f != test_f && need_epsilon != 0) {
        cout << "i   xn   vn   v2n   vn-v2n   hn   S   c1   c2" << endl;
        for (int i = 0; i < n_steps; i++) {
            cout << i << "   ";
            for (int j = 0; j < result[i].size(); ++j) {
                cout << result[i][j] << "   ";
            }
            cout << changes_of_the_step[i].first << "   " << changes_of_the_step[i].second;
            cout << endl;
        }
    }
    else if (f != test_f && need_epsilon == 0) {
        cout << "i   xn   vn   hn" << endl;
        for (int i = 0; i < n_steps; i++) {
            cout << i << "   ";
            for (int j = 0; j < result[i].size(); ++j) {
                cout << result[i][j] << "   ";
            }
            cout << endl;
        }
    }
    return 0;
}
