//наша первая лабораторная работа по численным методам
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <iomanip> // Для setw
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

pair<double, double> f_second_task(double u_1, double u_2, double x, double a, double b) {
    double du1_dx = u_2;
    double du2_dx = -pow(a, 2) * sin(u_1) - b * sin(x);

    return {du1_dx, du2_dx};
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

tuple<double, double,double> step_of_the_method_for_the_system(pair<double, double>(*f)(double, double, double, double, double), double xn, double v_1n, double v_2n, double h,double a,double b) {
    // Вычисление коэффициентов k1, k2, k3, k4 для u1 и u2
    double k1_u1 = f(xn, v_1n, v_2n,a,b).first;
    double k1_u2 = f(xn, v_1n, v_2n, a, b).second;

    double k2_u1 = f(xn + h / 2.0, v_1n + h / 2.0 * k1_u1, v_2n + h / 2.0 * k1_u2, a, b).first;
    double k2_u2 = f(xn + h / 2.0, v_1n + h / 2.0 * k1_u1, v_2n + h / 2.0 * k1_u2, a, b).second;

    double k3_u1 = f(xn + h / 2.0, v_1n + h / 2.0 * k2_u1, v_2n + h / 2.0 * k2_u2, a, b).first;
    double k3_u2 = f(xn + h / 2.0, v_1n + h / 2.0 * k2_u1, v_2n + h / 2.0 * k2_u2, a, b).second;

    double k4_u1 = f(xn + h, v_1n + h * k3_u1, v_2n + h * k3_u2, a, b).first;
    double k4_u2 = f(xn + h, v_1n + h * k3_u1, v_2n + h * k3_u2, a, b).second;

    // Обновление значений переменных u1 и u2
    xn = xn + h;  // Обновляем x на один шаг
    v_1n += (h / 6.0) * (k1_u1 + 2 * k2_u1 + 2 * k3_u1 + k4_u1);
    v_2n += (h / 6.0) * (k1_u2 + 2 * k2_u2 + 2 * k3_u2 + k4_u2);

    return {xn, v_1n, v_2n}; // Возвращаем обновленные значения
}

double calculate_S(double vn_plus_1, double v2n, double epsilon) {
    double S = abs(v2n - vn_plus_1) / 15.0; //2^4 - 1 = 15;
    return S;
}

double calculate_S_for_system(double v_1n_plus_1, double v_12n, double v_2n_plus_1, double v_22n, double epsilon) {
    double S_1 = abs(v_12n - v_1n_plus_1) / 15.0; //2^4 - 1 = 15;
    double S_2 = abs(v_22n - v_2n_plus_1) / 15.0;
    double S = sqrt(pow(S_1, 2) + pow(S_2, 2));
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

int check_the_point_for_system(double v_1n_plus_1, double v_12n, double v_2n_plus_1, double v_22n, double* h, double epsilon) {
    double S = calculate_S_for_system(v_1n_plus_1,v_12n,v_2n_plus_1,v_22n, epsilon);
    if ((epsilon / 15.0) <= S && S < epsilon) { //хорошая точка
        return 1; //значит, что точка хорошая и мы продолжаем счет
    }
    if (S < (epsilon / 15.0)) {
        *h = *h * 2;
        return 2; //аналогично с предыдущим случаем, только еще поменяли шаг
    }
    if (S >= epsilon) {//точка плохая, меняем шаг и начинаем той же точки (xn,v_1,v_2)
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

vector<vector<double>> runge_kutta_4th_order_for_system(pair<double, double>(*f)(double, double, double,double,double), double x0, double v_10, double v_20,double a,double b,//передаю функцию, которая отвечает за правую часть
                                                                    double h, int n_steps, double epsilon, int need_epsilon,
                                                                    vector< pair<double, double>>* changes_step) {
    double xn = x0;
    double v_1n = v_10;
    double v_2n = v_20;
    double hn = h;
    double v_12n, v_22n;
    int c1 = 0, c2 = 0, counter_for_c = 0;
    //создаю пары вне цикла, чтобы сто раз не выделялась память, а просто перезаписывалось значение
    tuple<double, double,double> half_step_point;
    tuple<double, double,double> new_point_with_half_step;
    tuple<double, double,double> new_point;

    // Вектор для хранения результатов
    vector<vector<double>> numerical_solution;

    if (need_epsilon) { // с контролем локальной погрешности
        vector<double> new_raw = { xn, v_1n,v_2n, NAN, NAN, NAN,hn };//положили первую строку сразу
        numerical_solution.push_back(new_raw);
        for (int i = 1; i < n_steps + 1; i++) {
            //добавить выход за границу
            while (1) {//пока не найдем хорошую точку
                c1 = 0, c2 = 0;
                new_point = step_of_the_method_for_the_system(f, xn, v_1n,v_2n, hn,a,b); //поcчитали новую точку с обычным шагом

                half_step_point = step_of_the_method_for_the_system(f, xn, v_1n,v_2n, hn / 2.0,a,b);//считаем точку с половинным шагом
                new_point_with_half_step = step_of_the_method_for_the_system(f, get<0>(half_step_point), get<1>(half_step_point), get<2>(half_step_point), hn / 2.0,a,b);
                v_12n = get<1>(new_point_with_half_step);
                v_22n = get<2>(new_point_with_half_step);

                int olp = check_the_point_for_system(get<1>(new_point), get<1>(new_point_with_half_step), get<2>(new_point), get<2>(new_point_with_half_step), &hn, epsilon); // 1 - точка хорошая, шаг тот же,
                //2 - точка хорошая, шаг в два раза больше,
                //3 - точка плохая, шаг в два раза меньше
                if (olp == 2) {//счетчик удвоения шага
                    c2 += 1;
                }
                else if (olp == 0) {//счетчик деления шага на два
                    c1 += 1;
                }

                if (olp) {//если точка хорошая
                    xn = get<0>(new_point);
                    v_1n = get<1>(new_point);//обновляем точку
                    v_2n = get<2>(new_point);//обновляем точку
                    vector<double> new_raw = { xn, v_1n,v_2n, v_12n,v_22n, (v_1n - v_12n),(v_2n - v_22n),calculate_S_for_system(v_1n,v_1n,v_12n,v_22n,epsilon),hn };//наш результат за этот шаг
                    numerical_solution.push_back(new_raw);
                    break;
                }
            }
            changes_step->push_back({ (*changes_step)[i - 1].first + c1,(*changes_step)[i - 1].second + c2 });//добавили счетчик изменения шага для этого шага
        }
    }
    else { //без контроля локальной погрешности
        vector<double> new_raw = { xn, v_1n,v_2n, hn };//положили первую строку сразу
        numerical_solution.push_back(new_raw);
        for (int i = 0; i < n_steps; i++) {
            //добавить выход за границу
            new_point = step_of_the_method_for_the_system(f, xn, v_1n, v_2n, hn,a,b);
            xn = get<0>(new_point);
            v_1n = get<1>(new_point);//обновляем точку
            v_2n = get<2>(new_point);//обновляем точку
            vector<double> new_raw = { xn, v_1n,v_2n, hn };//наш результат за этот шаг
            numerical_solution.push_back(new_raw);
        }
    }
    return numerical_solution;
}

int main() {
    setlocale(LC_ALL, "RUS");
    int p = 0;
    // Начальные условия
    double x0, right;  // [x(0), right]
    double v0,v_10,v_20;// v(0)
    double a, b;
    int C1 = 0, C2 = 0;
    vector<pair<double, double>> changes_of_the_step = { {C1,C2} };
    vector<pair<double, double>>for_test_task;
    // Шаг и количество шагов
    double h;  // Шаг по времени
    int n_steps;  // Количество шагов
    double epsilon = 0.1;
    int need_epsilon = 1;
    vector<vector<double>> result;
    double(*f0)(double, double);
    double(*f1)(double, double);
    pair<double, double>(*f2)(double, double, double, double, double);

    cout << "Какую задачу вы хотите решить? 0.Тестовая. 1.Первая 2.Вторая : ";
    cin >> p;
    cout << endl;
    switch (p) {
    case 0:
        cout << "Введите левую и правую границы:" << endl;
        cout << "left = ";
        cin >> x0;
        cout << "right = ";
        cin >> right;
        cout << endl;

        cout << "Введите начальные условия:" << endl;
        cout << "v0 = ";
        cin >> v0;
        cout << endl;

        cout << "Введите шаг и количество шагов: " << endl;
        cout << "h = ";
        cin >> h;
        h = abs(h);
        cout << "n = ";
        cin >> n_steps;
        n_steps = abs(n_steps);
        cout << endl;

        cout << "Нужно ли контролировать локальную погрешность? Если да - введите 1, если нет - введите 0: " << endl;
        cin >> need_epsilon;
        cout << endl;
        if (need_epsilon) {
            cout << "Введите параметр контроля погрешности: " << endl;
            cout << "e = ";
            cin >> epsilon;
            cout << endl;
            epsilon = abs(epsilon);
        }
        f0 = test_f;
        //Вызываем метод Рунге-Кутта
        result = runge_kutta_4th_order(f0, x0, v0, h, n_steps, epsilon, need_epsilon, &changes_of_the_step, &for_test_task);

        if (f0 == test_f && need_epsilon == 0) {
            cout << setw(5) << "i" << setw(10) << "xn" << setw(10) << "vn" << setw(10) << "hn" << " " << setw(15) << "un" << setw(15) << "un - vn" << endl;
            for (int i = 0; i < n_steps; i++) {
                cout << setw(5) << i;
                for (int j = 0; j < result[i].size(); ++j) {
                    cout << setw(10) << result[i][j];
                }
                cout << " " << setw(15) << for_test_task[i].first << setw(15) << for_test_task[i].second;
                cout << endl;
            }
        }
        else {
            cout << setw(5) << "i" << setw(15) << "xn" << setw(15) << "vn" << setw(15) << "v2n" << setw(15) << "vn - v2n" << setw(15) << "hn" << setw(15) << "S" << setw(5) << "c1" << setw(5) << "c2" << setw(15) << "un" << " " << setw(10) << "un - vn" << endl;
            for (int i = 0; i < n_steps; i++) {
                cout << setw(5) << i;
                for (int j = 0; j < result[i].size(); ++j) {
                    cout << setw(15) << result[i][j];
                }
                cout << setw(5) << changes_of_the_step[i].first << setw(5) << changes_of_the_step[i].second;
                cout << setw(15) << for_test_task[i].first << " " << setw(10) << for_test_task[i].second;
                cout << endl;
            }
        }
        break; // завершение выполнения блока case
    case 1:
        cout << "Введите левую и правую границы:" << endl;
        cout << "left = ";
        cin >> x0;
        cout << "right = ";
        cin >> right;
        cout << endl;

        cout << "Введите начальные условия:" << endl;
        cout << "v0 = ";
        cin >> v0;
        cout << endl;

        cout << "Введите шаг и количество шагов: " << endl;
        cout << "h = ";
        cin >> h;
        h = abs(h);
        cout << "n = ";
        cin >> n_steps;
        n_steps = abs(n_steps);
        cout << endl;

        cout << "Нужно ли контролировать локальную погрешность? Если да - введите 1, если нет - введите 0: " << endl;
        cin >> need_epsilon;
        cout << endl;
        if (need_epsilon) {
            cout << "Введите параметр контроля погрешности: " << endl;
            cout << "e = ";
            cin >> epsilon;
            cout << endl;
            epsilon = abs(epsilon);
        }
        f1 = f_first_task;
        //Вызываем метод Рунге-Кутта
        result = runge_kutta_4th_order(f1, x0, v0, h, n_steps, epsilon, need_epsilon, &changes_of_the_step, &for_test_task);
        if (f1 ==f_first_task && need_epsilon != 0) {
            cout << setw(5) << "i" << setw(15) << "xn" << setw(15) << "vn" << setw(15) << "v2n" << setw(15) << "vn - v2n" << setw(15) << "hn" << setw(15) << "S" << setw(7) << "c1" << setw(7) << "c2" << endl;
            for (int i = 0; i < n_steps; i++) {
                cout << setw(5) << i;
                for (int j = 0; j < result[i].size(); ++j) {
                    cout << setw(15) << result[i][j];
                }
                cout << setw(7) << changes_of_the_step[i].first << setw(7) << changes_of_the_step[i].second;
                cout << endl;
            }
        }
        else{
            cout << setw(5) << "i" << setw(15) << "xn" << setw(15) << "vn" << setw(15) << "hn" << endl;
            for (int i = 0; i < n_steps; i++) {
                cout << setw(5) << i;
                for (int j = 0; j < result[i].size(); ++j) {
                    cout << setw(15) << result[i][j];
                }
                cout << endl;
            }
        }
        break;
    case 2:
        cout << "Введите левую и правую границы:" << endl;
        cout << "left = ";
        cin >> x0;
        cout << "right = ";
        cin >> right;
        cout << endl;

        cout << "Введите начальные условия:" << endl;
        cout << "v_10 = ";
        cin >> v_10;
        cout << "v_20 = ";
        cin >> v_20;
        cout << "a = ";
        cin >> a;
        cout << "b = ";
        cin >> b;
        cout << endl;

        cout << "Введите шаг и количество шагов: " << endl;
        cout << "h = ";
        cin >> h;
        h = abs(h);
        cout << "n = ";
        cin >> n_steps;
        n_steps = abs(n_steps);
        cout << endl;

        cout << "Нужно ли контролировать локальную погрешность? Если да - введите 1, если нет - введите 0: " << endl;
        cin >> need_epsilon;
        cout << endl;
        if (need_epsilon) {
            cout << "Введите параметр контроля погрешности: " << endl;
            cout << "e = ";
            cin >> epsilon;
            cout << endl;
            epsilon = abs(epsilon);
        }
        f2 = f_second_task;
        //Вызываем метод Рунге-Кутта
        result = runge_kutta_4th_order_for_system(f2, x0, v_10,v_20,a,b, h, n_steps, epsilon, need_epsilon, &changes_of_the_step);
        if (f2 == f_second_task && need_epsilon != 0) {
            cout << setw(5) << "i" << setw(14) << "xn" << setw(14) << "v_1n" << setw(14)<< "v_2n" << setw(14) << "v_12n" << setw(14) << "v_22n" << setw(14) << "v_1n - v_12n" << setw(14) << "v_2n - v_22n" << setw(14) << "hn" << setw(14) << "S" << setw(7) << "c1" << setw(7) << "c2" << endl;
            for (int i = 0; i < n_steps; i++) {
                cout << setw(5) << i;
                for (int j = 0; j < result[i].size(); ++j) {
                    cout << setw(14) << result[i][j];
                }
                cout << setw(7) << changes_of_the_step[i].first << setw(7) << changes_of_the_step[i].second;
                cout << endl;
            }
        }
        else {
            cout << setw(5) << "i" << setw(15) << "xn" << setw(15) << "v_1n" << setw(15)<< "v_2n" << setw(15) << "hn" << endl;
            for (int i = 0; i < n_steps; i++) {
                cout << setw(5) << i;
                for (int j = 0; j < result[i].size(); ++j) {
                    cout << setw(15) << result[i][j];
                }
                cout << endl;
            }
        }
        break;
    default:
        cout<<"Вы не правильно ввели номер задачи!";
        break;
    }

    //вывод
    return 0;
}
