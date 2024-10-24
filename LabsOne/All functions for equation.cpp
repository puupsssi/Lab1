#include"Names of Functions.h"
double u_test(double x, double u0) {
    return u0 * exp(x);
}

double test_f(double x, double v) {
    return v;
}

double f_first_task(double x, double v) {
    return (x / (pow(x, 2) + 1)) * pow(v, 2) + v - pow(v, 3) * sin(10 * x);
}

pair<double, double> step_of_the_method_for_equation(double(*f)(double, double), double xn, double vn, double h) {
    // Вычисление коэффициентов k1, k2, k3, k4 - для метода рунге кутта 4го порядка
    double k1 = f(xn, vn);
    double k2 = f(xn + h / 2.0, vn + h / 2.0 * k1);
    double k3 = f(xn + h / 2.0, vn + h / 2.0 * k2);
    double k4 = f(xn + h, vn + h * k3);

    // Обновление значений переменных x и v
    xn = xn + h;  // Обновляем x на один шаг
    vn = vn + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);  // Обновляем v
    return { xn, vn };
}
double calculate_S(double vn_plus_1, double v2n, double epsilon) {
    //величина для контроля локальной погрешности
    double S = abs(v2n - vn_plus_1) / 15.0; //2^4 - 1 = 15;
    return S;
}
int check_the_point(double vn_plus_1, double v2n, double* h, double epsilon) {
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
vector<vector<double>> runge_kutta_4th_order(double(*f)(double, double), double x0, double v0,//передаю функцию, которая отвечает за правую часть
    double h, int n_steps, double epsilon, int need_epsilon, double right_border, double epsilon_border,
    vector<pair<double, double>>* changes_step, vector<pair<double, double>>* test_task) {
    double xn = x0; //текущая точка
    double vn = v0; //текущее значение
    double hn = h;  //текущий шаг
    double v2n = v0;//переменная для точек, вычисленных в половинным шагом
    int c1 = 0, c2 = 0, is_task_test = (f == test_f); //c1 - счетчик делений шага пополам, c2 - счетчик удвоений шага, 
                                                        //is_task_test решаем ли мы сейчас тестовую задачу 
    //создаю пары вне цикла, чтобы сто раз не выделялась память, а просто перезаписывалось значение
    pair<double, double> half_step_point;//промежуточная точка - половина шага
    pair<double, double> new_point_with_half_step;//основная точка, посчитанная с половинным шагом
    pair<double, double>new_point;//здесь будем хранить сосчитанную точку

    // Вектор для хранения результатов
    vector<vector<double>> numerical_solution;
    test_task->push_back({ u_test(xn, v0), abs(u_test(xn, v0) - vn) });
    if (need_epsilon) { // с контролем локальной погрешности
        vector<double> new_raw = { xn, vn, NAN, NAN, NAN,hn };//положили первую строку сразу
        numerical_solution.push_back(new_raw);
        for (int i = 1; xn < right_border && i < n_steps + 1; i++) {
            //добавить выход за границу
            while (1) {//пока не найдем хорошую точку
                c1 = 0, c2 = 0;
                new_point = step_of_the_method_for_equation(f, xn, vn, hn); //поcчитали новую точку с обычным шагом

                half_step_point = step_of_the_method_for_equation(f, xn, vn, hn / 2.0);//считаем точку с половинным шагом
                new_point_with_half_step = step_of_the_method_for_equation(f, half_step_point.first, half_step_point.second, hn / 2.0);
                v2n = new_point_with_half_step.second;

                int olp = check_the_point(new_point.second, new_point_with_half_step.second, &hn, epsilon); // 1 - точка хорошая, шаг тот же,
                //2 - точка хорошая, шаг в два раза больше,
                //0 - точка плохая, шаг в два раза меньше
                if (olp == 2) {//счетчик удвоения шага
                    c2 += 1;
                }
                else if (olp == 0) {//счетчик деления шага на два
                    c1 += 1;
                    changes_step->push_back({ (*changes_step)[i - 1].first + c1,(*changes_step)[i - 1].second + c2 });//добавили счетчик изменения шага для этого шага
                }

                if (olp) {//если точка хорошая
                    xn = new_point.first;
                    vn = new_point.second;//обновляем точку
                    vector<double> new_raw = { xn, vn, v2n, (vn - v2n), calculate_S(vn,v2n,epsilon),hn };//наш результат за этот шаг
                    numerical_solution.push_back(new_raw);
                    if (is_task_test) {//для тестовой задачи еще аналитическое решение
                        test_task->push_back({ u_test(xn,v0),abs(u_test(xn,v0) - vn) });
                    }
                    break;
                }
            }
            changes_step->push_back({ (*changes_step)[i - 1].first + c1,(*changes_step)[i - 1].second + c2 });//добавили счетчик изменения шага для этого шага 
            if (right_border - epsilon_border <= xn && xn <= right_border + epsilon_border)
            {
                cout << "Досчитали до границы, последняя граничная точка имеет координаты: " << xn << " , " << vn << endl << endl;
            }
            else if (xn + hn > right_border + epsilon_border)
            {
                hn = right_border + epsilon_border - xn;
                new_point = step_of_the_method_for_equation(f, xn, vn, hn);
                xn = new_point.first;
                vn = new_point.second;
                vector<double> new_raw = { xn, vn, v2n, (vn - v2n), calculate_S(vn,v2n,epsilon),hn };//наш результат за этот шаг
                numerical_solution.push_back(new_raw);
                if (is_task_test) {//для тестовой задачи еще аналитическое решение
                    test_task->push_back({ u_test(xn,v0),abs(u_test(xn,v0) - vn) });
                }
                break;
            }
            changes_step->push_back({ (*changes_step)[i - 1].first + c1,(*changes_step)[i - 1].second + c2 });//добавили счетчик изменения шага для этого шага 
        }
    }
    else { //без контроля локальной погрешности
        vector<double> new_raw = { xn, vn, hn };//положили первую строку сразу
        numerical_solution.push_back(new_raw);
        for (int i = 0; xn < right_border && i < n_steps; i++) {
            new_point = step_of_the_method_for_equation(f, xn, vn, hn);
            xn = new_point.first;
            vn = new_point.second;
            if (is_task_test) {
                new_raw = { xn, vn, hn };//наш результат за этот шаг
                numerical_solution.push_back(new_raw);
                test_task->push_back({ u_test(xn,v0),abs(u_test(xn,v0) - vn) });
            }
            else {
                vector<double> new_raw = { xn, vn, hn };//наш результат за этот шаг
                numerical_solution.push_back(new_raw);
            }
            if (right_border - epsilon_border <= xn && xn <= right_border + epsilon_border)
            {
                cout << "Досчитали до границы, последняя граничная точка имеет координаты: " << xn << " , " << vn << endl<<endl;
            }
            else if (xn + hn > right_border + epsilon_border)
            {
                hn = right_border + epsilon_border - xn;
                new_point = step_of_the_method_for_equation(f, xn, vn, hn);
                xn = new_point.first;
                vn = new_point.second;
                if (is_task_test) {
                    new_raw = { xn, vn, hn };//наш результат за этот шаг
                    numerical_solution.push_back(new_raw);
                    test_task->push_back({ u_test(xn,v0),abs(u_test(xn,v0) - vn) });
                }
                else {
                    vector<double> new_raw = { xn, vn, hn };//наш результат за этот шаг
                    numerical_solution.push_back(new_raw);
                }
                break;
            }
        }
    }
    return numerical_solution;
}