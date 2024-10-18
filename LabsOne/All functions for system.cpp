#include"Names of Functions.h"
// Функция f(x, v) для системы уравнений


pair<double, double> f_second_task(double u_1, double u_2, double x, double a, double b) {
    double du1_dx = u_2;
    double du2_dx = -pow(a, 2) * sin(u_1) - b * sin(x);

    return { du1_dx, du2_dx };
}


tuple<double, double, double> step_of_the_method_for_the_system(pair<double, double>(*f)(double, double, double, double, double), double xn, double v_1n, double v_2n, double h, double a, double b) {
    // Вычисление коэффициентов k1, k2, k3, k4 для u1 и u2 - для метода рунге кутта 4го порядка
    double k1_u1 = f(xn, v_1n, v_2n, a, b).first;
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

    return { xn, v_1n, v_2n }; // Возвращаем обновленные значения
}


double calculate_S_for_system(double v_1n_plus_1, double v_12n, double v_2n_plus_1, double v_22n, double epsilon) {
    double S_1 = abs(v_12n - v_1n_plus_1) / 15.0; //2^4 - 1 = 15;
    double S_2 = abs(v_22n - v_2n_plus_1) / 15.0;
    double S = sqrt(pow(S_1, 2) + pow(S_2, 2));//усредненные величины для контроля локальной погрешности
    return S;
}


int check_the_point_for_system(double v_1n_plus_1, double v_12n, double v_2n_plus_1, double v_22n, double* h, double epsilon) {
    double S = calculate_S_for_system(v_1n_plus_1, v_12n, v_2n_plus_1, v_22n, epsilon);
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


vector<vector<double>> runge_kutta_4th_order_for_system(pair<double, double>(*f)(double, double, double, double, double), double x0, double v_10, double v_20, double a, double b,//передаю функцию, которая отвечает за правую часть
    double h, int n_steps, double epsilon, int need_epsilon, double right_boarder, double epsilon_boarder,
    vector< pair<double, double>>* changes_step) {
    double xn = x0;//текущая точка
    double v_1n = v_10;//текущее значение первой компоненты
    double v_2n = v_20;//текущее значение второй компоненты
    double hn = h; //текущий шаг
    double v_12n, v_22n;//переменные для точек, вычисленных в половинным шагом
    int c1 = 0, c2 = 0; // c1 - счетчик делений шага пополам, c2 - счетчик удвоений шага
    //создаю точки вне цикла, чтобы сто раз не выделялась память, а просто перезаписывалось значение
    tuple<double, double, double> half_step_point;//промежуточная точка - половина шага
    tuple<double, double, double> new_point_with_half_step;//основная точка, посчитанная с половинным шагом
    tuple<double, double, double> new_point;//здесь будем хранить сосчитанную точку

    // Вектор для хранения результатов
    vector<vector<double>> numerical_solution;

    if (need_epsilon) { // с контролем локальной погрешности
        vector<double> new_raw = { xn, v_1n,v_2n, NAN, NAN, NAN,hn };//положили первую строку сразу
        numerical_solution.push_back(new_raw);
        for (int i = 1; i < n_steps + 1; i++) {
            //добавить выход за границу
            while (1) {//пока не найдем хорошую точку
                c1 = 0, c2 = 0;
                new_point = step_of_the_method_for_the_system(f, xn, v_1n, v_2n, hn, a, b); //поcчитали новую точку с обычным шагом

                half_step_point = step_of_the_method_for_the_system(f, xn, v_1n, v_2n, hn / 2.0, a, b);//считаем точку с половинным шагом
                new_point_with_half_step = step_of_the_method_for_the_system(f, get<0>(half_step_point), get<1>(half_step_point), get<2>(half_step_point), hn / 2.0, a, b);
                v_12n = get<1>(new_point_with_half_step);
                v_22n = get<2>(new_point_with_half_step);

                int olp = check_the_point_for_system(get<1>(new_point), get<1>(new_point_with_half_step),
                                                        get<2>(new_point), get<2>(new_point_with_half_step), &hn, epsilon); 
                                                                        // 1 - точка хорошая, шаг тот же,
                                                                        //2 - точка хорошая, шаг в два раза больше,
                                                                        //0 - точка плохая, шаг в два раза меньше
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
            new_point = step_of_the_method_for_the_system(f, xn, v_1n, v_2n, hn, a, b);
            xn = get<0>(new_point);
            v_1n = get<1>(new_point);//обновляем точку
            v_2n = get<2>(new_point);//обновляем точку
            vector<double> new_raw = { xn, v_1n,v_2n, hn };//наш результат за этот шаг
            numerical_solution.push_back(new_raw);
        }
    }
    return numerical_solution;
}