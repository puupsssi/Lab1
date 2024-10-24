//наша первая лабораторная работа по численным методам
#pragma once
#include"Names of Functions.h"

int main() {
    setlocale(LC_ALL, "RUS");
    int p = 0;
    // Начальные условия
    double x0, right_border;  // [x(0), right_border]
    double v0,v_10,v_20;// v(0)
    double A=0, B=0;//параметры для второй задачи
    int C1 = 0, C2 = 0;
    vector<pair<double, double>> changes_of_the_step = { {C1,C2} };//изменения шага
    vector<pair<double, double>>for_test_task;//для тестовой задачи: пары(точное значение функции u(x) тестовой задачи, разница значений точного и численного решений в данной точке)
    // Шаг и количество шагов
    double h;  // Шаг по времени
    int n_steps;  // Количество шагов
    double epsilon = 0.01;//параметр контроля локальной погрешности
    int need_epsilon = 1;//нужно ли контролировать локальную погрешность
    double epsilon_border = 0.01;//в какой окрестности от границы мы должны завершить счет
    vector<vector<double>> result;
    double(*f0)(double, double);
    double(*f1)(double, double);
    pair<double, double>(*f2)(double, double, double, double, double);

    cout << "Какую задачу вы хотите решить?\n0.Тестовая\n1.Первая\n2.Вторая\nВведите номер: ";
    cin >> p;
    cout << endl;
    switch (p) {
    case 0:
        cout << "Введите левую и правую границы:" << endl;
        cout << "left = ";
        cin >> x0;
        cout << "right_boarder = ";
        cin >> right_border;
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
        result = runge_kutta_4th_order(f0, x0, v0, h, n_steps, epsilon, need_epsilon,right_border,epsilon_border, &changes_of_the_step, &for_test_task);
        n_steps = size(result);
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
            cout << setw(5) << "i" << setw(15) << "xn" << setw(15) << "vn" << setw(15) << "v2n" << setw(15) << "vn - v2n" << setw(15) << "S" << setw(15) << "hn" << setw(5) << "c1" << setw(5) << "c2" << setw(15) << "un" << " " << setw(10) << "un - vn" << endl;
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
        cout << "right_border = ";
        cin >> right_border;
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
        result = runge_kutta_4th_order(f1, x0, v0, h, n_steps, epsilon, need_epsilon, right_border, epsilon_border, &changes_of_the_step, &for_test_task);
        n_steps = size(result);
        if (f1 ==f_first_task && need_epsilon != 0) {
            cout << setw(5) << "i" << setw(15) << "xn" << setw(15) << "vn" << setw(15) << "v2n" << setw(15) << "vn - v2n" << setw(15) << "S" << setw(15) << "hn" << setw(7) << "c1" << setw(7) << "c2" << endl;
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
        cout << "right_boarder = ";
        cin >> right_border;
        cout << endl;

        cout << "Введите начальные условия:" << endl;
        cout << "v_10 = ";
        cin >> v_10;
        cout << "v_20 = ";
        cin >> v_20;
        cout << "a = ";
        cin >> A;
        cout << "b = ";
        cin >> B;
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
        result = runge_kutta_4th_order_for_system(f2, x0, v_10,v_20,A,B, h, n_steps, epsilon, need_epsilon, right_border, epsilon_border, &changes_of_the_step);
        n_steps = size(result);
        if (f2 == f_second_task && need_epsilon != 0) {
            cout << setw(5) << "i" << setw(14) << "xn" << setw(14) << "v_1n" << setw(14)<< "v_2n" << setw(14) << "v_12n" << setw(14) << "v_22n" << setw(14) << "v_1n - v_12n" << setw(14) << "v_2n - v_22n" << setw(14) << "S" << setw(14) << "hn" << setw(7) << "c1" << setw(7) << "c2" << endl;
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
