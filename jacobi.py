import numpy as np
import matplotlib.pyplot as plt
# исходные данные
A = [[5, 2, -1], [-4, 8, 3], [2, -2, 5]]
b = [12, 24, 9]

from convergence import isDiagDominant
'''
В СЛАУ вида:
a11 * x1 + a12 * x2 + a13 * x3 = b1
a21 * x1 + a22 * x2 + a23 * x3 = b2
a31 * x1 + a32 * x2 + a33 * x3 = b3
....
решения уравнений можно выразить следующим образом:
x1 = (1 / a11) * (b1 - a12 * x2 - a13 * x3)
x2 = (1 / a22) * (b2 - a21 * x1 - a23 * x3)
x3 = (1 / a33) * (b3 - a31 * x1 - a32 * x2)
....

тогда формулы для метода якоби:
x1 ^(k + 1) = (1 / a11) * (b1 - a12 * x2 ^ (k) - a13 * x3 ^ (k))
x2 ^(k + 1) = (1 / a22) * (b2 - a21 * x1 ^ (k) - a23 * x3 ^ (k))
x3 ^(k + 1) = (1 / a33) * (b3 - a31 * x1 ^ (k) - a32 * x2 ^ (k))
.....
k - метод итерации
'''

def errorCalc(x, xn):
    max_err = 0
    # по каждому решению ищем максимальное различие между соответствующими решениями x между текущим и прошлым шагом, знак роли не играет
    for i in range(len(x)):
        err = abs(x[i] - xn[i])
        if err > max_err:
            max_err = err
    # максимальное различие и есть ошибка - возвращем её
    return max_err

# Пример действия 
# вычисляет всю сумму решения, которую нужно отнять от b
# a11 * x1 + a12 * x2 + a13 * x3
# поэтому в решении в конце отнимаем диагональную часть a11 * x1
def fullBracket(A, x, i):
    res = 0
    for k in range(len(A)):
        res += A[i][k] * x[k]
    return res

# eps - точность, которую хотим получить
def jacobi(A, b, eps):
    # первичные значения, с которых начинаем вычисления
    start = [0] * len(A)

    # решения x по всем шагам
    X_res = []
    # пихаем первичные значения (чтобы начать вычисления с них)
    X_res.append(start)
    # задаем максимальную ошибку
    error = float('inf')
    # решение текущего шага
    xn = [0] * len(A)
    # номер итерации
    k = 0
    # совершаем итерации пока не достигнем желаемой точности
    while error > eps:
        # инкрементируем итерацию (нулевая итерация уже сделана, это первичные значения)
        k += 1
        # получаем значения решений с последней итерации
        # используем slice, т.к. прямое присвоение копирует
        # не значение массива, а ссылку на него
        # тогда как slice генерирует новый объект
        x = X_res[-1][:]
        # вычисляем значения текущей итерации по формуле
        for i in range(len(A)):
            xn[i] = (1 / float(A[i][i])) * (b[i] - (fullBracket(A, x, i) - A[i][i] * x[i]))
        # вставляем текущий шаг в список шагов
        X_res.append(xn[:])
        # вычисляем текущую ошибку
        error = errorCalc(x, xn)
        # если никак не можем подобраться к решению, ошибка скачет и прошло много итераций
        if error > 10000 or k == 1000:
            raise(RuntimeError('Итерации не сходятся'))
    # после выполнения цикла возвращаем список решений
    return X_res

def printPlot(solution):
    solution = np.array(solution)
    xPlots = list()
    fig, ax = plt.subplots(3)
    for i in range(len(solution[0])):
        xPlots.append(solution[:, i])
        ax[i].set_xlabel('iter')
        ax[i].set_ylabel("x"+str(i+1))
        ax[i].plot( range(0, len(xPlots[i])), xPlots[i], marker='.', linestyle='-')
    print(xPlots)
    plt.show()

if(isDiagDominant(A)):
    accuracy = 0.001
    solution = jacobi(A, b, accuracy)
    print("Решено за " + str(len(solution) - 1) + " итераций с точностью " + str(accuracy))
    print("Решение: ")
    print(solution[-1])
    printPlot(solution)
else:
    print("Метод Якоби может не сходиться, матрица не имеет строчного диагонального преобладания")