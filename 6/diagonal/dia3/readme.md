# Точная диагонализация с использованием состояний с полуимпульсом

##  Исходный код
hchain_mkp
dsyev (процедура диагонализации LAPAC; вещественная симметричная матрица)
 
## Инструкции по запуску и примеры результатов

Точная диагонализация цепочки Гейзенберга S = 1/2
(с использованием полуимпульсных состояний)

Программа: hchan_mkp

Обратите внимание, что за исключением k = 0 и k = pi (целочисленная метка k = N / 2), симметричные (p = 1) и антисимметричные (p = -1) состояния являются вырожденными (так же, как состояния с регулярным импульсом с импульсом k и -k вырождены). Таким образом, программа выполняет расчет для p = -1 и p = + 1 только для k = 0 и pi, и только для p = + 1 для других импульсов.

Инструкции по запуску
Ввод: файл read.in, содержащий:
Столбец 1: Размер системы N (целое число)
Столбец 2: Расчетный максимальный размер блочного гильбертова пространства rm (целое число)

Пример read.in (N = 8, rm = 20)

    8 20
Вывод: Файл eig.dat, содержащий:
Для каждого сектора с фиксированным количеством вращений вверх nu и импульсом k (целое число, k = 0, ..., N / 2) и p (-1 и +1):
Строка 1: nu, k, p и размер блока nst, за которыми следуют nst строк с:
Столбец 1: Число собственных значений (целое число)
Столбец 2: Собственное значение энергии (действительное)
Столбец 3: Квантовое число спина (действительное)
Файл low.dat, содержащий наименьшее собственное состояние для каждого сектора:
Столбец 1: nu (целое число)
Столбец 2: k (целое число)
Столбец 3: p (целое число)
Столбец 4: Собственное значение энергии (действительное)
Столбец 5: Квантовое число спина (действительное)
Примеры и комментарии
Ниже приведены результаты для N = 6. Обратите внимание, что для этой небольшой системы многие секторы не имеют состояний. Файл eig.dat:
nu = 0 k = 0 p = -1 nst = 0
nu = 0 k = 0 p = 1 nst = 1
    0 0,0000000000 3,000000000000
nu = 0 k = 1 p = 1 nst = 0
nu = 0 k = 2 p = 1 nst = 0
nu = 0 k = 3 p = -1 nst = 0
nu = 0 k = 3 p = 1 nst = 0
nu = 1 k = 0 p = -1 nst = 0
nu = 1 k = 0 p = 1 nst = 1
    0 0,0000000000 3,000000000000
nu = 1 k = 1 p = 1 nst = 1
    0 0,0000000000 2,000000000000
nu = 1 k = 2 p = 1 nst = 1
    0 0,0000000000 2,000000000000
nu = 1 k = 3 p = -1 nst = 1
    0 0,0000000000 2,000000000000
nu = 1 k = 3 p = 1 nst = 0
nu = 2 k = 0 p = -1 nst = 0
nu = 2 k = 0 p = 1 nst = 3
    0 -2,1180339887 1,000000000000
    1 0.1180339887 1.0000000000
    2 0.0000000000 3.0000000000
nu = 2 k = 1 p = 1 nst = 2
    0 -1.0000000000 1.0000000000
    1 1.0000000000 2.0000000000
nu = 2 k = 2 p = 1 nst = 3
    0 -1.2807764064 1.0000000000
    1 0,0000000000 2,000000000000
    2 0,0000000000 1,000000000000
nu = 2 k = 3 p = -1 nst = 1
    0 0,0000000000 2,000000000000
nu = 2 k = 3 p = 1 nst = 1
    0 0,0000000000 1,000000000000
nu = 3 k = 0 p = -1 nst = 1
    0 0,0000000000 0,0000000000
nu = 3 k = 0 p = 1 nst = 3
    0 -2,1180339887 1,000000000000
    1 0.1180339887 1.0000000000
    2 0.0000000000 3.0000000000
nu = 3 k = 1 p = 1 nst = 3
    0 -1.0000000000 1.0000000000
    1 -0,5000000000 0,0000000000
    2 0,0000000000 2,000000000000
nu = 3 k = 2 p = 1 nst = 3
    0 -1.2807764064 1.0000000000
    1 -0.0000000000 2.0000000000
    2 0,0000000000 1,000000000000
nu = 3 k = 3 p = -1 nst = 3
    0 -2,8027756377 0,0000000000
    1 -0.5000000000 2.0000000000
    2 0,0000000000 0,0000000000
nu = 3 k = 3 p = 1 nst = 1
    0 0,0000000000 1,000000000000
Файл low.dat (состояние с наименьшей энергией в каждом секторе, содержащем состояния):

    0 0 1 0.0000000000 3.0000000000 1
    1 0 1 0.0000000000 3.0000000000 1
    1 1 1 0.0000000000 2.0000000000 1
    1 2 1 0.0000000000 2.0000000000 1
    1 3 -1 0,0000000000 2,0000000000 1
    2 0 1 -2,1180339887 1,000000000000 3
    2 1 1 -1.0000000000 1.0000000000 2
    2 2 1 -1.2807764064 1.0000000000 3
    2 3 -1 0,0000000000 2,0000000000 1
    2 3 1 0.0000000000 1.0000000000 1
    3 0-1 0,0000000000 0,0000000000 1
    3 0 1 -2,1180339887 1,000000000000 3
    3 1 1 -1.0000000000 1.0000000000 3
    3 2 1 -1.2807764064 1.0000000000 3
    3 3 -1 -2,8027756377 0,0000000000 3
    3 3 1 0.0000000000 1.0000000000 1
Интересно посмотреть на состояния с nu = 2 и k = pi для N = 8. Когда четность не используется, есть 3 вырожденных состояния с энергией E = 0 и различным полным спином. Затем диагонализация дает нецелое число S, как обсуждалось с программой hchain_mk;

 nu = 2 k = 4 nst = 4
    0 -0,0000000000 2,0009877794
    1 0,0000000000 2,9990995860
    2 0,0000000000 2,0002724281
    3 1.0000000000 2.0000000000
При использовании паритета часть этого вырождения снимается. Состояния в секторах p = -1 и p = + 1 равны

nu = 2, k = 4, p = -1, nst = 2
    0 0,0000000000 2,7015621187
    1 0,0000000000 2,3722813233
nu = 2, k = 4, p = 1, nst = 2
    0 0,0000000000 2,000000000000
    1 1.0000000000 2.0000000000
В секторе p = + 1 имеется одно состояние с E = 0, со спином S = 2. В секторе p = -1 все еще есть 2 вырожденных состояния с разными спинами, которые не разрешаются этой диагонализацией (хотя на основе полученной информации мы все еще можем в принципе вычислить спины, поскольку только одно конкретное смешивание двух целочисленных спинов может дать ожидаемые значения, предоставленные программой).