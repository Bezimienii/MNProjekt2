import os
import random

import pandas as pd
import matplotlib.pyplot as ptl

starting_point = 0
forrandom = 60

def lagrange_coefs(xlist: list, index: int, givenx: float):
    multres = 1.0

    for i in range(len(xlist)):
        if i != index:
            multres *= (givenx - xlist[i]) / (xlist[index] - xlist[i])

    return multres


def lagrangefunc(x: float, xlist: list, ylist: list):
    result = 0.0

    for i in range(len(xlist)):
        result += ylist[i] * lagrange_coefs(xlist, i, x)

    return result

# for random
# uneven = [0]
#         for j in range(0, 512, forrandom):
#             if j < 512-forrandom:
#                 ifAdded = True
#                 while ifAdded:
#                     valnew = random.randint(j, j+forrandom-1)
#                     if valnew not in uneven:
#                         uneven.append(valnew)
#                         ifAdded = False
#
#         xmeasured = [float(data.iat[j, 0]) for j in range(len(data))]
#         ymeasured = [float(data.iat[j, 1]) for j in range(len(data))]
#
#         xinterpolated = [float(data.iat[j, 0]) for j in range(0, uneven[len(uneven)-1]+1)]
#         yinterpolated = [float(0) for _ in range(0, uneven[len(uneven)-1]+1)]
#
#         xlagrange = []
#         ylagrange = []
#
#         for j in uneven:
#             xlagrange.append(float(data.iat[j, 0]))
#             ylagrange.append(float(data.iat[j, 1]))
#             yinterpolated[j] = float(data.iat[j, 1])
#
#         jump = uneven[0]
#         for j in range(len(yinterpolated)):
#             if j != jump or j != len(yinterpolated) - 1:
#                 yinterpolated[j] = lagrangefunc(xmeasured[j], xlagrange, ylagrange)
#             else:
#                 jump = uneven.index(jump)
#                 jump = uneven[jump+1]
#                 jump = ymeasured.index(jump)

def lagrange(which: int):
    dir_list = os.listdir("./data")
    for i in dir_list:

        data = pd.read_csv("./data/" + i, sep=r"\s[0-9]|\\t|,", engine='python')
        collist = data.columns.values.tolist()
        if collist[0] == '0':
            newrow = pd.DataFrame({collist[0]: [collist[0]], collist[1]: [collist[1]]})
            data = pd.concat([newrow, data], ignore_index=True)

        uneven = [0]
        for j in range(0, 512, forrandom):
            if j < 512-forrandom:
                ifAdded = True
                while ifAdded:
                    valnew = random.randint(j, j+forrandom-1)
                    if valnew not in uneven:
                        uneven.append(valnew)
                        ifAdded = False

        xmeasured = [float(data.iat[j, 0]) for j in range(len(data))]
        ymeasured = [float(data.iat[j, 1]) for j in range(len(data))]

        xinterpolated = [float(data.iat[j, 0]) for j in range(0, uneven[len(uneven)-1]+1)]
        yinterpolated = [float(0) for _ in range(0, uneven[len(uneven)-1]+1)]

        xlagrange = []
        ylagrange = []

        for j in uneven:
            xlagrange.append(float(data.iat[j, 0]))
            ylagrange.append(float(data.iat[j, 1]))
            yinterpolated[j] = float(data.iat[j, 1])

        jump = uneven[0]
        for j in range(len(yinterpolated)):
            if j != jump or j != len(yinterpolated) - 1:
                yinterpolated[j] = lagrangefunc(xmeasured[j], xlagrange, ylagrange)
            else:
                jump = uneven.index(jump)
                jump = uneven[jump+1]
                jump = ymeasured.index(jump)

        fname = os.path.splitext(i)[0]

        ptl.cla()
        # ptl.yscale('log')
        ptl.plot(xmeasured, ymeasured, 'ro')
        ptl.plot(xlagrange, ylagrange, 'bo')
        ptl.plot(xinterpolated, yinterpolated, 'g')
        ptl.legend(['Pomierzone', 'Wybrane', 'Interpolowane'])
        ptl.xlabel("Droga przebyta [m]")
        ptl.ylabel("Wysokosc [m]")
        ptl.savefig("./imageslagrange/" + fname + ".png")


def splinesforU(pointnum: int, xes: list):
    alen = 4*(pointnum-1)
    A = [[float(0) for _ in range(alen)] for _ in range(alen)]
    A[0][0] = 1.0

    # 2
    for i in range(0, 4):
        A[1][i] = float((xes[1]-xes[0]) ** i)

    # 7
    A[2][2] = 2.0

    # 6
    A[3][2] = 2.0
    A[3][3] = 6.0 * (xes[1]-xes[0])
    A[3][6] = -2.0

    for i in range(1, pointnum - 1):

        # 3
        A[4 * i][4 * i] = 1.0

        # 5
        A[4 * i + 1][4 * i - 3] = 1.0
        A[4 * i + 1][4 * i + 1] = -1.0
        A[4 * i + 1][4 * i - 2] = float(2 * (xes[i]-xes[i-1]))
        A[4 * i + 1][4 * i - 1] = float(3 * ((xes[i]-xes[i-1]) ** 2))

        # 4
        for j in range(0, 4):
            A[4 * i + 2][4 * i + j] = float((xes[i+1]-xes[i]) ** j)

        if i != (pointnum - 2):
            # 6
            A[4 * i + 3][4 * i + 2] = 2.0
            A[4 * i + 3][4 * i + 3] = float(6 * (xes[i+1]-xes[i]))
            A[4 * i + 3][4 * i + 6] = -2.0

    # 8
    A[alen - 1][alen - 2] = 2.0
    A[alen - 1][alen - 1] = float(6 * (xes[len(xes)-1]-xes[len(xes)-2]))

    return A

def LU(pointnum: int, xes: list):
    A = splinesforU(pointnum, xes)
    U = A
    L = [[0.0 if i != j else 1.0 for j in range(len(A[0]))] for i in range(len(A))]

    for i in range(len(A) - 1):

        for j in range(i + 1, len(A)):
            L[j][i] = U[j][i] / U[i][i]
            for k in range(i, len(A)):
                U[j][k] -= (U[i][k] * L[j][i])

    return L, U


def valofsplinefunc(leftend: float, coefs: list, point: float):
    valofdiff = point - leftend

    res = 0.0
    for i in range(len(coefs)):
        res += coefs[i] * (valofdiff ** i)

    return res

# for random splain points
# uneven = [0]
#         for j in range(0, 512, forrandom):
#             if j < 512-forrandom:
#                 ifAdded = True
#                 while ifAdded:
#                     valnew = random.randint(j, j+forrandom-1)
#                     if valnew not in uneven:
#                         uneven.append(valnew)
#                         ifAdded = False
#
#
#         xsplines = [float(data.iat[j, 0]) for j in uneven]
#         ysplines = [float(data.iat[j, 1]) for j in uneven]
#
#         lastindex = uneven[len(uneven)-1]
#
#         xinterpolated = [float(data.iat[j, 0]) for j in range(0, lastindex + 1)]
#         yinterpolated = [0.0 for _ in range(0, lastindex + 1)]
#
#         for j in uneven:
#             yinterpolated[j] = float(data.iat[j, 1])
#
#         points = [j for j in uneven]
#         L, U = LU(len(points), xsplines)
#
#         b = []
#         for j in range(4*(len(points)-1)):
#             b.append(0)
#
#         b[0] = ysplines[0]
#         b[1] = ysplines[1]
#         for j in range(1, len(ysplines) - 1):
#             b[4 * j] = ysplines[j]
#             b[4 * j + 2] = ysplines[j+1]
#
#         y = [j for j in b]
#
#         for j in range(len(L)):
#             for k in range(j):
#                 y[j] -= L[j][k] * y[k]
#             y[j] /= L[j][j]
#
#         x = [j for j in y]
#
#         for j in range(len(U) - 1, -1, -1):
#             for k in range(j + 1, len(U)):
#                 x[j] -= U[j][k] * x[k]
#             x[j] /= U[j][j]
#
#         indexofleft = uneven[0]
#         indexofright = uneven[1]
#         indexofcoefs = 0
#         coefsforfunc = [0.0, 0.0, 0.0, 0.0]
#         for j in range(4):
#             coefsforfunc[j] = x[4 * indexofcoefs + j]
#         anew = True
#
#         for j in range(uneven[0] + 1, lastindex):
#             if j == indexofright:
#                 indexofleft = uneven[uneven.index(indexofleft)+1]
#                 indexofright = uneven[uneven.index(indexofright)+1]
#                 indexofcoefs += 1
#                 anew = False
#             else:
#                 if not anew:
#                     for k in range(4):
#                         coefsforfunc[k] = x[4 * indexofcoefs + k]
#                 yinterpolated[j - starting_point] = valofsplinefunc(xmeasured[indexofleft], coefsforfunc,
#                                                                           xmeasured[j])

def splines(gap: int):
    dir_list = os.listdir("./data")
    for i in dir_list:

        data = pd.read_csv("./data/" + i, sep=r"\s[0-9]|\\t|,", engine='python')
        collist = data.columns.values.tolist()
        if collist[0] == '0':
            newrow = pd.DataFrame({collist[0]: [collist[0]], collist[1]: [collist[1]]})
            data = pd.concat([newrow, data], ignore_index=True)


        xmeasured = [float(data.iat[j, 0]) for j in range(len(data))]
        ymeasured = [float(data.iat[j, 1]) for j in range(len(data))]

        xsplines = [float(data.iat[j, 0]) for j in range(starting_point, 512, gap)]
        ysplines = [float(data.iat[j, 1]) for j in range(starting_point, 512, gap)]

        lastindex = starting_point
        while lastindex + gap <= 511:
            lastindex += gap

        xinterpolated = [float(data.iat[j, 0]) for j in range(starting_point, lastindex + 1)]
        yinterpolated = [0.0 for _ in range(starting_point, lastindex + 1)]

        points = [i for i in range(starting_point, 512, gap)]
        L, U = LU(len(points), xsplines)

        for j in range(starting_point, 512, gap):
            yinterpolated[j - starting_point] = ymeasured[j]

        b = []
        for j in range(4*(len(points)-1)):
            b.append(0)

        b[0] = ysplines[0]
        b[1] = ysplines[1]
        for j in range(1, len(ysplines) - 1):
            b[4 * j] = ysplines[j]
            b[4 * j + 2] = ysplines[j+1]

        y = [j for j in b]

        for j in range(len(L)):
            for k in range(j):
                y[j] -= L[j][k] * y[k]
            y[j] /= L[j][j]

        x = [j for j in y]

        for j in range(len(U) - 1, -1, -1):
            for k in range(j + 1, len(U)):
                x[j] -= U[j][k] * x[k]
            x[j] /= U[j][j]

        for j in x:
            print(j)

        indexofleft = starting_point
        indexofright = starting_point + gap
        indexofcoefs = 0
        coefsforfunc = [0.0, 0.0, 0.0, 0.0]
        for j in range(4):
            coefsforfunc[j] = x[4 * indexofcoefs + j]
        anew = True
        print()
        for j in range(starting_point + 1, lastindex):
            if j == indexofright:
                indexofleft += gap
                indexofright += gap
                indexofcoefs += 1
                anew = False
            else:
                if not anew:
                    for k in range(4):
                        coefsforfunc[k] = x[4 * indexofcoefs + k]
                yinterpolated[j - starting_point] = valofsplinefunc(xmeasured[indexofleft], coefsforfunc,
                                                                          xmeasured[j])

        fname = os.path.splitext(i)[0]

        ptl.cla()
        #ptl.yscale('log')
        ptl.plot(xmeasured, ymeasured, 'ro')
        ptl.plot(xsplines, ysplines, 'bo')
        ptl.plot(xinterpolated, yinterpolated, 'g')
        ptl.legend(['Pomierzone', 'Wybrane', 'Interpolowane'])
        ptl.xlabel("Droga przebyta [m]")
        ptl.ylabel("Wysokosc [m]")
        ptl.savefig("./imagessplains/" + fname + ".png")

lagrange(60)
#splines(14)

