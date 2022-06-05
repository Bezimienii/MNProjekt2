import os
import pandas as pd
import matplotlib.pyplot as ptl

starting_point = 2


def lagrange_coefs(xlist: list, index: int, givenx: float):

    multres = 1.0

    for i in range(len(xlist)):
        if i != index:
            multres *= (givenx - xlist[i])/(xlist[index] - xlist[i])

    return multres


def lagrangefunc(x: float, xlist: list, ylist: list):

    result = 0.0

    for i in range(len(xlist)):
        result += ylist[i]*lagrange_coefs(xlist, i, x)

    return result


def lagrange(which: int):

    dir_list = os.listdir("./data")
    for i in dir_list:

        data = pd.read_csv("./data/"+i, sep=r"\s[0-9]|\\t|,", engine='python')
        collist = data.columns.values.tolist()
        if collist[0] == '0':
            newrow = pd.DataFrame({collist[0]: [collist[0]], collist[1]: [collist[1]]})
            data = pd.concat([newrow, data], ignore_index=True)

        xmeasured = [float(data.iat[j, 0]) for j in range(len(data))]
        ymeasured = [float(data.iat[j, 1]) for j in range(len(data))]

        yinterpolated = [float(0) for _ in range(len(data))]

        xlagrange = []
        ylagrange = []

        for j in range(starting_point, (len(data)), which):
            xlagrange.append(float(data.iat[j, 0]))
            ylagrange.append(float(data.iat[j, 1]))
            yinterpolated[j] = float(data.iat[j, 1])
            if j+which > (len(data)-1):
                xlagrange.append(float(data.iat[(len(data)-1), 0]))
                ylagrange.append(float(data.iat[(len(data)-1), 1]))
                yinterpolated[j] = float(data.iat[(len(data)-1), 1])

        jump = starting_point
        for j in range(len(yinterpolated)):
            if j != jump or j != len(yinterpolated) - 1:
                yinterpolated[j] = lagrangefunc(xmeasured[j], xlagrange, ylagrange)
            else:
                jump += which

        fname = os.path.splitext(i)[0]

        ptl.cla()
        # ptl.yscale('log')
        ptl.plot(xmeasured, ymeasured, 'ro')
        ptl.plot(xlagrange, ylagrange, 'bo')
        ptl.plot(xmeasured, yinterpolated, 'g')
        ptl.legend(['Pomierzone', 'Wybrane', 'Interpolowane'])
        ptl.xlabel("Droga przebyta [m]")
        ptl.ylabel("Wysokosc [m]")
        ptl.savefig("./imageslagrange/"+fname+".png")


def splinesforU(pointnum: int, h: int):

    alen = 4*(pointnum-1)
    A = [[0 for _ in range(alen)] for _ in range(alen)]

    # wedlug wykladu 5

    # 1
    A[0][0] = 1

    # 2
    for i in range(0, 4):
        A[1][i] = h**i

    # 7
    A[2][2] = 2

    # 6
    A[3][2] = 2
    A[3][3] = 6*h
    A[3][6] = -2

    for i in range(1, pointnum):

        # 3
        A[4*i][4*i] = 1

        # 4
        for j in range(0, 4):
            A[4*i+1][4*i+j] = h**i

        # 5
        A[4*i+2][4*i-3] = 1
        A[4*i+2][4*i+1] = -1
        A[4*i+2][4*i-2] = 2*h
        A[4*i+2][4*i-1] = 3 * (h**2)

        if i != pointnum-1:
            # 6
            A[4*i+3][4*i+2] = 2
            A[4*i+3][4*i+3] = 6*h
            A[4*i+3][4*i+6] = -2

    # 8
    A[alen-1][alen-2] = 2
    A[alen-1][alen-1] = 6*h

    return A

def LUpivoting(pointnum: int, h: int):

    A = splinesforU(pointnum, h)
    U = A
    L = [[float(0) if i != j else float(1) for j in range(len(A[0]))] for i in range(len(A))]
    P = [[float(0) if i != j else float(1) for j in range(len(A[0]))] for i in range(len(A))]

    for i in range(len(A) - 1):

        ind = 0
        maxval = 0.0

        for j in range(i, len(A)):
            for k in range(0, i+1):
                if abs(U[j][k]) > maxval and ind != j:
                    ind = j
                    maxval = abs(U[j][k])

        rowkforU = [U[i][j] for j in range(i, len(A))]
        rowindforU = [U[ind][j] for j in range(i, len(A))]

        for j in range(i, len(A)):
            U[i][j] = rowindforU[j]
            U[ind][j] = rowkforU[j]

        rowkforL = [L[i][j] for j in range(0, i)]
        rowindforL = [L[ind][j] for j in range(0, i)]

        for j in range(0, i):
            L[i][j] = rowindforL[j]
            L[ind][j] = rowkforL[j]

        rowkforP = [P[i][j] for j in range(len(A))]
        rowindforP = [P[ind][j] for j in range(len(A))]

        for j in range(len(A)):
            P[i][j] = rowindforP[j]
            P[ind][j] = rowkforP[j]

        for j in range(i + 1, len(A)):
            L[j][i] = U[j][i] / U[i][i]
            for k in range(i, len(A)):
                U[j][k] -= (U[i][k]*L[j][i])

lagrange(40)
