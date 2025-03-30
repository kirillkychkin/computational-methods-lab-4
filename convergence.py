def isDiagDominant(A):
    for i in range(len(A)):
        row = A[i]
        diagEl = 0
        # если есть A[i][i]
        if(i < len(row)):
            diagEl = abs(row[i])
        sumEls = 0
        for j in range(len(row)):
            # не диагональные элементы
            if(j != i):
                sumEls += abs(row[j])
        if(diagEl <= sumEls):
            return False
    return True
