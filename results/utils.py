def xpos(i, j, nnodes):
    if i == j:
        print("i == j in xpos")
        return -1

    if i > j:
        return xpos(j, i, nnodes);

    # print(i * nnodes + j - ((i + 1) * (i + 2)) / 2)
    return i * nnodes + j - ((i + 1) * (i + 2)) / 2


def xxpos(i, j, nnodes):
    # print(i * nnodes + j)
    return i * nnodes + j

def upos(i, nnodes):
    return (nnodes * nnodes) + i - 1

def ypos(i, j, nnodes):

    if i == j: 
        print("variable y does not exist for same i and j")
        return -1

    n = nnodes * nnodes + i * (nnodes - 1) + j
    if i < j:
        n -= 1

    # print(n)
    return n

def rev(nnodes, fun, target):
    if fun == upos:
        for i in range(nnodes):
            if fun(i, nnodes) == target:
                print(i+1)
                return
    else:
        for i in range(nnodes):
            for j in range(nnodes):

                if fun(i, j, nnodes) == target:
                    print(i+1, j+1)
                    return

    print("not found")

