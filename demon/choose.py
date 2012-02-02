def choose(n, k):
    r = 1
    d = n - k
    if k > n:
        return 0
    if d > k:
        k = d
        d = n - k
    while n > k:
        r *= n
        n -= 1
        while d > 1 and not r % d:
            r /= d
            d -= 1
    return r

if __name__=='__main__':
    import sys

    if len(sys.argv) == 3:
        n, k = (int(arg) for arg in sys.argv[1:])
    else:
        n = input('n = ')
        k = input('k = ')

    print choose(n, k)

