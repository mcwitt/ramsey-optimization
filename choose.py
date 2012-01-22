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
