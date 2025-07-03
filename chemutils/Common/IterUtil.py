
import itertools


def izip_longest(*args, **kwds):
    """
    Lazily zip and fill with fillvalues.

    >>> [(a, b) for a, b in izip_longest('ABCD', 'xy', fillvalue='-')]
    [('A', 'x'), ('B', 'y'), ('C', '-'), ('D', '-')]

    # izip_longest('ABCD', 'xy', fillvalue='-') --> Ax By C- D-
    """

    fillvalue = kwds.get("fillvalue")

    def sentinel(counter=([fillvalue] * (len(args) - 1)).pop):
        yield counter()  # yields the fillvalue, or raises IndexError

    fillers = itertools.repeat(fillvalue)
    iters = [itertools.chain(it, sentinel(), fillers) for it in args]
    try:
        for tup in itertools.izip(*iters):
            yield tup
    except IndexError:
        pass


def combinations(iterable, r):
    """Pure python implementation of itertools.combinations(), which
    is included in Python 2.6, but not available in 2.5

    All two-combinations of a list with 3 elements
    >>> [c for c in combinations([1, 2, 3], 2)]
    [(1, 2), (1, 3), (2, 3)]

    Only one three-combination of a three element list
    >>> [c for c in combinations([1, 2, 3], 3)]
    [(1, 2, 3)]

    Can't pick four from a list of three
    >>> [c for c in combinations([1, 2, 3], 4)]
    []
    """
    pool = tuple(iterable)
    n = len(pool)
    if r > n:
        return
    indices = range(r)
    yield tuple(pool[i] for i in indices)
    while True:
        for i in reversed(range(r)):
            if indices[i] != i + n - r:
                break
        else:
            return
        indices[i] += 1
        for j in range(i + 1, r):
            indices[j] = indices[j - 1] + 1
        yield tuple(pool[i] for i in indices)


def combinations_with_replacement(iterable, r):
    """Pure python implementation of itertools.combinations_with_replacement(), which
    is included in Python 2.6, but not available in 2.5
    """

    # combinations_with_replacement('ABC', 2) --> AA AB AC BB BC CC
    pool = tuple(iterable)
    n = len(pool)
    if not n and r:
        return
    indices = [0] * r
    yield tuple(pool[i] for i in indices)
    while True:
        for i in reversed(range(r)):
            if indices[i] != n - 1:
                break
        else:
            return
        indices[i:] = [indices[i] + 1] * (r - i)
        yield tuple(pool[i] for i in indices)


def permutations(iterable, r=None):
    """Pure python implementation of itertools.permutations(), which
    is included in Python 2.6, but not available in 2.5

    >>> [p for p in permutations([1, 2, 3], 2)]
    [(1, 2), (1, 3), (2, 1), (2, 3), (3, 1), (3, 2)]
    """
    pool = tuple(iterable)
    n = len(pool)
    r = n if r is None else r
    if r > n:
        return
    indices = range(n)
    cycles = range(n, n - r, -1)
    yield tuple(pool[i] for i in indices[:r])
    while n:
        for i in reversed(range(r)):
            cycles[i] -= 1
            if cycles[i] == 0:
                indices[i:] = indices[i + 1 :] + indices[i : i + 1]
                cycles[i] = n - i
            else:
                j = cycles[i]
                indices[i], indices[-j] = indices[-j], indices[i]
                yield tuple(pool[i] for i in indices[:r])
                break
        else:
            return


def _test():
    import doctest

    # print 'Actually about to call doctest'
    doctest.testmod()


if __name__ == "__main__":
    _test()
