def fft(vals, modulus, domain):
    if len(vals) == 1:
        return vals
    L = fft(vals[::2], modulus, domain[::2])
    R = fft(vals[1::2], modulus, domain[::2])
    o = [0 for i in vals]
    for i, (x, y) in enumerate(zip(L, R)):
        y_times_root = y * domain[i]
        o[i] = (x + y_times_root) % modulus
        o[i + len(L)] = (x - y_times_root) % modulus
    return o


def fft_test():
    modulus = 337
    values = [3, 1, 4, 1, 5, 9, 2, 6]
    domain = [1, 85, 148, 111, 336, 252, 189, 226]
    assert fft(values[::4], modulus, domain[::4]) == [8, 335]
    assert fft(values, modulus, domain) == [31, 70, 109, 74, 334, 181, 232, 4]


def modular_inverse(x, n):
    return pow(x, n - 2, n)


def inverse_fft(vals, modulus, domain):
    vals = fft(vals, modulus, domain)
    return [x * modular_inverse(len(vals), modulus) % modulus for x in [vals[0]] + vals[1:][::-1]]
