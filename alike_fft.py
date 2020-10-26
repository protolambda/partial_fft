def alike_fft(evens, odds, modulus, domain) -> (list, list):
    print("alike_fft evens ->>>>>>>>>>>>>>>>>>>>>>", evens, "odds", odds)
    # assert is_power_of_2(a_vals) and is_power_of_2(b_vals)
    if len(evens) == 1:
        x = evens[0]
        y = odds[0]
        y_times_root = y * domain[0]
        a = [(x + y_times_root) % modulus]
        b = [(x - y_times_root) % modulus]

        print("---")
        print("y_times_root", y_times_root)
        print("x", x)
        print("y", y)
        print("a", a)
        print("b", b)
        print("---")
        return a, b

    half = len(evens)
    halfhalf = half // 2

    L0, L1 = alike_fft(evens[::2], evens[1::2], modulus, domain[::2])
    print("alike_fft L", L0 + L1)
    R0, R1 = alike_fft(odds[::2], odds[1::2], modulus, domain[::2])
    print("alike_fft R", R0 + R1)
    a = [0] * half
    b = [0] * half
    # for i, (x, y) in enumerate(zip(L0+L1, R0+R1)):
    #     y_times_root = y*domain[i]
    #     a[i] = (x+y_times_root) % modulus
    #     b[i] = (x-y_times_root) % modulus

    print("evens len", len(evens))
    print("odds len", len(evens))
    print("halfhalf", halfhalf)
    print("a len", len(a))
    print("b len", len(b))

    print("L0: ", L0)
    print("R0: ", R0)
    for i, (x, y) in enumerate(zip(L0, R0)):
        print(f"alike_fft evens out: ", i, x, y)
        y_times_root = y * domain[i * 2]
        a[i] = (x + y_times_root) % modulus
        a[halfhalf + i] = (x - y_times_root) % modulus

    print("L1: ", L1)
    print("R1: ", R1)
    for i, (x, y) in enumerate(zip(L1, R1)):
        print(f"alike_fft odds out: ", i, x, y)
        y_times_root = y * domain[1 + i * 2]
        b[i] = (x + y_times_root) % modulus
        b[halfhalf + i] = (x - y_times_root) % modulus
    print("alike_fft a", a)
    print("alike_fft b", b)
    return a, b


modulus = 337
#
# values = [3, 5]
# domain = [1, 336]
# even_values = values[::2]
# odd_values = values[1::2]
# print("even_values", even_values)
# print("odd_values", odd_values)
#
# assert alike_fft(even_values, odd_values, modulus, domain) == ([8], [335])

values = [3, 1, 4, 1, 5, 9, 2, 6]
domain = [1, 85, 148, 111, 336, 252, 189, 226]
even_values = values[::2]
odd_values = values[1::2]
print("even_values", even_values)
print("odd_values", odd_values)
assert alike_fft(even_values, odd_values, modulus, domain) == ([31, 109, 334, 232], [70, 74, 181, 4])
