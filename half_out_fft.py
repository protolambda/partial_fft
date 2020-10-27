def half_out_fft(values, modulus, domain) -> list:
    if len(values) == 2:
        x = values[0]
        y = values[1]
        y_times_root = y * domain[0]
        a = [(x + y_times_root) % modulus]
        return a

    half = len(values) // 2
    halfhalf = half // 2

    L0 = half_out_fft(values[::2], modulus, domain[::2])
    print("half_fft L", L0)
    R0 = half_out_fft(values[1::2], modulus, domain[::2])
    print("half_fft R", R0)
    a = [0] * half

    print("values len", len(values))
    print("half", half)
    print("halfhalf", halfhalf)
    print("a len", len(a))

    for i, (x, y) in enumerate(zip(L0, R0)):
        print(f"half_fft evens out: ", i, x, y)
        y_times_root = y * domain[i * 2]
        a[i] = (x + y_times_root) % modulus
        a[halfhalf + i] = (x - y_times_root) % modulus

    print("half_fft a", a)
    return a


modulus = 337

inputs = [3, 5]
domain = [1, 336]
# Full expected outputs: [8, 335]
assert half_out_fft(inputs, modulus, domain[::4]) == [8]

inputs = [3, 1, 4, 1, 5, 9, 2, 6]
domain = [1, 85, 148, 111, 336, 252, 189, 226]
# Full expected outputs: [31, 70, 109, 74, 334, 181, 232, 4]
assert half_out_fft(inputs, modulus, domain) == [31, 109, 334, 232]
