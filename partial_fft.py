from classic_fft import fft, inverse_fft, modular_inverse


def partial_fft_base_case(even_val: int, a: int, modulus: int, inverse_domain: int) -> (int, int):
    # In a regular FFT:
    # x = evens[0]
    # y = odds[0]
    # y_times_root = y*domain[0]
    # a = [(x+y_times_root) % modulus]
    # b = [(x-y_times_root) % modulus]
    # So:
    # a-x = y_times_root
    # y = y_times_root/domain[0] = y_times_root * inverse_domain[0]
    x = even_val
    y_times_root = (a - x) % modulus
    y = (y_times_root * inverse_domain) % modulus
    odd_val = y
    b = (x - y_times_root) % modulus
    print("---")
    print("y_times_root", y_times_root)
    print("x", x)
    print("y", y)
    print("a", a)
    print("b", b)
    return odd_val, b


def partial_fft(evens: list, a: list, modulus: int, domain: list, inverse_domain: list) -> (list, list):
    if len(evens) == 1:
        assert len(a) == 1
        odd_val, b = partial_fft_base_case(evens[0], a[0], modulus, inverse_domain[0])
        return [odd_val], [b]

    half = len(evens)
    halfhalf = half // 2

    inverse_of_2 = modular_inverse(2, modulus)
    L0 = [0] * halfhalf
    R0 = [0] * halfhalf
    for i, (a_half0, a_half1) in enumerate(zip(a[:halfhalf], a[halfhalf:])):
        # x = L0[i]
        # y = R0[i]
        # y_times_root = y*domain[i*2] = R0[i]*domain[i*2]
        # a_half0 = x+y_times_root = L0[i] + R0[i]*domain[i*2]
        # a_half1 = x-y_times_root = L0[i] - R0[i]*domain[i*2]
        # L0[i] = (a_half0+a_half1)//2
        L0[i] = (((a_half0 + a_half1) % modulus) * inverse_of_2) % modulus
        # Either: R0[i]*domain[i*2] =  a_half0 - L0[i]
        #     Or: R0[i]*domain[i*2] = -a_half1 + L0[i]
        # Either: R0[i] = (a_half0 - L0[i]) / domain[i*2]
        #     Or: R0[i] = (-a_half1 + L0[i]) / domain[i*2]
        R0[i] = (((a_half0 - L0[i]) % modulus) * inverse_domain[i * 2]) % modulus
        # R0[i] = (((-a_half1 + L0[i]) % modulus) * inverse_domain[i*2]) % modulus
    print("L0: ", L0)
    print("R0: ", R0)

    b = [0] * half
    odds0, L1 = partial_fft(evens[::2], L0, modulus, domain[::2], inverse_domain[::2])
    odds1, R1 = partial_fft(evens[1::2], R0, modulus, domain[::2], inverse_domain[::2])

    print("L1: ", L1)
    print("R1: ", R1)

    for i, (x, y) in enumerate(zip(L1, R1)):
        print(f"alike_fft odds out: ", i, x, y)
        y_times_root = y * domain[1 + i * 2]
        b[i] = (x + y_times_root) % modulus
        b[halfhalf + i] = (x - y_times_root) % modulus

    print("a: ", a)
    print("b: ", b)

    return [odds0[i // 2] if i % 2 == 0 else odds1[i // 2] for i in range(half)], b


modulus = 337

# For 2 values: [8, 335]
values = [3, 5]
domain = [1, 336]
inverse_domain = [modular_inverse(d, modulus) for d in domain]
coeffs = [8, 335]

even_values = values[::2]
even_coeffs = coeffs[::2]
resolved_odd_values, resolved_odd_coeffs = partial_fft(even_values, even_coeffs, modulus, domain[::4], inverse_domain[::4])
print("resolved_odd_values", resolved_odd_values)
print("resolved_odd_coeffs", resolved_odd_coeffs)

reconstructed_values = [even_values[i // 2] if i % 2 == 0 else resolved_odd_values[i // 2] for i in range(2)]
reconstructed_coeffs = [even_coeffs[i // 2] if i % 2 == 0 else resolved_odd_coeffs[i // 2] for i in range(2)]
assert fft(reconstructed_values, modulus, domain) == reconstructed_coeffs
assert inverse_fft(reconstructed_coeffs, modulus, domain) == reconstructed_values

# For 8 values:
values = [3, 1, 4, 1, 5, 9, 2, 6]
domain = [1, 85, 148, 111, 336, 252, 189, 226]
inverse_domain = [modular_inverse(d, modulus) for d in domain]
coeffs = [31, 70, 109, 74, 334, 181, 232, 4]

even_values = values[::2]
even_coeffs = coeffs[::2]
resolved_odd_values, resolved_odd_coeffs = partial_fft(even_values, even_coeffs, modulus, domain, inverse_domain)
print("resolved_odd_values", resolved_odd_values)
print("resolved_odd_coeffs", resolved_odd_coeffs)

reconstructed_values = [even_values[i // 2] if i % 2 == 0 else resolved_odd_values[i // 2] for i in range(8)]
reconstructed_coeffs = [even_coeffs[i // 2] if i % 2 == 0 else resolved_odd_coeffs[i // 2] for i in range(8)]
assert fft(reconstructed_values, modulus, domain) == reconstructed_coeffs
assert inverse_fft(reconstructed_coeffs, modulus, domain) == reconstructed_values


values = [3, 1, 4, 1]
domain = [1, 85, 148, 111, 336, 252, 189, 226]
inverse_domain = [modular_inverse(d, modulus) for d in domain]
coeffs = [0, 0, 0, 0]

resolved_odd_values, resolved_odd_coeffs = partial_fft(values, coeffs, modulus, domain, inverse_domain)
print("resolved_odd_values", resolved_odd_values)
print("resolved_odd_coeffs", resolved_odd_coeffs)

reconstructed_values = [even_values[i // 2] if i % 2 == 0 else resolved_odd_values[i // 2] for i in range(8)]
reconstructed_coeffs = [even_coeffs[i // 2] if i % 2 == 0 else resolved_odd_coeffs[i // 2] for i in range(8)]
assert fft(reconstructed_values, modulus, domain) == reconstructed_coeffs
assert inverse_fft(reconstructed_coeffs, modulus, domain) == reconstructed_values

