from classic_fft import fft, modular_inverse


def input_extension_fft_base_case(even_val: int, a: int, modulus: int, inverse_domain: int) -> int:
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
    return odd_val


def input_extension_fft(evens: list, odds: list, a: list, modulus: int, domain: list, inverse_domain: list, inv2: int) -> list:
    if len(evens) == 1:
        assert len(a) == 1
        assert odds[0] == None
        odd_val = input_extension_fft_base_case(evens[0], a[0], modulus, inverse_domain[0])
        return [odd_val]

    half = len(evens)
    halfhalf = half // 2

    L0 = [0] * halfhalf
    R0 = [0] * halfhalf
    for i, (a_half0, a_half1) in enumerate(zip(a[:halfhalf], a[halfhalf:])):
        # x = L0[i]
        # y = R0[i]
        # y_times_root = y*domain[i*2] = R0[i]*domain[i*2]
        # a_half0 = x+y_times_root = L0[i] + R0[i]*domain[i*2]
        # a_half1 = x-y_times_root = L0[i] - R0[i]*domain[i*2]
        # L0[i] = (a_half0+a_half1)//2
        L0[i] = (((a_half0 + a_half1) % modulus) * inv2) % modulus
        # Either: R0[i]*domain[i*2] =  a_half0 - L0[i]
        #     Or: R0[i]*domain[i*2] = -a_half1 + L0[i]
        # Either: R0[i] = (a_half0 - L0[i]) / domain[i*2]
        #     Or: R0[i] = (-a_half1 + L0[i]) / domain[i*2]
        R0[i] = (((a_half0 - L0[i]) % modulus) * inverse_domain[i * 2]) % modulus
        # R0[i] = (((-a_half1 + L0[i]) % modulus) * inverse_domain[i*2]) % modulus

    odds0 = input_extension_fft(evens[::2], evens[1::2], L0, modulus, domain[::2], inverse_domain[::2], inv2)
    odds1 = input_extension_fft(odds[::2], odds[1::2], R0, modulus, domain[::2], inverse_domain[::2], inv2)

    odds = [odds0[i // 2] if i % 2 == 0 else odds1[i // 2] for i in range(half)]
    return odds


def input_extension_fft_test(half_values, domain, even_coeffs):
    modulus = 337

    inverse_of_2 = modular_inverse(2, modulus)

    assert len(half_values) * 2 == len(even_coeffs) * 2 == len(domain)
    inverse_domain = [modular_inverse(d, modulus) for d in domain]

    half_full_values = half_values + [None] * len(half_values)
    resolved_odd_values = input_extension_fft(half_full_values[::2], half_full_values[1::2], even_coeffs, modulus, domain, inverse_domain, inverse_of_2)
    print("resolved_odd_values", resolved_odd_values)

    reconstructed_values = half_values + resolved_odd_values
    print("reconstructed_values", reconstructed_values)
    assert fft(reconstructed_values, modulus, domain)[::2] == even_coeffs

# input_extension_fft_test([3], [1, 336], [8])
# input_extension_fft_test([3, 1, 4, 1], [1, 85, 148, 111, 336, 252, 189, 226], [31, 109, 334,  232])
# input_extension_fft_test([3, 1, 4, 1], [1, 85, 148, 111, 336, 252, 189, 226], [0, 0, 0, 0])

