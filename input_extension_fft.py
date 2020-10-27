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


# values here is only the first half of the original input
# "a" are the even-indexed expected outputs

# The output is just the odds list
# "odds" is called odds since it's split that way on the deepest level, and merged back up.
# However, as an outside caller, it will match the 2nd half that's missing after all input half_values.
def input_extension_fft(half_values: list, a: list, modulus: int, domain: list, inverse_domain: list, inv2: int) -> list:
    if len(a) == 1:
        assert len(a) == 1
        assert len(half_values) == 1
        odd_val = input_extension_fft_base_case(half_values[0], a[0], modulus, inverse_domain[0])
        return [odd_val]

    half = len(a)
    assert len(half_values) == half
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

    odds0 = input_extension_fft(half_values[::2], L0, modulus, domain[::2], inverse_domain[::2], inv2)
    odds1 = input_extension_fft(half_values[1::2], R0, modulus, domain[::2], inverse_domain[::2], inv2)

    odds = [odds0[i // 2] if i % 2 == 0 else odds1[i // 2] for i in range(half)]
    return odds


def input_extension_fft_test(half_inputs, domain, even_outputs):
    modulus = 337

    inverse_of_2 = modular_inverse(2, modulus)

    assert len(half_inputs) * 2 == len(even_outputs) * 2 == len(domain)
    inverse_domain = [modular_inverse(d, modulus) for d in domain]

    resolved_second_half_inputs = input_extension_fft(half_inputs, even_outputs, modulus, domain, inverse_domain, inverse_of_2)
    print("resolved_second_half_inputs", resolved_second_half_inputs)

    reconstructed_inputs = half_inputs + resolved_second_half_inputs
    print("reconstructed_inputs", reconstructed_inputs)
    assert fft(reconstructed_inputs, modulus, domain)[::2] == even_outputs

# input_extension_fft_test([3], [1, 336], [8])
# input_extension_fft_test([3, 1, 4, 1], [1, 85, 148, 111, 336, 252, 189, 226], [31, 109, 334,  232])
# input_extension_fft_test([3, 1, 4, 1], [1, 85, 148, 111, 336, 252, 189, 226], [0, 0, 0, 0])

