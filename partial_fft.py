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
    # print("---")
    # print("y_times_root", y_times_root)
    # print("x", x)
    # print("y", y)
    # print("a", a)
    # print("b", b)
    # print("---")
    return odd_val, b

# values here is only the first half of the original input
# "a" are the even-indexed expected outputs
# The output is (odds, b)
# "odds" is called odds since it's split that way on the deepest level, and merged back up.
# However, as an outside caller, it will match the 2nd half that's missing after all input half_values.
# And "b" are the derived odd-index outputs.
def partial_fft(half_values: list, a: list, modulus: int, domain: list, inverse_domain: list, inv2: int) -> (list, list):
    # print("alike_fft evens ->>>>>>>>>>>>>>>>>>>>>>", evens, "odds", odds)
    if len(a) == 1:
        assert len(a) == 1
        assert len(half_values) == 1
        odd_val, b = partial_fft_base_case(half_values[0], a[0], modulus, inverse_domain[0])
        return [odd_val], [b]

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
    # print("L0: ", L0)
    # print("R0: ", R0)

    b = [0] * half
    odds0, L1 = partial_fft(half_values[::2], L0, modulus, domain[::2], inverse_domain[::2], inv2)
    odds1, R1 = partial_fft(half_values[1::2], R0, modulus, domain[::2], inverse_domain[::2], inv2)

    # print("L1: ", L1)
    # print("R1: ", R1)

    for i, (x, y) in enumerate(zip(L1, R1)):
        # print(f"alike_fft odds out: ", i, x, y)
        y_times_root = y * domain[1 + i * 2]
        b[i] = (x + y_times_root) % modulus
        b[halfhalf + i] = (x - y_times_root) % modulus

    # print("a: ", a)
    # print("b: ", b)

    odds = [odds0[i // 2] if i % 2 == 0 else odds1[i // 2] for i in range(half)]
    return odds, b


def partial_fft_test(half_inputs, domain, even_outputs):
    modulus = 337

    inverse_of_2 = modular_inverse(2, modulus)

    assert len(half_inputs) * 2 == len(even_outputs) * 2 == len(domain)
    inverse_domain = [modular_inverse(d, modulus) for d in domain]

    resolved_second_half_inputs, resolved_odd_outputs = partial_fft(half_inputs, even_outputs, modulus, domain, inverse_domain, inverse_of_2)
    print("resolved_second_half_inputs", resolved_second_half_inputs)
    print("resolved_odd_outputs", resolved_odd_outputs)

    reconstructed_inputs = half_inputs + resolved_second_half_inputs
    reconstructed_outputs = [even_outputs[i // 2] if i % 2 == 0 else resolved_odd_outputs[i // 2] for i in range(len(even_outputs)+len(resolved_odd_outputs))]
    print("reconstructed_inputs", reconstructed_inputs)
    print("reconstructed_outputs", reconstructed_outputs)
    assert fft(reconstructed_inputs, modulus, domain) == reconstructed_outputs
    assert inverse_fft(reconstructed_outputs, modulus, domain) == reconstructed_inputs


# partial_fft_test([3], [1, 336], [8])
# partial_fft_test([3, 1, 4, 1], [1, 85, 148, 111, 336, 252, 189, 226], [31, 109, 334,  232])
# partial_fft_test([3, 1, 4, 1], [1, 85, 148, 111, 336, 252, 189, 226], [0, 0, 0, 0])

