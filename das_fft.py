from classic_fft import fft, inverse_fft, modular_inverse

# def das_fft_base_case(a: int) -> int:
#     # In a regular FFT:
#     # x = evens[0]
#     # y = odds[0]
#     # y_times_root = y*domain[0]
#     # a = [(x+y_times_root) % modulus]
#     # b = [(x-y_times_root) % modulus]
#     # So:
#     # x = b+y_times_root
#     # y = 0  # Second half of the inputs is 0
#     # y_times_root = (y * domain) % modulus = 0
#     # x = (a - y_times_root) % modulus
#     # x = a
#     # b = (x - y_times_root) % modulus
#     # b = x
#     # a, b, x all the same!
#     return a

# values here is only the first half of the original input
# "a" are the even-indexed expected outputs
# The returned list is "b", the derived odd-index outputs.
def das_fft(a: list, modulus: int, domain: list, inverse_domain: list, inv2: int) -> list:
    if len(a) == 2:
        a_half0 = a[0]
        a_half1 = a[1]
        x = (((a_half0 + a_half1) % modulus) * inv2) % modulus
        # y = (((a_half0 - x) % modulus) * inverse_domain[0]) % modulus     # inverse_domain[0] will always be 1
        y = (a_half0 - x) % modulus

        y_times_root = y * domain[1]
        return [
            (x + y_times_root) % modulus,
            (x - y_times_root) % modulus
        ]

    if len(a) == 1:  # for illustration purposes, neat simplification. Inputs are always a power of two, dead code.
        return a

    half = len(a)
    halfhalf = half // 2

    L0 = [0] * halfhalf
    R0 = [0] * halfhalf
    for i, (a_half0, a_half1) in enumerate(zip(a[:halfhalf], a[halfhalf:])):
        L0[i] = (((a_half0 + a_half1) % modulus) * inv2) % modulus
        R0[i] = (((a_half0 - L0[i]) % modulus) * inverse_domain[i * 2]) % modulus

    L1 = das_fft(L0, modulus, domain[::2], inverse_domain[::2], inv2)
    R1 = das_fft(R0, modulus, domain[::2], inverse_domain[::2], inv2)

    b = [0] * half
    for i, (x, y) in enumerate(zip(L1, R1)):
        y_times_root = y * domain[1 + i * 2]
        b[i] = (x + y_times_root) % modulus
        b[halfhalf + i] = (x - y_times_root) % modulus

    return b


def das_fft_test(domain, even_outputs):
    modulus = 337

    inverse_of_2 = modular_inverse(2, modulus)

    assert len(even_outputs) * 2 == len(domain)
    inverse_domain = [modular_inverse(d, modulus) for d in domain]

    half = len(even_outputs)

    resolved_odd_outputs = das_fft(even_outputs, modulus, domain, inverse_domain, inverse_of_2)
    print("resolved_odd_outputs", resolved_odd_outputs)

    reconstructed_outputs = [even_outputs[i // 2] if i % 2 == 0 else resolved_odd_outputs[i // 2] for i in range(2*half)]
    print("reconstructed_outputs", reconstructed_outputs)
    reconstructed_inputs = inverse_fft(reconstructed_outputs, modulus, domain)
    print("reconstructed_inputs", reconstructed_inputs)

    assert reconstructed_inputs[half:] == [0] * half
    assert fft(reconstructed_inputs, modulus, domain) == reconstructed_outputs

#
# das_fft_test([1, 336], [8])
# das_fft_test([1, 85, 148, 111, 336, 252, 189, 226], [31, 109, 334,  232])
# das_fft_test([1, 85, 148, 111, 336, 252, 189, 226], [0, 0, 0, 0])

