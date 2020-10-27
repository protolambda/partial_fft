# FFT and recovery code, copy-pasta from research repo.
# Original by Vitalik & Dankrad here: https://github.com/ethereum/research/blob/master/verkle/fft.py
# Some minor modifications were made for debugging purposes, mostly `debug_bigs` for debugging.

#from py_ecc import optimized_bls12_381 as b
b = object()  # unused here


def _simple_ft(vals, modulus, roots_of_unity):
    L = len(roots_of_unity)
    o = []
    for i in range(L):
        last = b.Z1 if type(vals[0]) == tuple else 0
        for j in range(L):
            if type(vals[0]) == tuple:
                last = b.add(last, b.multiply(vals[j], roots_of_unity[(i*j)%L]))
            else:
                last += vals[j] * roots_of_unity[(i*j)%L]
        o.append(last if type(last) == tuple else last % modulus)
    # debug_bigs("simple ft out %d" % len(o), o)
    return o

def _fft(vals, modulus, roots_of_unity):
    # debug_bigs("_fft vals %d" % len(vals), vals)
    # debug_bigs("_fft rootz %d" % len(roots_of_unity), roots_of_unity)
    if len(vals) <= 4 and type(vals[0]) != tuple:
        #return vals
        return _simple_ft(vals, modulus, roots_of_unity)
    elif len(vals) == 1 and type(vals[0]) == tuple:
        return vals

    L = _fft(vals[::2], modulus, roots_of_unity[::2])
    R = _fft(vals[1::2], modulus, roots_of_unity[::2])
    o = [0 for i in vals]
    for i, (x, y) in enumerate(zip(L, R)):
        # print("i: %d, x: %d\ny: %d\nroot: %d" % (i, x, y, roots_of_unity[i]))
        # print(f"i {i}, x {x} y {y}")
        y_times_root = (b.multiply(y, roots_of_unity[i]) if type(y) == tuple else y*roots_of_unity[i]) % modulus
        # print("ytimesroot: %d" % y_times_root)
        o[i] = b.add(x, y_times_root) if type(x) == tuple else (x+y_times_root) % modulus
        # print("out i: %d" % o[i])
        o[i+len(L)] = b.add(x, b.neg(y_times_root)) if type(x) == tuple else (x-y_times_root) % modulus
        # print("out i+half: %d" % o[i+len(L)])

    # debug_bigs("_fft out %d" % len(o), o)
    return o

def expand_root_of_unity(root_of_unity, modulus):
    # Build up roots of unity
    rootz = [1, root_of_unity]
    while rootz[-1] != 1:
        rootz.append((rootz[-1] * root_of_unity) % modulus)
    return rootz

def fft(vals, modulus, root_of_unity, inv=False):
    rootz = expand_root_of_unity(root_of_unity, modulus)
    # Fill in vals with zeroes if needed
    if len(rootz) > len(vals) + 1:
        vals = vals + [0] * (len(rootz) - len(vals) - 1)
    if inv:
        # debug_bigs("expanded rootz", rootz)
        # Inverse FFT
        invlen = pow(len(vals), modulus-2, modulus)
        if type(vals[0]) == tuple:
            return [b.multiply(x, invlen) for x in
                    _fft(vals, modulus, rootz[:0:-1])]
        else:
            return [(x*invlen) % modulus for x in
                    _fft(vals, modulus, rootz[:0:-1])]
    else:
        # Regular FFT
        return _fft(vals, modulus, rootz[:-1])

# Evaluates f(x) for f in evaluation form
def inv_fft_at_point(vals, modulus, root_of_unity, x):
    if len(vals) == 1:
        return vals[0]
    # 1/2 in the field
    half = (modulus + 1)//2
    # 1/w
    inv_root = pow(root_of_unity, len(vals)-1, modulus)
    # f(-x) in evaluation form
    f_of_minus_x_vals = vals[len(vals)//2:] + vals[:len(vals)//2]
    # e(x) = (f(x) + f(-x)) / 2 in evaluation form
    evens = [(f+g) * half % modulus for f,g in zip(vals, f_of_minus_x_vals)]
    # o(x) = (f(x) - f(-x)) / 2 in evaluation form
    odds = [(f-g) * half % modulus for f,g in zip(vals, f_of_minus_x_vals)]
    # e(x^2) + coordinate * x * o(x^2) in evaluation form
    comb = [(o * x * inv_root**i + e) % modulus for i, (o, e) in enumerate(zip(odds, evens))]
    return inv_fft_at_point(comb[:len(comb)//2], modulus, root_of_unity ** 2 % modulus, x**2 % modulus)

def shift_domain(vals, modulus, root_of_unity, factor):
    if len(vals) == 1:
        return vals
    # 1/2 in the field
    half = (modulus + 1)//2
    # 1/w
    inv_factor = pow(factor, modulus - 2, modulus)
    half_length = len(vals)//2
    # f(-x) in evaluation form
    f_of_minus_x_vals = vals[half_length:] + vals[:half_length]
    # e(x) = (f(x) + f(-x)) / 2 in evaluation form
    evens = [(f+g) * half % modulus for f,g in zip(vals, f_of_minus_x_vals)]
    print('e', evens)
    # o(x) = (f(x) - f(-x)) / 2 in evaluation form
    odds = [(f-g) * half % modulus for f,g in zip(vals, f_of_minus_x_vals)]
    print('o', odds)
    shifted_evens = shift_domain(evens[:half_length], modulus, root_of_unity ** 2 % modulus, factor ** 2 % modulus)
    print('se', shifted_evens)
    shifted_odds = shift_domain(odds[:half_length], modulus, root_of_unity ** 2 % modulus, factor ** 2 % modulus)
    print('so', shifted_odds)
    return (
            [(e + inv_factor * o) % modulus for e, o in zip(shifted_evens, shifted_odds)] +
            [(e - inv_factor * o) % modulus for e, o in zip(shifted_evens, shifted_odds)]
    )

def shift_poly(poly, modulus, factor):
    factor_power = 1
    inv_factor = pow(factor, modulus - 2, modulus)
    o = []
    for p in poly:
        o.append(p * factor_power % modulus)
        factor_power = factor_power * inv_factor % modulus
    return o

def mul_polys(a, b, modulus, root_of_unity):
    rootz = [1, root_of_unity]
    while rootz[-1] != 1:
        rootz.append((rootz[-1] * root_of_unity) % modulus)
    # debug_bigs("mul polys rootz", rootz)
    if len(rootz) > len(a) + 1:
        a = a + [0] * (len(rootz) - len(a) - 1)
    if len(rootz) > len(b) + 1:
        b = b + [0] * (len(rootz) - len(b) - 1)

    # debug_bigs("mul polys a", a)
    # debug_bigs("mul polys b", b)
    x1 = _fft(a, modulus, rootz[:-1])
    # debug_bigs("mul polys x1", x1)
    x2 = _fft(b, modulus, rootz[:-1])
    # debug_bigs("mul polys x2", x2)
    return _fft([(v1*v2)%modulus for v1,v2 in zip(x1,x2)],
                modulus, rootz[:0:-1])


# Evaluate a polynomial at a point
def eval_poly_at(coeffs, x, modulus):
    y = 0
    power_of_x = 1
    for i, p_coeff in enumerate(coeffs):
        y = (y + (power_of_x * p_coeff)) % modulus
        power_of_x = (power_of_x * x) % modulus
    return y

def debug_bigs(msg, bigs):
    print(f"---{msg}---")
    for i, v in enumerate(bigs):
        print("#%4d: %s" % (i, "<nil>" if v is None else "%d" % v))
    print("")

# Calculates modular inverses [1/values[0], 1/values[1] ...]
def multi_inv(values, modulus):
    partials = [1]
    for i in range(len(values)):
        partials.append(partials[-1] * values[i] % modulus)
    inv = pow(partials[-1], modulus - 2, modulus)
    outputs = [0] * len(values)
    for i in range(len(values), 0, -1):
        outputs[i-1] = partials[i-1] * inv % modulus
        inv = inv * values[i-1] % modulus
    return outputs

# Generates q(x) = poly(k * x)
def p_of_kx(poly, modulus, k):
    o = []
    power_of_k = 1
    for x in poly:
        o.append(x * power_of_k % modulus)
        power_of_k = (power_of_k * k) % modulus
    # debug_bigs("pOfKX", o)
    return o

# Return (x - root**positions[0]) * (x - root**positions[1]) * ...
# possibly with a constant factor offset
def _zpoly(positions, modulus, roots_of_unity):
    # If there are not more than 4 positions, use the naive
    # O(n^2) algorithm as it is faster
    if len(positions) <= 4:
        root = [1]
        for pos in positions:
            x = roots_of_unity[pos]
            root.insert(0, 0)
            for j in range(len(root)-1):
                root[j] -= root[j+1] * x
        o = [x % modulus for x in root]
        # debug_bigs("_zpoly small out", o)
        return o
    else:
        # Recursively find the zpoly for even indices and odd
        # indices, operating over a half-size subgroup in each
        # case
        left = _zpoly([x//2 for x in positions if x%2 == 0],
                      modulus, roots_of_unity[::2])
        right = _zpoly([x//2 for x in positions if x%2 == 1],
                       modulus, roots_of_unity[::2])
        invroot = roots_of_unity[-1]
        # Offset the result for the odd indices, and combine
        # the two
        o = mul_polys(left, p_of_kx(right, modulus, invroot),
                      modulus, roots_of_unity[1])
    # Deal with the special case where mul_polys returns zero
    # when it should return x ^ (2 ** k) - 1
    if o == [0] * len(o):
        # debug_bigs("_zpoly zero out", o)
        return [1] + [0] * (len(o) - 1) + [modulus - 1]
    else:
        # debug_bigs("_zpoly out", o)
        return o

def zpoly(positions, modulus, root_of_unity):
    # Precompute roots of unity
    rootz = [1, root_of_unity]
    while rootz[-1] != 1:
        rootz.append((rootz[-1] * root_of_unity) % modulus)
    return _zpoly(positions, modulus, rootz[:-1])

def erasure_code_recover(vals, modulus, root_of_unity):
    # Generate the polynomial that is zero at the roots of unity
    # corresponding to the indices where vals[i] is None
    z = zpoly([i for i in range(len(vals)) if vals[i] is None],
              modulus, root_of_unity)
    # debug_bigs("z", z)
    zvals = fft(z, modulus, root_of_unity)
    # debug_bigs("zvals", zvals)

    # Pointwise-multiply (vals filling in zero at missing spots) * z
    # By construction, this equals vals * z
    vals_with_zeroes = [x or 0 for x in vals]
    p_times_z_vals = [x*y % modulus for x,y in zip(vals_with_zeroes, zvals)]
    # debug_bigs("p_times_z_vals", p_times_z_vals)
    p_times_z = fft(p_times_z_vals, modulus, root_of_unity, inv=True)
    # debug_bigs("p_times_z", p_times_z)

    # Keep choosing k values until the algorithm does not fail
    # Check only with primitive roots of unity
    for k in range(2, modulus):
        if pow(k, (modulus - 1) // 2, modulus) == 1:
            continue
        invk = pow(k, modulus - 2, modulus)
        # Convert p_times_z(x) and z(x) into new polynomials
        # q1(x) = p_times_z(k*x) and q2(x) = z(k*x)
        # These are likely to not be 0 at any of the evaluation points.
        p_times_z_of_kx = [x * pow(k, i, modulus) % modulus
                           for i, x in enumerate(p_times_z)]
        # debug_bigs("p_times_z_of_kx", p_times_z_of_kx)

        p_times_z_of_kx_vals = fft(p_times_z_of_kx, modulus, root_of_unity)
        # debug_bigs("p_times_z_of_kx_vals", p_times_z_of_kx_vals)

        z_of_kx = [x * pow(k, i, modulus) for i, x in enumerate(z)]
        # debug_bigs("z_of_kx", z_of_kx)
        z_of_kx_vals = fft(z_of_kx, modulus, root_of_unity)
        # debug_bigs("z_of_kx_vals", z_of_kx_vals)

        # Compute q1(x) / q2(x) = p(k*x)
        inv_z_of_kv_vals = multi_inv(z_of_kx_vals, modulus)
        # debug_bigs("inv_z_of_kv_vals", inv_z_of_kv_vals)
        p_of_kx_vals = [x*y % modulus for x,y in
                        zip(p_times_z_of_kx_vals, inv_z_of_kv_vals)]
        # debug_bigs("p_of_kx_vals", p_of_kx_vals)
        p_of_kx = fft(p_of_kx_vals, modulus, root_of_unity, inv=True)
        # debug_bigs("p_of_kx", p_of_kx)

        # print("invk: %d" % invk)
        # Given q3(x) = p(k*x), recover p(x)
        p_of_x = [x * pow(invk, i, modulus) % modulus
                  for i, x in enumerate(p_of_kx)]
        # debug_bigs("p_of_x", p_of_x)
        output = fft(p_of_x, modulus, root_of_unity)
        # debug_bigs("output", output)

        # Check that the output matches the input
        success = True
        for inpd, outd in zip(vals, output):
            success *= (inpd is None or inpd == outd)
        if not success:
            continue

        # Output the evaluations if all good
        return output


# ====================================================================================================================


rootOfUnityCandidates = {
    # 512: 12531186154666751577774347439625638674013361494693625348921624593362229945844,
    # 256: 21071158244812412064791010377580296085971058123779034548857891862303448703672,
    # 128: 3535074550574477753284711575859241084625659976293648650204577841347885064712,
    # 64:  6460039226971164073848821215333189185736442942708452192605981749202491651199,
    # 32:  32311457133713125762627935188100354218453688428796477340173861531654182464166,
    # 16:  35811073542294463015946892559272836998938171743018714161809767624935956676211,
    2         : 52435875175126190479447740508185965837690552500527637822603658699938581184512,
    4         : 3465144826073652318776269530687742778270252468765361963008,
    8         : 23674694431658770659612952115660802947967373701506253797663184111817857449850,
    16        : 14788168760825820622209131888203028446852016562542525606630160374691593895118,
    32        : 36581797046584068049060372878520385032448812009597153775348195406694427778894,
    64        : 31519469946562159605140591558550197856588417350474800936898404023113662197331,
    128       : 47309214877430199588914062438791732591241783999377560080318349803002842391998,
    256       : 36007022166693598376559747923784822035233416720563672082740011604939309541707,
    512       : 4214636447306890335450803789410475782380792963881561516561680164772024173390,
    1024      : 22781213702924172180523978385542388841346373992886390990881355510284839737428,
    2048      : 49307615728544765012166121802278658070711169839041683575071795236746050763237,
    4096      : 39033254847818212395286706435128746857159659164139250548781411570340225835782,
    8192      : 32731401973776920074999878620293785439674386180695720638377027142500196583783,
    16384     : 39072540533732477250409069030641316533649120504872707460480262653418090977761,
    32768     : 22872204467218851938836547481240843888453165451755431061227190987689039608686,
    65536     : 15076889834420168339092859836519192632846122361203618639585008852351569017005,
    131072    : 15495926509001846844474268026226183818445427694968626800913907911890390421264,
    262144    : 20439484849038267462774237595151440867617792718791690563928621375157525968123,
    524288    : 37115000097562964541269718788523040559386243094666416358585267518228781043101,
}

WIDTH = 256

from das_fft import das_fft
from classic_fft import modular_inverse

MODULUS = 52435875175126190479447740508185965837690552500527637822603658699938581184513
ROOT_OF_UNITY = rootOfUnityCandidates[WIDTH]
domain = expand_root_of_unity(ROOT_OF_UNITY, MODULUS)[:-1]

inverse_domain = [modular_inverse(d, MODULUS) for d in domain]
inverse_of_2 = modular_inverse(2, MODULUS)

# even data will be the original data of interest
even_data = [i*42 % MODULUS for i in range(WIDTH//2)]
debug_bigs("even_data", even_data)

odd_data = das_fft(even_data, MODULUS, domain, inverse_domain, inverse_of_2)
debug_bigs("odd data", odd_data)

extended_data = [even_data[i // 2] if i % 2 == 0 else odd_data[i // 2] for i in range(WIDTH)]
debug_bigs("extended data", extended_data)

# The odd_data was constructed in such a way that the second half of coefficients are all zero.
coeffs_poc = fft(extended_data, MODULUS, ROOT_OF_UNITY, inv=True)
assert coeffs_poc[WIDTH//2:] == [0]*(WIDTH//2)

import random

subset_data = [i for i in extended_data]
for i in random.sample(range(WIDTH), WIDTH//2):
    subset_data[i] = None

debug_bigs("subset", subset_data)

recovered = erasure_code_recover(subset_data, MODULUS, ROOT_OF_UNITY)
debug_bigs("recovered", recovered)

# Recovered our data!
assert recovered[::2] == even_data
