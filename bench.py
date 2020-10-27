# bench
import time
from classic_fft import fft, modular_inverse
from partial_fft import partial_fft
from other_partial_fft import other_partial_fft
from input_extension_fft import input_extension_fft
from output_extension_fft import output_extension_fft

MODULUS = 52435875175126190479447740508185965837690552500527637822603658699938581184513

rootOfUnityCandidates = {
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

def expand_root_of_unity(root_of_unity, modulus):
    # Build up roots of unity
    rootz = [1, root_of_unity]
    while rootz[-1] != 1:
        rootz.append((rootz[-1] * root_of_unity) % modulus)
    return rootz

def bench(scale: int):
    width = 2**scale

    modulus = MODULUS
    root_of_unity = rootOfUnityCandidates[width]
    domain = expand_root_of_unity(root_of_unity, modulus)[:-1]
    assert len(domain) == width
    inverse_domain = [modular_inverse(d, modulus) for d in domain]
    inverse_of_2 = modular_inverse(2, modulus)

    data = [i for i in range(width)]
    N = 200
    start = time.time()
    for i in range(N):
        output = fft(data, modulus, domain)
        assert len(output) == len(data)
    end = time.time()
    diff = end-start
    ns = diff*1e9
    print("          FFT: scale_%-5d %10d ops %15.0f ns/op" % (scale, N, ns/N))

    partial_data = data[:width//2]
    even_outputs = [0] * (width//2)
    start = time.time()
    for i in range(N):
        other_inputs, other_outputs = partial_fft(partial_data, even_outputs, modulus, domain, inverse_domain, inverse_of_2)
        assert len(other_inputs) == width//2
        assert len(other_outputs) == width//2
    end = time.time()
    diff = end-start
    ns = diff*1e9
    print("  Partial FFT: scale_%-5d %10d ops %15.0f ns/op" % (scale, N, ns/N))

    partial_data = data[width//2:]
    even_outputs = [0] * (width//2)
    start = time.time()
    for i in range(N):
        other_inputs, other_outputs = other_partial_fft(partial_data, even_outputs, modulus, domain, inverse_domain, inverse_of_2)
        assert len(other_inputs) == width//2
        assert len(other_outputs) == width//2
    end = time.time()
    diff = end-start
    ns = diff*1e9
    print("OtherPart FFT: scale_%-5d %10d ops %15.0f ns/op" % (scale, N, ns/N))

    partial_data = data[:width//2]
    even_outputs = [0] * (width//2)
    start = time.time()
    for i in range(N):
        other_inputs = input_extension_fft(partial_data, even_outputs, modulus, domain, inverse_domain, inverse_of_2)
        assert len(other_inputs) == width//2
    end = time.time()
    diff = end-start
    ns = diff*1e9
    print(" InputExt FFT: scale_%-5d %10d ops %15.0f ns/op" % (scale, N, ns/N))

    partial_data = data[:width//2]
    even_outputs = [0] * (width//2)
    start = time.time()
    for i in range(N):
        other_outputs = output_extension_fft(partial_data, even_outputs, modulus, domain, inverse_domain, inverse_of_2)
        assert len(other_outputs) == width//2
    end = time.time()
    diff = end-start
    ns = diff*1e9
    print("OutputExt FFT: scale_%-5d %10d ops %15.0f ns/op" % (scale, N, ns/N))


bench(8)

