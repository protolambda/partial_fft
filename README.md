# Partial FFT

This repo shows how to build a special kind of FFT:
- it takes N/2 values and N/2 coefficients
- it outputs some valid set of N/2 values and N/2 coefficients

First, have a look at the classic FFT: [`classic_fft.py`](classic_fft.py)

Then note that we can:
 - separate the even and odd outputs (`a` and `b`)
 - separate the even and odd inputs (`evens` and `odds`)
 
It's still an FFT, and does the same thing.
See [`alike_fft.py`](alike_fft.py).

Then, can we drop some of the work? Sure,
we can not generate any of the odd-index coefficients.
See [`half_out.py`](half_out_fft.py).

And now a magic trick: can we drop the 2nd half of the inputs,
 and deduce what they could have been, along with the missing outputs?

Should be possible, see [`partial_fft.py`](partial_fft.py) (**testing in progress, there may be _bugs_**)

The idea is then to extend `N` values into `2N`, suiting an FFT with `2N` values, all even coefficients being `0`.

TODO: error-correction compatibility. Currently it works by padding N zeros to the coefficients, instead of making all evens zero.
Maybe it still works, or alternatively, refactor the `partial_fft`.

Alternatively, use the `partial FFT` to compute the inverse FFT. Provide N consecutive zeroes, and get the odd missing values.

And can it be even faster? Yes, kind of. For data-extension we are only interested in doubling the values, not the coefficients.
So we can run a specialized form of the partial FFT which only does half of the work: [`input_extension_fft.py`](input_extension_fft.py).

## Benchmarks

```
         FFT: scale_8            200 ops         1222398 ns/op
 Partial FFT: scale_8            200 ops         1106892 ns/op
InputExt FFT: scale_8            200 ops          544612 ns/op
```

