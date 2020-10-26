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

And now a magic trick: can we drop the odd half of the inputs,
 and deduce what they could have been, along with the missing outputs?

Should be possible, see [`partial_fft.py`](partial_fft.py) (**work in progress, _bugs_ being fixed**)

The idea is then to extend `N` values into `2N`, suiting an FFT with `2N` values, all even coefficients being `0`.

TODO: error-correction compatibility. Currently it works by padding N zero values to the coefficients, instead of making all evens zero.
Maybe it still works, or alternatively, refactor the `partial_fft`.

