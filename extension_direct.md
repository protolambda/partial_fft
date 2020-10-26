Can we retrieve values to extend the evaluation of a polynomial with, given some coefficients? With no more work than a single FFT? Yes!

```
Known: half the values (the unextended data), half the coefficients (the zero padding)

Inputs:    v0, v1, v2, v3     (values)
Outputs:   c0, c1, c2, c3     (coeffs)
Domain:    d0, d1, d2, d3

v0, v1 = known
v2, v2 = unknown

c0, c1 = known
c1, c2 = unknown


What we want:
fft([v0, v1, v2, v3], M, [d0, d1, d2, d3]) == [c0, c1, c2, c3]

fft([v0, v1, v2, v3], M, [d0, d1, d2, d3])
  L = fft([v0, v2], M, [d0, d2])
    L = fft([v0], M, [d0]) = v0
    R = fft([v2], M, [d2]) = v2
    oL0 = v0 + v2*d2
    oL1 = v0 - v2*d2
  R = fft([v1, v3], M, [d0, d2])
    L = fft([v1], M, [d0]) = v1
    R = fft([v3], M, [d2]) = v3
    oR0 = v1 + v3*d2
    oR1 = v1 - v3*d2
  o0 = oL0 + oR0*d0 = v0 + v2*d2 + (v1 + v3*d2)*d0 = v0 + v2*d2 + v1*d0 + v3*d2*d0
  o2 = oL0 - oR0*d0 = v0 + v2*d2 - (v1 + v3*d2)*d0 = v0 + v2*d2 - v1*d0 - v3*d2*d0
  o1 = oL1 + oR1*d1 = v0 - v2*d2 + (v1 - v3*d2)*d1 = v0 - v2*d2 + v1*d1 - v3*d2*d1
  o3 = oL1 - oR1*d1 = v0 - v2*d2 - (v1 - v3*d2)*d1 = v0 - v2*d2 - v1*d1 + v3*d2*d1

o0 = v0 + v1*d0 + v2*d2 + v3*d2*d0 = c0
o2 = v0 - v1*d0 + v2*d2 - v3*d2*d0 = c1
o1 = v0 + v1*d1 - v2*d2 - v3*d2*d1 = c2
o3 = v0 - v1*d1 - v2*d2 + v3*d2*d1 = c3

------------

Example solution (there may be more, and this is only just a small FFT)

c0 = 0
c1 = 0
c2 = unknown coefficient
c3 = unknown coefficient

----------
v0 + v1*d0 + v2*d2 + v3*d2*d0 = 0
v0 - v1*d0 + v2*d2 - v3*d2*d0 = 0

2*(v0 + v2*d2) = 0        (add the two equations)
v2*d2 = -v0
v2=-v0/d2

v0 - v1*d0 + v2*d2 = v3*d2*d0
v3 = (v0 - v1*d0 + v2*d2)/(d2*d0)
```
