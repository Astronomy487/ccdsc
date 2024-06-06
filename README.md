This is a JavaScript library to handle arithmetic with the [complex numbers](https://en.wikipedia.org/wiki/Complex_number), the [dual numbers](https://en.wikipedia.org/wiki/Dual_number), and the [split-complex numbers)](https://en.wikipedia.org/wiki/Split-complex_number) all together. The complex numbers introduce a unit $i$ such that $i^2=-1$. Likewise, the dual numbers introduce $\epsilon$ such that $\epsilon^2=0$, and the split-complex numbers introduce $j$ such that $j^2=1$.

You can include this in your own projects by downloading and using [ccdsc.js](https://github.com/Astronomy487/ccdsc/blob/main/ccdsc.js). Please share with me if you do anything with this library! You can also reference the version hosted on GitHub Pages:

```html
<script src="https://astronomy487.github.io/ccdsc/ccdsc.js">
```

Everything is housed in the global object `ccdsc`, which has the following methods:

- `ccdsc.make(x)`
- `ccdsc.add(a, b, ...)`
- `ccdsc.divide(a, b)`
- `ccdsc.subtract(a, b)`
- `ccdsc.abs(x)`
- `ccdsc.abs(x)`
- `ccdsc.sign(x)`
- `ccdsc.exp(x)`
- `ccdsc.expm1(x)`
- `ccdsc.log(x)`
- `ccdsc.log10(x)`
- `ccdsc.log2(x)`
- `ccdsc.log1p(x)`
- `ccdsc.pow(a, b)`
- `ccdsc.sqrt(x)`
- `ccdsc.cbrt(x)`

It also has 24 trigonometric/hyperbolic functions and their inverses, named according to the pattern {∅|a}{sin|cos|tan|csc|sec|cot}{∅|h}. **The inverse trigonometric/hyperbolic functions are not working correctly right now.**

Once you have a CCDSC number, you can get it to text using `.toString()` (or `.toTex()` or `.toHtml()`, if you prefer).

More information about this library is on [this webpage](https://astronomy487.github.io/ccdsc).