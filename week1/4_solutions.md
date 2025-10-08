# Detailed Mathematical Derivations 

Given
- Hazard (failure rate): $z(t)=h(t)=\dfrac{0.4}{0.2t+1}$
- Cumulative hazard: $Z(t)=H(t)=\displaystyle\int_0^t z(u)\,du$

---

## 1) Cumulative hazard Z(t)

Substitution: $w=0.2u+1,\ dw=0.2\,du,\ du=5\,dw$, limits $u:0\to t \Rightarrow w:1\to 0.2t+1$:
$$
Z(t)=\int_0^t \frac{0.4}{0.2u+1}\,du
=\int_{1}^{0.2t+1}\frac{0.4}{w}\cdot 5\,dw
=2\int_{1}^{0.2t+1}\frac{1}{w}\,dw
=2\ln(0.2t+1).
$$

Result:
$$\boxed{Z(t)=2\ln(0.2t+1)}$$

---

## 2) Reliability R(t)

Using $R(t)=\exp(-Z(t))$:
$$
R(t)=\exp\!\big(-2\ln(0.2t+1)\big)
=(0.2t+1)^{-2}
=\frac{1}{(0.2t+1)^2}.
$$

---

## 3) PDF and CDF

$$
f(t)=z(t)\,R(t)=\frac{0.4}{0.2t+1}\cdot\frac{1}{(0.2t+1)^2}
=\frac{0.4}{(0.2t+1)^3},
\qquad
F(t)=1-R(t)=1-\frac{1}{(0.2t+1)^2}.
$$

---

## 4) Mean Time to Failure (MTTF)

Substitution $u=0.2t+1,\ du=0.2\,dt,\ dt=5\,du$, limits $t:0\to\infty \Rightarrow u:1\to\infty$:
$$
\mathrm{MTTF}=\int_0^\infty R(t)\,dt
=\int_0^\infty \frac{1}{(0.2t+1)^2}\,dt
=5\int_{1}^{\infty}u^{-2}\,du
=5\left[-\frac{1}{u}\right]_{1}^{\infty}
=5.
$$

Result:
$$\boxed{\mathrm{MTTF}=5}$$

---

## 5) Second moment and variance

Direct form with $u=0.2t+1,\ t=5(u-1),\ dt=5\,du$:
$$
E[T^2]=\int_0^\infty t^2 f(t)\,dt
=\int_0^\infty t^2 \frac{0.4}{(0.2t+1)^3}\,dt
=50\int_{1}^{\infty}\frac{(u-1)^2}{u^3}\,du
=50\int_{1}^{\infty}\left(\frac{1}{u}-\frac{2}{u^2}+\frac{1}{u^3}\right)du.
$$
Since $\int_{1}^{\infty}\frac{1}{u}\,du=\infty$, we obtain
$$
\boxed{E[T^2]=\infty\ \Rightarrow\ \mathrm{Var}(T)=\infty,\ \sigma=\infty}
$$

(Equivalent via parts: $E[T^2]=-\int_0^\infty t^2\,dR(t)=[-t^2R(t)]_0^\infty+2\int_0^\infty tR(t)\,dt$, and $t^2R(t)\not\to0$.)

---

## 6) Median and MRLT

Median from $R(t_{0.5})=0.5$:
$$
\frac{1}{(0.2t_{0.5}+1)^2}=0.5
\Rightarrow (0.2t_{0.5}+1)^2=2
\Rightarrow 0.2t_{0.5}+1=\sqrt{2}
\Rightarrow t_{0.5}=5(\sqrt{2}-1)\approx 2.0711.
$$

Mean residual life:
$$
m(t)=\frac{\int_t^\infty R(u)\,du}{R(t)}
=\frac{\dfrac{5}{0.2t+1}}{\dfrac{1}{(0.2t+1)^2}}
=5(0.2t+1)=t+5.
$$

---

## Mode

PDF: $f(t)=0.4(0.2t+1)^{-3},\ t\ge 0$.

Derivative:
$$
f'(t)=0.4(-3)(0.2t+1)^{-4}\cdot 0.2
= -0.24(0.2t+1)^{-4} < 0\ \ (t\ge 0).
$$

Hence $f(t)$ is strictly decreasing on $[0,\infty)$, maximum at boundary:
$$
\boxed{\operatorname{mode}=0,\ \ f(0)=0.4.}
$$

---

## 7) Hazard type

$$
z(t)=\frac{0.4}{0.2t+1}=0.4(0.2t+1)^{-1},\quad
z'(t)=0.4(-1)(0.2t+1)^{-2}\cdot 0.2
=-\frac{0.08}{(0.2t+1)^2}<0,
$$
so decreasing hazard rate (DHR).

---

## Summary

- $Z(t)=2\ln(0.2t+1)$
- $R(t)=1/(0.2t+1)^2$, $f(t)=0.4/(0.2t+1)^3$, $F(t)=1-1/(0.2t+1)^2$
- $\mathrm{MTTF}=5$
- $E[T^2]=\infty \Rightarrow \mathrm{Var}(T)=\infty,\ \sigma=\infty,\ \mathrm{CV}=\infty$
- $t_{0.5}=5(\sqrt{2}-1)\approx 2.0711$, $m(t)=t+5$, mode $=0$, DHR