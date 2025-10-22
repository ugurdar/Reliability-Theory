---
header-includes:
  - \usepackage{amsmath}
  - \usepackage{amssymb}
---

# Üstel Dağılım

## Reliability
$$ R(t) = Pr(X > t) = \int_t^\infty \frac{1}{\theta} e^{-x/\theta} dx = [-e^{-x/\theta}]_t^\infty = e^{-t/\theta} $$

## Arıza Hızı (Hazard)
$$ z(t) = \frac{f(t)}{R(t)} = \frac{\frac{1}{\theta} e^{-t/\theta}}{e^{-t/\theta}} = \frac{1}{\theta} \quad \text{(sabit)} $$

## Ortalama Arızaya Kadar Geçen Süre (MTTF)
$$ E[X] = \int_0^\infty x \frac{1}{\theta} e^{-x/\theta} dx $$

Parçalı integral ile:

$u = x$
$dv = \frac{1}{\theta} e^{-x/\theta} dx \implies v = -e^{-x/\theta}$
$$ E[X] = [-x e^{-x/\theta}]_0^\infty + \int_0^\infty e^{-x/\theta} dx = 0 + \theta = \theta $$

## İkinci Moment ve Varyans
Değişken değişimi: $y = x/\theta$, $x = \theta y$, $dx = \theta dy$

$$ E[X^2] = \int_0^\infty x^2 \frac{1}{\theta} e^{-x/\theta} dx = \theta^2 \int_0^\infty y^2 e^{-y} dy = \theta^2 \cdot 2! = 2\theta^2 $$

$$ Var(X) = E[X^2] - (E[X])^2 = 2\theta^2 - \theta^2 = \theta^2 $$

## Mod
$$ f'(x) = -\frac{1}{\theta^2} e^{-x/\theta} < 0 \quad (x > 0) \implies f \text{ tekdüze azalır} \implies t_{mod} = 0 $$

## Medyan
$$ \int_0^{t_{med}} \frac{1}{\theta} e^{-x/\theta} dx = [-e^{-x/\theta}]0^{t{med}} = 1 - e^{-t_{med}/\theta} = \frac{1}{2} $$

$$ e^{-t_{med}/\theta} = \frac{1}{2} \implies t_{med} = \theta \ln 2 $$

## Ortalama Kalan Ömür (MRL)
$$ MRL(x) = \frac{1}{R(x)} \int_x^\infty R(t) dt = \frac{1}{e^{-x/\theta}} \int_x^\infty e^{-t/\theta} dt $$

$$ = e^{x/\theta} \left[ -\theta e^{-t/\theta} \right]_x^\infty = e^{x/\theta} (0 + \theta e^{-x/\theta}) = \theta $$

### Özet
$$ R(t) = e^{-t/\theta}, \quad z(t) = \frac{1}{\theta}, \quad MTTF = \theta, \quad Var(X) = \theta^2, \quad t_{mod} = 0, \quad t_{med} = \theta \ln 2, \quad MRL(x) = \theta $$


# Gamma Dağılımı (şekil $\alpha$, ölçek $\beta$) 

**Olasılık yoğunluğu (pdf):**
$$
f(x)=\frac{1}{\Gamma(\alpha)\,\beta^{\alpha}}\,x^{\alpha-1} e^{-x/\beta},\qquad x\ge 0,\ \alpha>0,\ \beta>0.
$$

---

## Sağkalım (Reliability) $R(t)$
Tanım:
$$
R(t)=P(X>t)=\int_{t}^{\infty} f(x)\,dx.
$$
Pdf’yi yerine koyup değişken dönüşümü $y=x/\beta$ ($x=\beta y,\,dx=\beta\,dy$) ile:
$$
R(t)=\frac{1}{\Gamma(\alpha)\,\beta^{\alpha}}\int_{t}^{\infty} x^{\alpha-1}e^{-x/\beta}\,dx
= \frac{1}{\Gamma(\alpha)}\int_{t/\beta}^{\infty} y^{\alpha-1}e^{-y}\,dy
= \frac{\Gamma(\alpha,\,t/\beta)}{\Gamma(\alpha)}.
$$
Burada $\Gamma(\alpha, x)=\int_{x}^{\infty} y^{\alpha-1}e^{-y}\,dy$ **üst (upper) eksik gamma** fonksiyonudur.

## Arıza hızı (Hazard) $z(t)$
$$
z(t)=\frac{f(t)}{R(t)}=
\frac{\frac{1}{\Gamma(\alpha)\,\beta^{\alpha}}\,t^{\alpha-1}e^{-t/\beta}}{\frac{\Gamma(\alpha,\,t/\beta)}{\Gamma(\alpha)}}
= \frac{t^{\alpha-1}e^{-t/\beta}}{\beta^{\alpha}\,\Gamma(\alpha,\,t/\beta)}.
$$

## MTTF $=\mathbb E[X]$
$$
\mathbb E[X]=\int_{0}^{\infty} x\,f(x)\,dx
=\frac{1}{\Gamma(\alpha)\,\beta^{\alpha}}\int_{0}^{\infty} x^{\alpha}e^{-x/\beta}\,dx
$$
$y=x/\beta\Rightarrow x=\beta y,\,dx=\beta dy$:
$$
\mathbb E[X]=\frac{1}{\Gamma(\alpha)\,\beta^{\alpha}}\int_{0}^{\infty} (\beta y)^{\alpha}e^{-y}\,\beta\,dy
=\frac{\beta^{\alpha+1}}{\Gamma(\alpha)\,\beta^{\alpha}}\int_{0}^{\infty} y^{\alpha}e^{-y}\,dy
=\beta\frac{\Gamma(\alpha+1)}{\Gamma(\alpha)}=\alpha\,\beta.
$$

## İkinci moment ve **Varyans**
$$
\mathbb E[X^2]=\int_{0}^{\infty} x^2 f(x)\,dx
=\frac{1}{\Gamma(\alpha)\,\beta^{\alpha}}\int_{0}^{\infty} x^{\alpha+1}e^{-x/\beta}\,dx
\overset{y=x/\beta}{=}\frac{1}{\Gamma(\alpha)}\beta^{2}\int_{0}^{\infty} y^{\alpha+1}e^{-y}\,dy
=\beta^{2}\frac{\Gamma(\alpha+2)}{\Gamma(\alpha)}
=\beta^{2}\alpha(\alpha+1).
$$
$$
\operatorname{Var}(X)=\mathbb E[X^2]-(\mathbb E[X])^2
=\beta^{2}\alpha(\alpha+1)-(\alpha\beta)^2
=\alpha\,\beta^{2}.
$$

## Mod $t_{\text{mod}}$
$$
\frac{d}{dx}\log f(x)=\frac{d}{dx}[(\alpha-1)\ln x - x/\beta + C]
=\frac{\alpha-1}{x}-\frac{1}{\beta}.
$$
Sıfıra eşitleyince $x=(\alpha-1)\beta$. Dolayısıyla

$$
t_{\text{mod}}=
\begin{cases}
(\alpha-1)\beta, & \alpha>1,\\[2pt]
0, & \alpha\le 1\ \text{(tepe 0’da).}
\end{cases}
$$

## Medyan $t_{\text{med}}$
Medyan, $F(t_{\text{med}})=\tfrac12$ koşulundan gelir:
$$
F(t)=\int_{0}^{t} f(x)\,dx
=\frac{1}{\Gamma(\alpha)}\int_{0}^{t/\beta} y^{\alpha-1}e^{-y}\,dy
=\frac{\gamma(\alpha,\,t/\beta)}{\Gamma(\alpha)}.
$$
Burada $\gamma(\alpha,x)=\int_{0}^{x} y^{\alpha-1}e^{-y}\,dy$ **alt (lower) eksik gamma**. 
Kapalı form **genelde yoktur**; $t_{\text{med}}$ şu denklemin çözümüdür:
$$
\frac{\gamma(\alpha,\,t_{\text{med}}/\beta)}{\Gamma(\alpha)}=\tfrac12
\quad\Longleftrightarrow\quad
 t_{\text{med}}=\beta\cdot \gamma^{-1}\!\Big(\alpha,\ \tfrac12\,\Gamma(\alpha)\Big).
$$
(Özel durum: $\alpha=1$ \Rightarrow üstel; $t_{\text{med}}=\beta\ln 2$.)

## Ortalama kalan ömür (MRL) $\text{MRL}(x)$
Tanım:
$$
\text{MRL}(x)=\frac{\int_{x}^{\infty} R(t)\,dt}{R(x)}
=\mathbb E[X-x\mid X>x].
$$
Önce payı hesaplayalım:
$$
\int_{x}^{\infty} R(t)\,dt
=\int_{x}^{\infty}\Big(\int_{t}^{\infty} f(u)\,du\Big)dt
=\int_{x}^{\infty}(u-x)f(u)\,du
=\underbrace{\int_{x}^{\infty} u f(u)\,du}_{A}-x\,R(x).
$$
$A$ için:
$$
A=\frac{1}{\Gamma(\alpha)\,\beta^{\alpha}}\int_{x}^{\infty} u^{\alpha}e^{-u/\beta}\,du
=\frac{\beta}{\Gamma(\alpha)}\int_{x/\beta}^{\infty} y^{\alpha}e^{-y}\,dy
=\frac{\beta\,\Gamma(\alpha+1,\,x/\beta)}{\Gamma(\alpha)}.
$$
Dolayısıyla $R(x)=\dfrac{\Gamma(\alpha,\,x/\beta)}{\Gamma(\alpha)}$ ile birlikte
$$
\text{MRL}(x)=\frac{A-xR(x)}{R(x)}
=\beta\,\frac{\Gamma(\alpha+1,\,x/\beta)}{\Gamma(\alpha,\,x/\beta)}-x.
$$

---

# Weibull Dağılımı (şekil $k$, ölçek $\theta$) — $z(t)$ notasyonu ile

**Olasılık yoğunluğu (pdf):**
$$
f(x)=\frac{k}{\theta}\left(\frac{x}{\theta}\right)^{k-1} \exp\!\left[-\left(\frac{x}{\theta}\right)^{k}\right],\qquad x\ge 0,\ k>0,\ \theta>0.
$$

**Sağkalım (Reliability) $R(t)$:**
$$
R(t)=Pr(X>t)=\int_{t}^{\infty} f(x)\,dx
=\int_{t}^{\infty}\frac{k}{\theta}\left(\frac{x}{\theta}\right)^{k-1}\exp\!\left[-\left(\frac{x}{\theta}\right)^{k}\right]dx.
$$
Değişken dönüşümü $y=(x/\theta)^k\Rightarrow dy=(k/\theta)(x/\theta)^{k-1}dx$ ile
$$
R(t)=\int_{(t/\theta)^{k}}^{\infty} e^{-y}\,dy
=\Big[-e^{-y}\Big]_{(t/\theta)^{k}}^{\infty}
=e^{-(t/\theta)^{k}}.
$$

**Arıza hızı (Hazard) $z(t)$:**
$$
z(t)=\frac{f(t)}{R(t)}
=\frac{\tfrac{k}{\theta}\left(\tfrac{t}{\theta}\right)^{k-1}e^{-(t/\theta)^k}}{e^{-(t/\theta)^k}}
=\frac{k}{\theta}\left(\frac{t}{\theta}\right)^{k-1}.
$$

---

## MTTF $(=\mathbb E[X])$ — integral ile
$$
\mathbb E[X]=\int_{0}^{\infty} R(t)\,dt
=\int_{0}^{\infty} \exp\!\left[-\left(\frac{t}{\theta}\right)^{k}\right] dt.
$$
$y=(t/\theta)^k\Rightarrow dt=\frac{\theta}{k} y^{\frac{1}{k}-1}dy$:
$$
\mathbb E[X]=\frac{\theta}{k}\int_{0}^{\infty} y^{\frac{1}{k}-1} e^{-y}\,dy
=\frac{\theta}{k}\,\Gamma\!\left(\frac{1}{k}\right)
=\theta\,\Gamma\!\left(1+\frac{1}{k}\right),
$$
çünkü $\Gamma(1+u)=u\,\Gamma(u)$.

## İkinci moment ve **Varyans**
$$
\mathbb E[X^2]=\int_{0}^{\infty} x^2 f(x)\,dx
=\int_{0}^{\infty} x^2 \frac{k}{\theta}\left(\frac{x}{\theta}\right)^{k-1} e^{-(x/\theta)^k} dx.
$$
$y=(x/\theta)^k\Rightarrow x=\theta y^{1/k},\ dx=\frac{\theta}{k} y^{\frac{1}{k}-1}dy$:
$$
\mathbb E[X^2]=\int_{0}^{\infty} (\theta y^{1/k})^2 \frac{k}{\theta} y^{\frac{k-1}{k}} e^{-y}\,\frac{\theta}{k} y^{\frac{1}{k}-1}dy
=\theta^{2}\int_{0}^{\infty} y^{\frac{2}{k}} e^{-y}\,dy
=\theta^{2}\,\Gamma\!\left(1+\frac{2}{k}\right).
$$
Dolayısıyla
$$
\operatorname{Var}(X)=\mathbb E[X^2]-(\mathbb E[X])^{2}
=\theta^{2}\Big[\Gamma\!\left(1+\frac{2}{k}\right)-\Gamma\!\left(1+\frac{1}{k}\right)^{2}\Big].
$$

## Mod $t_{\text{mod}}$
$$
\frac{d}{dt}\log f(t)=\frac{d}{dt}\Big[\log k-\log\theta+(k-1)\log t-(k-1)\log\theta-\left(\frac{t}{\theta}\right)^{k}\Big]
=\frac{k-1}{t}-\frac{k}{\theta}\left(\frac{t}{\theta}\right)^{k-1}.
$$
Sıfıra eşitleyince $(k-1)/t = (k/\theta)(t/\theta)^{k-1} \Rightarrow (t/\theta)^{k}=(k-1)/k$.  

$$
t_{\text{mod}}=
\begin{cases}
\theta\left(\dfrac{k-1}{k}\right)^{1/k}, & k>1,\\[6pt]
0, & k\le 1\ \text{(tepe 0’da).}
\end{cases}
$$

## Medyan $t_{\text{med}}$ — integralden
$$
\int_{0}^{t_{\text{med}}} f(x)\,dx=\frac{1}{2}
\quad\Longleftrightarrow\quad
1-e^{-(t_{\text{med}}/\theta)^{k}}=\frac{1}{2}
\ \Rightarrow\ t_{\text{med}}=\theta(\ln 2)^{1/k}.
$$

## Ortalama kalan ömür (MRL) $\text{MRL}(x)$
Tanım (integralle):
$$
\text{MRL}(x)=\frac{\int_{x}^{\infty} R(t)\,dt}{R(x)}
=\frac{\int_{x}^{\infty} e^{-(t/\theta)^{k}}\,dt}{e^{-(x/\theta)^{k}}}.
$$
$y=(t/\theta)^k\Rightarrow dt=\frac{\theta}{k} y^{\frac{1}{k}-1}dy$ ve alt sınır $a=(x/\theta)^k$:
$$
\int_{x}^{\infty} e^{-(t/\theta)^k} dt
=\frac{\theta}{k}\int_{a}^{\infty} y^{\frac{1}{k}-1} e^{-y}\,dy
=\frac{\theta}{k}\,\Gamma\!\left(\frac{1}{k},\,a\right),
$$
burada $\Gamma(s,a)=\int_{a}^{\infty} y^{s-1}e^{-y}dy$ **üst eksik gamma**. Dolayısıyla
$$
\text{MRL}(x)=\frac{\frac{\theta}{k}\Gamma(1/k,\,a)}{e^{-a}}
=\frac{\theta}{k}e^{a}\Gamma\!\left(\frac{1}{k},\,a\right)
=\theta\,e^{a}\Gamma\!\left(1+\frac{1}{k},\,a\right)-x,
$$
son eşitlik $\Gamma(s+1,a)=s\,\Gamma(s,a)+a^{s}e^{-a}$ özdeşliğinden gelir.

---

## Parametrizasyon notu
- Weibull’un iki parametresi vardır: şekil $k$ ve ölçek $\theta$ (bazı kaynaklarda $\lambda$ veya $\eta$).
- Üstel dağılımda olduğu gibi “ölçek $\theta$” yerine hız (rate) $r=1/\theta$ yazılabilir, fakat Weibull’da bu daha az kullanılır. Bu durumda
  $$
  f(x)=k\,r\,(r x)^{k-1}\,e^{-(r x)^{k}},\quad
  R(t)=e^{-(r t)^{k}},\quad
  z(t)=k\,r\,(r t)^{k-1}.
  $$
- Özel durum: $k=1$ olursa Weibull $\text{Exp}(\theta)$’ye iner; $z(t)=1/\theta$ sabit olur.

---

### Özet
- $R(t)=e^{-(t/\theta)^{k}}$, $z(t)=\dfrac{k}{\theta}\left(\dfrac{t}{\theta}\right)^{k-1}$.
- $\text{MTTF}=\theta\,\Gamma\!\left(1+\dfrac{1}{k}\right)$, $\operatorname{Var}(X)=\theta^{2}\!\left[\Gamma\!\left(1+\dfrac{2}{k}\right)-\Gamma\!\left(1+\dfrac{1}{k}\right)^{2}\right]$.
- $t_{\text{mod}}=\theta\left(\dfrac{k-1}{k}\right)^{1/k}$ ($k>1$; yoksa $0$).
- $t_{\text{med}}=\theta(\ln 2)^{1/k}$.
- $\text{MRL}(x)=\dfrac{\theta}{k}e^{(x/\theta)^{k}}\Gamma\!\left(\dfrac{1}{k},\,(x/\theta)^{k}\right)\ =\ \theta\,e^{(x/\theta)^{k}}\Gamma\!\left(1+\dfrac{1}{k},\,(x/\theta)^{k}\right)-x$.

---
