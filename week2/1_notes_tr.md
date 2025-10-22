## Üstel Dağılım
Parametre: $\theta > 0$

Olasılık yoğunluk fonksiyonu (PDF):
$$
f(t) = \frac{1}{\theta} e^{-t/\theta}, \quad t \geq 0
$$

Kümülatif dağılım fonksiyonu (CDF):
$$
F(t) = 1 - e^{-t/\theta}
$$

**MTTF (Ortalama Arızaya Kadar Geçen Süre) aşağıdaki şekilde tanımlanır:**
$$
\mathrm{MTTF} = E(T) = \int_0^\infty t\,f(t)\,dt
$$

$f(t) = -R'(t)$ olduğundan:
$$
E(T) = \int_0^\infty t\,(-R'(t))\,dt = -\int_0^\infty t\,R'(t)\,dt
$$

Parçalı integrasyon (integration by parts) uygula:
- $u = t \implies du = dt$
- $dv = R'(t)\,dt \implies v = R(t)$

Parçalı integrasyon formülü:
$$
\int u\,dv = u\,v - \int v\,du
$$

Buna göre:
$$
-\int_0^\infty t\,R'(t)\,dt = -\left[ t\,R(t) \right]_0^\infty + \int_0^\infty R(t)\,dt
$$

Sınır terimi:
- $t \to \infty$ iken $R(t) \to 0$ ve $t\,R(t) \to 0$
- $t = 0$ iken $t\,R(0) = 0$

Sonuç olarak:
$$
-\left[ t\,R(t) \right]_0^\infty = 0
$$

Yani:
$$
E(T) = \int_0^\infty R(t)\,dt
$$

Bu, ortalama arıza süresinin güvenilirlik fonksiyonunun altında kalan alan olduğunu gösterir.

---
