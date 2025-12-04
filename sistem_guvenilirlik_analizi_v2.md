# Seri, Paralel ve Kompleks Sistemlerin Güvenilirlik Analizi

## Kapsamlı Teorik ve Pratik Rehber

---

# KISIM I: OLASILIK TEORİSİ TEMELLERİ

## 1. Olasılık Aksiyomları (Kolmogorov)

Bir $\Omega$ örneklem uzayı ve $\mathcal{F}$ sigma-cebiri üzerinde tanımlı $P$ olasılık ölçüsü için:

1. **Negatif olmama**: $P(A) \geq 0$, her $A \in \mathcal{F}$ için
2. **Normalleştirme**: $P(\Omega) = 1$
3. **Sayılabilir toplamsal**: Ayrık $A_1, A_2, \ldots$ olayları için:
   $$P\left(\bigcup_{i=1}^{\infty} A_i\right) = \sum_{i=1}^{\infty} P(A_i)$$

## 2. Temel Olasılık Kuralları

### Tümleyen Kuralı
$$P(A') = 1 - P(A)$$

### Toplama Kuralı (Genel)
$$P(A \cup B) = P(A) + P(B) - P(A \cap B)$$

### Koşullu Olasılık
$$P(A|B) = \frac{P(A \cap B)}{P(B)}, \quad P(B) > 0$$

### Çarpım Kuralı (Genel)
$$P(A \cap B) = P(A|B) \cdot P(B) = P(B|A) \cdot P(A)$$

### Toplam Olasılık Teoremi
$B_1, B_2, \ldots, B_n$ tam bir bölüntü ise:
$$P(A) = \sum_{i=1}^{n} P(A|B_i) \cdot P(B_i)$$

### Bayes Teoremi
$$P(B_i|A) = \frac{P(A|B_i) \cdot P(B_i)}{\sum_{j=1}^{n} P(A|B_j) \cdot P(B_j)}$$

## 3. Bağımsızlık Kavramı

### İki Olayın Bağımsızlığı

$A$ ve $B$ olayları **bağımsızdır** ancak ve ancak:
$$P(A \cap B) = P(A) \cdot P(B)$$

Eşdeğer ifadeler:
- $P(A|B) = P(A)$ (B'nin olması A'yı etkilemez)
- $P(B|A) = P(B)$ (A'nın olması B'yi etkilemez)

### n Olayın Karşılıklı Bağımsızlığı

$A_1, A_2, \ldots, A_n$ olayları **karşılıklı bağımsızdır** ancak ve ancak her $k \leq n$ ve her $\{i_1, i_2, \ldots, i_k\} \subseteq \{1, 2, \ldots, n\}$ alt kümesi için:

$$P(A_{i_1} \cap A_{i_2} \cap \cdots \cap A_{i_k}) = P(A_{i_1}) \cdot P(A_{i_2}) \cdots P(A_{i_k})$$

**Önemli**: İkişerli bağımsızlık, karşılıklı bağımsızlığı garanti ETMEZ!

### Örnek: İkişerli Bağımsız ama Karşılıklı Bağımsız Değil

İki adil zar atılsın:
- $A$: İlk zar çift
- $B$: İkinci zar çift  
- $C$: Toplam çift

$P(A) = P(B) = P(C) = 1/2$

İkişerli kontrol:
- $P(A \cap B) = 1/4 = P(A) \cdot P(B)$ ✓
- $P(A \cap C) = 1/4 = P(A) \cdot P(C)$ ✓
- $P(B \cap C) = 1/4 = P(B) \cdot P(C)$ ✓

Üçlü kontrol:
- $P(A \cap B \cap C) = 1/4$ (her ikisi çift ise toplam çift)
- $P(A) \cdot P(B) \cdot P(C) = 1/8$
- $1/4 \neq 1/8$ ✗

**Sonuç**: İkişerli bağımsız ama karşılıklı bağımsız değil!

---

# KISIM II: GÜVENİLİRLİK TEORİSİ TEMELLERİ

## 1. Temel Tanımlar

### Güvenilirlik Fonksiyonu (Survival Function)

$$R(t) = P(T > t) = 1 - F(t)$$

- $T$: Arızaya kadar geçen süre (rastgele değişken)
- $F(t)$: Birikimli dağılım fonksiyonu (CDF)
- $R(t)$: Sağkalım fonksiyonu

**Özellikler**:
- $R(0) = 1$ (başlangıçta sistem çalışıyor)
- $R(\infty) = 0$ (sonunda her sistem arızalanır)
- $R(t)$ azalan bir fonksiyondur

### Güvensizlik Fonksiyonu (Unreliability)

$$Q(t) = P(T \leq t) = F(t) = 1 - R(t)$$

### Olasılık Yoğunluk Fonksiyonu (PDF)

$$f(t) = \frac{dF(t)}{dt} = -\frac{dR(t)}{dt}$$

### Tehlike Oranı Fonksiyonu (Hazard Rate)

$$Z(t) = \frac{f(t)}{R(t)} = -\frac{d}{dt}\ln R(t)$$

**Fiziksel yorum**: $Z(t) \cdot \Delta t \approx P(t < T \leq t + \Delta t | T > t)$

Yani, $t$ anına kadar sağ kalmış bir sistemin, sonraki $\Delta t$ sürede arızalanma olasılığı.

### Birikimli Tehlike Fonksiyonu

$$H(t) = \int_0^t Z(u) \, du = -\ln R(t)$$

$$R(t) = e^{-H(t)} = \exp\left(-\int_0^t Z(u) \, du\right)$$

### Ortalama Arızaya Kadar Süre (MTTF)

$$MTTF = E[T] = \int_0^{\infty} t \cdot f(t) \, dt = \int_0^{\infty} R(t) \, dt$$

## 2. Üstel Dağılım (Sabit Tehlike Oranı)

Parametrelendirme: $\theta$ = ortalama ömür (MTTF)

### Fonksiyonlar

| Fonksiyon | Formül |
|-----------|--------|
| PDF | $f(t) = \frac{1}{\theta} e^{-t/\theta}$ |
| CDF | $F(t) = 1 - e^{-t/\theta}$ |
| Güvenilirlik | $R(t) = e^{-t/\theta}$ |
| Tehlike Oranı | $Z(t) = \frac{1}{\theta}$ (SABİT) |
| MTTF | $E[T] = \theta$ |
| Varyans | $Var(T) = \theta^2$ |

### Hafızasızlık Özelliği

Üstel dağılımın en önemli özelliği **hafızasızlık**tır:

$$P(T > t + s | T > t) = P(T > s)$$

**İspat**:
$$P(T > t + s | T > t) = \frac{P(T > t + s)}{P(T > t)} = \frac{e^{-(t+s)/\theta}}{e^{-t/\theta}} = e^{-s/\theta} = P(T > s)$$

**Yorum**: Sistem $t$ süre çalıştıysa, bundan sonra ne kadar çalışacağı $t$'den bağımsızdır. Sistem "eskimez" - yeni gibidir!

### Örnek: Elektronik Bileşen

Bir elektronik bileşenin ortalama ömrü $\theta = 10000$ saat.

**a)** 1000 saat çalışma olasılığı:
$$R(1000) = e^{-1000/10000} = e^{-0.1} = 0.9048$$

**b)** 5000 saat çalıştıktan sonra 1000 saat daha çalışma olasılığı:
$$P(T > 6000 | T > 5000) = P(T > 1000) = 0.9048$$

(Hafızasızlık nedeniyle aynı!)

**c)** İlk 1000 saatte arızalanma olasılığı:
$$Q(1000) = 1 - R(1000) = 1 - 0.9048 = 0.0952$$

## 3. Weibull Dağılımı (Değişken Tehlike Oranı)

3-parametreli Weibull: $T \sim Weibull(\beta, \theta, \gamma)$

- $\beta$: Şekil parametresi (shape)
- $\theta$: Ölçek parametresi (scale)
- $\gamma$: Konum parametresi (location/threshold)

### Fonksiyonlar

$$R(t) = \exp\left[-\left(\frac{t-\gamma}{\theta}\right)^\beta\right], \quad t > \gamma$$

$$Z(t) = \frac{\beta}{\theta}\left(\frac{t-\gamma}{\theta}\right)^{\beta-1}$$

### Tehlike Oranı Davranışı

| $\beta$ değeri | $Z(t)$ davranışı | Yorum |
|----------------|------------------|-------|
| $\beta < 1$ | Azalan | Erken arızalar (infant mortality) |
| $\beta = 1$ | Sabit | Üstel dağılım (rastgele arızalar) |
| $\beta > 1$ | Artan | Yaşlanma/aşınma (wear-out) |

### Küvet Eğrisi (Bathtub Curve)

Gerçek sistemler genellikle üç fazdan geçer:

```
Z(t)
  │
  │╲                              ╱
  │ ╲                            ╱
  │  ╲                          ╱
  │   ╲________________________╱
  │    
  └──────────────────────────────── t
     Faz I      Faz II      Faz III
   (Erken)   (Rastgele)  (Yaşlanma)
   β < 1      β = 1        β > 1
```

---

# KISIM III: SERİ SİSTEMLER

## 1. Tanım ve Yapı

Seri sistemde **tüm bileşenlerin çalışması** gerekir. Tek bir bileşenin arızası sistemi durdurur.

### Blok Diyagramı

```
    ┌───┐   ┌───┐   ┌───┐       ┌───┐
  ──┤ 1 ├───┤ 2 ├───┤ 3 ├─ ... ─┤ n ├──
    └───┘   └───┘   └───┘       └───┘
```

### Yapı Fonksiyonu

$$\phi(\mathbf{x}) = \min(x_1, x_2, \ldots, x_n) = \prod_{i=1}^{n} x_i$$

Burada $x_i \in \{0, 1\}$ (0: arızalı, 1: çalışıyor)

## 2. Bağımsız Bileşenler Durumu

### Temel Formül

Bileşenler bağımsız ise, kesişim olasılığı çarpıma eşittir:

$$R_{sys}(t) = P(T_1 > t \cap T_2 > t \cap \cdots \cap T_n > t)$$

Bağımsızlık varsayımı ile:

$$R_{sys}(t) = P(T_1 > t) \cdot P(T_2 > t) \cdots P(T_n > t) = \prod_{i=1}^{n} R_i(t)$$

### Tehlike Oranı (Bağımsız Bileşenler)

$$Z_{sys}(t) = -\frac{d}{dt}\ln R_{sys}(t) = -\frac{d}{dt}\sum_{i=1}^{n}\ln R_i(t) = \sum_{i=1}^{n} Z_i(t)$$

**Önemli sonuç**: Seri sistemde tehlike oranları TOPLANIR!

### Üstel Dağılım için

Her bileşen $\theta_i$ ortalama ömre sahipse:

$$R_{sys}(t) = \prod_{i=1}^{n} e^{-t/\theta_i} = \exp\left(-t\sum_{i=1}^{n}\frac{1}{\theta_i}\right)$$

Sistem de üstel dağılımlıdır:
$$\theta_{sys} = \frac{1}{\sum_{i=1}^{n} 1/\theta_i}$$

$$MTTF_{sys} = \theta_{sys} = \left(\sum_{i=1}^{n} \frac{1}{\theta_i}\right)^{-1}$$

### Örnek 1: Basit Seri Sistem

Üç bileşenli seri sistem: $R_1 = 0.95$, $R_2 = 0.90$, $R_3 = 0.85$

$$R_{sys} = 0.95 \times 0.90 \times 0.85 = 0.7268$$

**Yorum**: En düşük güvenilirlik 0.85 iken, sistem güvenilirliği 0.7268'e düştü!

### Örnek 2: n Özdeş Bileşen

Her bileşenin güvenilirliği $R = 0.99$:

| n | $R_{sys} = R^n$ | Güvenilirlik Kaybı |
|---|-----------------|-------------------|
| 1 | 0.9900 | - |
| 5 | 0.9510 | %4.9 |
| 10 | 0.9044 | %9.5 |
| 50 | 0.6050 | %39.5 |
| 100 | 0.3660 | %63.4 |
| 500 | 0.0066 | %99.3 |

### Örnek 3: Üstel Dağılım

Üç bileşen: $\theta_1 = 1000$ saat, $\theta_2 = 2000$ saat, $\theta_3 = 5000$ saat

$$\frac{1}{\theta_{sys}} = \frac{1}{1000} + \frac{1}{2000} + \frac{1}{5000} = 0.001 + 0.0005 + 0.0002 = 0.0017$$

$$\theta_{sys} = \frac{1}{0.0017} = 588.2 \text{ saat}$$

$t = 100$ saat için:
$$R_{sys}(100) = e^{-100/588.2} = e^{-0.17} = 0.8437$$

### Örnek 4: Seri Sisteme Bileşen Eklemenin Etkisi

Mevcut sistem: $R_{sys} = 0.95$

Yeni bileşen eklenirse ($R_{new} = 0.98$):
$$R_{sys,new} = 0.95 \times 0.98 = 0.931$$

**Sonuç**: Her yeni bileşen güvenilirliği AZALTIR!

## 3. Bağımlı Bileşenler Durumu

### Neden Bağımlılık Oluşur?

1. **Ortak stres**: Aynı ortam koşulları (sıcaklık, titreşim, nem)
2. **Ortak neden arızaları (CCF)**: Tek bir olay birden fazla bileşeni etkiler
3. **Yük paylaşımı**: Bir bileşen arızalandığında diğerlerine yük biner
4. **Kaskad arızalar**: Bir arıza diğerlerini tetikler
5. **Ortak tedarikçi**: Aynı üreticiden gelen kusurlu parti

### Matematiksel İfade

Bağımlı durumda:
$$P(A \cap B) \neq P(A) \cdot P(B)$$

Bunun yerine:
$$P(A \cap B) = P(A) \cdot P(B|A)$$

### Korelasyon ile Modelleme

İki bileşen arasındaki korelasyon $\rho$ ile:

$$P(A \cap B) = P(A) \cdot P(B) + \rho \cdot \sqrt{P(A)(1-P(A))P(B)(1-P(B))}$$

### Örnek 5: Pozitif Korelasyonlu Seri Sistem

İki bileşen: $R_1 = R_2 = 0.90$, korelasyon $\rho = 0.3$

**Bağımsız varsayım**:
$$R_{sys,ind} = 0.90 \times 0.90 = 0.81$$

**Bağımlı hesap**:
$Q_1 = Q_2 = 0.10$

Arıza olasılıkları için:
$$P(Q_1 \cap Q_2) = Q_1 \cdot Q_2 + \rho \cdot \sqrt{Q_1(1-Q_1)Q_2(1-Q_2)}$$
$$= 0.01 + 0.3 \times \sqrt{0.09 \times 0.09} = 0.01 + 0.3 \times 0.09 = 0.037$$

Sistem güvensizliği (inclusion-exclusion):
$$Q_{sys} = Q_1 + Q_2 - P(Q_1 \cap Q_2) = 0.10 + 0.10 - 0.037 = 0.163$$

$$R_{sys,dep} = 1 - 0.163 = 0.837$$

**Karşılaştırma**: $R_{sys,dep} = 0.837 > R_{sys,ind} = 0.81$

**Yorum**: Pozitif korelasyon seri sistemde güvenilirliği ARTIRIR! Çünkü bileşenler birlikte arızalanma eğiliminde, yani biri çalışıyorsa diğeri de muhtemelen çalışıyordur.

### Ortak Neden Arızaları (CCF) - Beta Faktör Modeli

Toplam arıza oranı iki kısma ayrılır:
- Bağımsız arızalar: $(1-\beta) \cdot Z_i(t)$
- Ortak neden arızaları: $\beta \cdot Z_{CCF}(t)$

$\beta$ tipik olarak 0.01 - 0.10 arasındadır.

### Örnek 6: CCF ile Seri Sistem

İki özdeş bileşen: $\theta = 10000$ saat, $\beta = 0.05$

Bağımsız arıza oranı: $(1-0.05)/10000 = 0.000095$ /saat
CCF oranı: $0.05/10000 = 0.000005$ /saat

Sistem arıza oranı:
$$Z_{sys} = 2 \times 0.000095 + 0.000005 = 0.000195 \text{ /saat}$$

(Bağımsız varsayımda $Z_{sys} = 2/10000 = 0.0002$ /saat olurdu)

---

# KISIM IV: PARALEL SİSTEMLER

## 1. Tanım ve Yapı

Paralel sistemde **en az bir bileşenin çalışması** yeterlidir.

### Blok Diyagramı

```
       ┌───┐
    ┌──┤ 1 ├──┐
    │  └───┘  │
    │  ┌───┐  │
  ──┼──┤ 2 ├──┼──
    │  └───┘  │
    │  ┌───┐  │
    └──┤ n ├──┘
       └───┘
```

### Yapı Fonksiyonu

$$\phi(\mathbf{x}) = \max(x_1, x_2, \ldots, x_n) = 1 - \prod_{i=1}^{n}(1 - x_i)$$

## 2. Bağımsız Bileşenler Durumu

### Temel Formül

Sistem ancak TÜM bileşenler arızalandığında arızalanır:

$$Q_{sys}(t) = P(T_1 \leq t \cap T_2 \leq t \cap \cdots \cap T_n \leq t)$$

Bağımsızlık ile:
$$Q_{sys}(t) = \prod_{i=1}^{n} Q_i(t) = \prod_{i=1}^{n} [1 - R_i(t)]$$

$$R_{sys}(t) = 1 - \prod_{i=1}^{n} [1 - R_i(t)]$$

### Örnek 7: Basit Paralel Sistem

Üç bileşenli paralel sistem: $R_1 = 0.90$, $R_2 = 0.85$, $R_3 = 0.80$

$$Q_{sys} = (1-0.90)(1-0.85)(1-0.80) = 0.10 \times 0.15 \times 0.20 = 0.003$$

$$R_{sys} = 1 - 0.003 = 0.997$$

**Yorum**: En yüksek güvenilirlik 0.90 iken, sistem 0.997'ye çıktı!

### Örnek 8: n Özdeş Bileşen

Her bileşenin güvenilirliği $R = 0.80$:

| n | $R_{sys} = 1-(1-R)^n$ | İyileşme |
|---|----------------------|----------|
| 1 | 0.8000 | - |
| 2 | 0.9600 | +20% |
| 3 | 0.9920 | +24% |
| 4 | 0.9984 | +25% |
| 5 | 0.9997 | +25% |

### Tehlike Oranı (Paralel Sistem)

Üstel bileşenler için bile sistem tehlike oranı sabit DEĞİLDİR:

$$Z_{sys}(t) = \frac{\sum_{i=1}^{n} Z_i(t) \cdot \prod_{j \neq i} Q_j(t)}{\prod_{i=1}^{n} Q_i(t) - \prod_{i=1}^{n} Q_i(t) + R_{sys}(t)}$$

Bu karmaşık ifade, paralel sistemlerin analizinin seri sistemlerden daha zor olduğunu gösterir.

### MTTF (Özdeş Üstel Bileşenler)

$$MTTF_{sys} = \theta \sum_{i=1}^{n} \frac{1}{i} = \theta \left(1 + \frac{1}{2} + \frac{1}{3} + \cdots + \frac{1}{n}\right)$$

### Örnek 9: MTTF Karşılaştırması

$\theta = 1000$ saat olan bileşenler:

| n | MTTF (saat) | Çarpan |
|---|-------------|--------|
| 1 | 1000 | 1.00 |
| 2 | 1500 | 1.50 |
| 3 | 1833 | 1.83 |
| 4 | 2083 | 2.08 |
| 5 | 2283 | 2.28 |
| 10 | 2929 | 2.93 |

## 3. Bağımlı Bileşenler Durumu

### Pozitif Korelasyonun Etkisi

Paralel sistemde pozitif korelasyon güvenilirliği AZALTIR!

### Örnek 10: Pozitif Korelasyonlu Paralel Sistem

İki bileşen: $R_1 = R_2 = 0.90$, korelasyon $\rho = 0.3$

**Bağımsız varsayım**:
$$R_{sys,ind} = 1 - (0.10)(0.10) = 0.99$$

**Bağımlı hesap**:
$$P(Q_1 \cap Q_2) = 0.01 + 0.3 \times 0.09 = 0.037$$

$$R_{sys,dep} = 1 - 0.037 = 0.963$$

**Karşılaştırma**: $R_{sys,dep} = 0.963 < R_{sys,ind} = 0.99$

**Yorum**: Pozitif korelasyon paralel sistemde güvenilirliği AZALTIR! Çünkü biri arızalandığında diğeri de muhtemelen arızalanır - yedeklilik amacı boşa çıkar.

### Korelasyonun Seri vs Paralel Etkisi

| Korelasyon | Seri Sistem | Paralel Sistem |
|------------|-------------|----------------|
| $\rho > 0$ (pozitif) | Güvenilirlik ARTAR | Güvenilirlik AZALIR |
| $\rho = 0$ (bağımsız) | Temel formül | Temel formül |
| $\rho < 0$ (negatif) | Güvenilirlik AZALIR | Güvenilirlik ARTAR |

### Sezgisel Açıklama

**Seri sistemde** (hepsinin çalışması lazım):
- Pozitif korelasyon: "Biri çalışıyorsa hepsi çalışıyor" → İYİ
- Negatif korelasyon: "Biri çalışıyorsa diğeri arızalı" → KÖTÜ

**Paralel sistemde** (birinin çalışması yeterli):
- Pozitif korelasyon: "Biri arızalıysa hepsi arızalı" → KÖTÜ
- Negatif korelasyon: "Biri arızalıysa diğeri çalışıyor" → İYİ

### Örnek 11: CCF ile Paralel Sistem

İki özdeş bileşen: $\theta = 10000$ saat, $\beta = 0.05$ (CCF faktörü)

**Bağımsız model**:
$$R_{sys,ind}(t) = 1 - (1 - e^{-t/10000})^2$$

$t = 1000$ saat için: $R_{sys,ind} = 0.9909$

**CCF modeli** (basitleştirilmiş):
Ortak neden arızası olasılığı: $Q_{CCF} = 1 - e^{-\beta t/\theta} = 1 - e^{-0.05 \times 1000/10000} = 0.005$

$$R_{sys,CCF} \approx R_{sys,ind} \times (1 - Q_{CCF}) = 0.9909 \times 0.995 = 0.986$$

CCF, paralel sistemin avantajını önemli ölçüde azaltır!

---

# KISIM V: k-out-of-n SİSTEMLER

## 1. Tanım

n bileşenden **en az k tanesinin çalışması** gereken sistem.

Gösterim: k/n veya k-out-of-n:G (G = Good, çalışan)

### Özel Durumlar

| Sistem | k değeri | Eşdeğer |
|--------|----------|---------|
| Paralel | k = 1 | 1-out-of-n |
| Seri | k = n | n-out-of-n |
| Oylama | k = ⌈n/2⌉ | Çoğunluk |

## 2. Bağımsız Özdeş Bileşenler

### Binom Dağılımı ile Hesap

$$R_{sys} = P(X \geq k) = \sum_{j=k}^{n} \binom{n}{j} R^j (1-R)^{n-j}$$

Burada $X$ = çalışan bileşen sayısı, $X \sim Binomial(n, R)$

### Örnek 12: 2-out-of-3 Sistem

$R = 0.90$ için:

$$R_{sys} = \binom{3}{2}(0.9)^2(0.1)^1 + \binom{3}{3}(0.9)^3(0.1)^0$$
$$= 3 \times 0.81 \times 0.1 + 1 \times 0.729 \times 1$$
$$= 0.243 + 0.729 = 0.972$$

### Örnek 13: Uçak Motor Sistemi

4 motorlu uçak, farklı gereksinimler için:

| Gereksinim | k | $R_{sys}$ (R=0.95) |
|------------|---|-------------------|
| En az 1 motor | 1 | 0.99999 |
| En az 2 motor | 2 | 0.99945 |
| En az 3 motor | 3 | 0.98598 |
| Tüm motorlar | 4 | 0.81451 |

### Örnek 14: Oylama Sistemi (TMR - Triple Modular Redundancy)

3 bilgisayar, çoğunluk oyu ile karar (2-out-of-3):

$R_{computer} = 0.99$ için:
$$R_{TMR} = 3(0.99)^2(0.01) + (0.99)^3 = 0.0297 + 0.9703 = 0.999973$$

Tek bilgisayara göre iyileşme: $0.999973 / 0.99 = 1.0101$ (küçük ama önemli!)

## 3. Farklı Güvenilirlikli Bileşenler

Bileşenler özdeş değilse, durum uzayı (state enumeration) metodu kullanılır.

### Örnek 15: 2-out-of-3, Farklı Güvenilirlikler

$R_1 = 0.95$, $R_2 = 0.90$, $R_3 = 0.85$

Başarılı durumlar (en az 2 çalışıyor):

| Durum | 1 | 2 | 3 | Olasılık |
|-------|---|---|---|----------|
| 111 | ✓ | ✓ | ✓ | $0.95 \times 0.90 \times 0.85 = 0.72675$ |
| 110 | ✓ | ✓ | ✗ | $0.95 \times 0.90 \times 0.15 = 0.12825$ |
| 101 | ✓ | ✗ | ✓ | $0.95 \times 0.10 \times 0.85 = 0.08075$ |
| 011 | ✗ | ✓ | ✓ | $0.05 \times 0.90 \times 0.85 = 0.03825$ |

$$R_{sys} = 0.72675 + 0.12825 + 0.08075 + 0.03825 = 0.974$$

## 4. Bağımlı Bileşenlerle k-out-of-n

### Örnek 16: Pozitif Korelasyon Etkisi

2-out-of-3 sistem, $R_i = 0.90$, çiftler arası $\rho = 0.2$

Bu durumda multinomial/copula modelleri gerekir. Basitleştirilmiş yaklaşım:

Bağımsız: $R_{sys,ind} = 0.972$

Pozitif korelasyon genellikle k-out-of-n sistemlerde güvenilirliği AZALTIR (paralel sisteme benzer etki).

---

# KISIM VI: SERİ-PARALEL KOMBİNE SİSTEMLER

## 1. Çözüm Stratejisi

1. Sistemi içten dışa doğru analiz et
2. Önce en iç paralel/seri grupları tek bileşene indirge
3. Adım adım dışa doğru ilerle

## 2. Seri-Paralel Yapı

```
       ┌───┐
    ┌──┤ 1 ├──┐
    │  └───┘  │     ┌───┐
  ──┤         ├─────┤ 3 ├──
    │  ┌───┐  │     └───┘
    └──┤ 2 ├──┘
       └───┘
    [Paralel]      [Seri]
```

### Örnek 17

$R_1 = 0.90$, $R_2 = 0.85$, $R_3 = 0.95$

**Adım 1**: Paralel grup (1,2):
$$R_{12} = 1 - (1-0.90)(1-0.85) = 1 - 0.015 = 0.985$$

**Adım 2**: Seri bağlantı:
$$R_{sys} = R_{12} \times R_3 = 0.985 \times 0.95 = 0.936$$

## 3. Paralel-Seri Yapı

```
    ┌───┐   ┌───┐
 ┌──┤ 1 ├───┤ 2 ├──┐
 │  └───┘   └───┘  │
─┤                 ├─
 │  ┌───┐   ┌───┐  │
 └──┤ 3 ├───┤ 4 ├──┘
    └───┘   └───┘
```

### Örnek 18

Tüm bileşenler $R = 0.90$

**Adım 1**: Seri gruplar:
- Üst dal: $R_{12} = 0.90 \times 0.90 = 0.81$
- Alt dal: $R_{34} = 0.90 \times 0.90 = 0.81$

**Adım 2**: Paralel bağlantı:
$$R_{sys} = 1 - (1-0.81)(1-0.81) = 1 - 0.0361 = 0.964$$

## 4. Çok Aşamalı Sistem

### Örnek 19: Beş Aşamalı Üretim Hattı

```
    ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐
    │   │   │ B │   │ C │   │ E │   │   │
 ───┤ A ├───┤ B ├───┤ C ├───┤ E ├───┤ G ├───
    │   │   │ B │   │ C │   │ E │   │   │
    └───┘   └───┘   └───┘   └───┘   └───┘
    n=1     n=2     n=3     n=2     n=1
```

$R_A = 0.95$, $R_B = 0.85$, $R_C = 0.80$, $R_E = 0.90$, $R_G = 0.95$

**Çözüm**:
- Aşama A: $R_A = 0.95$
- Aşama B (2 paralel): $R_{BB} = 1 - (0.15)^2 = 0.9775$
- Aşama C (3 paralel): $R_{CCC} = 1 - (0.20)^3 = 0.992$
- Aşama E (2 paralel): $R_{EE} = 1 - (0.10)^2 = 0.99$
- Aşama G: $R_G = 0.95$

$$R_{sys} = 0.95 \times 0.9775 \times 0.992 \times 0.99 \times 0.95 = 0.866$$

---

# KISIM VII: KOMPLEKS SİSTEMLER

## 1. Tanım

Basit seri-paralel indirgeme ile çözülemeyen sistemler **kompleks sistem** olarak adlandırılır.

### Klasik Örnek: Köprü (Bridge) Yapısı

```
         1
    A●───────●C
     │╲     ╱│
     │ ╲   ╱ │
   2 │  ╲5╱  │ 4
     │   ╳   │
     │  ╱ ╲  │
     │ ╱   ╲ │
    B●───────●D
         3
```

Bu yapı, bileşen 5 nedeniyle basit seri-paralel değildir.

## 2. Çözüm Yöntemleri

### Yöntem 1: Decomposition (Ayrıştırma)

Toplam olasılık teoremine dayanır:

$$R_{sys} = R_k \cdot R_{sys|k \text{ çalışıyor}} + Q_k \cdot R_{sys|k \text{ arızalı}}$$

### Örnek 20: Köprü Sistemi (Decomposition)

Tüm bileşenler $R = 0.90$

**Bileşen 5'i anahtar seç:**

**Durum 1**: 5 çalışıyor ($P = 0.90$)
```
    A●───────●C
     │       │
     │   5   │
     │   ●   │
     │       │
    B●───────●D
```
Sol taraf (1 veya 2): $R_L = 1 - (0.1)(0.1) = 0.99$
Sağ taraf (3 veya 4): $R_R = 1 - (0.1)(0.1) = 0.99$
$R_{sys|5 \text{ çalışıyor}} = 0.99 \times 0.99 = 0.9801$

**Durum 2**: 5 arızalı ($P = 0.10$)
```
    A●───────●C
     │       │
     │       │
     │       │
    B●───────●D
```
Yol 1-4 (seri): $0.90 \times 0.90 = 0.81$
Yol 2-3 (seri): $0.90 \times 0.90 = 0.81$
Paralel: $R_{sys|5 \text{ arızalı}} = 1 - (0.19)(0.19) = 0.9639$

$$R_{sys} = 0.90 \times 0.9801 + 0.10 \times 0.9639 = 0.8821 + 0.0964 = 0.9785$$

### Yöntem 2: Event Space (Durum Uzayı)

Tüm $2^n$ durumu enumerate et, başarılı olanları topla.

### Örnek 21: Köprü Sistemi (Event Space)

5 bileşen → $2^5 = 32$ durum

A'dan D'ye yol var mı? (1-4, 2-3, 1-5-3, 2-5-4 yolları)

| Durum | 1 | 2 | 3 | 4 | 5 | Yol Var? | Olasılık |
|-------|---|---|---|---|---|----------|----------|
| 11111 | 1 | 1 | 1 | 1 | 1 | Evet | $R^5$ |
| 11110 | 1 | 1 | 1 | 1 | 0 | Evet | $R^4 Q$ |
| 11101 | 1 | 1 | 1 | 0 | 1 | Evet | $R^4 Q$ |
| ... | | | | | | | |
| 00000 | 0 | 0 | 0 | 0 | 0 | Hayır | $Q^5$ |

Tüm başarılı durumların olasılıkları toplanır.

### Yöntem 3: Minimal Yol / Minimal Kesim

**Minimal Yol (Minimal Path)**: Sistemin çalışması için gerekli en küçük bileşen kümesi

Köprü için minimal yollar: $\{1,4\}, \{2,3\}, \{1,5,3\}, \{2,5,4\}$

**Minimal Kesim (Minimal Cut)**: Sistemin arızalanması için yeterli en küçük bileşen kümesi

Köprü için minimal kesimler: $\{1,2\}, \{3,4\}, \{1,5,4\}, \{2,5,3\}, \{1,3\}, \{2,4\}$

### Sınırlar (Bounds)

**Alt sınır** (minimal kesimlerden):
$$R_{sys} \geq \prod_{j} \left(1 - \prod_{i \in C_j} Q_i\right)$$

**Üst sınır** (minimal yollardan):
$$R_{sys} \leq 1 - \prod_{j} \left(1 - \prod_{i \in P_j} R_i\right)$$

---

# KISIM VIII: STANDBY (BEKLEME) SİSTEMLERİ

## 1. Türleri

| Tür | Bekleme Durumu | Arıza Oranı | Avantaj/Dezavantaj |
|-----|----------------|-------------|-------------------|
| **Cold Standby** | Tamamen kapalı | $Z_s = 0$ | En yüksek güvenilirlik, anahtarlama gecikmesi |
| **Warm Standby** | Düşük güç | $0 < Z_s < Z$ | Orta seviye |
| **Hot Standby** | Tam güç | $Z_s = Z$ | Hızlı geçiş, paralel sisteme eşdeğer |

## 2. Cold Standby (Mükemmel Anahtarlama)

### Güvenilirlik Fonksiyonu

n özdeş üstel bileşen ($\theta$ ortalama ömür):

$$R_{sys}(t) = e^{-t/\theta} \sum_{i=0}^{n-1} \frac{(t/\theta)^i}{i!}$$

Bu, Poisson dağılımının CDF'idir!

### MTTF

$$MTTF_{sys} = n \cdot \theta$$

### Örnek 22: Cold Standby vs Paralel

$\theta = 1000$ saat, $n = 3$ bileşen

**Paralel sistem**:
$$MTTF_{paralel} = 1000 \times (1 + 1/2 + 1/3) = 1833 \text{ saat}$$

**Cold standby**:
$$MTTF_{cold} = 3 \times 1000 = 3000 \text{ saat}$$

$t = 2000$ saat için güvenilirlik:

**Paralel**: $R_{par}(2000) = 1 - (1 - e^{-2})^3 = 0.753$

**Cold standby**: $R_{cold}(2000) = e^{-2}(1 + 2 + 2) = e^{-2} \times 5 = 0.677$

**Yorum**: Cold standby daha yüksek MTTF verir ama belirli zamanlarda paralel daha iyi olabilir!

## 3. Kusurlu Anahtarlama

Anahtarlama başarı olasılığı $p$ ise:

$$R_{sys}(t) = R_1(t) + p \cdot \int_0^t f_1(u) R_2(t-u) \, du$$

### Örnek 23: %90 Anahtarlama Güvenilirliği

2 bileşenli cold standby, $\theta = 1000$ saat, $p = 0.90$

$$MTTF_{sys} = \theta + p \cdot \theta = 1000 + 900 = 1900 \text{ saat}$$

(Mükemmel anahtarlamada 2000 saat olurdu)

---

# KISIM IX: BAĞIMLILIK ANALİZİ - İLERİ KONULAR

## 1. Copula Modelleri

Bileşen ömürleri arasındaki bağımlılığı modellemek için copula fonksiyonları kullanılır.

### Frank Copula

$$C(u, v) = -\frac{1}{\alpha} \ln\left(1 + \frac{(e^{-\alpha u} - 1)(e^{-\alpha v} - 1)}{e^{-\alpha} - 1}\right)$$

$\alpha > 0$: pozitif bağımlılık
$\alpha < 0$: negatif bağımlılık
$\alpha \to 0$: bağımsızlık

### Örnek 24: Copula ile Paralel Sistem

İki bileşen, üstel marjinaller, Frank copula ($\alpha = 2$)

$$R_{sys}(t) = R_1(t) + R_2(t) - C(R_1(t), R_2(t))$$

## 2. Ortak Yük (Common Load) Modeli

Bileşenler aynı rastgele yüke maruz kalırsa:

$$R_{sys}(t) = \int_0^{\infty} \prod_{i=1}^{n} P(C_i > s) \cdot f_S(s) \, ds$$

Burada $C_i$ bileşen kapasitesi, $S$ ortak yük.

## 3. Marshall-Olkin Modeli

Üç tür şok:
- Şok 1: Sadece bileşen 1'i etkiler (oran $\lambda_1$)
- Şok 2: Sadece bileşen 2'yi etkiler (oran $\lambda_2$)
- Şok 12: Her ikisini de etkiler (oran $\lambda_{12}$)

$$R(t_1, t_2) = \exp(-\lambda_1 t_1 - \lambda_2 t_2 - \lambda_{12} \max(t_1, t_2))$$

---

# KISIM X: ÖZET VE KARŞILAŞTIRMA

## Ana Formüller

| Sistem | Bağımsız Güvenilirlik | Bağımlılık Etkisi |
|--------|----------------------|-------------------|
| **Seri** | $R_{sys} = \prod_i R_i$ | $\rho > 0 \Rightarrow R \uparrow$ |
| **Paralel** | $R_{sys} = 1 - \prod_i Q_i$ | $\rho > 0 \Rightarrow R \downarrow$ |
| **k-out-of-n** | Binom formülü | $\rho > 0 \Rightarrow R \downarrow$ (genelde) |
| **Cold Standby** | Poisson toplamı | Anahtarlama hatası $\Rightarrow R \downarrow$ |

## Tasarım Önerileri

1. **Kritik bileşenleri belirle**: Seri sistemlerde en zayıf halka
2. **Yedeklilik stratejisi**: 
   - Bağımsız bileşenler → Paralel tercih et
   - Bağımlı bileşenler → Cold standby düşün
3. **CCF'yi minimize et**: Farklı üreticiler, farklı teknolojiler
4. **Korelasyonu hesaba kat**: Bağımsızlık varsayımı tehlikeli olabilir

## Kaynaklar

1. Ebeling, C.E. (2010). *An Introduction to Reliability and Maintainability Engineering*
2. Modarres, M. (2006). *Risk Analysis in Engineering*
3. Barlow, R.E. & Proschan, F. (1975). *Statistical Theory of Reliability and Life Testing*
4. Nelson, W. (1982). *Applied Life Data Analysis*
5. Hoyland, A. & Rausand, M. (1994). *System Reliability Theory*
