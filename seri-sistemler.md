# KISIM III: SERİ SİSTEMLER

## 1. Mantık ve Yapı

**Seri sistem:** Bütün bileşenlerin çalışması gerekir; **tek bir bileşen bile arızalanırsa sistem durur.**

- \(X_i(t)\): \(i\). bileşen, \(t\) anında **çalışıyorsa** 1, **arızalıysa** 0.
- Sistem ancak tümü 1 ise çalışır:

\[
\phi(\mathbf{x}) = \min(x_1,\dots,x_n) = \prod_{i=1}^n x_i.
\]

Bu, “en zayıf halka” sezgisinin matematik karşılığıdır.

## 2. Bağımsız Bileşenler: Güvenilirlik

\(T_i\): \(i\). bileşenin ömrü.  

Sistemin \(t\) anına kadar çalışıyor olması demek:

\[
\{T_{\text{sys}} > t\} = \{T_1 > t, \dots, T_n > t\} = \bigcap_{i=1}^n \{T_i > t\}.
\]

Bileşenler **bağımsız** ise:

\[
R_{\text{sys}}(t) = P(T_{\text{sys}} > t) 
= P\Big(\bigcap_{i=1}^n \{T_i > t\}\Big) 
= \prod_{i=1}^n P(T_i > t)
= \prod_{i=1}^n R_i(t).
\]

Yani **seri sistem güvenilirliği = bileşen güvenilirliklerinin çarpımı**.

## 3. Tehlike Oranı (Hazard) – Neden Toplanıyor?

Tanım:

\[
Z_i(t) = \frac{f_i(t)}{R_i(t)}, \quad
R_i(t) = P(T_i>t).
\]

Seri sistem için:

\[
R_{\text{sys}}(t)=\prod_{i=1}^n R_i(t)
\quad\Rightarrow\quad
\ln R_{\text{sys}}(t) = \sum_{i=1}^n \ln R_i(t).
\]

Türevi al:

\[
\frac{d}{dt} \ln R_{\text{sys}}(t) 
= \sum_{i=1}^n \frac{R_i'(t)}{R_i(t)}
\quad\Rightarrow\quad
- Z_{\text{sys}}(t)=\sum_{i=1}^n (-Z_i(t)).
\]

Dolayısıyla

\[
Z_{\text{sys}}(t) = \sum_{i=1}^n Z_i(t).
\]

> Seri sistemde **risk oranları toplanır**: “ne kadar çok bileşen, o kadar toplam risk”.

## 4. Üstel Ömürler – Temiz Örnek

Her bileşen üstel:

- \(T_i \sim \text{Exp}(\theta_i)\), yani
  \[
  R_i(t) = e^{-t/\theta_i},\quad Z_i(t) = 1/\theta_i.
  \]

Seri sistem için:

\[
R_{\text{sys}}(t) = \prod_{i=1}^n e^{-t/\theta_i}
= \exp\Big(-t \sum_{i=1}^n \tfrac{1}{\theta_i}\Big).
\]

Bu da yine **üstel** bir dağılım:

\[
T_{\text{sys}} \sim \text{Exp}(\theta_{\text{sys}}), \quad
\frac{1}{\theta_{\text{sys}}} = \sum_{i=1}^n \frac{1}{\theta_i}.
\]

Dolayısıyla

\[
MTTF_{\text{sys}} = \theta_{\text{sys}} = \Big(\sum_{i=1}^n \frac{1}{\theta_i}\Big)^{-1}.
\]

**Örnek:**  
\(\theta_1=1000\), \(\theta_2=2000\), \(\theta_3=5000\):

\[
\frac{1}{\theta_{\text{sys}}} = 0.001 + 0.0005 + 0.0002 = 0.0017
\Rightarrow \theta_{\text{sys}} \approx 588.2\ \text{saat}.
\]

Tek tek ortalamalar 1000–5000 arasında, ama sistemin ortalama ömrü ≈ 588 saat. Seri bağlantının “sistemi zayıflattığı” açık.

## 5. Bağımlı Bileşenler – Korelasyon Neden Önemli?

Bağımsız değilsek:

\[
R_{\text{sys}}(t)=P\Big(\bigcap_{i=1}^n \{T_i>t\}\Big)
\neq \prod_i P(T_i>t).
\]

İki bileşen için güvenilirlik \(R_1, R_2\), arıza olasılıkları \(Q_1, Q_2\) olsun.  
Çok kaba bir korelasyon yaklaşımı (olay düzeyinde):

\[
P(Q_1 \cap Q_2)
= Q_1 Q_2 + \rho\,\sqrt{Q_1(1-Q_1)Q_2(1-Q_2)}.
\]

- **Seri sistem:** Sistem arızası = “en az bir bileşen arızalı”.  
  Pozitif korelasyon → “ya ikisi birden arızalı ya ikisi de sağlam” →  
  \(\Rightarrow\) **seri sistem güvenilirliği artabilir** (karma durumlar azalıyor).

Paralelde bunun tam tersi etkiyi göreceğiz.

---

# KISIM IV: PARALEL SİSTEMLER

## 1. Mantık ve Yapı

Paralel sistemde **en az bir bileşenin çalışması yeterlidir.**  

- Sistem ancak **tüm bileşenler arızalıysa** çöker.

Bu yüzden paralel sistem “yedeklilik” sağlar.

## 2. Bağımsız Bileşenler: Temel Formül

\[
Q_i(t) = 1 - R_i(t).
\]

Sistem ancak

\[
\{T_1 \le t, \dots, T_n \le t\}
\]

olduğunda arızalanır:

\[
Q_{\text{sys}}(t) = P\Big(\bigcap_{i=1}^n\{T_i \le t\}\Big).
\]

Bağımsızlık varsayarsak:

\[
Q_{\text{sys}}(t) = \prod_{i=1}^n Q_i(t).
\]

Dolayısıyla güvenilirlik:

\[
R_{\text{sys}}(t) = 1 - Q_{\text{sys}}(t) 
= 1 - \prod_{i=1}^n [1-R_i(t)].
\]

Bu formül tamamen **tümleyen olay** fikrinden geliyor.

### Örnek: Basit Paralel Sistem

3 bileşen:  
\(R_1 = 0.90,\ R_2 = 0.85,\ R_3 = 0.80\).

\[
Q_{\text{sys}} = 0.10 \times 0.15 \times 0.20 = 0.003
\Rightarrow R_{\text{sys}} = 0.997.
\]

En iyi bileşen 0.90 ama sistem **0.997** → paralel yedekliliğin gücü.

## 3. Özdeş Bileşenler

Her bileşen için \(R_i(t) = R(t)\):

\[
R_{\text{sys}}(t) = 1 - [1-R(t)]^n.
\]

**Örnek:** Her biri \(R(t)=0.8\):

- \(n=1\): \(R_{\text{sys}}=0.8\)  
- \(n=2\): \(1-(0.2)^2 = 0.96\)  
- \(n=3\): \(1-(0.2)^3 = 0.992\)

n büyüdükçe belirli bir \(t\) için sistem güvenilirliği 1’e yaklaşır.

## 4. Üstel Ömür ve MTTF

Bileşenler üstel: \(T_i\sim \text{Exp}(\theta)\), **özdeş ve bağımsız**.

- Seri sistem → minimum \(T_{\min} = \min T_i\).
- Paralel sistem → maksimum \(T_{\max} = \max T_i\).

Paralel sistem ömrü = **maksimum** olduğu için \(T_{\max}\)’in beklenen değeri:

\[
MTTF_{\text{sys}} = E[T_{\max}] 
= \theta\sum_{i=1}^n \frac{1}{i}.
\]

Bu, sıralı istatistikler (n adet üstelin sıralanması) teorisinden gelir.

**Örnek:** \(\theta=1000\) saat

- \(n=1\): \(MTTF=1000\)
- \(n=2\): \(1000(1+1/2)=1500\)
- \(n=3\): \(1000(1+1/2+1/3)\approx 1833\)

Artış var ama **azalan marjinal kazanç**: İlk yedek çok faydalı, sonrakilerin katkısı daha küçük.

## 5. Bağımlı Bileşenler – Korelasyon Paraleli Neden Bozar?

Paralel sistemde sistem arızası = “hepsi arızalı”:

\[
Q_{\text{sys}}(t) = P(Q_1 \cap \cdots \cap Q_n).
\]

Pozitif korelasyon → “biri bozulursa diğerleri de büyük ihtimalle bozuluyor”.  
Dolayısıyla \(\,P(Q_1 \cap Q_2)\) **bağımsızlığa göre büyüyor**, bu da \(Q_{\text{sys}}\)’i büyütüyor → **güvenilirlik düşüyor**.

İki bileşen örneği (her biri \(R=0.9, Q=0.1\)):

- Bağımsız:  
  \(Q_{\text{sys}} = 0.1\times 0.1 = 0.01\), \(R_{\text{sys}} = 0.99\).
- Pozitif korelasyonlu (örneğin \(\rho = 0.3\)):  
  \(P(Q_1\cap Q_2)=0.01+0.3\times 0.09=0.037\),  
  \(R_{\text{sys}} = 1-0.037 = 0.963 < 0.99\).

**Özet:**  
- Seri + pozitif korelasyon ⇒ **güvenilirlik genelde artma eğiliminde**  
- Paralel + pozitif korelasyon ⇒ **güvenilirlik azalma eğiliminde**

---

# KISIM V: k-out-of-n SİSTEMLER

## 1. Tanım

- \(n\) bileşenden **en az \(k\) tanesinin çalışması** gereken sistem.
- Notasyon: **k-out-of-n:G** (G = Good).

Özel durumlar:

- Seri sistem: \(k=n\)
- Paralel sistem: \(k=1\)
- Çoğunluk oylama: \(k = \lceil n/2 \rceil\)

## 2. Bağımsız ve Özdeş Bileşenler

Her bileşen için \(R_i(t)=R(t)\), \(Q(t)=1-R(t)\).  

Çalışan bileşen sayısı \(X\):

\[
X \sim \text{Binomial}(n, R(t)).
\]

Sistem çalışıyor ⇔ \(X \ge k\):

\[
R_{\text{sys}}(t) = P(X \ge k)
= \sum_{j=k}^{n} \binom{n}{j} [R(t)]^{j}[1-R(t)]^{n-j}.
\]

**Formülün kökeni:** “kaç tanesi çalışıyor?” sorusunu ikili (başarı/başarısızlık) denemeler gibi görüp Binom dağılımı ile yazmak.

### Örnek: 2-out-of-3 Sistem

Her bileşen için \(R(t)=0.9\):

\[
R_{\text{sys}} 
= \binom{3}{2}(0.9)^2(0.1)^1 + \binom{3}{3}(0.9)^3(0.1)^0
= 3\cdot 0.81\cdot 0.1 + 0.729
= 0.243 + 0.729
= 0.972.
\]

Bu yapı hem seriden (3/3) hem paralelden (1/3) farklı, arada bir tasarım.

## 3. Farklı Güvenilirlikli Bileşenler

Bileşen güvenilirlikleri \(R_1,\dots,R_n\) farklıysa, Binom formülü yok.  

En temiz yöntem: **durum uzayı (state enumeration)**.

Örneğin 2-out-of-3 sistem, \(R_1=0.95, R_2=0.9, R_3=0.85\).  
Çalışan senaryolar (en az 2 çalışan) tek tek yazılır, her birinin olasılığı çarpılarla bulunup toplanır.

---

# KISIM VI: SERİ-PARALEL KOMBİNE SİSTEMLER

## 1. Çözüm Stratejisi

Kompleks görünse bile yöntem basit:

1. **Blok diyagramını çiz.**
2. “**İçten dışa**” git:
   - Önce en içteki seri/paralel grubu **tek bir eşdeğer blok** haline indir.
   - Sonra bir üst kademeye geç.
3. Bunu sistem tek blok kalana kadar tekrarla.

Hepsi bağımsız ise, her blok için:

- Seri: \(R = \prod R_i\)
- Paralel: \(R = 1-\prod(1-R_i)\)

## 2. Seri-Paralel Örneği

Önce paralel, sonra seri:

- Paralel blok (1,2)
- Sonra 3 ile seri

Adımlar:

1. \(\displaystyle R_{12} = 1 - (1-R_1)(1-R_2)\)
2. \(\displaystyle R_{\text{sys}} = R_{12}\cdot R_3.\)

## 3. Paralel-Seri Örneği

Önce seri, sonra paralel:

- Üst dal: 1 ve 2 seri
- Alt dal: 3 ve 4 seri
- Sonra bu iki dal paralel

Adımlar:

1. \(R_{12}=R_1R_2,\quad R_{34}=R_3R_4\)
2. Paralel:
   \[
   R_{\text{sys}} = 1 - (1-R_{12})(1-R_{34}).
   \]

Bu örneklerde yaptığımız tek şey:  
**seri blok → çarp**, **paralel blok → 1-(1-R)’lerin çarpımı**.

---

# KISIM VII: KOMPLEKS SİSTEMLER

Basit seri/paralel indirgeme ile çözülemeyen yapılara **kompleks sistem** denir (köprü/bridge gibi).

## 1. Köprü Yapısı ve Neden Zor?

Klasik bridge:

- Dört kenar (1–4), ortada çapraz bir 5. eleman.
- 5 numara yüzünden sistem ne tam seri ne tam paralel bloklara ayrılabiliyor.

Bu durumda üç temel teknik kullanılır:

1. **Decomposition (ayrıştırma)**  
2. **Durum uzayı (2ⁿ durum) enumerate etme**  
3. **Minimal yollar / minimal kesimler**

## 2. Decomposition Fikri

Bir bileşeni “anahtar” seç → üzerinde koşullu olasılık kullan:

\[
R_{\text{sys}} 
= P(\text{çalışıyor})
= P(\text{sys çalışıyor, 5 çalışıyor}) + P(\text{sys çalışıyor, 5 arızalı}).
\]

Bu da:

\[
R_{\text{sys}}
= R_5 \cdot R_{\text{sys}|5\text{ çalışıyor}} 
+ Q_5 \cdot R_{\text{sys}|5\text{ arızalı}}.
\]

Geri kalan iki koşullu sistem artık **seri+paralel kombinasyonuna** indirgenebilir, orada Kısım VI’daki kuralları kullanırsın.

## 3. Minimal Yol / Minimal Kesim

- **Minimal yol**: Sistemin çalışmasını garanti eden, gereksiz eleman içermeyen en küçük bileşen kümesi.
- **Minimal kesim**: Sistemi mutlaka durduran, gereksiz eleman içermeyen en küçük kümeler.

Bunları kullanarak **üst ve alt sınırlar (bounds)** elde edilir; tam analitik çözüm zor olduğunda iş görür.

---

# KISIM VIII: STANDBY (BEKLEME) SİSTEMLERİ

Paralel sistemde **hepsi aynı anda aktif**. Standby sistemde yedek elemanlar **boşta bekler**, sadece gerektiğinde devreye girer.

## 1. Türler (Sezgi)

| Tür           | Beklerken durumu      | Arıza oranı        | Paralel ile ilişkisi        |
|---------------|------------------------|--------------------|-----------------------------|
| Cold standby  | Tam kapalı             | neredeyse 0        | “Sırayla kullanılan yedek” |
| Warm standby  | Düşük güç              | 0 ile Z arasında   | Arada                      |
| Hot standby   | Tam güçte aktif        | Z                  | Paralel sisteme yakın      |

Hot standby ≈ paralel.  
Cold standby ≈ “biri biter diğerine geçilir”.

## 2. Cold Standby – Üstel Bileşenler

n adet üstel bileşen, ortalama ömür \(\theta\), hepsi **cold standby** (mükemmel anahtarlama):

- Sistem ömrü ≈ **toplam ömürlerin toplamı** → n tane üstelin toplamı \(\Rightarrow \text{Gamma}(n,\theta)\).
- Bu dağılımın MTTF’i:

\[
MTTF_{\text{sys}} = n\cdot \theta.
\]

Bu, “sırayla kullanıyoruz, her birinin ortalaması \(\theta\), toplamda \(n\theta\)” sezgisine tam uyar.

Güvenilirlik fonksiyonu:

\[
R_{\text{sys}}(t) 
= e^{-t/\theta} \sum_{i=0}^{n-1} \frac{(t/\theta)^i}{i!}.
\]

Bu ifade, Poisson CDF formunun bire bir aynısıdır.

---

# KISIM IX: BAĞIMLILIK İÇİN İLERİ MODELLER (ÇOK KISA)

Gerçek sistemlerde bileşen ömürleri çoğu zaman bağımsız değildir. Örnek ileri modeller:

- Ortak çevresel yük (sıcaklık, titreşim)
- Ortak neden arızaları (CCF)
- Şok modelleri (Marshall–Olkin)
- Genel bağımlılık için **copula**’lar (Frank, Clayton, Gumbel…)

Bu araçlar, “marjinler” (tek tek ömür dağılımları) sabitken, aralarındaki **bağımlılık yapısını** parametrize eder.

---

# KISIM X: KISA ÖZET

- **Seri sistem:**
  - Çalışma olayı = **hepsi çalışıyor**.
  - \(R_{\text{sys}} = \prod R_i\).
  - Hazard’lar **toplanır**: \(Z_{\text{sys}} = \sum Z_i\).
  - Pozitif korelasyon çoğu durumda sistem lehine.

- **Paralel sistem:**
  - Çalışma olayı = **en az biri çalışıyor**.
  - \(R_{\text{sys}} = 1 - \prod(1-R_i)\).
  - Pozitif korelasyon yedeklerin “birlikte çökmesine” yol açıp güvenilirliği düşürür.

- **k-out-of-n:**
  - “En az k çalışan” ifadesini Binom ile yazarız (özdeş ise).
  - Seri/paralelin genellemesi.

- **Seri-paralel kombine sistemler:**
  - Blok diyagramı çiz, içten dışa indirgeme yap.
  - Seri blok → çarp, paralel blok → 1-(1-R)’lerin çarpımı.

- **Kompleks sistemler (bridge vb.):**
  - Decomposition (koşullu olasılık), durum uzayı, minimal yol/kesim, bounds gibi araçlar.

- **Standby sistemler:**
  - Cold standby: sıralı kullanım, MTTF ≈ \(n\theta\).
  - Hot standby: paralel sisteme benzer davranış.
