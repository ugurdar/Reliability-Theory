# Week 1: Reliability Theory - Temel Kavramlar ve Hazard Fonksiyonu

## 📚 Ders Notları - Güvenilirlik Teorisi

### 🎯 Hafta 1 Konuları
- Temel güvenilirlik fonksiyonları
- Hazard rate (bozulma oranı) fonksiyonu
- Fonksiyonlar arası matematiksel ilişkiler

---

## 📊 Temel Tanımlar

### 🔹 Bozulma Süresi
**T**: Bir sistemin çalışmaya başladığı andan bozulana kadar geçen süre
- Rastgele değişken (random variable)
- t ≥ 0 (negatif olamaz)

### 🔹 Güvenilirlik Fonksiyonu - R(t)
**Tanım:** Bir birimin t anına kadar bozulmama (sağlam kalma) olasılığı

$$R(t) = P(T > t)$$

**Özellikler:**
- Azalmayan (non-increasing) fonksiyon
- $\lim_{t \to -\infty} R(t) = 1$ (başlangıçta mükemmel güvenilirlik)
- $\lim_{t \to \infty} R(t) = 0$ (sonsuza kadar çalışamaz)
- $R(0) = 1$ (başlangıçta çalışır durumda)

### 🔹 Kümülatif Dağılım Fonksiyonu - F(t)
**Tanım:** Bir birimin t anından önce veya t anında bozulma olasılığı

$$F(t) = P(T \leq t) = \int_{-\infty}^{t} f(u)du$$

**Özellikler:**
- Azalmayan (non-decreasing) fonksiyon
- Sağdan sürekli
- $\lim_{t \to -\infty} F(t) = 0$
- $\lim_{t \to \infty} F(t) = 1$

### 🔹 Olasılık Yoğunluk Fonksiyonu - f(t)
**Tanım:** Anlık bozulma olasılığı oranı

$$f(t) = \frac{dF(t)}{dt} = \lim_{\Delta t \to 0} \frac{F(t+\Delta t) - F(t)}{\Delta t}$$

**Alternatif tanım:**
$$f(t) = \lim_{\Delta t \to 0} \frac{P(t < T < t+\Delta t)}{\Delta t}$$

**Özellikler:**
- $t > 0$ için $f(t) \geq 0$
- $\int_{-\infty}^{\infty} f(t)dt = 1$ (normalizasyon koşulu)

---

## 🔗 Fonksiyonlar Arası Temel İlişkiler

### 📐 R(t) ve F(t) İlişkisi
$$R(t) + F(t) = 1$$
$$R(t) = 1 - F(t)$$

### 📐 f(t) ve R(t) İlişkisi
$$f(t) = -\frac{dR(t)}{dt}$$

### 📐 Belirli Aralıkta Bozulma Olasılığı
(t₁, t₂) aralığında bozulma olasılığı:
$$P(t_1 < T < t_2) = F(t_2) - F(t_1) = R(t_1) - R(t_2)$$

---

## ⚡ Hazard Fonksiyonu (Bozulma Oranı) - h(t) veya z(t)

### 🎯 Kavramsal Tanım
Bir birimin **t zamanına kadar çalışır durumda olduğu bilindiğinde**, bir sonraki anlık zaman aralığında $(t, t+\Delta t)$ bozulma olasılığının oranı.

### 📊 Matematiksel Tanım
$$h(t) = \lim_{\Delta t \to 0} \frac{P(t < T < t+\Delta t \mid T > t)}{\Delta t}$$

### 🔄 Temel Hazard Formülleri

#### 1. Temel Oran Formülü
$$h(t) = \frac{f(t)}{R(t)}$$

#### 2. Logaritmik Türev Formülü
$$h(t) = -\frac{d}{dt} \ln(R(t)) = -\frac{1}{R(t)} \cdot \frac{dR(t)}{dt}$$

#### 3. Kümülatif Hazard Fonksiyonu
$$H(t) = \int_0^t h(u)du = -\ln(R(t))$$

### 🔄 Ters Dönüşümler

#### Hazard'dan Güvenilirlik
$$R(t) = e^{-H(t)} = e^{-\int_0^t h(u)du}$$

#### Hazard'dan PDF
$$f(t) = h(t) \cdot R(t) = h(t) \cdot e^{-\int_0^t h(u)du}$$

### ⚖️ Önemli Eşitsizlik
$$f(t) \leq h(t)$$ 
(PDF her zaman hazard rate'den küçük veya eşittir)