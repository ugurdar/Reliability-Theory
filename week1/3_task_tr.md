# Hafta 1: Güvenilirlik Teorisi Analizi

## Problem Tanımı

Verilen hazard rate (bozulma oranı) fonksiyonu:
$$h(t) = \frac{0.4}{0.2t + 1}$$

## Görevler

### Görev 1: Analitik Türetme
Aşağıdaki fonksiyonları analitik olarak türetin:

1.1. **Kümülatif Hazard Fonksiyonu** $H(t)$
1.2. **Güvenilirlik Fonksiyonu** $R(t)$ 
1.3. **Olasılık Yoğunluk Fonksiyonu** $f(t)$
1.4. **Kümülatif Dağılım Fonksiyonu** $F(t)$

### Görev 2: İstatistiksel Metrikler
Aşağıdaki güvenilirlik metriklerini hesaplayın:

2.1. **Ortalama Arıza Süresi (MTTF)**
2.2. Yaşam süresi dağılımının **Varyansı**
2.3. **Standart Sapma**
2.4. **Medyan Yaşam Süresi** $t_{0.5}$
2.5. Dağılımın **Modu**
2.6. **Değişkenlik Katsayısı**

### Görev 3: İleri Analiz
3.1. **Ortalama Kalan Yaşam Süresi (MRLT)** fonksiyonu $m(t)$'yi türetin
3.2. MRLT'yi şu zaman noktalarında değerlendirin: t = 0, 1, 2, 5, 10
3.3. **Hazard rate tipini** belirleyin (artıcı/azalıcı/sabit)

### Görev 4: Görselleştirme
Aşağıdaki grafikleri oluşturun:

4.1. Hazard rate fonksiyonu $h(t)$
4.2. Güvenilirlik fonksiyonu $R(t)$
4.3. Olasılık yoğunluk fonksiyonu $f(t)$
4.4. Kümülatif dağılım fonksiyonu $F(t)$
4.5. Ortalama kalan yaşam süresi $m(t)$
4.6. Tüm fonksiyonları içeren birleşik dashboard

### Görev 5: Simülasyon Çalışması
5.1. Ters dönüşüm yöntemi kullanarak verilen dağılımdan **50 birimlik örneklem** üretin
5.2. Örneklemden **ampirik istatistikleri** hesaplayın:
   - Örneklem ortalaması
   - Örneklem varyansı
   - Örneklem medyanı
   - Örneklem modu (yaklaşık)
   - Örneklem değişkenlik katsayısı

5.3. Teorik ve ampirik değerler arasında **karşılaştırma tablosu** oluşturun
5.4. **İstatistiksel doğrulama** yapın:
   - Kolmogorov-Smirnov uyum iyiliği testi
   - Görsel karşılaştırma grafikleri

### Görev 6: Dokümantasyon ve Yorumlama
6.1. Tüm **matematiksel türevlemeleri** adım adım belgeleyin
6.2. Sonuçların **mühendislik yorumunu** yapın
6.3. Hazard rate deseninin **pratik etkilerini** tartışın
6.4. **Teorik vs ampirik** bulguları karşılaştırın
6.5. Kapsamlı **analiz raporu** yazın

## Çıktılar

Aşağıdaki dosyaları oluşturun:

### Scriptler ve Kod
- `solutions.R` - Ana analitik hesaplamalar
- `visualization.R` - Grafik fonksiyonları ve çizimler
- `simulation_study.R` - Örneklem üretimi ve ampirik analiz
- `app.R` - İnteraktif Shiny uygulaması

### Dokümantasyon
- `solutions.md` - Tam matematiksel çözümler ve türevlemeler
- `results_report.md` - Analiz sonuçları ve yorumlar  
- `simulation_report.md` - Ampirik çalışma bulguları
- `notes.md` - Çalışma notları ve temel kavramlar

### Üretilen Dosyalar
- `sample_data.csv` - Üretilen 50 birimlik örneklem
- `comparison_table.csv` - Teorik vs ampirik istatistikler
- Tekil grafik dosyaları (PNG formatı)
- Birleşik dashboard görselleştirmesi

## Matematiksel İlişkiler

Çalışılacak temel formüller:
- Hazard ve güvenilirlik ilişkisi: $R(t) = \exp(-H(t))$
- Hazard ve güvenilirlikten PDF: $f(t) = h(t) \cdot R(t)$
- MTTF tanımı: $\mu = \int_0^{\infty} R(t) dt$
- MRLT tanımı: $m(t) = \frac{\int_t^{\infty} R(u) du}{R(t)}$

## Uygulama Notları

- Mümkün olduğunca analitik çözümler kullanın
- Tekrarlanabilir rastgele örnekleme için seed ayarlayın
- Tüm kodlarda açık yorumlar ekleyin
- Profesyonel kalitede görselleştirmeler oluşturun
- Teorik sonuçları simülasyonla doğrulayın