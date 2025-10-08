
# Week 1: Reliability Theory - Fundamental Concepts and Hazard Function

## 📚 Lecture Notes - Reliability Theory

### 🎯 Week 1 Topics
- Basic reliability functions
- Hazard rate (failure rate) function
- Mathematical relationships between functions

---

## 📊 Fundamental Definitions

### 🔹 Failure Time
**T**: Time from when a system starts operating until it fails
- Random variable
- t ≥ 0 (cannot be negative)

### 🔹 Reliability Function - R(t)
**Definition:** The probability that a unit will not fail (survive) up to time t

$$R(t) = P(T > t)$$

**Properties:**
- Non-increasing function
- $\lim_{t \to -\infty} R(t) = 1$ (perfect reliability at start)
- $\lim_{t \to \infty} R(t) = 0$ (cannot operate forever)
- $R(0) = 1$ (operational at start)

### 🔹 Cumulative Distribution Function - F(t)
**Definition:** The probability that a unit fails before or at time t

$$F(t) = P(T \leq t) = \int_{-\infty}^{t} f(u)du$$

**Properties:**
- Non-decreasing function
- Right-continuous
- $\lim_{t \to -\infty} F(t) = 0$
- $\lim_{t \to \infty} F(t) = 1$

### 🔹 Probability Density Function - f(t)
**Definition:** Instantaneous failure probability rate

$$f(t) = \frac{dF(t)}{dt} = \lim_{\Delta t \to 0} \frac{F(t+\Delta t) - F(t)}{\Delta t}$$

**Alternative definition:**
$$f(t) = \lim_{\Delta t \to 0} \frac{P(t < T < t+\Delta t)}{\Delta t}$$

**Properties:**
- $f(t) \geq 0$ for $t > 0$
- $\int_{-\infty}^{\infty} f(t)dt = 1$ (normalization condition)

---

## 🔗 Fundamental Relationships Between Functions

### 📐 R(t) and F(t) Relationship
$$R(t) + F(t) = 1$$
$$R(t) = 1 - F(t)$$

### 📐 f(t) and R(t) Relationship
$$f(t) = -\frac{dR(t)}{dt}$$

### 📐 Failure Probability in Specific Interval
Failure probability in interval $(t_1, t_2)$:
$$P(t_1 < T < t_2) = F(t_2) - F(t_1) = R(t_1) - R(t_2)$$

---

## ⚡ Hazard Function (Failure Rate) - h(t) or z(t)

### 🎯 Conceptual Definition
The rate of failure probability in the next instantaneous time interval $(t, t+\Delta t)$, **given that a unit has operated successfully up to time t**.

### 📊 Mathematical Definition
$$h(t) = \lim_{\Delta t \to 0} \frac{P(t < T < t+\Delta t \mid T > t)}{\Delta t}$$

### 🔄 Basic Hazard Formulas

#### 1. Basic Rate Formula
$$h(t) = \frac{f(t)}{R(t)}$$

#### 2. Logarithmic Derivative Formula
$$h(t) = -\frac{d}{dt} \ln(R(t)) = -\frac{1}{R(t)} \cdot \frac{dR(t)}{dt}$$

#### 3. Cumulative Hazard Function
$$H(t) = \int_0^t h(u)du = -\ln(R(t))$$

### 🔄 Inverse Transformations

#### Reliability from Hazard
$$R(t) = e^{-H(t)} = e^{-\int_0^t h(u)du}$$

#### PDF from Hazard
$$f(t) = h(t) \cdot R(t) = h(t) \cdot e^{-\int_0^t h(u)du}$$

### ⚖️ Important Inequality
$$f(t) \leq h(t)$$ 
(PDF is always less than or equal to hazard rate)


