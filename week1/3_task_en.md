# Week 1: Reliability Theory Analysis

## Problem Statement

Given the hazard rate function:
$$h(t) = \frac{0.4}{0.2t + 1}$$

## Tasks

### Task 1: Analytical Derivation
Derive the following functions analytically:

1.1. **Cumulative Hazard Function** $H(t)$
1.2. **Reliability Function** $R(t)$ 
1.3. **Probability Density Function** $f(t)$
1.4. **Cumulative Distribution Function** $F(t)$

### Task 2: Statistical Metrics
Calculate the following reliability metrics:

2.1. **Mean Time to Failure (MTTF)**
2.2. **Variance** of the lifetime distribution
2.3. **Standard Deviation**
2.4. **Median Life Time** $t_{0.5}$
2.5. **Mode** of the distribution
2.6. **Coefficient of Variation**

### Task 3: Advanced Analysis
3.1. Derive the **Mean Residual Life Time (MRLT)** function $m(t)$
3.2. Evaluate MRLT at time points: t = 0, 1, 2, 5, 10
3.3. Determine the **hazard rate type** (increasing/decreasing/constant)

### Task 4: Visualization
Create the following plots:

4.1. Hazard rate function $h(t)$
4.2. Reliability function $R(t)$
4.3. Probability density function $f(t)$
4.4. Cumulative distribution function $F(t)$
4.5. Mean residual life time $m(t)$
4.6. Combined dashboard with all functions

### Task 5: Simulation Study
5.1. Generate a **50-unit sample** from the given distribution using inverse transform method
5.2. Calculate **empirical statistics** from the sample:
   - Sample mean
   - Sample variance
   - Sample median
   - Sample mode (approximate)
   - Sample coefficient of variation

5.3. Create **comparison table** between theoretical and empirical values
5.4. Perform **statistical validation**:
   - Kolmogorov-Smirnov goodness-of-fit test
   - Visual comparison plots

### Task 6: Documentation and Interpretation
6.1. Document all **mathematical derivations** step by step
6.2. Provide **engineering interpretation** of results
6.3. Discuss **practical implications** of the hazard rate pattern
6.4. Compare **theoretical vs empirical** findings
6.5. Write comprehensive **analysis report**

## Deliverables

Create the following files:

### Scripts and Code
- `solutions.R` - Main analytical calculations
- `visualization.R` - Plotting functions and graphs
- `simulation_study.R` - Sample generation and empirical analysis
- `app.R` - Interactive Shiny application

### Documentation
- `solutions.md` - Complete mathematical solutions and derivations
- `results_report.md` - Analysis results and interpretations  
- `simulation_report.md` - Empirical study findings
- `notes.md` - Study notes and key concepts

### Generated Files
- `sample_data.csv` - Generated 50-unit sample
- `comparison_table.csv` - Theoretical vs empirical statistics
- Individual plot files (PNG format)
- Combined dashboard visualization

## Mathematical Relationships

Key formulas to work with:
- Relationship between hazard and reliability: $R(t) = \exp(-H(t))$
- PDF from hazard and reliability: $f(t) = h(t) \cdot R(t)$
- MTTF definition: $\mu = \int_0^{\infty} R(t) dt$
- MRLT definition: $m(t) = \frac{\int_t^{\infty} R(u) du}{R(t)}$

## Implementation Notes

- Use analytical solutions wherever possible
- Set seeds for reproducible random sampling
- Include clear comments in all code
- Create professional-quality visualizations
- Validate theoretical results with simulation