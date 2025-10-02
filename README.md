# AMS 595 / DCS 525 — Project 2  
**Mandelbrot Fractal Analysis**

Author: *Isaac Odoom Amanor*  

---

## 📖 Overview
This project explores the **Mandelbrot set** using numerical and computational methods.  
The workflow is divided into four tasks:

1. **Visualization** — Generate Mandelbrot fractal images using the escape-time algorithm.  
2. **Boundary Extraction** — Apply automatic bracketing + bisection to find the fractal boundary.  
3. **Polynomial Fit** — Approximate the top boundary with a degree–15 polynomial.  
4. **Boundary Length** — Compute the length of the polynomial curve via arc-length integration.  

The project demonstrates how fractal geometry can be visualized, approximated, and quantified numerically.

---

## ⚙️ Methods
- Implemented in **MATLAB**.  
- Escape-time iteration for visualization.  
- Bracketing and bisection for robust boundary detection.  
- Polynomial fitting using `polyfit` (order 15).  
- Arc-length computed with numerical integration.  

---

## 📊 Results
- Escape-time visualization of Mandelbrot set.  
- Boundary points (top and bottom) extracted accurately.  
- Smooth polynomial fit to the top boundary.  
- Boundary length on interval $[-1.759850, \, 0.367730]$:  
  **L ≈ 3.9963**

---

## 🖼️ Figures

| Task 1: Visualization | Task 2: Boundary Extraction |
|-----------------------|-----------------------------|
| ![Figure 1](Figure_1.pdf) | ![Figure 2](Figure_2.pdf) |

| Task 3: Polynomial Fit | Task 4: Arc-Length Integrand |
|------------------------|------------------------------|
| ![Figure 3](Figure_3.pdf) | ![Figure 4](Figure_4.pdf) |

---

## 📂 Repository Structure
