    # How-Long-is-the-Coast-of-Britain-
MATLAB project for AMS 595/DCS 525 exploring Mandelbrot’s coastline paradox. Implements fractal iteration, bisection boundary detection, polynomial fitting, and arc-length integration to approximate the Mandelbrot set boundary length.
# AMS 595 / DCS 525 Project 2: How Long is the Coast of Britain?

This repository contains MATLAB code, figures, and report files for Project 2.

## Overview
We approximate the length of the Mandelbrot set boundary as an analogue to the coastline paradox:
1. **Task 1:** Implemented `fractal(c)` to test divergence of points in the Mandelbrot set.
2. **Task 2:** Used the bisection method with an indicator function to locate boundary points.
3. **Task 3:** Fitted a 15th-order polynomial to the boundary.
4. **Task 4:** Computed curve length numerically via arc-length integration.

## Results
- **Top boundary length:** ~3.9963  
- **Top + bottom boundary length:** ~7.99  

## Repository Layout
- `code/` – MATLAB scripts and functions for each task  
- `figures/` – Generated plots used in the report  
- `report/` – Final LaTeX report  

## Instructions
Run each task script in order (`task1_fractal.m`, `task2_bisection.m`, etc.), or run `run_all.m` to execute the entire pipeline.

## Author
Isaac Odoom Amanor – AMS 595/DCS 525, Stony Brook University

