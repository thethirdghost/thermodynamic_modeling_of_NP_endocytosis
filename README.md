# thermodynamic\_modeling\_of\_NP\_endocytosis

This repository contains the Python implementation of a thermodynamic model to evaluate how nanoparticle \(NP\) size influences membrane wrapping during endocytosis, as described in our research\.

## Overview

This work adopts a  energetic model developed by Agudo Canalejo and Lipowsky (ACS Nano 2015) to characterize the nanoparticle endocytosis rate dependent on particle diameter, for both clathrin\-independent \(CIE\) and clathrin\-mediated \(CME\) endocytosis pathways\. The model focuses on size\-related membrane wrapping effects, and directly computes normalized endocytosis rate using a unified physical formula\.

The core formula for endocytosis rate is defined as:
$y = \frac{X^2 - (1 + a \cdot X)^2}{X^3}$
Where negative values are set to 0 to conform to physical constraints, and the final results are normalized to a maximum value of 1\.

### Variable Definition

- $D$: Nanoparticle diameter \(nm\), ranging from 1 to 200 nm

- $R_w$: Characteristic adhesion length scale = 48\.6 nm

- $D_w = 2R_w$

- $X = D / D_w$

- $a = R_w \cdot m_{bo}$

- $m_{bo}$: Membrane local spontaneous curature \(different for CIE and CME\)


## Requirements

Python 3\.12\+
numpy \>= 1\.21\.0
matplotlib \>= 3\.4\.0

Install dependencies via pip:

```bash
pip install numpy matplotlib
```

## How to Run

1. Clone this repository

```bash
git clone https://github.com/[Your-Username]/[Repo-Name].git
cd [Repo-Name]
```

2. Run the main script

```bash
python englufment_rate.py
```

The script will:

- Calculate endocytosis rate over NP diameter from 1 nm to 200 nm with high\-resolution sampling

- Generate independent figures for CIE and CME pathways

- Mark the optimal diameter with maximum endocytosis rate on plots

- Print peak diameter values for two pathways in the terminal

## Key Parameters

| Parameter      | Description                          | Unit | CIE Value | CME Value |
| -------------- | ------------------------------------ | ---- | --------- | --------- |
| $R_w$          | Characteristic adhesion length       | nm   | 48\.6     | 48\.6     |
| $m_{bo}$       | Local spontaneous curature           | nm-1   | \-0\.1    | \-0\.025  |
| Diameter range | Simulated NP diameter range          | nm   | 1 \~ 200  | 1 \~ 200  |



## Output

### Figures

Two standalone figures \(8 × 5 inches\):

- X\-axis: NP Diameter \(nm\), range: 0 \~ 200 nm

- Y\-axis: Normalized Endocytosis rate, range: 0 \~ 1\.05 \(y\-axis ticks hidden\)

- Dashed vertical line: Optimal diameter with maximum endocytosis rate

- Text label: Exact value of peak diameter

- Gray dashed horizontal line: Baseline at y = 0

### Terminal Output

Print the optimal diameter for CIE and CME respectively:

```Plain Text
CIE 峰值直径 = XX.X nm
CME 峰值直径 = XX.X nm
```

> （注：文档部分内容可能由 AI 生成）
