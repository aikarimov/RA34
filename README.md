# RA34 - Implementation of the Rational Approximation Method for Stiff Initial Value Problems
## A supplementary material for the paper *A. Karimov et al. "Rational Approximation Method for Stiff Initial Value Problems"*

This repository contains all necessary files for computational experiments presented in the article. 
To draw stability regions, run
```matlab
Stabreg_RA;
```
![Image 1](https://github.com/aikarimov/RA34/blob/main/StabReg.png)

To draw time-domain solution of the VDPL test problem, run
```matlab
Test_varstep_VDPL;
```
![Image 2](https://github.com/aikarimov/RA34/blob/main/VDPL-time.jpg)

To draw time-domain solution of the VDPL test problem with element-wise division, run
```matlab
Test_varstep_VDPL_elwise;
```
![Image 3](https://github.com/aikarimov/RA34/blob/main/VDPL1000-elwise.jpg)

To obtain time-domain solutions for the Riccati, HIRES and FNR test problems, respectively, run 
```matlab
Test_varstep_RICCATI;
Test_varstep_HIRES;
Test_varstep_FNR;
```

To draw the performance plot for the VDPL problem, run
```matlab
Performance_varstep_VDPL;
```
![Image 4](https://github.com/aikarimov/RA34/blob/main/VDPL-perf.jpg)

For the Riccati, HIRES and FNR test problems, respectively, run
```matlab
Performance_varstep_RICCATI;
Performance_varstep_HIRES;
Performance_varstep_FNR;
```

Figures in the paper are saved in `.jpg` format using `export_fig` add-on by the command
```matlab
export_fig(gcf,'FNR-time-zoom.jpg','-r300','-q100','-transparent');
```

The structure of the files is the following, as exemplified by the VDPL test problems. 
The function `RA34_VDPL` solves the VDPL problem with the variable stepsize controller using RA4(3) method. On each step, it executes `RA34_step_VDPL` which contains analytical derivatives for the particular problem. Similarly, the Taylor series method for the VDPL problem is implemented in `Taylor34filtVDPL` function, which executes `Taylor34_step_1step_VDPL` on each step.

Functions `ERK34` and `Lobatto34` implement Runge-Kutta 4(3) and LobattoIIIC 4(3) methods, which do not need analytical differentiation and therefore are universal solvers.

Use `syms_dF` for finding high-order function derivatives in matrix form.

