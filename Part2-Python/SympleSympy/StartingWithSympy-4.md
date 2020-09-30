---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.12
    jupytext_version: 1.6.0
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# The biharmonic equation for Stokes flow

$$\nabla^4 \psi = -{\rm Ra} \frac{\partial T}{\partial x_1}$$
   

$$\nabla^4 \equiv \nabla^2 ( \nabla ^2) \equiv 
                    \left( \frac{\partial ^4}{\partial x_1^4} + 
                    \frac{\partial ^2}{\partial x_1^2} \frac{\partial ^2}{\partial x_2^2} +
                    \frac{\partial ^4}{\partial x_2^4} \right)$$ 
                    
                    

Solutions take the form:

$ \psi_0 = \sin 2πx \cosh 2πy $ and  $\psi_1 = y \sin 2πx  \cosh 2πy$ 

```{code-cell} ipython3
import sympy
import math
import numpy as np
```

## Symbolic approach

```{code-cell} ipython3
from sympy.core.symbol import Symbol
from sympy.core.function import Function

x       = Symbol('x')
y       = Symbol('y')
k       = Symbol('k', positive=True)

A       = Symbol('A')
B       = Symbol('B')
C       = Symbol('C')
D       = Symbol('D')

Ra      = Symbol('Ra', positive=True)
psi     = Function('psi')
T       = Function('T')

# Potential solution

PSI1 = A * sympy.sin(k*x) * sympy.cosh(k*y) + B * sympy.sin(k*x) * sympy.sinh(k*y) 
PSI1
```

Define the biharmonic operator and the biharmonic equation

```{code-cell} ipython3
BH = psi(x,y).diff(x,4) + psi(x,y).diff((x,2),(y,2)) + psi(x,y).diff(y,4) 
eq1=sympy.Eq(BH, 0)
eq2=sympy.Eq(BH, Ra*T(x,y).diff(x))
```

Substitute the potential solution and verify that it is an eigenfunction of $\nabla^4$

```{code-cell} ipython3
BH1 = BH.replace(psi(x,y), PSI1)
# This is too complicated to simplify in this form ... (try it !)
# BH1.simplify()
```

```{code-cell} ipython3
BH2 = BH1.doit().simplify()
(BH2 / PSI1).simplify()
```

There is another pair of solutions, $y * \psi_0$

```{code-cell} ipython3
PSI2 = y * PSI1

BH1 = BH.replace(psi(x,y), PSI2).doit().simplify()
BH2 = BH1.doit().simplify()
(BH2 / PSI2).simplify()
```

These don't quite satisfy the equation on their own but it is possible to combine the four potential solutions to find a complete eigenfunction

```{code-cell} ipython3

```
