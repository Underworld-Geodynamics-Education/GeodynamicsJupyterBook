# Introduction to the Finite Element Method

$$
\newcommand{\dGamma}{\mathbf{d}\boldsymbol{\Gamma}}
\newcommand{\erfc}{\mbox{\rm erfc}}
\newcommand{\curly}{\sf }
\newcommand{\Red     }[1]{\textcolor[rgb]{0.7,0.0,0.0}{ #1}} 
\newcommand{\Green   }[1]{\textcolor[rgb]{0.0,0.7,0.0}{ #1}} 
\newcommand{\Blue    }[1]{\textcolor[rgb]{0.0,0.0,0.7}{ #1}} 
\newcommand{\Emerald }[1]{\textcolor[rgb]{0.0,0.7,0.3}{ #1}} 
$$

## Prelude to Finite Elements: The Variational Calculus

This is something of an aside but it is absolutely necessary to understand the variational method in order to follow how Finite Element Methods work

### Example: How Short is a Straight Line ?

It's intuitively obvious that the shortest distance between two points on a plane is just a straight line. To demonstrate this mathematically is more tricky.

In words, the procedure goes like this: of all possible curves between the two points, find  one (if it exists) which  always becomes longer if it is altered in any way.   (Thinking physically, if the points were linked by a rubber band, then to disturb it from the shortest curve would require additional energy, no matter what that disturbance looked like).

<!-- %% FIGURE: variation for straight line 
\begin{figure}[h]           
    \begin{center}
             \includegraphics[width=0.66\linewidth]{Diagrams/varn_strt.png}
            \caption[]{What curve gives the shortest distance between two
            points in a plane ?  And how to prove it !}
            \label{fig:varn2}
        \end{center}
\end{figure} -->

![What curve gives the shortest distance between two points in a plane ?  And how can we prove it !][variational]

[variational]: ./Diagrams/varn_strt.png

We decide, arbitrarily, that we will make $x$ the independent variable and find a curve $y(x)$ which satisfies the minimum distance requirement.

The distance along a curve is a path integral:

$$ 
    S = \int_{x_1}^{x_2} \frac{ds}{dx} dx   
$$

where

$$ 
        \frac{ds}{dx} = \sqrt{1 + \left( \frac{dy}{dx} \right)^2}
$$

We seek to find $y(x)$ which minimizes $S$. First consider the function

$$ 
    \begin{split}
    Y(x,\alpha) &= y(x) + \alpha \eta(x) \\
    \frac{\partial Y}{\partial x} \equiv Y' &= y'(x) + \alpha \eta'(x)
    \end{split}
$$

where $\eta(x)$ is an arbitrary function which is differentiable and vanishes at $x=x_1,x_2$. This is the variation which we apply to some curve $y(x)$ to see if it gets shorter or longer. Note the notation for the derivative which will be useful as we procede. The optimal path will minimize  $S(\alpha)$, when $\alpha=0$, i.e.

$$ 
    \left. \frac{\partial S(\alpha)}{\partial \alpha} \right|_{\alpha = 0} = 0
$$

In our current notation

 $$ 
     S(\alpha) = \int_{x_1}^{x_2} \sqrt{1 + (Y')^2}
 $$

and 

 $$ 
     \frac{\partial S(\alpha)}{\partial \alpha} = 
         \int_{x_1}^{x_2}  \frac{1}{2} \frac{1}{\sqrt{1 + (Y')^2}} . 2 Y' \frac{\partial Y'}{\partial \alpha} dx
$$

from the definition of $Y'(x,\alpha)$

$$ 
    \frac{\partial Y'}{\partial \alpha} = \eta'(x)
$$

So we now must solve

$$ 
    \left. \frac{\partial S(\alpha)}{\partial \alpha} \right|_{\alpha = 0} =
        \int_{x_1}^{x_2} \frac{y'(x)\eta'(x)}{\sqrt{1 + (y')^2}} dx = 0
$$

We integrate by parts to give

$$ 
    \left. \frac{\partial S(\alpha)}{\partial \alpha} \right|_{\alpha = 0} =
        \left[   \frac{y'(x)\eta(x)}{\sqrt{1 + (y')^2}} \right]_{x_1}^{x_2} - 
        \int_{x_1}^{x_2} \eta{x} \frac{d}{dx} \frac{y'(x)}{\sqrt{1 + (y')^2}} dx = 0
$$

The first term on the RHS vanishes because $\eta$ vanishes at the boundaries. The second term is valid for arbitrary $\eta(x)$ which implies

$$    
    \frac{d}{dx} \frac{y'(x)}{\sqrt{1 + (y')^2}} dx = 0
$$

which in turn implies $y'=$constant, i.e. the equation of a straight line. 

The important things to note here are that an integral method can be used to solve a simple geometrical problem and that the method itself includes the boundary conditions as a natural consequence of the way it is set up.

The variational method can be generalized to solve more important problems. In particular, instead of solving each problem as we have for the straight line / distance question above, we solve a generic problem whose solutions we can apply immediately.

### Generalisation

The general form works like this. To find the function $y(x)$ which produces a stationary value of the functional 

$$
    J=\int_{x_1}^{x_2} F(x,y,y') dx
$$

we work through the same procedure as above, and use the same arguments concerning the arbitrary nature of the variation to obtain

$$ 
\frac{d}{dx}\frac{\partial F}{\partial y'} - \frac{\partial F}{\partial y} = 0
$$

This is known as the Euler equation. Hamilton's principle states that mechanical systems evolve such that the integral 

$$ 
J=\int_{t_1}^{t_2} L dt
$$

is stationary. Here $L$ is the Lagrangian of the system which is identified with a combination of  the work done on the system and the kinetic energy of the system, e.g. potential energy - kinetic energy.

Application of the Euler equation to each direction independently recovers Newton's law ($F=ma$). In complex geometries and with difficult boundary conditions, the variational form may be easier to solve than the differential or "strong" form. _This shows us that there are equivalent integral representations for the standard mechanical equations we are accustomed to using --- variational or weak forms versus differential or strong forms_.


Although this may seem complicated and of rather theoretical interest, in fact it runs throughout finite element methods, and the concept must be familiar in order to follow how FEM works. Advantages of using variational forms of the equations include:

- The simplification of the construction of the governing equations in the sense that scalar quantities --- energies, potentials --- are considered in place of forces, displacements etc. There is also the possibility that such formulations can be derived more-or-less automatically for previously unknown systems.

- Governing equations may be more directly accessible since "unimportant" variables such as internal forces doing no net work do not appear in the variational form.

- When dealing with approximate solutions, the variational form often allows a broader range of trial functions than for the standard differential form. This happens because some boundary conditions are implicit in the formulation and hence are not imposed on the trial functions themselves.


### Example of Variational Forms for FEM

See Klaus-Jürgen Bathe's book for a more complete outline of this approach to FEM. Although this is not the best text to explain how to _implement_ a finite element code, he does a great job of explaining the link between variational methods and weak forms of various equations of progressively increasing complexity.

<!-- %% FIGURE: heated bar
    \begin{figure}[h]           
        \begin{center}
             \includegraphics[width=0.66\linewidth]{Diagrams/heatcond.png}
            \caption[]{A Slab of Material Subjected to an sudden onset of heating at $Q$ on one side 
            at time $t=0$}
            \label{fig:heat1}
        \end{center}
    \end{figure} -->

![A Slab of Material Subjected to an sudden onset of heating at $Q$ on one side at time $t=0$][heated-bar]

[heated-bar]: ./Diagrams/heatcond.png


The functional governing the temperature in the [block of material](#heated-bar) is

$$ 
    \Pi = \int_0^L \frac{1}{2} k \left(\frac{\partial \theta}{\partial x} \right)^2 dx -  
                    \int_0^L \theta q^B dx - \theta(0,t) Q
$$

where $q^B$ is an internal heat generation rate. The fixed boundary condition is $\theta(L,t) =\theta_i$.

This is our generalized problem with 

$$ 
        F = \frac{1}{2}k  \left(\frac{\partial \theta}{\partial x} \right)^2 - \theta q^B
            = \frac{1}{2}k {\theta'}^2 - \theta q^B 
$$

which produces a stationary functional if

$$ 
    \frac{d}{dx}\frac{\partial F}{\partial \theta'} - \frac{\partial F}{\partial \theta} = 0    
$$

or, in other words

$$ 
k\frac{d^2 \theta}{d x^2} = q^B 
$$

Which we recognize to be the governing differential equation. The variational
statement also contains the natural boundary condition

$$ 
    k \left. \frac{\partial \theta}{\partial x}    \right|_{x=0} + Q = 0
$$
    
This is a clear demonstration that the standard form of the equations plus certain boundary conditions can be fully wrapped up in integral form and are exactly equivalent to the standard form.
The major difficulty is in how we produce the correct functional in the first place, especially if we want to avoid first deriving the differential form of the equations and back-calculating as we done above.
            
### Extension to Approximate Methods}   

The problem above is simple enough that the integral or differential forms of the equations can be solved directly. In general, however, we anticipate dealing with problems where no closed form of solution exists. Under these circumstances approximate solutions are desirable. In particular, there is a class of approximation methods which use families of trial functions to obtain a best fit approximation to the solution. These naturally develop into finite element algorithms as we shall soon see.

### Formulation of a General Problem

We consider a steady-state problem characterized by the following 
strong form

$$ 
{\cal{L}} (\phi) = f
$$     

where ${\cal{L}}$ is a linear differential operator acting on the 
(unknown) state variable $\phi$ in responce to a forcing function $f$.
Boundary conditions are

$$ 
    {\cal{B}}_i [\phi] = \left. q_i \right|_{\mbox{\small at boundary } S_i}    \;\;\; i=1,2,\ldots
$$

The operator should be symmetric 

$$ 
\int_\Omega v {\cal{L}}(u) d\Omega =  \int_\Omega u {\cal{L}}(v) d\Omega 
$$

and positive definite

$$ 
\int_\Omega u {\cal{L}}(u) d\Omega > 0
$$

$\Omega$ is the domain of the operator and $u$ and $v$ are any functions which satisfy the boundary conditions.  

![A rod subject to end load --- Young's modulus, $E$, density, $\rho$, cross sectional area, $A$][elastic-rod]

[elastic-rod]: ./Diagrams/rodpull.png


Consider the 1D example of a bar subject to a steady end load. The response is the solution to 

$$ 
    -EA\frac{\partial^2 u}{\partial x^2} = 0
$$ 

subject to the boundary conditions

$$ 
    \begin{split}
        \left. u \right|_{x=0} & = 0 \\
        \left. EA\frac{\partial u}{\partial x} \right|_{x=L} = R
    \end{split}
    
$$

We therefore identify

$$
 \begin{eqnarray*}
     {\cal{L}}= -EA\frac{\partial^2 u}{\partial x^2} & \phi = u & f = 0 \\
     & B_1 = 1 ; q_1 = 0 & \\
     & B_2 =  EA\frac{\partial }{\partial x}; q_2 = R & 
 \end{eqnarray*}
$$

To check symmetry and positive definiteness of the operator we consider $R=0$ since the operator properties are independent of the actual load. Integrating by parts gives

$$ 
    \begin{split}
    \int_0^L-EA\frac{\partial^2 u}{\partial x^2} v dx & =   
        - \left. EA\frac{\partial u}{\partial x} v \right|_0^L + \int_0^L EA \frac{\partial u}{\partial x} \frac{\partial v}{\partial x} dx \\
         &= - \left. EA\frac{\partial u}{\partial x} v \right|_0^L  - \left. + EA u \frac{\partial v}{\partial x}  \right|_0^L
         - \int_0^L EA\frac{\partial^2 v}{\partial x^2} u dx
    \end{split}  
$$

Application of boundary conditions demonstrates that the operator is symmetric by our definition and it
is positive definite as well because:

$$ 
\int_0^L-EA\frac{\partial^2 u}{\partial x^2} u dx =     
    - \left. EA\frac{\partial u}{\partial x} u \right|_0^L + \int_0^L EA \frac{\partial u}{\partial x} \frac{\partial u}{\partial x} dx         =
    0 + \int_0^L EA \left(\frac{\partial u}{\partial x}\right)^2 dx         
$$

Suppose we now search for approximate solutions of the form

$$ 
    \bar{\phi} = \sum_{i=1}^{n} a_i \Phi_i  
$$
where $\Phi_i$ are linearly independent trial functions and the $a_i$ are the unknown weights for each of the functions.

In weighted residuals methods, the expansion is used directly on the strong form of the equations. $\Phi_i$ are chosen so as to satisfy all boundary conditions and then we seek to minimize a residual

$$ 
    R = f - {\cal{L}}( \sum_{i=1}^{n} a_i \Phi_i    )
$$     
    
#### Least Squares Method

Minimize the square of the residual with respect to $a_i$

$$ 
    \frac{\partial}{\partial a_i} \int_\Omega R^2 d\Omega = 0   \;\;\; i = 1,2,\ldots
$$

This method produces a symmetric coefficient matrix regardless of the
properties of the operator. 
    
#### Galerkin Method
To determine $a_i$, solve the n equation system
$$ 
    \int_\Omega N_i R d\Omega = 0    \;\;\; i = 1,2,\ldots
 $$     
over the solution domain $\Omega$. This method produces a symmetric, positive definite coefficient matrix if the operator is symmetric and positive definite.
    
#### Ritz Method     
The Ritz method does not operate on the residual of the strong problem, but minimizes the weak form of the problem with respect to each of the unknown parameters $a_i$ in the usual variational manner. The trial functions no longer need satisfy the natural boundary conditions of the problem as these are wrapped up in the variational form. (Again, see Bathe for discussion and examples). 
    
The Galerkin method can be extended to include a term which minimizes the violation of natural boundary conditions, and thus permits the use of a wider range of trial functions

$$ 
    \int_\Omega N_i R d\Omega  + \int_\Gamma N_i R_B d\Gamma = 0     \;\;\; i = 1,2,\ldots
    \label{eq:galerk1}
$$     

However, it does not now necessarily produce a symmetric matrix even for a symmetric operator. However, if the equation (\ref{eq:galerk1} ) is integrated once by parts, it yields a symmetric form and also reduces the order of derivatives inside the integral. This means that the trial functions need be of lower order, and it makes the Galerkin formulation equivalent to the Ritz formulation.

Weighted residual formulations have one advantage: they can be used whether or not there exists a functional corresponding to the particular problem.  As a result is used this method is used extensively in constructing finite element methods.

## Towards a Matrix Method
    
One of the dominant features of the Finite Element literature is that it is filled with matrix algebra. The fact that differential equations can be rendered into matrices seems at first to be mysterious. However, it results quite naturally from the discretization of the problem, and the parameterization of the discrete equations through a limited set of unknown parameters.

Thus, before getting deeply involved in the arcane lore of Finite Elements, We give one example of a genuinely discrete system and show how it generates a matrix problem quite naturally.
    

![(a) A system of three carts interconnected by springs of different stiffnesses and in turn connected to an end wall. (b) The element equilibrium diagram for the spring $k_2$][mine-carts]

[mine-carts]: ./Diagrams/carts.png

    
Consider the system illustrated in the [figure (a) above](#mine-carts) --- three freely rolling carts attached by springs. There are three loads applied $R_1,R_2,R_3$, one to each cart, and we wish to determine the equilibrium displacements $U_1,U_2,U_3$.  We can [illustrate graphically](#mine-carts) the equilibrium condition for one of the springs based on its internal degrees of freedom and effective external load. This equilibrium is:

$$ 
\begin{split}
        k_2 (U_1 - U_2) &= {F_1}^{(2)} \\
        k_2 (U_2 - U_1) &= {F_2}^{(2)}
\end{split}
$$

or

$$ 

k_2 \left[  \begin{array}{cc}  1 & -1 \\ -1 & 1 \end{array} \right]
        \left[      \begin{array}{c} U_1 \\ U_2 \end{array} \right] = 
        \left[      \begin{array}{c} {F_1}^{(2)} \\ {F_2}^{(2)}\end{array} \right] 

$$

This is essentially the same for all five spring elements except for $k_1$ which is anchored at one end and so satisfies

$$ 
    k_1 U_1 = {F_1}^{(1)}
$$
 
The equilibrium relations for the system as a whole are 

$$ 
\begin{split}
    & {F_1}^{(1)} +             {F_1}^{(2)} +           {F_1}^{(3)} +           {F_1}^{(4)} = R_1 \\
    & {F_2}^{(2)} +             {F_2}^{(3)} +           {F_2}^{(5)} = R_2 \\
    & {F_3}^{(4)} +             {F_3}^{(5)} = R_3
\end{split}
\label{eq:globeq}
$$

If we now write all five equilibrium relations in terms of all available degrees of freedom we obtain a form like this

$$ 
k_2 \left[  \begin{array}{ccc}  k_2 & -k_2 & 0  \\ -k_2 & k_2 & 0 \\ 0 & 0 & 0 \end{array} \right]
        \left[      \begin{array}{c} U_1 \\ U_2 \\ U_3 \end{array} \right] = 
        \left[      \begin{array}{c} {F_1}^{(2)} \\ {F_2}^{(2)} \\ 0\end{array} \right] 
$$

which can also be written

$$ 
    \mathbf{K}^{(2)}\mathbf{U} = \mathbf{F}^{(2)}
$$

in each case.

Thus the global equilibrium requirement of \ref{eq:globeq} become

$$ 
    \mathbf{K}\mathbf{U} = \mathbf{R}
 $$

where

$$ 
    \mathbf{K} = \left[
    \begin{array}{ccc}
        (k_1 + k_2 +k_3 + k_4) & -(k_2+k_3) & -k_4 \\
        -(k_2+k_3) & (k_2 + k_3+k_4) & -k_5 \\
        -k_4 & -k_5 & (k_4 + k_5)
    \end{array}
    \right]
 $$ 

Note, by the way, that $ \mathbf{K}$ is symmetric. An important observation is that 

$$ 
    \mathbf{K} = \sum_{i=1}^{5} \mathbf{K}^{(i)}
 $$

The individual element stiffnesses can be summed to form a global stiffness matrix. A symmetric, positive definite matrix problem can be solved in numerous different ways, many of which are easy to look up in textbooks !
A very simple problem like this has captured much of what we need to do in arbitrarily complex finite element computations. The construction of local element equilbrium problems based on the interaction of each available degree of freedom with every other is followed by the assembly into a global problem by summing the contributions of the individual elements.

Note that the formulation of the local equilibrium conditions are done in a symmetric manner (if we change this degree of freedom how does the balance change, and then what if we change this degree of freedom ?) rather than trying to simplify the system.

The degrees of freedom are related by elastic spring constants here. In our more abstract formulations we will replace spring constants by coefficients obtained from the variational method but the form is {\em exactly the same}.


%% --- 

