

# Navier-Stokes Equation

In the previous section ({doc}`./ContinuumMechanics`), we looked a very general way use conservation principles to construct equations describing how sources, fluxes and material transport of varous quantities balance one another through time. This approach can also be used to generate conservation equations for additional physical variables
(for example, angular momentum, electric current). 

One thing is still missing, though, and that is how we express the flux term for a quantity in terms of gradients in the quantity itself. We did this without any comment in the conservation of heat energy because we used Fourier's law  $\mathbf{F} = -k \nabla T$ to express the heat flux as a linear function of the temperature gradient. This is an empirical relationship that contains a physical constant that we are not able to derive from first principles.
This missing relationship is known as a **constitutive law**.

## Viscosity

For a fluid, the constitutive relationship relating momentum flux to momentum gradients introduces a new physical parameter, **viscosity** which relates stress to rates of change in velocity gradients and a second viscosity relating to rates of volume change. 


If we make the assumption that the fluid is incompressible ($\nabla \cdot \mathbf{v}$) which we already mentioned in deriving the temperature equation, we only have to worry about the viscosity that relates to velocity gradients 

$$
\sigma_{ij} = \eta \left( \frac{\partial v_i}{\partial x_j} + \frac{\partial v_j}{\partial x_i}\right) - p\delta_{ij}
$$

(An expanded derivation for
this is in Landau & Lifschitz[REF] and they also talk extensively about rate of volume change effects). 

%TODO - Landau / Lifshitz reference etc from below

If the viscosity is constant, then we can substitute the constitutive
law into the stress-divergence term of the momentum conservation
equation. In index notation once again, 

$$
\begin{split}
                \nabla \cdot \boldsymbol{\sigma} & =
                    \frac{\partial}{\partial x_j} \eta 
                        \left( \frac{\partial v_i}{\partial x_j} +  \frac{\partial v_j}{\partial x_i} \right)
                        - \frac{\partial p}{\partial x_j} \delta_{ij} \\
                & = \eta \frac{\partial^2 v_i}{\partial x_j \partial x_j} +
                            \eta \frac{\partial^2 v_j}{\partial x_i \partial x_j} -
                                \frac{\partial p}{\partial x_i}\\
                & = \eta \nabla^2 \mathbf{v} +
                             \eta \textcolor[rgb]{0.0,0.7,0.0}{ \nabla (\nabla \cdot \mathbf{v})} - \nabla p       
\end{split}
$$

the second (green) term in this final form must vanish because of the
incompressibility assumption, so the momentum conservation equation
becomes

```{math}
:label: eq:navstokes
\textcolor[rgb]{0.7,0.0,0.0}{   \rho \left( \frac{\partial \mathbf{v}}{\partial t}
                            + (\mathbf{v} \cdot \nabla) \mathbf{v} \right) =
                            \eta \nabla^2 \mathbf{v} - \nabla p 
                            - g\rho\hat{\mathbf{z}}      }
```

This is the Navier-Stokes equation. Once again some new notation has
shown up uninvited. The Laplacian operator $\nabla^2$ is defined (in a
scalar context) as 

```{math}
\begin{split}
            \nabla^2 \phi & = \nabla \cdot \nabla \phi \\
                                    & = \frac{\partial^2 \phi}{\partial x^2} + 
                                            \frac{\partial^2 \phi}{\partial y^2} +  
                                            \frac{\partial^2 \phi}{\partial z^2}  
                                            \mathrm{\hspace{1cm} (Cartesian)}
\end{split}
```
and in a vector context as 

``` {math}
\begin{split}
                \nabla^2 \mathbf{u} & = \nabla \nabla \cdot \mathbf{u} - \nabla \times (\nabla \times \mathbf{u}) \\
                                                    & = \mathbf{i} \nabla \cdot \nabla u_x + 
                                                            \mathbf{j} \nabla \cdot \nabla u_y + 
                                                            \mathbf{k} \nabla \cdot \nabla u_z   
                                                            \mathrm{\hspace{1cm} (Cartesian)}
\end{split}
```

In Cartesian coordinates, the Laplacian operator has a simple form, and
the vector Laplacian is simply the scalar operator applied in each
direction. In other coordinate systems this operator becomes
substantially more complicated but this notation is handy because the
equations always look the same no matter the coordinate system.

## Boussinesq Approximation, Equation of State, Density Variations

The equation which relates pressure, temperature and density is known as
the equation of state. For the equations derived so far, we have
specified an incompressible fluid and that means the density 
can't change at all. 

Can you see the problem ?  We have also included a source term for momentum
which relies on gravity acting on *density variations*. 
This conflict is typical in fluid mechanics: 
simplifying assumptions if taken to their logical limit imply
only trivial solutions (no motions are possible, that sort of thing)

We have to do something to get around this problem and the simplest thing
to do is assume density *changes* are typically
small relative to the magnitude of the density.
Any terms in the equations 
which are multiplied by density can therefore assume that density is a large,
constant value. Terms which contain gradients of density or density
variations (where the constant terms cancel out) need to consider 
the equation of state. 

We call this the *Boussinesq
approximation* and it is only good for nearly-
incompressible fluids. In the Navier-Stokes equation, the hydrostatic
pressure does not influence the velocity field at all. Only
*differences* in density drive fluid flow, and so the only term
where we need to worry about density not being constant is 
in the gravitational body forces.

In the case of density variations due to temperature, the equation of
state is this:

```{math}
:label: eq:boussinesq-state
\rho = \rho_0 \left(1 - \alpha ( T-T_0 )\right)
            \label{eq:state}
```
where $\rho_0$ is the density at a reference temperature $T_0$. $\alpha$
is the coefficient of thermal expansion. $\alpha$ is generally *much* smaller
than one which justifies making the approximation in the first place.


The energy and momentum conservation equations are now coupled
through the buoyancy forces,

$$
        g\rho\hat{\mathbf{z}} = g \rho_0 \left(1 - \alpha(T-T_0)\right)
$$

Density variations due to pressure produce a perfectly vertical,
isotropic forcing term on the momentum conservation equation. In the
steady state case, this is balanced by the hydrostatic pressure
gradient. (The isotropic term does not contribute at all the the
deviatoric part of the stress equation and thus cannot induce steady
flow). This is the reason we can ignore the vertical density gradient due to the
weight of the fluid.

Another density variation is the one that results from variation in
chemical composition from one fluid element to another. An example
would be 
where two immiscible fluids live in the same region. Now density
variations might be large -- the fluid domains must be considered
separately if we still want to assume small density changes due to 
temperature.

## Fluid Motion: Eulerian v. Lagrangian

We have seen the $\mathbf{v} \cdot \nabla$ operator a number of times
now. The presence of this term causes major complications when we
want to find solutions to the equations of continuum
mechanics since it introduces a strong non-linearity.
It is this term which is behind turbulence in high speed flows.

This term describes the **advection** of a quantity which accounts for the passive
transport of something (temperature, momentum), by the motion of the
fluid. Advection also presents some serious headaches in numerical
methods and has spawned a vast  literature devoted to efficient and
accurate solution methods.

One way to side-step  
the complexity of advection is to consider an
elemental volume of space which *moves with the fluid*. The surface flux
term from equation {eq}`eq:cons1`  vanishes at once if there is no
motion of material in or out of the boundary. In reality, we 
are just sweeping the problem under the rug -- it doesn't actually
go away altogether -- but it can be useful to do this
because every once in a while that term will cancel out before we need to 
actually do anything with it ! And, sometimes, when it comes to numerical methods,
different views of the equations can lead to better approximations.

Mathematically, we introduce a new notation:

$$
\frac{D \phi}{D t} = \frac{d}{dt} \phi\left[ x_1(t),x_2(t),x_3(t),t \right]
$$

where the change in the reference position $(x_1(t),x_2(t),x_3(t))$ is
governed by the local flow velocity:

$$
\frac{d x_1}{d t} = v_1 \;\;\;\;
                    \frac{d x_2}{d t} = v_2  \;\;\;\;
                    \frac{d x_3}{d t} = v_3
$$

and this keeps the reference point moving with the fluid. 
Differentiating gives 

$$
        \frac{D \phi}{D t} =    \frac{\partial \phi}{\partial t} +
                                \frac{\partial \phi}{\partial x_1}\frac{d x_1}{d t} + 
                                \frac{\partial \phi}{\partial x_2}\frac{d x_2}{d t} +
                                \frac{\partial \phi}{\partial x_3}\frac{d x_3}{d t} 
$$

which leads to 

$$
                \frac{D \phi}{D t} =   \frac{\partial \phi}{\partial t} +
                                    v_1 \frac{\partial \phi}{\partial x_1} + 
                                    v_2 \frac{\partial \phi}{\partial x_2} +
                                    v_3 \frac{\partial \phi}{\partial x_3}  
$$

which is equivalent to

$$
                \frac{D \phi}{D t} = \frac{\partial \phi}{\partial t} + 
                (\mathbf{v} \cdot \nabla) \phi                         
$$

If we think of $\phi$ as the concentration of a dye in the fluid, then
the above is a conservation equation assuming the dye does not diffuse
and has no sources or sinks.

Viewed from a reference frame locked to *a particular fluid element*,
the energy conservation equation becomes 

$$
        \frac{D T}{Dt} =
        \kappa \nabla^2 T + \frac{H}{C_p}
$$

and the momentum conservation equation now becomes 

$$
        \frac{D \mathbf{v} }{D t} =
                \eta \nabla^2 \mathbf{v} - \nabla P     
                - g\rho\hat{\mathbf{z}}
$$

This is a considerably more compact way of writing the equations, but
it is now necessary to use a coordinate system which is locked into the
fluid and rapidly deforms as the fluid flows. Before long, the
coordinate system is unimaginably complex -- the advection problem
returns in another guise. This formulation is known as the Lagrangian
formulation and contrasts with the Eulerian viewpoint which is fixed in
space.

From the numerical point of view, however, this approach can have some
distinct advantages. The computer can often track the distorted
coordinate system far better than it can handle successive applications
of the $\mathbf{v} \cdot \nabla$ operator at a fixed point in space. We
will return to this point later.

## Reading Material

G. F. Davies. *Dynamic Earth*. Cambridge University Press, New York,
1999.

J. Grotzinger, T. H. Jordan, F. Press, and R. Siever. *Understanding
Earth*. W. H. Freeman & Co, 5 edition, 2006. URL
<http://bcs.whfreeman.com/understandingearth5e/>.

L. D. Landau and E. M. Lifshitz. *Fluid Mechanics*, volume 6 of *Course
of Theoretical Physics*. Pergamon Press, 1959.

G. Schubert, D. L. Turcotte, and P. Olson. *Mantle Convection in Earth
and Planets*. Cambridge University Press, UK, 2001.

D.L. Turcotte and G. Schubert. *Geodynamics*. John Wiley and Sons, New
York, 1982.
