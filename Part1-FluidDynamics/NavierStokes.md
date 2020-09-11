

# Navier-Stokes Equation

In Chapter ??. we ... 


<!-- Constitutive Laws
================= -->

The formulation above is quite general, and can be extended where
necessary to include conservation laws for additional physical variables
(for example, angular momentum, electric current). Specific to the type
of material which is deforming is the constitutive law which describes
the stress. In the case of an incompressible fluid, the stress is
related to strain rate through a viscosity, $\eta$, and to the pressure,
$p$. Incompressibility is expressed as $$\nabla \cdot \mathbf{u} = 0$$

This is a tighter constraint than mass conservation and emerges from
equation ([\[eq:masscons\]](#eq:masscons){reference-type="ref"
reference="eq:masscons"}) when $\rho$ is assumed to be constant. With
this assumption, the constitutive law is written [ (A derivation for
this is in Landau & Lifschitz). ]{style="color: 0.0,0.7,0.3"}
$$\sigma_{ij} = \eta \left( \frac{\partial v_i}{\partial x_j} + \frac{\partial v_j}{\partial x_i}\right) - p\delta_{ij}$$

If the viscosity is constant, then we can substitute the constitutive
law into the stress-divergence term of the momentum conservation
equation. In index notation once again, $$\begin{split}
                \nabla \cdot \boldsymbol{\sigma} & =
                    \frac{\partial}{\partial x_j} \eta 
                        \left( \frac{\partial v_i}{\partial x_j} +  \frac{\partial v_j}{\partial x_i} \right)
                        - \frac{\partial p}{\partial x_j} \delta_{ij} \\
                & = \eta \frac{\partial^2 v_i}{\partial x_j \partial x_j} +
                            \eta \frac{\partial^2 v_j}{\partial x_i \partial x_j} -
                                \frac{\partial p}{\partial x_i}\\
                & = \eta \nabla^2 \mathbf{v} +
                             \eta \textcolor[rgb]{0.0,0.7,0.0}{ \nabla (\nabla \cdot \mathbf{v})} - \nabla p       
            \end{split}$$

the second term in this final form must vanish because of the
incompressibility assumption, so the momentum conservation equation
becomes

$$\textcolor[rgb]{0.7,0.0,0.0}{   \rho \left( \frac{\partial \mathbf{v}}{\partial t}
                            + (\mathbf{v} \cdot \nabla) \mathbf{v} \right) =
                            \eta \nabla^2 \mathbf{v} - \nabla p 
                            - g\rho\hat{\mathbf{z}}      }
                \label{eq:navstokes}$$

This is the Navier-Stokes equation. Once again some new notation has
shown up uninvited. The Laplacian operator $\nabla^2$ is defined (in a
scalar context) as $$\begin{split}
            \nabla^2 \phi & = \nabla \cdot \nabla \phi \\
                                    & = \frac{\partial^2 \phi}{\partial x^2} + 
                                            \frac{\partial^2 \phi}{\partial y^2} +  
                                            \frac{\partial^2 \phi}{\partial z^2}  \text{\hspace{1cm} (Cartesian)}
            \end{split}$$

and in a vector context as $$\begin{split}
                \nabla^2 \mathbf{u} & = \nabla \nabla \cdot \mathbf{u} - \nabla \times (\nabla \times \mathbf{u}) \\
                                                    & = \mathbf{i} \nabla \cdot \nabla u_x + 
                                                            \mathbf{j} \nabla \cdot \nabla u_y + 
                                                            \mathbf{k} \nabla \cdot \nabla u_z   \text{\hspace{1cm} (Cartesian)}
            \end{split}$$

In Cartesian coordinates, the Laplacian operator has a simple form, and
the vector Laplacian is simply the scalar operator applied in each
direction. In other coordinate systems this operator becomes
substantially more elaborate.

Boussinesq Approximation, Equation of State, Density Variations
===============================================================

The equation which relates pressure, temperature and density is known as
the equation of state. For the equations derived so far, we have
specified an incompressible fluid, for which no density variations are
possible. However, we have also included a source term for momentum
which relies on gravity acting on *density variations*.

This conflict is typical of fluid mechanics: simplifying assumptions if
taken to their logical limit imply no motion or some other trivial
solution to the equations.

In this case, we make the assumption that density changes are typically
small relative to the overall magnitude of the density itself. Terms
which are scaled by density can therefore assume that it is a large
constant value. Terms which contain gradients of density or density
variations should consider the equation of state. This is the Boussinesq
approximation and is only a suitable appropriation for nearly-
incompressible fluids. In the Navier-Stokes equation, the hydrostatic
pressure does not influence the velocity field at all. Only
*differences* in density drive fluid flow, and so the sole term in which
density needs to be considered variable is that of the gravitational
body forces.

In the case of density variations due to temperature, the equation of
state is simply $$\rho = \rho_0 \left(1 - \alpha ( T-T_0 )\right)
            \label{eq:state}$$

where $\rho_0$ is the density at a reference temperature $T_0$. $\alpha$
is the coefficient of thermal expansion. It is generally much smaller
than one, making the Boussinesq approximation a reasonable choice.

The energy and momentum conservation equations thus become coupled
through the term
$$g\rho\hat{\mathbf{z}} = g \rho_0 \left(1 - \alpha(T-T_0)\right)$$

Density variations due to pressure produce a perfectly vertical,
isotropic forcing term on the momentum conservation equation. In the
steady state case, this is balanced by the hydrostatic pressure
gradient. (The isotropic term does not contribute at all the the
deviatoric part of the stress equation and thus cannot induce steady
flow). We therefore ignore the vertical density gradient due to the
fluid overburden.

Another density variation is that which results from variation in
chemical composition from one fluid element to another. This is the case
where two immiscible fluids live in the same region. Now density
variations might be large -- the fluid domains must be considered
separately.

Advection and the Lagrangian Formulation
----------------------------------------

We have seen the $\mathbf{v} \cdot \nabla$ operator a number of times
now. The presence of this term causes major difficulties in continuum
mechanics since it introduces a strong non-linearity into the momentum
equation. It is this term which produces turbulence in high speed flows
etc.

This term is the 'advection' term which accounts for the passive
transport of information (temperature, momentum, by the motion of the
fluid. Advection also presents some serious headaches in numerical
methods and has spawned entire literatures devoted to efficient and
accurate solution methods.

One obvious way to avoid the problem of advection is to consider an
elemental volume of space which *moves with the fluid*. The surface flux
term from equation ([\[eq:cons1\]](#eq:cons1){reference-type="ref"
reference="eq:cons1"}) vanishes immediately

Mathematically, we introduce a new notation (of course), as follows:
$$\frac{D \phi}{D t} = \frac{d}{dt} \phi[x_1(t),x_2(t),x_3(t),t]$$

where the change in the reference position $(x_1(t),x_2(t),x_3(t))$ is
governed by the local flow velocity:
$$\frac{d x_1}{d t} = v_1 \mbox{\hspace{1cm}}
                    \frac{d x_2}{d t} = v_2 \mbox{\hspace{1cm}}
                    \frac{d x_3}{d t} = v_3$$

which keeps the reference point moving with the fluid. Differentiating
gives $$\begin{aligned}
                \frac{D \phi}{D t} &= \frac{\partial \phi}{\partial t}
                                    \frac{\partial \phi}{\partial x_1}\frac{d x_1}{d t} + 
                                    \frac{\partial \phi}{\partial x_2}\frac{d x_2}{d t} +
                                    \frac{\partial \phi}{\partial x_3}\frac{d x_3}{d t} \nonumber \\
        \intertext{and so leads to}
                \frac{D \phi}{D t} &=   \frac{\partial \phi}{\partial t}
                                    v_1 \frac{\partial \phi}{\partial x_1} + 
                                    v_2 \frac{\partial \phi}{\partial x_2} +
                                    v_3 \frac{\partial \phi}{\partial x_3}  \nonumber \\
        \intertext{which is equivalent to}
                \frac{D \phi}{D t} &= \frac{\partial \phi}{\partial t} + (\mathbf{v} \cdot \nabla) \phi                         
        \end{aligned}$$

If we think of $\phi$ as the concentration of a dye in the fluid, then
the above is a conservation equation assuming the dye does not diffuse
and has no sources or sinks.

Viewed from a reference frame locked to *a particular fluid element*,
the energy conservation equation becomes $$% \rho  ??
                \frac{D T}{Dt} =
                        \kappa \nabla^2 T + \frac{H}{C_p}$$

and the momentum conservation equation now becomes $$\rho %% ?
            \frac{D \mathbf{v} }{D t} =
                            \eta \nabla^2 \mathbf{v} - \nabla P     
                            - g\rho\hat{\mathbf{z}}$$

This is a considerably more compact way of writing the equations, but we
have only really succeeded in hiding the nasty term under the rug, since
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



Reading Material
================

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
