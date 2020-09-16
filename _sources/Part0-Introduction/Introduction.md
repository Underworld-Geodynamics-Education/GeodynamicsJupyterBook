# Introduction 

Our goal is to develop an quantitative understanding of the dynamic
processes within the Earth and sister planets as I have sketched in the
Figure to the right. These dynamic processes are largely driven by the
internal heat of the planet escaping to the surface through whatever
mechanisms are available. Some of the heat is left over from the
original formation of the planet, and the rest originates in the decay
of radioactive elements. In the Earth's early history and elsewhere in
the solar system, tidal heating, chemical segregation, and impacts have
all played a role in supplying the interior heat budget. 


```{figure} Diagrams/EarthProcessesPlume.png
---
width: 100mm
name: mantle-slice
---
A schematic of the
Earth's interior and something similar for Venus and Mars would be, on
the face of it, much simpler
```


This is because the Earth's dynamics is completely dominated by Plate
Tectonics --- a unique manifestation of interior heat loss as far as we
are aware. Part of our task is to understand why plate tectonics is a
possible outcome of a hot planet, but also why it is not the only
possible outcome. If we can also understand how the different modes are
selected for planets of different size, composition, and heat budget,
then we have a powerful way to predict the geological behaviour of
extrasolar planets. Plate tectonics creates a number of very efficient
cycling mechanisms which link the interior of the Earth and the
Atmosphere and Oceans; it may prove to be an essential ingredient for
the kind of friendly world we expect to be needed to nurture
(intelligent) life.

## What is a model ?

Global scale geodynamics is a discipline where we cannot do controlled
experiments on the basic processes we are studying. We rely on observing
the Earth and the other terrestrial planets and moons and looking for
multiple manifestations of the same processes under different conditions
to give us control on certain parameters.

While it is not possible to do experiments at the planetary scale over
geological time, it is possible to perform experiments at a physically
manageable size and, by careful scaling, to generalise the results to
geologically relevant space and time-scales. If these processes of
interest can be understood through a mathematical description, then the
equations are automatically applicable at geological time and space
scales provided the assumptions which go into developing the
mathematical model are still valid.

We will frequently be talking about "modeling" --- people mean many
different things by this and all of the following fall under the general
concept of modeling:

-   Laboratory based physical models which can be scaled to give
    meaningful, quantitative insight into deformation at geological
    scales.

-   The building of mathematical descriptions of the world and their use
    to approximate physical \"reality\".

-   Computational solution of these descriptions (where needed) and the
    concept of a numerical model.

-   How to go about constructing a model

-   How to go about using a model (these are quite different things !).

-   Exploring parameter variation to understand the dominant effects.

## The mathematics you need

These notes contains material at different levels. There is broadly
descriptive content which is intended to introduce the subject and lead
up to the more advanced mathematical content. It should be possible to
follow these notes without detailed knowledge of how the mathematical
results are obtained, but it is expected that you can use the results in
exercises for to solve real problems in geophysics.

Familiarity with vector & tensor notation & index notation is required
for understanding the kinds of equations we will be dealing with and I
interchange them a little. On the assumption that these things are
disconcerting (at best) the first few times, I tend to write everything
out in full -- at least for the Cartesian case. In other geometries,
vector notation generally still holds, but the definition of the
operators can be very different, and it is always worth checking before
using them.

## Programming skills

All of the examples that we use in this book are based around `python` code
which is easy to pick up and also very powerful. Many python codes are actually
`fortran` or `C` codes in disguise with python providing a more convenient
way to operate and coordinate than stacks of punched cards or text-based input files.

Python is a very flexible language which is constantly improving (some would say).
The flexibility can make it hard to know where to start or how to find 
documentation on the function or module you are using. Jupyter notebooks add
interactivity to the language and make it much more straightforward to 
dig into the details of the code you are running. We use Jupyter notebooks extensively
for examples in this book and these can be run live.
