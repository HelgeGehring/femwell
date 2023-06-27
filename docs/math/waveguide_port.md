# Waveguide port boundary conditions

Here we introduce the math necessary for waveguide port boundary conditions, which launch a certain mode into the waveguide and absorb the reflections{cite}`Jin2015`

Let's start with the simplest case: a parallel-plate waveguide. In this case, the modes can simply be described by

$$
h_m(y) = \sqrt{\frac{v_m}{b}}\cos\frac{m\pi y}{b},
\quad
v_m = 
\begin{cases}
    1, &m=0 \\
    2, &m \neq 0
\end{cases}
$$

where $b$ is the width of the waveguide. 

## Bibliography

```{bibliography}
:style: unsrt
:filter: docname in docnames
```