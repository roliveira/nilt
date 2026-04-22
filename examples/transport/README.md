# Transport - Diffusion Examples

## 1. Sphere diffusion (`sphere_diffusion.cpp`)

Diffusion from a solid sphere of radius *a* into an infinite surrounding
medium. Tracks the average concentration inside the sphere.

**Laplace domain:**

$$\bar{C}_{\text{avg}}(s) = \frac{C_0}{s}\,\frac{3}{q}\left[\coth q - \frac{1}{q}\right], \qquad q = a\sqrt{\frac{s}{D}}$$

**Analytical (series):**

$$C_{\text{avg}}(t) = C_0\,\frac{6}{\pi^2}\sum_{n=1}^{\infty}\frac{1}{n^2}\exp\!\left(-\frac{Dn^2\pi^2 t}{a^2}\right)$$

**Parameters:** $D = 10^{-9}$ m$^2$/s, $a = 1$ mm, $C_0 = 1$ mol/m$^3$.

**Reference:** Crank (1975), §6.3.

## 2. Cylinder diffusion - axisymmetric (`cylinder_diffusion.cpp`)

Diffusion from a long solid cylinder of radius *a* into an infinite medium.
This is a 2D (radial) problem; axial symmetry reduces it to one spatial
variable.

**Laplace domain:**

$$\bar{C}_{\text{avg}}(s) = \frac{C_0}{s}\,\frac{2}{q}\,\frac{I_1(q)}{I_0(q)}, \qquad q = a\sqrt{\frac{s}{D}}$$

where $I_0$, $I_1$ are modified Bessel functions of the first kind.

**Analytical (series):**

$$C_{\text{avg}}(t) = 4\,C_0 \sum_{n=1}^{\infty}\frac{1}{\alpha_n^2}\exp\!\left(-\frac{D\alpha_n^2 t}{a^2}\right)$$

where $\alpha_n$ are the positive zeros of $J_0(\alpha) = 0$.

**Parameters:** $D = 5 \times 10^{-10}$ m$^2$/s, $a = 0.5$ mm, $C_0 = 1$ mol/m$^3$.

**Reference:** Crank (1975), §5.3.

## 3. Advection-diffusion plume - **asymmetric** (`advection_plume_2d.cpp`)

Instantaneous point release of mass M in a uniform flow field (velocity v in
the x-direction) with isotropic diffusion D. The plume is elongated
downstream - no radial symmetry.

**Laplace domain** (via Galilean substitution):

$$\bar{C}(x,y,s) = \frac{M}{2\pi D}\,\exp\!\left(\frac{vx}{2D}\right)K_0\!\left(r\sqrt{\frac{v^2}{4D^2}+\frac{s}{D}}\right)$$

**Analytical solution:**

$$C(x,y,t) = \frac{M}{4\pi Dt}\,\exp\!\left(-\frac{(x-vt)^2 + y^2}{4Dt}\right)$$

**Parameters:** $D = 0.1$ m$^2$/s, $v = 0.5$ m/s, $M = 1$ kg/m.

**Reference:** Bear (1972), §10.6.

## Running

```bash
# from the examples/transport/build directory
./transport_sphere_diffusion
./transport_cylinder_diffusion
./transport_advection_plume_2d

# plot results (reads/writes from build/ by default)
python ../plot_sphere_diffusion.py
python ../plot_cylinder_diffusion.py
python ../plot_advection_plume_2d.py
```
