In order to test the correct calculation of radiative transfer, the theory must first be covered. The radiative transfer (RT) equation,

(1)
$$
dI_\nu = - I_\nu \kappa_\nu ds + \epsilon_\nu ds
$$

must be integrated to calculate the observed intensity. Simple applications of RT assume constant $\epsilon_\nu$ and $\kappa_\nu$ between voxels to obtain the observed intensity. The benefit of `kosmatau3d` is the improved calculation of the RT equation. Here the emissivity and absorption terms, $\epsilon_\nu$ and $\kappa_\nu$, are still used, but the first-order linear terms, $\Delta\epsilon_\nu = \frac{\epsilon_\nu}{\Delta s}$ and $\Delta\kappa_\nu = \frac{\kappa_\nu}{\Delta s}$, also appear. Here $\Delta s$ is the width of each voxel. This leads to,

(2)
$$
I_\nu = e^{-\kappa_\nu \Delta s - \frac{1}{2} \Delta \kappa_\nu \Delta s^2} \left[ \int_0^{\Delta s}\left( \epsilon_\nu + \Delta \epsilon_\nu s' \right) e^{\kappa_\nu s' + \frac{1}{2} \Delta \kappa_\nu s'^2} ds' + I_{bg} \right],
$$

and allows us to define specific cases for the way we integrate the RT equation:

 - *Basic*: $\Delta\kappa = 0$ and $| \kappa \Delta s | < 10^{-10}$
 - *Simple*: $\kappa > 10^3 \Delta \kappa \Delta s$
 - *Complex*: otherwise
   - *method 1*: $\Delta \kappa > 0$
   - *method 2*: $\Delta \kappa < 0$
   

The complex methods of calculating the RT equation are the most accurate, and should be used by default. The other methods are simplifications that can be made. The *Basic* case occurs when there is (mostly) no absorption. Since $\kappa_\nu = 0$ and $\Delta\kappa_\nu = 0$, the RT equation becomes, 

$$
\int_0^{I_\nu} dI_\nu = \int_0^{\Delta s} \left( \epsilon_\nu + \Delta\epsilon_\nu s' \right) ds', \\
I_\nu = \epsilon_\nu \Delta s + \frac{1}{2} \Delta\epsilon_\nu \left( \Delta s \right)^2 + I_{\nu,0}.
$$

The *Simple* case occurs when there is constant absorption, so $\Delta\kappa_\nu = 0$ and $\kappa_\nu \neq 0$. This has a more-complex equation:

$$
\int_0^{I_\nu} dI_\nu = e^{-\kappa_\nu \Delta s} \left( \int_0^{\Delta s} (\epsilon_\nu + \Delta \epsilon_\nu s') e^{\kappa_\nu s'} ds' + I_0 \right), \\
I_\nu = \left( \frac{\epsilon_\nu \kappa_\nu + \Delta\epsilon_\nu (\kappa_\nu \Delta s - 1)}{\kappa_\nu^2} \right) - e^{-\kappa_\nu \Delta s} \left( \frac{\epsilon_\nu \kappa_\nu - \Delta\epsilon_\nu}{\kappa_\nu^2} \right) + e^{-\kappa_\nu \Delta s} I_{\nu,0}.
$$

Finally, the *Complex* integration is the full integration of the RT equation:

$$
I_\nu = \tilde{I_\nu} + I_{bg}, \\
\tilde{I_\nu} = \frac{\Delta \epsilon_\nu}{\Delta \kappa_\nu} \left[ 1 - e^{-\kappa_\nu \Delta s - \frac{1}{2} \Delta \kappa_\nu \Delta s^2} \right] - \frac{\epsilon_\nu \Delta \kappa_\nu - \kappa_\nu \Delta \epsilon_\nu}{\sqrt{\Delta \kappa_\nu}^3} \sqrt{\frac{\pi}{2}} e^{-\frac{(\kappa_\nu + \Delta \kappa_\nu \Delta s)^2}{2 \Delta \kappa_\nu}} \left[ \mathrm{erfi} \left( \frac{\kappa_\nu}{\sqrt{2 \Delta \kappa_\nu}} \right) - \mathrm{erfi} \left( \frac{\kappa_\nu + \Delta \kappa_\nu \Delta s}{\sqrt{2 \Delta \kappa_\nu}} \right) \right].
$$

Here erfi() is the imaginary error function. The only distinction made is if $\Delta\kappa \lessgtr 0$.

$$
I_\nu = \frac{\Delta\epsilon_\nu}{\Delta\kappa_\nu} \left( 1 - e^{-\kappa_\nu \Delta s - \frac{1}{2} \Delta\kappa_\nu (\Delta s)^2} \right) - \frac{\epsilon_\nu \Delta\kappa_\nu - \kappa_\nu \Delta\epsilon_\nu}{\Delta\kappa_\nu} \sqrt{\frac{\pi}{2|\Delta\kappa_\nu|}} \left( e^{a_\nu^2 - b_\nu^2}\tilde{E}(a_\nu) - \tilde{E}(b_\nu) \right) + I_{\nu,0} e^{-\kappa_\nu \Delta s - \frac{1}{2} \Delta\kappa_\nu (\Delta s)^2}
$$

Here $a_\nu=\frac{\kappa_\nu}{\sqrt{2\Delta\kappa_\nu}}$ and $b_\nu=\frac{\kappa_\nu + \Delta\kappa_\nu \Delta s}{\sqrt{2 \Delta\kappa_\nu}}$, while $\tilde{E}$ is an approximation of the imaginary error function. The distinction we make is how we treat the error function if $a,b$ are real or imaginary values. This is what I believe to be the primary issue in the calculation of the RT equation: the handling of complex values. I need to perform the integration to properly understand how to treat the complex numbers, namely in the exponential function of $a$ and $b$. As seen below, this might be causing the erroneous negative observed intensities.