# maximum_drift calculator

Calculate the maximum drift of a host star from a companion with specified separation, using some simple orbital assumptions.

``maximum_drift`` takes a csv file of nearby Gaia objects, and determines whether each is bound, unclear or unbound based on the relative parallax and color of the objects. The ``maximum_drift`` that can be caused by each companion is then calculated, and compared to the observed linear drift of the host star.

### drift calculation

Based on the method described in Montagnier 2008 (PhD thesis), see https://www.theses.fr/2008GRE10289 - using formula (10) on page 78 to reproduce figure (11) on page 79. Briefly, the RV drift of a companion is ``dz``. For a companion observed with a projected separation ``x``, the true separation is:
$\rho = \sqrt{(x^2+z^2)}$.
The acceleration magnitude is a function of the true separation: for any true separation >= x, we can calculate the expected acceleration, and the component in the dz direction, assuming a circular orbit:

$\frac{d^2z}{dt^2} = \frac{Gm_2}{\rho^2}\cos(\arcsin(\frac{x}{r}))$.

See Montagnier + for full derivation.

### mass calculation

This code uses the Pecaut & Mamajek (2013, ApJS, 208, 9; http://adsabs.harvard.edu/abs/2013ApJS..208....9P) isochrones to convert Gaia magnitudes into masses (assuming distance == distance of the host star). See https://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt
