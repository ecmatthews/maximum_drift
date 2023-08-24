###############
maximum_drift calculator
###############

Calculate the maximum drift of a host star from a companion with specified separation, using some simple orbital assumptions.

Based on the method described in Montagnier 2008 (PhD thesis), see https://www.theses.fr/2008GRE10289

`maximum_drift` takes a csv file of nearby Gaia objects, and determines whether each is bound, unclear or unbound based on the relative parallax and color of the objects. The maximum_drift that can be caused by each companion is then calculated, and compared to the observed linear drift of the host star.
