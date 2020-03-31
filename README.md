# SELAM
Simulation of Epistasis Local adaptation, during Admxiture with Mate choice

￼￼# Recent additions: 3-31-2020

SELAM now includes a simple gene conversion model. In the current version this is specified as:

--gc [float, fraction] [float, length]

That is, we specify the fraction of crossovers that resolve as gene conversion events, and the mean length (in Morgans) of the GC tract. As currently specified, this means that the GC landscape is fixed with respect to crossovers. This may not be biologically realistic so please be careful. This current version is exploratory and should be used cautiously. 
