# W(35,25)
This repository contains weighing matrices which were discovered by Christopher Munro between August 12th 2026 and 5th September 2026. The construction of these matrices will be presented in a forthcoming paper.

Weighing matrices $W(n,k)$ are matrices of the form:

$$WW^\top=k I$$

where $I$ is an $n\times n$ identity matrix. They are special objects studied in the theory of combinatorial designs. The existence of a weighing matrix $W(35,25)$ has been listed as open in recent catalogues. This repository contains a number of different weighing matrices from distinct H-classes and various ranks taken over $\mathbb{F}_5$.

The classes I and II were the first to be found by an OR-Tools CP-SAT solve for the upper left $30\times 30$ block matrix, imposing a symmetry satisfying $W=PWP^\top$, fixing the border, and choosing the lower right constant block to be all 1, or all 0. The remaining classes were found through various tweaks to the original I and II.
