## Numerical Modeling of Heat Exchange in a Fluid Flowing Through a Conduit

This project focuses on the finite element numerical modeling of heat exchange in a flue gas conduit, simulating convective-diffusive thermal transfer. Based on a previously designed air-based residential boiler (you can find it @ gianmarco-san.github.io) and examines how heat dissipates from the combustion gases flowing through the chimney, both inside and outside the building. The objective is to quantify temperature distribution and optimize energy recovery.

The governing equation accounts for convective heat transfer within the flue gases and conductive transfer through the metal walls, using a cylindrical coordinate system for computational efficiency. The model incorporates boundary conditions for internal and external environments, including Dirichlet, Neumann, and Robin conditions, ensuring realistic heat exchange representation.

A finite element discretization approach was applied to solve the stationary and transient heat transfer problem. The mesh resolution was optimized to balance computational cost and accuracy, using Gaussian elimination for matrix resolution. The transient simulation models the system’s thermal evolution, validating that steady-state conditions match real operational scenarios.

A sensitivity analysis was conducted by varying the gas velocity, demonstrating its impact on heat dissipation and flue gas temperature at the outlet. Results confirmed that lower velocities increase heat loss, reducing the final temperature of exhaust gases.

