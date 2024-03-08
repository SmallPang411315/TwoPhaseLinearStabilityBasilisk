# TwoPhaseLinearStabilityBasilisk

This is a $Basilisk$ version that includes two-phase linear blobal stability analysis with general eigenvalue solver Arpack.

Firstly, we have to get a steady two-phase flow, and the interface is described by Level-Set function. In the folder [navier-stokes](navier-stokes/), [centered_SFD.h](navier-stokes/centered_SFD.h) helps you get a steady flow state through SFD method. During this process, , [level_set_reinitialization.h](level_set_reinitialization.h) deals with Hamilton–Jacobi-typed equations, which is ready for the computation of derivatives of interface. 

Secondly, [centered-LSA.h](navier-stokes/centered-LSA.h) is about the time integration of linearised perturbation equations. [two-phase-LS.h](two-phase-LS.h) reconstruct the steady flow state. [bcg-LSA.h](bcg-LSA.h) and [viscosity-LSA.h](viscosity-LSA.h) deal with linearised advection and viscosity terms, respectively. The most important term is linearised surface tension force that is exerted by [iforce-LSA.h](iforce-LSA.h), [tracer-LSA.h](tracer-LSA.h), and [tension-LSA.h](tension-LSA.h).

Thirdly, we link the general eigenvalue solver $Arpack$ with $Basilisk$. Given a appropriate time interval and random noise as an initial guess, calling the second process and we will have a new result. Then you just put this new result into $Arpack$ function. After some iterations, You can have some pairs of eigenslutions.

In [Validation_Schmidt_2021](Validation_Schmidt_2021/), some .c codes that are used to read base flow states, make an time integration of linearised eqautions and give the flow information to eigen-problem solver. And we push the case for Weber number $We=5$ as an example. In this folder [Initial](Validation_Schmidt_2021/Initial), some videos about perturbed velocities and perturbed interface are made in time for a vivid view. The most dangerous mode emerges as the integration time is long enough. Thus, we can have a quick convergence by using the linear flow field as initial guess instead of random noise.

In [JetLSA](JetLSA/), the file EigenVec comprises of all eigenvectors, and top eigenvalues are logged in file [out](JetLSA/out) for each iteration. You can run [EigenPlot.py](JetLSA/EigenPlot.py) to plot the eigenvectors. 

There are still many resutls, and we do not upload all because of storage. Sorry for that! And we appreciate the help from Simon Schmidt. Thanks! 
