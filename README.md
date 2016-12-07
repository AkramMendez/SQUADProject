# Qualitative modeling of regulatory networks using SQUAD.

This project exemplifies the implementation of the _Standardized Qualitative Analysis of Dynamics (SQUAD)_ method as explained in the Chapter ___"The SQUAD method for the qualitative analysis of regulatory networks"___.

See the [SQUAD_implementation.ipynb](https://github.com/caramirezal/SQUADBookChapter/blob/master/SQUAD_implementation.ipynb) tutorial for a basic explanation of how to implement the SQUAD method in R programming language. This tutorial shows how to simulate the dynamics of small regulatory networks from a toy model and from a known biological system as presented in Fig. 4 and Fig. 5 from the Chapter. 

## Translate a Boolean network to a continuous dynamical system
The [CartoonNetwork.R](https://github.com/caramirezal/SQUADBookChapter/blob/master/cartoonNetwork.R) file includes the regulatory network comprising three nodes _A_,_B_ and _C_ shown in Fig. 1 from the book Chapter.

[SQUADImplementation.R](https://github.com/caramirezal/SQUADBookChapter/blob/master/SQUADImplementation.R) file exemplifies the simulation of the regulatory network model defined in the [CartoonNetwork.R](https://github.com/caramirezal/SQUADBookChapter/blob/master/cartoonNetwork.R) file, first as a Boolean regulatory network model, and then as a continuous regulatory network using the SQUAD method. 

Additionally, [animation.gif](https://github.com/caramirezal/SQUADBookChapter/blob/master/animation.gif) shows an animation of the trajectories of the states of the network at different rotation angles that converge into two different steady states as shown in Fig. 3 from the Chapter.

