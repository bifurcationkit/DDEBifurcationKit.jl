# Periodic orbits computation

Consider the DDE problem written

$$\frac{du}{dt}=F(u(t), u_t;p).\tag{E}$$

A periodic solution $u^*$ with period $T$ satisfies $u^*(t+T)=u^*(t)$ for all time $t$.

We provide 1 methods for computing periodic orbits (PO):

2. one (Collocation) based on orthogonal collocation to discretize the above problem (E),


It is important to understand the pro and cons of each method to compute PO in large dimensions.


### Collocation method

The Collocation method is by far the most precise one. Additionally, the mesh can be automatically adapted during the continuation. The implementation will be improved for large dimensional systems.

## Important notes

We regroup here some important notes which are valid for all methods above. 


### 1. Finaliser
If you pass a `finalize` function to [continuation](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/library/#BifurcationKit.continuation), the following occurs:

1. If the newton solve was successfull, we update the phase condition every `updateSectionEveryStep`
2. we call the user defined finalizer `finalize`
