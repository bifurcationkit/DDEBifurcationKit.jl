# Periodic orbits computation

Consider the DDE problem written

$$\frac{du}{dt}=F(u(t), u_t;p).\tag{E}$$

A periodic solution $u^*$ with period $T$ satisfies $u^*(t+T)=u^*(t)$ for all time $t$.

We provide 1 method for computing periodic orbits (PO):

2. one (Collocation) based on orthogonal collocation to discretize the above problem (E),


## Important notes

We regroup here some important notes which are valid for all methods above. 


### 1. Finaliser
If you pass a `finalize` function to [continuation](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/library/#BifurcationKit.continuation), the following occurs:

1. If the newton solve was successful, we update the phase condition every `update_section_every_step`
2. we call the user defined finalizer `finalize`.
