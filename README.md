# Equivariant Graphic Statics
This project investigates graphic statics of a symmetric framework and its irreducible components.


```python
import CODE.cosheaves
import CODE.samples
import CODE.output
```

The module *samples* contains many sample framework geometries. Calling a function in sample returns a pair **(cell_complex, action)** where **cell_complex** encodes cell combinatorics and **action** encodes the specified group action on the vertices of the **cell_complex**.


```python
cell_complex, action = CODE.samples.boxed()
print("Plotting Framework:")
cell_complex.plot()
action.print_char_table()
```

    Box wheel with spokes sample
    Dihedral symmetry order 4
    Plotting Framework:
    


    
![png](README_files/README_3_1.png)
    


    Dihedral Group of order: 8
    [1. 1. 1. 1. 1. 1. 1. 1.]
    [ 1.  1.  1.  1. -1. -1. -1. -1.]
    [ 1. -1.  1. -1.  1. -1.  1. -1.]
    [ 1. -1.  1. -1. -1.  1. -1.  1.]
    [ 2.  0. -2. -0.  0.  0.  0.  0.]
    

The module *cosheaves* contains the force $\mathcal{F}$, constant $\overline{\mathbb{R}^2}$, and position $\mathcal{P}$ cosheaf classes. We check  that these are well-defined and the cosheaf maps that pass between are valid by running **tests.run_tests(F,J,P)**. These ensure that all the commutitive squares between the boundary maps, cosheaf maps, and cosheaf representations (for each group generators) are satisfied.


```python
J = CODE.cosheaves.constant_cosheaf(cell_complex, action)
F = CODE.cosheaves.force_cosheaf(cell_complex, action)
P = CODE.cosheaves.position_cosheaf(cell_complex, action)

CODE.output.run_tests(F,J,P)
```

    Checking commutativity of cosheaf representation and boundary matrices
    
    Group generator index: 1
    Force Cosheaf,    dim 1 error: 0.0
    Force Cosheaf,    dim 2 error: 0
    Constant Cosheaf, dim 1 error: 0.0
    Constant Cosheaf, dim 2 error: 0.0
    Position Cosheaf, dim 1 error: 0
    Position Cosheaf, dim 2 error: 0.0
    Group generator index: 4
    Force Cosheaf,    dim 1 error: 0.0
    Force Cosheaf,    dim 2 error: 0
    Constant Cosheaf, dim 1 error: 0.0
    Constant Cosheaf, dim 2 error: 0.0
    Position Cosheaf, dim 1 error: 0
    Position Cosheaf, dim 2 error: 0.0
    
    ---------------------------------------
    
    Check commutativity of cosheaf representation and cosheaf map matrices
    
    Group generator index: 1
    Force Cosheaf,    dim 0 error: 0.0
    Force Cosheaf,    dim 1 error: 0.0
    Force Cosheaf,    dim 2 error: 0
    Constant Cosheaf, dim 0 error: 0
    Constant Cosheaf, dim 1 error: 0.0
    Constant Cosheaf, dim 2 error: 0.0
    Group generator index: 4
    Force Cosheaf,    dim 0 error: 0.0
    Force Cosheaf,    dim 1 error: 0.0
    Force Cosheaf,    dim 2 error: 0
    Constant Cosheaf, dim 0 error: 0
    Constant Cosheaf, dim 1 error: 0.0
    Constant Cosheaf, dim 2 error: 0.0
    
    ---------------------------------------
    
    Check commutativity of boundary and cosheaf map matrices
    
    Force Cosheaf,    dim 1 error: 0.0
    Force Cosheaf,    dim 2 error: 0
    Constant Cosheaf, dim 1 error: 0
    Constant Cosheaf, dim 2 error: 0.0
    

Everything looking good! The first homology $H_1 \mathcal{F}$ encodes the self-stresses of the framework and is computed using the method **F.homology(dim)**. There is one degree of self-stress in this example. By classical graphic statics, there are three degrees of reciprocal diagrams for this example, the dimension of $H_2 \mathcal{P}$. We can plot the reciprocal diagram and the corresponding self-stress with the method **P.plot_both**. 


```python
CODE.output.print_homologies(F, J, P)

print("\n")
print("Basis for self-stresses")
print(F.homology(1))

#There is one degree of self-stress, so we choose to plot this first dimension only
print("\nSelf stress from the force cosheaf")
F.plot_self_stress([1])
#If there were two degrees of self-stress, we could imput:  F.plot_self_stress([0.5, 0.5]) for the average of both generators.

print("Reciprocal diagram from the position cosheaf and matched self-stress")
P.plot_both([1, 0, 0])
#The vector [1,0,0] corresponds to choosing the first basis generator of H_2 P
```

    Homology space dimensions:
       Force  Constant  Position
    H_2: 0 ------ 2 ------ 3
    H_1: 1 ------ 0 ------ 1
    H_0: 3 ------ 2 ------ 0
    
    
    Basis for self-stresses
    [[ 0.28867513]
     [ 0.28867513]
     [ 0.28867513]
     [ 0.28867513]
     [-0.40824829]
     [-0.40824829]
     [-0.40824829]
     [-0.40824829]]
    
    Self stress from the force cosheaf
    Blue edges in tension, red edges in compression.
    


    
![png](README_files/README_7_1.png)
    


    Reciprocal diagram from the position cosheaf and matched self-stress
    


    
![png](README_files/README_7_3.png)
    


We can distinguish the self-stress - reciprocal diagram pairs by irreducible representation of the underlying group. We observe that the only non-trivial self-stress - reciprocal diagram pair above corresponds to the first (trivial) irreducible representation of $D_8$.


```python
#Output the irreducible self-stress - reciprocal diagram pairs
CODE.output.irreducible_pairs(F,J,P)
```

    Irreducible number: 0 Dimension: 1.0
    Dimensions of the irreducible cosheaf homology:
         F --R^2-- P
    H_2: 0 -- 0 -- 1
    H_1: 1 -- 0 -- 0
    H_0: 0 -- 0 -- 0
    


    
![png](README_files/README_9_1.png)
    


    -------------------------------------------------------------------------------------------
    Irreducible number: 1 Dimension: 1.0
    Dimensions of the irreducible cosheaf homology:
         F --R^2-- P
    H_2: 0 -- 0 -- 0
    H_1: 0 -- 0 -- 1
    H_0: 1 -- 0 -- 0
    -------------------------------------------------------------------------------------------
    Irreducible number: 2 Dimension: 1.0
    Dimensions of the irreducible cosheaf homology:
         F --R^2-- P
    H_2: 0 -- 0 -- 0
    H_1: 0 -- 0 -- 0
    H_0: 0 -- 0 -- 0
    -------------------------------------------------------------------------------------------
    Irreducible number: 3 Dimension: 1.0
    Dimensions of the irreducible cosheaf homology:
         F --R^2-- P
    H_2: 0 -- 0 -- 0
    H_1: 0 -- 0 -- 0
    H_0: 0 -- 0 -- 0
    -------------------------------------------------------------------------------------------
    Irreducible number: 4 Dimension: 2.0
    Dimensions of the irreducible cosheaf homology:
         F --R^2-- P
    H_2: 0 -- 1 -- 1
    H_1: 0 -- 0 -- 0
    H_0: 1 -- 1 -- 0
    Degenerate picture, stress is trivial
    Degenerate picture, stress is trivial
    -------------------------------------------------------------------------------------------
    

We observe that the two degrees of translation correspond to the last irreducible representation of $D_8$. This irreducible has dimension 2. This framework isn't that interesting, so we consider the **klein_four** sample with $D_4$ symmetry.


```python
cell_complex2, action2 = CODE.samples.klein_four()
J2 = CODE.cosheaves.constant_cosheaf(cell_complex2, action2)
F2 = CODE.cosheaves.force_cosheaf(cell_complex2, action2)
P2 = CODE.cosheaves.position_cosheaf(cell_complex2, action2)

print("\n-----------------------------------------------")
CODE.output.irreducible_pairs(F2, J2, P2)

```

    Klein four sample 
    Dihedral symmetry order 4
    
    -----------------------------------------------
    Irreducible number: 0 Dimension: 1.0
    Dimensions of the irreducible cosheaf homology:
         F --R^2-- P
    H_2: 0 -- 0 -- 3
    H_1: 3 -- 0 -- 0
    H_0: 0 -- 0 -- 0
    


    
![png](README_files/README_11_1.png)
    



    
![png](README_files/README_11_2.png)
    



    
![png](README_files/README_11_3.png)
    


    -------------------------------------------------------------------------------------------
    Irreducible number: 1 Dimension: 1.0
    Dimensions of the irreducible cosheaf homology:
         F --R^2-- P
    H_2: 0 -- 0 -- 0
    H_1: 0 -- 0 -- 2
    H_0: 2 -- 0 -- 0
    -------------------------------------------------------------------------------------------
    Irreducible number: 2 Dimension: 1.0
    Dimensions of the irreducible cosheaf homology:
         F --R^2-- P
    H_2: 0 -- 1 -- 2
    H_1: 1 -- 0 -- 0
    H_0: 1 -- 1 -- 0
    


    
![png](README_files/README_11_5.png)
    



    
![png](README_files/README_11_6.png)
    


    -------------------------------------------------------------------------------------------
    Irreducible number: 3 Dimension: 1.0
    Dimensions of the irreducible cosheaf homology:
         F --R^2-- P
    H_2: 0 -- 1 -- 1
    H_1: 0 -- 0 -- 0
    H_0: 1 -- 1 -- 0
    Degenerate picture, stress is trivial
    -------------------------------------------------------------------------------------------
    

We can gain a better understanding of the homology spaces by understanding the extended euler characteristic for each cosheaf. Take the force cosheaf chain complex $C\mathcal{F}$, decomposed as a direct sum of irreducible chain complexes $C_1^{(j)} \mathcal{F} \xrightarrow{\partial} C_0^{(j)} \mathcal{F}$ for every irreducible representation $\mu^{(j)}$. Each chain complex has an *Euler characteristic equation*, relating the dimensions of $C^{(j)} \mathcal{F}$ and $H^{(j)} \mathcal{F}$ (known as *Maxwell's rule* over the force cosheaf $\mathcal{F}$). This Euler characteristic extends to a character counting rule, relating the number of copies of $\mu^{(j)}$ in the representation for $C^{(j)}\mathcal{F}$ and sub-representation for $H^{(j)} \mathcal{F}$.


```python
CODE.output.irreducible_euler_char(F2)
```

    Irreducible number: 0 Dimension: 1.0
    
    Chain Characters    # copies of irriducible    |    Homology Characters    # copies of irriducible
    ------------------  -------------------------  ---  ---------------------  -------------------------
    char over C_0:      9                          |    char over H_0:         0
    char over C_1:      12                         |    char over H_1:         3
    char over C_2:      0                          |    char over H_2:         0
                        --------                        --------
    alternating sum:    -3                         |                           -3
    
    ----------------------------------------------------------------------------------------
    
    Irreducible number: 1 Dimension: 1.0
    
    Chain Characters    # copies of irriducible    |    Homology Characters    # copies of irriducible
    ------------------  -------------------------  ---  ---------------------  -------------------------
    char over C_0:      9                          |    char over H_0:         2
    char over C_1:      7                          |    char over H_1:         0
    char over C_2:      0                          |    char over H_2:         0
                        --------                        --------
    alternating sum:    2                          |                           2
    
    ----------------------------------------------------------------------------------------
    
    Irreducible number: 2 Dimension: 1.0
    
    Chain Characters    # copies of irriducible    |    Homology Characters    # copies of irriducible
    ------------------  -------------------------  ---  ---------------------  -------------------------
    char over C_0:      10                         |    char over H_0:         1
    char over C_1:      10                         |    char over H_1:         1
    char over C_2:      0                          |    char over H_2:         0
                        --------                        --------
    alternating sum:    0                          |                           0
    
    ----------------------------------------------------------------------------------------
    
    Irreducible number: 3 Dimension: 1.0
    
    Chain Characters    # copies of irriducible    |    Homology Characters    # copies of irriducible
    ------------------  -------------------------  ---  ---------------------  -------------------------
    char over C_0:      10                         |    char over H_0:         1
    char over C_1:      9                          |    char over H_1:         0
    char over C_2:      0                          |    char over H_2:         0
                        --------                        --------
    alternating sum:    1                          |                           1
    
    ----------------------------------------------------------------------------------------
    
    

The characters of the chain complex representation $\chi(C^{(j)}\mathcal{F})$ above count the number of force asignment chains fixed by the $\mu^{(j)}$ symmetry. For instance, take the trivial irreducible representation $\mu^{(0)}$ with complete $D_4$ symmetry. The number of 1-chains fixed by every group element of $D_4$ is equal to the number of edge orbits of the edges under the $D_4$ permutation action. We see that $\chi(C_1^{(0)} \mathcal{F})$ has twelve copies of $\chi(\mu^{(0)})$ meaning that there are twelve edge orbits. The number of edge orbits can be determined by counting the representative edges in one quadrant of the figure.
