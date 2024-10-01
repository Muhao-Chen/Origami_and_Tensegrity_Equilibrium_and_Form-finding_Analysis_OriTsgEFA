# Origami and Tensegrity Equilibrium and Form-finding Analysis (OriTsgEFA)

### **Welcome to OriTsgEFA Software!**

#### General Information

This research focuses on integrated origami and tensegrity structures in equilibrium and form-finding analysis. Using an analytical approach, the software enables modeling and evaluating both origami and tensegrity paradigms as a cohesive system. Parameters such as nodal coordinates and hinge angles among origami panels characterize tensegrity and origami structures. This software is released as open-source to assist researchers interested in this multidisciplinary field.

This user guide provides details on all aspects of the software to improve its usability. Questions and suggestions for enhancement are highly encouraged. The software specializes in performance modeling, structural design, and nonlinear static simulations of origami and tensegrity systems. The primary contributions of this software are in the following areas:
#### Modeling:
1. Enable modeling of any origami and tensegrity structures via nodal coordinates, node connectivity, and hinge angles between panels.
2. Provide the capability to specify constraints on nodal coordinates, such as limiting movements in the X-, Y-, or Z-axes for particular nodes.
3. Allow for the grouping of structural members within the software.

#### Statics:
1. Execute prestress and mechanism mode evaluations through singular value decomposition of the equilibrium matrix.
2. Conduct stiffness and stability analyses.
3. Resolve equilibrium equations considering any applied external forces and the nonlinearity of geometry and material properties.
4. Facilitate simulations of forced motions by inputting nodal sequences or adjusting the rest lengths of specific structure members.

The name OriTsgEFA is suggested to be pronounced as "Ori Tenseg EFA." The software provides a platform for static analysis of origami and tensegrity structures, encompassing:

1. Load analysis that accounts for either elastic or plastic deformations in bars and strings.
2. Both infinitesimal and large deformation analyses for a better understanding of structural stress and actuation strategies.
3. The ability to adapt to various boundary conditions, such as fixing or applying static loads at any nodes in any direction (e.g., gravitational force, specified forces, or moments).
4. Stiffness analysis features that include calculations of eigenvalues and their modes.

A foundational understanding of undergraduate-level linear algebra, material or continuum mechanics, finite element methods, and basic MATLAB knowledge is advised for effective software utilization. The software is developed for:

* 64-bit Windows  
* MATLAB  

**Note:** While the software is compatible with MATLAB versions later than 2009a on platforms like Win7/Win10/Mac OS/Linux/Win XP/Win Vista, using the most up-to-date release of MATLAB is strongly recommended. Further information regarding MATLAB versions can be accessed [here](https://en.wikipedia.org/wiki/MATLAB).


#### LICENSE

    /* This Source Code Form is subject to the terms of the Mozilla Public
     * License, v. 2.0. If a copy of the MPL was not distributed with this
     * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 
---

***<font size=4>Origami and Tensegrity Equilibrium and Form-finding Analysis (OriTsgEFA) folder contains the following parts:</font>***

---


#### setup.m

To utilize the OriTsgEFA software, initially execute the `setup.m` script. Open MATLAB and run `setup.m`, which performs the following actions:

- Incorporate all requisite library functions into the MATLAB path,
- Include the Software Verification and Examples folders,
- Add the User Guide,
- Add the Videos directory,
- Add the JOSS paper.

**Note:** It is imperative to run `setup.m` each time MATLAB is initiated before executing other files.

#### JossPaper

This directory houses the academic paper associated with the software, complete with source files and reference materials. The paper offers an introductory background, a project summary, its applications, citations, and more.

#### Function Library

Stored in this folder are the following:

All the functions are essential for origami-tensegrity statics analysis. To perform said analysis, follow the guidelines outlined in `User_Guide.pdf`.

#### Software Verification and Examples

This folder incorporates:

1. Five Statics Examples  
Herein, examples are presented to validate and showcase the statics capabilities of this software.

2. `auto all static test. m`  File  
A `.m` file is provided to execute all five statics examples sequentially.

#### User Guide

This directory contains `User_Guide.pdf`, which elaborates on software usage, application variety, and developer onboarding.

#### Videos

This folder showcases intriguing animations related to origami-tensegrity systems.

---

### **Help Desk**

We are eager to assist with any inquiries you may have. Clearly articulate your questions and email Muhao Chen: <muhaochen@uky.edu> Shuo Ma: <mashuo@zjut.edu.cn>. Your correspondence is appreciated.

---

### Acknowledgment

The authors extend their heartfelt gratitude to Mr. Xueshi Wang and Dr. Hongying Zhang for their invaluable assistance. Indeed, many thanks!

----

### **Become Part of the OriTsgEFA Community and Contribute**

#### How to Contribute

Your feedback and contributions are valuable. For seamless communication, adhere to consistent terminology.

1. Fork the repository,
2. Either submit a pull request or send your queries and contributions to the help desk.

Replies will be promptly provided.

#### Coding Standards

- MATLAB (version 2009a or later)
- Comments delineating function inputs and outputs
- Consistent nomenclature as specified



#### Nomenclature

##### Geometry: 
    N: initial node positions
    n: nodal coordinate vector
    C_b: bar connectivity
    C_s: string connectivity
    C: connectivity matrix of the structure
    C_h: hinge connectivity
    C_rh: rigid hinge connectivity
    Ca: triangle element connectivity
    E_n: transformation matrix from hinge element to total nodes
    S: clustering matrix
    ne: amount of members
    nn: amount of nodes
    n_h: amount of hinges
    n_p: amount of panels
    a: a vector containing the index of free nodal coordinate
    b: a vector containing the index of fixed nodal coordinate
    Ia: a matrix locating the free nodal coordinate
    Ib: a matrix locating the fixed nodal coordinate
    Gp: group matrix
    l: length of members
    l0_c: rest length of truss members  
    theta_0: initial angle of hinges  
##### Statics
    A_1: equilibrium matrix with no constraints, force density as the variable
    A_1c: Equilibrium matrix with group, force density as the variable.
    A_1a: equilibrium matrix with constraints, force density as the variable
    A_1ac: equilibrium matrix with constraints and group, force density as the variable
    A_2: equilibrium matrix with no constraints, force as variable
    A_2c: Equilibrium matrix with group, force as the variable.
    A_2a: equilibrium matrix with constraints, force as variable
    A_2ac: equilibrium matrix with constraints and group, force as variable
    phTpn: Jacobian matrix of origami
    G: Hessian matrix of origami
    A_o: the equilibrium matrix of the whole structure
    k_h: the hinge rotational stiffness
    w0: external force
    q_gp: force density vector in group
    q: force density vector
    t_gp: force vector in group 
    t: force vector
    index_b: number of bars
    index_s: number of strings
    A_b: cross-sectional area of bars
    A_s: cross-sectional area of strings
    A_gp: cross-sectional area of all members in a group
    A_c: cross-sectional area of all members
    E_c: Young's modulus of members
    Kt_aa: tangent stiffness of the tensegrity
    Kg_aa: geometry tangent stiffness of the tensegrity
    Ke_aa: material tangent stiffness of the tensegrity
    K_t_oa: tangent stiffness of the whole structure with constraints
    rho: material density of members
    mass: mass vector of members
    substep: substeps in static analysis
    w_t: external force in substeps
    dnb_t: forced motion of fixed nodal coordinate in substeps
    l0_ct: the changing rest length difference of truss elements
    theta_0_t: the changing angle difference of hinge elements 
    Fhis: the load factor
    t_t: the force of members in substeps
    n_t: nodal coordinate vector in substeps
    M_out: hinge moment in every step
    l_out: truss member length in every step
