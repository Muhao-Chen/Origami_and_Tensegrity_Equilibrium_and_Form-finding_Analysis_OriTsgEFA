---

title: 'OriTsgEFA: Origami and Tensegrity Equilibrium and Form-finding Analysis'


tags:
  - Origami and Tensegrity systems
  - Flexible structures
  - Prestressable structures
  - Elastic and plastic deformation
  - Stiffness of origami and tensegrity structures 

authors:
  - name: Shuo Ma
    orcid: 0000-0003-3789-2893
    affiliation: 1

  - name: Muhao Chen^[corresponding author]
    orcid: 0000-0003-1812-6835
    affiliation: 2

  - name: Robert E. Skelton
    orcid: 0000-0001-6503-9115
    affiliation: 3
    
affiliations:
 - name: College of Civil Engineering, Zhejiang University of Technology, Hangzhou, Zhejiang, China
   index: 1

 - name: Department of Mechanical and Aerospace Engineering, University of Kentucky, Lexington, KY, USA
   index: 2

 - name: Department of Aerospace Engineering, Texas A&M University, College Station, TX, USA
   index: 3

date: 31 August 2023
bibliography: paper.bib
---

# Summary

To fully grasp the foundational concept of integrating origami and tensegrity frameworks, it is imperative to first understand the intrinsic characteristics, advantages, and limitations of each individual paradigm. Originating from the ancient traditions of Chinese paper invention and Japanese artistry, the art form known as origami combines the Japanese terms "oru" (to fold) and "kami" (paper) [@meloni2021engineering]. Initially celebrated in Asian ceremonies like weddings and festivals, origami gained global recognition in the 1950s through renowned artist Akira Yoshizawa. Yoshizawa authored 18 books and crafted more than 50,000 models, dramatically elevating the art form's status [@meloni2021engineering]. Over the years, rigorous academic research has transitioned origami from a purely aesthetic discipline to an engineering asset used in applications like space solar panels [@chen2019autonomous], energy absorption systems [@li2019origami], and robotics [@ze2022soft]. The advantages of origami systems are manifold: (1). Complexity of Form: Origami allows the creation of intricate shapes through complex folding techniques [@hernandez2019active;@lang2012origami]. (2). Compactability: Two-dimensional sheets can be easily compressed into compact forms and later expanded into voluminous structures [@huang2020complex;@zhang2019fold]. (3). Cost-Efficient Manufacturing: The transformation of 2D sheets into 3D or 4D structures simplifies the assembly process, making it highly cost-effective [@ge2014active;@zhao20183d]. (4). Metamaterial Characteristics: Many origami forms feature tunable properties like negative Poisson's ratios, structural bistability, and self-locking capabilities [@kamrava2017origami;@yasuda2015reentrant]. However, the system is not without its limitations: (1). Low Stiffness: While actuators or locking mechanisms can be added to increase stiffness, these modifications necessitate additional design complexity. (2). Control Challenges: The transient dynamics during shape morphing are difficult to manage. While actuation at the hinges can partially mitigate this, it adds further intricacy to the design.

The concept of tensegrity was initially formulated by artists Ioganson in 1921 and Snelson in 1948, showcasing a structural framework comprising compressible elements like bars or struts and tensile components such as strings or cables [@snelson1965continuous]. The term "tensegrity" itself was coined by Buckminster Fuller in 1959, melding 'tensile' and 'integrity' [@fuller1959us]. Subsequent decades of research have demonstrated tensegrity's immense utility in crafting lightweight structures, deployable platforms, and robotics, with applications spanning space platforms [@roffman2021morphing], cable domes [@ma2020design], and space landers [@kim2020rolling]. Strengths of Tensegrity Systems include: (1). Mass Efficiency: Due to axial loading on all one-dimensional structural elements, the system offers lightweight design solutions [@chen2023minimal;@fraddosio2019minimal]. (2). Model Accuracy: The absence of bar bending allows for highly accurate structural models, enabling precise control [@kan2021novel;@ma2022dynamics]. (3). Adjustable Parameters: Structure parameters such as prestress levels and string lengths can be adjusted to achieve varying equilibrium states, stiffness, and shapes [@feng2020optimal;@trinh2022force]. (4). Metamaterial Properties: The configuration of bar-string patterns can yield materials with unique properties [@amendola2019optimal;@bauer2021tensegrity]. (5). Morphological Control: Actuation of the strings allows for substantial shape morphing [@rhodes2019compact;@zhou2021distributed]. Limitations of Tensegrity Systems: (1). Lack of Enclosure: The structure is typically open and can be augmented with membranes or shells, though research in this domain remains limited. (2). Complex Joints: Given that a single node may connect multiple bars, each with varying degrees of freedom, the creation of joints tends to involve non-standard types.

The analysis of static properties in origami and tensegrity structures is critical for understanding their inherent characteristics. The OriTsg framework formulates equations for static analysis based on the following key assumptions: (1). Bars and string elements are subjected only to axial loading and are joined by frictionless ball joints. (2). Both bar and string elements can undergo elastic or plastic deformation. (3). Any rotation along the longitudinal axis of a structural member is disregarded. (4). Each structural member is assumed to be homogeneous in composition and has a uniform cross-section, resulting in a uniform distribution of mass along its length. (5). A string is only capable of pulling; if it tries to push, its tension is set to zero. (6). The origami structure is modeled based on a bar-and-hinge approach, where origami panels act as bars and their junctions as hinges. (7). The stress in the hinges is considered to be independent of the axial stress induced by either tension or compression [@liu2017nonlinear;@filipov2017bar]. Building upon Finite Element Method (FEM) and Lagrangian approaches, with nodal vectors serving as generalized coordinates [@chen2020general;@ma2022tensegrity], the OriTsg statics framework has been developed.

The software, referred to as OriTsgEFA, is recommended to be pronounced as "Ori Tenseg EFA." It offers static analysis capabilities for both origami and tensegrity structures. In terms of static analysis, OriTsgEFA enables the calculation of equilibrium configurations under various conditions, including external forces, nodal displacement, and changes in string rest lengths. The software accounts for nonlinear material and geometrical properties and employs the Modified Newton Method to ensure convergence to a stable equilibrium state.


# Statement of need

In the realm of structural engineering and design, there exists a divide over which paradigm—tensegrity or origami—is more fundamental. Advocates for tensegrity point to its widespread application in natural systems such as cellular structures [@wang2001mechanical], elbow mechanics [@scarr2012consideration], and spider fibers [@fraternali2020tensegrity]. Conversely, supporters of the origami paradigm also present compelling evidence from nature, including the folding mechanisms observed in flowers [@caratelli2022paper], horse-chestnut leaves [@kresling2012origami], and earwig wings [@faber2018bioinspired]. Given the inherent strengths of both paradigms, it is posited that their integration could offer synergistic advantages. Yet, the literature reflecting the amalgamation of these two systems remains scant. For instance, Rohmer et al. explored shape memory alloy-based tensegrity cylinders integrated with origami structures [@rohmer2015experimental]. Miranda et al. demonstrated the use of tensegrity D-Bars in actuating origami sunscreens for energy-efficient buildings [@miranda2020mechanics]. Similarly, Park proposed a geodesic dome designed for Martian agriculture that leverages both tensegrity and origami elements [@hong2021tensegami]. Despite these isolated studies, analytical research focusing on the fusion of these paradigms remains underdeveloped. This software aims to bridge this gap by providing a unified static analysis framework for both origami and tensegrity systems. It accomplishes this by formulating governing equations that allow the origami and tensegrity paradigms to be modeled and analyzed as an integrated system.

The software developed offers a comprehensive framework for modeling and analyzing both origami and tensegrity paradigms as a unified system. The capabilities of the software extend to various types of static analyses, including (1). Load analysis, accommodating both elastic and plastic deformations in bars and strings. (2). Analysis of both infinitesimal and large-scale deformations. (3). Versatility in handling diverse boundary conditions, such as fixed points or nodes subjected to various types of static loads in any given direction—be it gravitational forces, specified forces, or moments. (4). Stiffness analysis featuring eigenvalue computations and mode shapes. By providing nuanced insights into structural integrity, material properties, and overall performance, this software enhances our capacity to design and deploy large-scale, multifaceted structures.



# References

