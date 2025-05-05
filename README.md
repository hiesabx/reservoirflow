# ReservoirFlow: Reservoir Simulation and Engineering Library in Python

<!-- 
![five_spot_single_phase](/docs\source\user_guide\tutorials\tutorial_five_spot_single_phase\grid_animated.gif)\ 
![five_spot_single_phase](https://drive.google.com/uc?id=11NhTbAU_lA768yiEAsoA18SshMjDtRqZ)\
![five_spot_single_phase](https://drive.google.com/thumbnail?id=1mQ276IokIJBUQMZN2BOcGiV6pLwihguT&sz=w1000)\
*Example: Pressure Distribution of Single Phase Flow in Five Spot Wells Patterns.*
-->

**Table of Content:**

- [ReservoirFlow: Reservoir Simulation and Engineering Library in Python](#reservoirflow-reservoir-simulation-and-engineering-library-in-python)
  - [Introduction](#introduction)
  - [Installation](#installation)
  - [Import Convention](#import-convention)
  - [Version](#version)
  - [License](#license)
  - [Disclaimer](#disclaimer)
  - [Contact](#contact)

## Introduction

*ReservoirFlow* is a modern open-source Python library developed by [Hiesab](https://www.hiesab.com/en/); a startup company specialized in advanced analytics, computing, and automation founded in 2024 with a mission is to accelerate R&D for AI applications and data solutions in different fields including Science, Engineering, and Education, see [About Us](https://reservoirflow.hiesab.com/about_us.html).

<!-- <p align="center">
<iframe src="https://drive.google.com/file/d/1JCD7W_5vJsUqNWf99NBTYjsLSMeVPEEa/preview"
allow="autoplay"  width="100%" height="500" frameborder="0" scrolling="auto" class="iframe-full-height" allowfullscreen></iframe>
</p> -->

*ReservoirFlow* is designed to study and model the process of fluid flow in porous media related to subsurface energy storage systems, reservoir simulation and engineering. *ReservoirFlow* is the first reservoir simulator based on physics-informed neural network models and one of its kind in a sense that it allows comparing and combining analytical solutions, numerical solutions, and neurical solutions (i.e. solutions based on artificial neural networks). *ReservoirFlow* is planned to be a central platform between education and industry where scientific papers are implemented and distributed in a standard and accessible format with coding examples, tutorials, and trainings.

<!--
تدفق المكامن هي مكتبة حديثة مفتوحة المصدر تم تطويرها بواسطة شركة حساب وهي مصممة لدراسة ومحاكاة ظاهرة تدفق الموائع في الوسط المسامي المعروفة باسم محاكاة وهندسة المكامن.
-->

*ReservoirFlow* is designed based on the modern Python stack for data science, scientific computing, machine learning, and deep learning with the objective to support high-performance computing including multithreading, parallelism, GPU, and TPU. Throughout our computing problems based on large simulation models, intensive benchmarking well be carefully designed and carried out to evaluate the performance of computing software (i.e. frameworks) and hardware (e.g. GPUs). The outcome of these tests will be used to further improve the performance of *ReservoirFlow* and to provide materials with recommendations about available computing tools, techniques and frameworks. *ReservoirFlow* is planned to support different backends including [NumPy](https://numpy.org/doc/stable/index.html), [SciPy](https://scipy.org/), [JAX](https://jax.readthedocs.io/en/latest/index.html), [PyTorch](https://pytorch.org/), [TensorFlow](https://www.tensorflow.org/), and more.

<p align="center">
  <img src="https://drive.google.com/thumbnail?id=1mQ276IokIJBUQMZN2BOcGiV6pLwihguT&sz=w1000" 
  alt="GIF image"
  title="Example: Pressure Distribution of Single Phase Flow in Five Spot Wells Patterns." width="100%" height="100%"
  />
  <em>
  Example: Pressure Distribution of Single Phase Flow in Five Spot Wells Patterns, see <a href="https://reservoirflow.hiesab.com/user_guide/tutorials/tutorial_five_spot_single_phase.html">Tutorials</a>.
  </em>
</p>

*ReservoirFlow* brings reservoir simulation and engineering to the Python ecosystem to empower automation in intelligent fields where engineers and specialists can deploy their models in containers that will be ready to make real-time optimization for any well in the field. In contrast to commercial black-box software where reservoir simulation studies are relatively isolated, important actions can be immediately predicted and made available for the field hardware to execute. A special attention well be given to provide solutions for environmentally friendly projects with a clear objective to reduce emissions. We are committed to extend our tools to cover the topic of Carbon Capture and Storage (CCS) especially $CO_2$ Underground Storage. In addition, we are looking forward to covering a wider range of topics from Reservoir Engineering including: Pressure Transient Analysis (PTA) and Rate Transient Analysis (RTA), Enhanced Oil Recovery (EOR), Improved Oil Recovery (IOR), Pressure-Volume-Temperature (PVT), Equation-of-State (EOS), etc.

*ReservoirFlow* aims to achieving a high quality open research for reservoir simulation and engineering to provide solutions that combine the strength of scientific computing with the power of deep learning for different applications such as: reverse computing, interpolation or extrapolation, etc. Below are few examples of the problems that will be tackled in the future:

- Real-time reservoir management and production optimization using Cloud Computing and IoT.
- Reinforcement learning to achieve better production strategies for specific goals (e.g. maximize recovery, accelerate production).
- History matching using machine learning.
- Advanced computing such as GPU, TPU and Quantum Computing.
- Scientific Machine learning using Physics-informed neural networks (PINNs) or DeepONets.

An open-source reservoir simulation and engineering library within the Python ecosystem is also very important to students, universities, researchers, engineers, and practitioners. Unlike the common monopolistic approach in the Oil and Gas industry where software is usually offered as a closed black-box at a high cost, we plan to make our tools accessible and freely available to everyone except for commercial-use where an explicit authorization will be required. We aim to offer our sponsors the commercial-use license with other benefits including trainings, custom features, studies, and more. On the other hand, our license allows universities, students, academics, and researchers to use our tools directly for teaching or publication just with a proper referencing. Therefore, the growth of this tool can only be taken as a positive growth for a new community that we try to create. However, this requires a huge support to meet the upcoming challenges that we are looking for, see [Support Us](https://reservoirflow.hiesab.com/support_us.html).

## Installation

Install `reservoirflow` directly from [PyPi](https://pypi.org/project/reservoirflow/):

```console
$ pip install reservoirflow
```

For more information about the installation process, see: [Getting Started](https://reservoirflow.hiesab.com/user_guide/getting_started/getting_started.html) in the documentation.

## Import Convention

The following convention is used to import `reservoirflow` after installation:

```python
>>> import reservoirflow as rf
```

The abbreviation `rf` refers to `reservoirflow` where all modules under this library can be accessed. `rf` is also used throughout the documentation. We recommend our users to stick with this convention.

## Version

[Semantic Versioning](https://semver.org/) is used for the version numbers. Since this library is still under active development, `major=0` is used until the first stable version is reached. The first version `v0.1.0b1` was released on May 1, 2025. The current version is `v0.1.0b1`. To know about which features are currently supported, check [Capabilities](https://reservoirflow.hiesab.com/capabilities.html).

**Version History:**

| |**Version**  |**Status**         |**Release Date (dd.mm.yyyy)**  |
|-|-            |-                  |-                              |
|1|`v0.1.0b1`     |current version    |01.05.2025                     |
|2|`v0.1.0rc1`     |*under development*|*ongoing*                      |

## License

This work is licensed under a [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License](https://creativecommons.org/licenses/by-nc-sa/4.0/). Detailed license can be found [here](https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode). The current license does not allow commercial use without an explicit authorization from the author. For commercial applications, please [contact](#contact) us.

## Disclaimer

*ReservoirFlow* is in a Beta state and still under active development where most modules are not fully functional. We expect our users and partners to collaborate with us to improve this library to reach a production state withing the following years.

*ReservoirFlow* is developed and copyrighted by Hiesab. Third-party components are copyrighted by their respective authors.

## Contact

Check our website: [hiesab.com](https://www.hiesab.com/en/). Check Contact page.
