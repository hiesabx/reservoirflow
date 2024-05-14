# ReservoirFlow: a Reservoir Simulation and Engineering Library in Python

<!-- ![five_spot_single_phase](/docs\source\user_guide\tutorials\tutorial_five_spot_single_phase\grid_animated.gif)\ -->
![five_spot_single_phase](https://drive.google.com/uc?id=11NhTbAU_lA768yiEAsoA18SshMjDtRqZ)\
*Example: Pressure Distribution of Single Phase Flow in Five Spot Wells Patterns.*

**Table of Content:**

- [Introduction](#introduction)
- [Installation](#installation)
- [Import Convention](#import-convention)
- [Version](#version)
- [License](#license)
- [Contact](#contact)

## Introduction

*ReservoirFlow* is a modern open-source Python library developed by [Zakariya Abugrin](https://github.com/zakgrin). *ReservoirFlow* is designed to study and model the phenomenon of fluid flow in porous media known as: Reservoir Simulation and Engineering. Interestingly, *ReservoirFlow* is the first reservoir simulator based on physics-informed neural network models and one of its kind in a sense that it allows using and combining different solutions including analytical, numerical, and neurical (i.e. solutions based on artificial neural networks from machine learning). *ReservoirFlow* is planned to be a central platform between education and industry where scientific papers are implemented and distributed in a standard and accessible format with hands-on tutorials.

<!--
تدفق المكامن هي مكتبة حديثة مفتوحة المصدر تم تطويرها بواسطة زكريا أبوقرين وهي مصممة لدراسة ومحاكاة ظاهرة تدفق الموائع في الوسط المسامي المعروفة باسم محاكاة وهندسة المكامن.
-->

*ReservoirFlow* is designed based on the modern Python stack for data science, scientific computing, and machine learning with the objective to support high-performance computing including multithreading, parallelism, GPU, and TPU. Through out our computing problems (e.g. large simulation models), intensive benchmarking well be carefully designed and carried out to evaluate the performance of computing software (i.e. frameworks) and hardware (e.g. GPU). The outcome of this benchmarking will be used to further improve the performance of *ReservoirFlow* and to provide materials with recommendations about available computing tools, techniques and frameworks. *ReservoirFlow* is planned to support different backends including Python frameworks such as [NumPy](https://numpy.org/doc/stable/index.html), [SciPy](https://scipy.org/), [JAX](https://jax.readthedocs.io/en/latest/index.html), [PyTorch](https://pytorch.org/), [TensorFlow](https://www.tensorflow.org/), and more.

*ReservoirFlow* brings reservoir simulation and engineering to the Python ecosystem including data analytics and machine learning tools to empower automation in intelligent fields where engineers and specialists can deploy their models in containers that will be ready to make real-time optimization for any well in the field. In contrast to commercial black-box software where reservoir simulation studies are relatively isolated, important actions can be immediately predicted and made available for the field hardware to execute. A special attention well be given to provide solutions for environmentally friendly projects with a clear objective to reduce emissions. We are committed to extend our tools to cover the topic of Carbon Capture and Storage (CCS) especially $CO_2$ Underground Storage. In addition, we are looking forward to cover a wider range of topics from Reservoir Engineering including: Pressure Transient Analysis (PTA) and Rate Transient Analysis (RTA), Pressure-Volume-Temperature (PVT), Equation-of-State (EOS), etc.

*ReservoirFlow* aims to achieving a high quality open research and science for reservoir simulation and engineering. Solutions that can combine the strength of scientific computing with the power of machine learning (e.g. deep learning) for different objectives (e.g. reverse computing, interpolation or extrapolation, etc.). Below are few examples of the problems that will be tackled in the future:

- Real-time reservoir management and production optimization using Cloud Computing and IoT.
- Reinforcement learning to achieve better production strategies for specific goals (e.g. maximize recovery, accelerate production).
- History matching using machine learning.
- Advanced computing such as GPU, TPU and Quantum Computing.
- Scientific Machine learning using Physics-informed neural networks (PINNs) or DeepONets.

An open-source reservoir simulation and engineering library within the Python ecosystem is also very important to students, universities, researchers, engineers, and practitioners. Therefore, the growth of this tool can only be taken as a positive growth for a new community that we try to create. However, this requires a huge support to meet the upcoming challenges that we are looking for, see [Support Us](/support_us.html).

## Installation

Install `reservoirflow` directly from [PyPi](https://pypi.org/):

```bash
pip install reservoirflow
```

For more information about the installation process, see: [Getting Started](/user_guide/getting_started/getting_started.html) in the documentation.

## Import Convention

The following convention is used to import `reservoirflow` after installation:

```python
import reservoirflow as rf
```

The abbreviation `rf` refers to `reservoirflow` where all modules under this library can be accessed. `rf` is also used throughout the documentation. We recommend our users to stick with this convention.

## Version

[Semantic Versioning](https://semver.org/) is used for the version numbers. Since this library is still under active development, `major=0` is used until the first stable version is reached. The first version `v0.1.0` was released on January 1, 2024. The current version is `v0.1.0`. To know about which features are currently supported, check [Capabilities](capabilities.html).

**Version History:**

| |**Version**  |**Status**         |**Release Date (dd.mm.yyyy)**  |
|-|-            |-                  |-                              |
|1|`v0.1.0`     |current version    |01.01.2024                     |
|2|`v0.1.1`     |*under development*|*ongoing*                      |

## License

This work is licensed under a [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License](https://creativecommons.org/licenses/by-nc-sa/4.0/). Detailed license can be found [here](https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode). The current license does not allow commercial use without an explicit authorization from the author. For commercial applications, please [contact](#contact) us.

## Contact

You can contact us directly on the following email: <reservoirflow.contact@gmail.com>.
