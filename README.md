# ReservoirFlow: a Petroleum Reservoir Simulation Library in Python

![](/images/logo.png)
![](/images/five_spot_single_phase.gif)
*Example: Pressure Distribution of Single Phase Flow in Five Spot Wells Patterns*

**Table of Content:**

- [Introduction](#introduction)
- [Installation](#installation)
- [Import Convention](#import-convention)
- [Version](#version)
- [Vision](#vision)
- [Motivation](#motivation)
- [License](#license)

## Introduction

_ReservoirFlow_ <`reservoirflow`>, is a modern open-source Python library created by [Zakariya Abugrin](https://github.com/zakgrin). _ReservoirFlow_ helps mathematicians, researchers, and engineers to study and model the phenomenon of fluid flow in porous media using different approaches including analytical solutions, numerical solutions, and machine learning solutions. Studying different solution techniques for solving Partial Differential Equations (PDEs) gives great insight into developing more efficient solutions.

This library is designed based on the modern Python stack for data science and scientific computing. Hence, this library is directly integrated with AI and ML Python frameworks and libraries.

_ReservoirFlow_ is designed to open the door towards achieving the highest quality use cases and research for reservoir engineering and simulation. Solutions that can combine the strength of scientific computing with the power of machine learning including the state-of-the-art deep learning models. Below are some examples of the problems that will be tackled in the future:

- History matching using machine learning.
- Reinforcement learning to achieve better production strategies for specific goals (e.g. maximize recovery, accelerate production).
- Scientific Machine learning: using Physics-informed neural networks.
- Advanced computing such as GPU Computing, and Quantum Computing.
- Real-time reservoir management and production optimization using Cloud Computing and IoT.

The author aims to open a start-up that can provide accessible cutting-edge software and engineering solutions with hands-on training or consultation. For more information, please contact us at `reservoirflow@gmail.com`.

## Installation

- Install `reservoirflow` directly from [PyPi](https://pypi.org/):

    ```bash
    pip install reservoirflow
    ```

    You need to install python (see next step). In addition, You may need to create a project directory with a dedicated virtual environment.

## Import Convention

The following convention is used to import `reservoirflow` after installation:

```python
import reservoirflow as rf
```

The abbreviation `rf` refers to `reservoirflow` where all modules under this library can be accessed. `rf` is also used throughout the documentation. We recommend our users to stick with this convention.

## Version

[Semantic Versioning](https://semver.org/) is used for the version numbers. Since this library is still under development, `major=0` is used until the first stable version is reached. The first version `v0.1.0` was released on January 1, 2024. The current version is `v0.1.0`.

**Version History:**

||**Version**|**Status**|**Release Date (dd.mm.yyyy)**|
|-|-|-|-|
|1|`v0.1.0`|current version|01.01.2024|
|2|`v0.1.1`|*under development*|*ongoing*|

## Vision

This library brings reservoir simulation to the Python ecosystem including data analytics and machine learning tools which is in itself a huge benefit. Having this option allows for better integration and support with intelligent fields (e.g. software containers based on this library can be built and deployed to improve production and reservoir management such as real-time well production optimization). In contrast to commercial software where reservoir simulation studies are relatively isolated, important actions can be immediately predicted and made available for the field hardware to execute.

## Motivation

Accessible and reproducible research is very scarce in the oil & gas industry. This unfortunate reality led many reservoir engineers to be limited and dependent on commercial tools that are closed black-boxes which led to severe limitations when it comes to skills and innovation. For example, most reservoir engineers today can not easily express their scientific knowledge by coding and normally are not trained to do so. This perhaps is the main reason why the oil & gas industry is still behind in the AI revolution.

The author believes that this sad reality has to be changed as soon as possible. Reservoir engineers face difficult challenges in getting careers in or outside the oil & gas industry mainly due to the lack of soft skills. Otherwise, a reservoir engineer, who has a great mix of science (i.e. mathematics, statistics, scientific computing, programming, chemistry, physics, fluid mechanics and thermodynamics, geology, etc.), should have no reason to stay unemployed.

This sad reality can be changed if reservoir engineers are trained to express their scientific and engineering knowledge with coding. In addition, this will bring more excitement and job satisfaction. Those who learn to do reservoir simulation in Python should be familiar with the scientific computing and machine learning stack in Python. Therefore, it will be easy for them to change careers and work as data scientists, data analysts, or ML engineers in other industries.

On a larger scale, an open-source reservoir simulation library within the Python ecosystem will be very important to students, universities, and researchers. Therefore, the growth of this tool can only be taken as a positive growth for a new and healthy oil & gas community that we try to create. However, this brings us to the next topic of [Sponsorship](#sponsorship).

Unfortunately, this work is not sponsored. As a result, the current license does not allow for commercial applications. Companies who are interested to use and support this tool for commercial applications must get in touch to get a proper license.

## License

This work is licensed under a [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License](https://creativecommons.org/licenses/by-nc-sa/4.0/). Detailed license can be found [here](https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode).

The current license does not allow commercial use. The author would like to change the license to [BSD 3-Clause License](https://opensource.org/licenses/BSD-3-Clause) which allows for commercial use. However, this depends on the level of sponsorship and support that is offered for this project especially from oil and gas operation companies.
