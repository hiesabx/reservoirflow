# ReservoirFlow: a Petroleum Reservoir Simulation Library in Python

![logo](/images/logo.png)

![five_spot_single_phase](/images/five_spot_single_phase.gif)\\
*Example: Pressure Distribution of Single Phase Flow in Five Spot Wells Patterns*

**Table of Content:**

- [Introduction](#introduction)
- [Installation](#installation)
- [Import Convention](#import-convention)
- [Version](#version)
- [Mission](#mission)
- [Vision](#vision)
- [Philosophy](#philosophy)
- [Are You a Reservoir Engineer?](#are-you-a-reservoir-engineer)
- [License](#license)

## Introduction

*ReservoirFlow* <`reservoirflow`>, is a modern open-source Python library created by [Zakariya Abugrin](https://github.com/zakgrin). *ReservoirFlow* helps mathematicians, researchers, students and engineers to study and model the phenomenon of fluid flow in porous media using different approaches including analytical solutions, numerical solutions, and machine learning solutions. This library is designed based on the modern Python stack for data science and scientific computing. Hence, this library is directly integrated with AI and ML Python frameworks and libraries.

## Installation

Install `reservoirflow` directly from [PyPi](https://pypi.org/):

```bash
pip install reservoirflow
```

For more information about the installation process, see: [Getting Started](/user_guide/getting_started.html) in the documentation.

## Import Convention

The following convention is used to import `reservoirflow` after installation:

```python
import reservoirflow as rf
```

The abbreviation `rf` refers to `reservoirflow` where all modules under this library can be accessed. `rf` is also used throughout the documentation. We recommend our users to stick with this convention.

## Version

[Semantic Versioning](https://semver.org/) is used for the version numbers. Since this library is still under active development, `major=0` is used until the first stable version is reached. The first version `v0.1.0` was released on January 1, 2024. The current version is `v0.1.0`. To know about which features are currently supported, check [Capabilities](capabilities.html).

**Version History:**

||**Version**|**Status**|**Release Date (dd.mm.yyyy)**|
|-|-|-|-|
|1|`v0.1.0`|current version|01.01.2024|
|2|`v0.1.1`|*under development*|*ongoing*|

## Mission

*Bayanatics GmbH* is formed with a mission to provide accessible cutting-edge software, AI, and engineering solutions with hands-on training, consultation, and collaboration with both companies and universities. *Bayanatics* is not limited to oil & gas industry. Our aim is to be a central platform between education and industry where technological progress is shared, standardized, cleansed, and redistributed in a form of high quality open-source tools supported with an educational content.

A large enough pool of collaboration with universities and companies throughout our platform will be utilized to provide not only a better education but also a unified and shared technological progress. Companies who choose to sponsor our work, will get additional benefits including commercial use, shape the development of our tools, and make requests for custom features.

Collaboration with universities is essential ideally with some funding but not necessary as long as we have enough resources to collaborate. Research outcomes are used to update our tools which then are made available to universities to continue the development cycle. These tools then can be used in lectures to give students hands-on training and also to use them in their research or theses.

If this mission means a lot to you, then please see [Support Us](/support_us.html).

## Vision

*ReservoirFlow* represents our first trial to build an open-source tool with high quality educational content for the topic of fluid flow in porous media mainly related to subsurface reservoirs (e.g., petroleum reservoirs). *ReservoirFlow* brings reservoir simulation to the Python ecosystem including data analytics and machine learning tools which is in itself a huge benefit. Having this option allows for better integration and support with intelligent fields (e.g., software containers based on this library can be built and deployed to improve production and reservoir management such as real-time well production optimization). In contrast to commercial software where reservoir simulation studies are relatively isolated, important actions can be immediately predicted and made available for the field hardware to execute.

*ReservoirFlow* is designed to open the door towards achieving a high quality research and use cases for reservoir engineering and simulation. Solutions that can combine the strength of scientific computing with the power of machine learning including the state-of-the-art deep learning models. Below are few examples of the problems that will be tackled in the future:

- History matching using machine learning.
- Reinforcement learning to achieve better production strategies for specific goals (e.g. maximize recovery, accelerate production).
- Scientific Machine learning using Physics-informed neural networks (PINNs).
- Advanced computing such as GPU Computing, and Quantum Computing.
- Real-time reservoir management and production optimization using Cloud Computing and IoT.

An open-source reservoir simulation library within the Python ecosystem will be very important to students, universities, and researchers. Therefore, the growth of this tool can only be taken as a positive growth for a new and healthy oil & gas community that we try to create. However, this requires a huge support to meet the upcoming challenges that we are looking for, see [Support Us](/support_us.html).

## Philosophy

Our Philosophy is to try to attack the same problem from multiple angles by combining different solutions including analytical solutions, numerical solutions, and machine learning solutions. Studying different solution techniques to solve Partial Differential Equations (PDEs) gives great insight into developing more efficient solutions. However, this requires a combination of wide range of topics especially Mathematics, Physics, and Engineering, or what we shortly refer to as *Mathephysineering*. With the objective to improve the practical applications, *Mathephysineering* can improve our overall understanding especially with the power of other tools available from computer science and machine learning. *Mathephysineering* inspires us to revisit the basics of mathematics required to better describe and solve physical problems according to best engineering practices. This is the philosophy behind developing this tool as will be demonstrated within this project.

## Are You a Reservoir Engineer?

Accessible and reproducible research is very scarce in the oil & gas industry. This unfortunate reality led many reservoir engineers to be limited and dependent on commercial tools that are closed black-boxes which led to severe limitations when it comes to skills and innovation. For example, most reservoir engineers today can not easily express their scientific knowledge by coding and normally are not trained to do so. This perhaps is the main reason why the oil & gas industry is still behind in the AI revolution.

The author believes that this sad reality has to be changed as soon as possible. Reservoir engineers face difficult challenges in getting careers in or outside the oil & gas industry mainly due to the lack of soft skills. Otherwise, a reservoir engineer, who has a great mix of science (i.e. mathematics, statistics, scientific computing, programming, chemistry, physics, fluid mechanics and thermodynamics, geology, etc.), should have no reason to stay unemployed.

This sad reality can be changed if reservoir engineers are trained to express their scientific and engineering knowledge with coding. In addition, this will bring more excitement and job satisfaction. Those who learn to do reservoir simulation in Python should be familiar with the scientific computing and machine learning stack in Python. Therefore, it will be easier for them to change careers and work as data scientists, data analysts, or ML engineers in other industries.

## License

This work is licensed under a [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License](https://creativecommons.org/licenses/by-nc-sa/4.0/). Detailed license can be found [here](https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode). The current license does not allow commercial use.
