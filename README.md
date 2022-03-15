# Open Reservoir Simulation Library <`openresim`>

<!--- 
Petroleum Reservoir Simulation using Scientific Computing and Machine Learning With Python developed by [Zakariya ABUGRIN](https://github.com/zakgrin).
--->

## Introduction

Open Reservoir Simulation Library, shortly `openresim`, is a modern open-source python library created by [Zakariya ABUGRIN](https://github.com/zakgrin). `openresim` helps petroleum reservoir engineers to be proficient scientists and engineers who are able to express their practical and scientific knowledge by coding.

This library is designed based on the modern Python stack for data science and scientific computing. Hence, this library is directly integrated with AI and ML Python frameworks and libraries.

`openresim` is designed to achieve the highest quality use cases and research for reservoir engineering and simulation. Solutions that can combine the strength scientific computing with the power of machine learning including the state-of-the-art deep learning models. Below are some examples for the problems that will be tackled in future:

- Reinforcement Learning to achieve better production strategies.
- Deep Learning to perform history matching.

## Version

[Semantic Versioning](https://semver.org/) is used for the version numbers. Since this library is still under development, `major=0` is used until a stable version is reached.

The first version `v0.1.0` was released on April 1, 2022. Supported capabilities are shown below:

| **Feature**         | **Type**       | **Support** | **Starting From** |
| --------------------- | ---------------- | ----------- | ------------- |
| **Grid Type**       | Cartesian      | Yes         | `v0.1.0`    |
|                     | Radial         | No          | `-`         |
| **Dimension**       | 1D             | Yes         | `v0.1.0`    |
|                     | 2D             | No          | `v0.1.0`    |
|                     | 3D             | No          | `v0.1.0`    |
| **Phases**          | Single Phase   | Yes         | `v0.1.0`    |
|                     | Two Phases     | No          | `-`         |
|                     | Three Phases   | No          | `-`         |
|                     | Compositional  | No          | `-`         |
| **Compressibility** | Incompressible | Yes         | `v0.1.0`    |
|                     | Compressible   | No          | `v0.1.0`    |
| **Experiments**     | Core-Flooding  | No          | `-`         |
|                     | Slim-Tube      | No          | `-`         |
| **History Matching**| Conventional   | No          | `-`         |
|                     | Machine Learning | No        | `-`         |
| **Optimization**    |Reinforcement Learning | No   | `-`         |

## Convention

An excellent library needs two things: a clean code and a clean documentation. The author tried to follow most of Python enhancement proposals mentioned in [PEP0](https://peps.python.org/pep-0000/). Mainly:

- [PEP8](https://peps.python.org/pep-0008/) to produce a clean code that can be easily understood by others (PEP8 is given the highest priority).
- [PEP257](https://peps.python.org/pep-0257/) to produce a clean documentation.

The author follows the intuitive pythonic way ([PEP20](https://peps.python.org/pep-0020/)) to design this library. This should allow to build complex reservoir simulation models that are both innovative and efficient.

Future developers and engineers who intend to contribute to this library should always keep in mind this convention to keep the source code clean and easy to follow.

### Import Convention

The following convention is used to import `openresim` after installation:

```python
import openresim as rs
```

By this, engineers should remember that `rs` refers to __Reservoir Simulation__ and it is the official shortname for `openresim`.

## Motivation

Open-source and reproducible research is very scarce in the oil&gas industry. This unfortunate reality has lead many reservoir engineers to be limited and dependant to commercial tools that are closed black-boxes. This has led to severe limitations when it comes to skills and innovation. For example, engineers can not easily express their scientific knowledge by coding and normally are not trained to do so. 

The author believes that this sad reality is the main reason why reservoir engineers face difficult challenges to get careers outside of the oil&gas industry. Otherwise, a reservoir engineer, who has a great mix of science (i.e. mathematics, statistics, scientific computing, chemistry, physics, fluid mechanics and thermodynamics, geology, etc.), should have no reason to stay unemployed.

This sad reality can be changed if reservoir engineers are trained to express their scientific and engineering knowledge with coding. Those who learn to do reservoir simulation in python should be familiar with the python stack for scientific computing and machine learning. Therefore, it will be easy for them to change careers and work as data scientists, data analysts, or ML engineer in other industries.
## Installation

### Install python

### Install the python library

```bash
git clone https://github.com/zakgrin/openresim.git
cd openresim
```


### Setup a python environment

```bash
python -m venv .venv
```

# Choose one option:

# 1. For windows

```bash
source .venv/Scripts/activate
```

# 2. For Linux or Mac

```bash
source .venv/Scripts/activate
pip install -r requirements.txt
pip install .
```

## Tutorials

### Example 1

![](images/example_1_code.png)
![](images/example_1_3d.png)
