# Open Reservoir Simulation Library <`openresim`>

<!--- 
Petroleum Reservoir Simulation using Scientific Computing and Machine Learning With Python developed by [Zakariya ABUGRIN](https://github.com/zakgrin).
--->

## Introduction

Open Reservoir Simulation Library, shortly `openresim`, is a modern open-source python library created by [Zakariya ABUGRIN](https://github.com/zakgrin). `openresim` helps reservoir engineers to be proficient scientists and engineers who are able to express their practical and scientific knowledge by coding.

This library is designed based on the modern Python stack for data science and scientific computing. Hence, this library is directly integrated with AI and ML Python frameworks and libraries.

`openresim` is designed to open the door towards achieving the highest quality use cases and research for reservoir engineering and simulation. Solutions that can combine the strength scientific computing with the power of machine learning including the state-of-the-art deep learning models. Below are some examples for the problems that will be tackled in future:

- History matching using machine learning.
- Reinforcement learning to achieve better production strategies for specific goals (e.g. maximize recovery, accelerate production).
- Physics-informed neural networks and deep learning models.
- Advanced computing such as: GPU Computing, Quantum Computing, and Cluod Computing. 

## Version

[Semantic Versioning](https://semver.org/) is used for the version numbers. Since this library is still under development, `major=0` is used until a stable version is reached. The first version `v0.1.0` was released on April 1, 2022. Current version is `v0.1.0`. Supported capabilities are shown below:

| **Feature**          | **Type**       | **Support** | **Starting From** |
| --------------------- | ---------------- | ----------- | ------------- |
| **Grid Type**        | Cartesian      | Yes         | `v0.1.0`    |
|                      | Radial         | No          | `-`         |
| **Dimension**        | 1D             | Yes         | `v0.1.0`    |
|                      | 2D             | Yes         | `v0.1.0`    |
|                      | 3D             | Yes         | `v0.1.0`    |
| **Phases**           | Single Phase   | Yes         | `v0.1.0`    |
|                      | Two Phases     | No          | `-`         |
|                      | Three Phases   | No          | `-`         |
|                      | Compositional  | No          | `-`         |
| **Compressibility**  | Incompressible | Yes         | `v0.1.0`    |
|                      | Compressible   | No          | `v0.1.0`    |
| **Experiments**      | Core-Flooding  | No          | `-`         |
|                      | Slim-Tube      | No          | `-`         |
| **History Matching** | Conventional   | No          | `-`         |
|                      | Machine Learning | No        | `-`         |
| **Optimization**     |Reinforcement Learning | No   | `-`         |
| **Quantum Computing**|                | No          | `-`         |

## Motivation

Open-source and reproducible research is very scarce in the oil&gas industry. This unfortunate reality has lead many reservoir engineers to be limited and dependant on commercial tools that are closed black-boxes. This has led to severe limitations when it comes to skills and innovation. For example, reservoir engineers can not easily express their scientific knowledge by coding and normally are not trained to do so. This perhaps is the main reason why oil&gas industry is still behind in the AI revolution.

The author believes that this sad reality has to be changed as soon as possible. Reservoir engineers face difficult challenges to get careers in or outside of the oil&gas industry mainly due to the lack of soft-skills. Otherwise, a reservoir engineer, who has a great mix of science (i.e. mathematics, statistics, scientific computing, programming, chemistry, physics, fluid mechanics and thermodynamics, geology, etc.), should have no reason to stay unemployed.

This sad reality can be changed if reservoir engineers are trained to express their scientific and engineering knowledge with coding. In addition, this will bring more excitements and job satisfaction. Those who learn to do reservoir simulation in Python should be familiar with the scientific computing and machine learning stack in Python. Therefore, it will be easy for them to change careers and work as data scientists, data analysts, or ML engineer in other industries.

On a larger scale, an open-source reservoir simulation library within Python eco-system will be very important to students, universities, researches or even service companies who can use this tool for commercial studies and applications. Therefore, the growth of this tool can only be taken as a positive growth for a new and healthy oil&gas community that we try to create. However, this brings us to the next topic of [Sponsorship](#sponsorship).

## Sponsorship

Any open-source tool requires resources that cost both time and money. Although much of open-source work is driven by the good intension of having an accessible and reproducible research, open-source have also proved to be a very successful business model especially when the tool grows to a production level. Compare for example Python with other commercial tools (e.g. Matlab). Python is open-source but also much more successful!

To reach a stable and production level, we will need enough investment to develop, upgrade, and maintain this tool. Here comes the role of oil&gas companies (i.e. especially operation companies) who may benefit a lot from this tool to [Sponsorship on Patreon](https://www.patreon.com/zakgrin) this project. This will help us as a community to dedicate more resources and to add new features (i.e. let it be Quantum Computing!) or fix some bugs.

If you are an individual who may benefit from this work, consider supporting this project by [Sponsorship on Patreon](https://www.patreon.com/zakgrin) or by any means. Even by just starting using this tool and reporting some errors or insights. If you can't then that is fine. Our main duty is to always keep this tool 100% free and open-source.

## How Can I support This Project?

There are two ways to support this project. One way is by [Sponsorship](#sponsorship). The other way is to offer some technical support as explained below.

**Normal Users:**

As normal user, you can support this project in many ways:

- Give a star in GitHub for [`openresim`](https://github.com/zakgrin/openresim) (you need a GitHub Account).
- Start using `openresim` for your projects or thesis.
- Compare `openresim` with other commercial tools and openly share your feedback.
- Share your progress and experience in [LinkedIn](https://www.linkedin.com/feed/) or any other platfrom. Use `#openresim` hashtag so the community can interact with your work. Mention [Zakariya ABUGRIN](https://www.linkedin.com/in/zakariya-abugrin-45306987/) the author can give you a quick feedback.
- Report some issues you face under [Issues](https://github.com/zakgrin/openresim/issues).

**Developers:**

If you are a developer or want to be one, you need additionally to:

- [Fork](https://docs.github.com/en/get-started/quickstart/fork-a-repo) this repository in your machine and start learning about the source-code.
- [Contribute](https://docs.github.com/en/get-started/quickstart/contributing-to-projects) by creating a [Pull Request](https://github.com/zakgrin/openresim/pulls) to add features or solve bugs under [issues](https://github.com/zakgrin/openresim/issues). Please keep in mind that you need stick with the project [Convention and Rules](#convention-and-rules).

**Companies:**

If you are a company who may use this tool for commercial application, then the best way is to [Sponsorship](#sponsorship) this project. Additionally you can do the following:

- Use this tool for your internal projects.
- In case new features are needed, then open a new request under [issues](https://github.com/zakgrin/openresim/issues). In contrast, we also expect some [Sponsorship](#sponsorship) so we have the required resources to develop this feature specifically for you. 
- Encourage your RE teams to participate in this project as developers and not only as normal users. If you do not have any, hire one who can do so.
- Use our future offers that we will announce latter once the community grow enough (i.e. trainings, consultation, special use cases etc).

## Installation

### Install Python

### Install `openresim`

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

## Convention and Rules

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
