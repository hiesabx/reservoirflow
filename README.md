# ReadMe

`ReservoirFlow: a Petroleum Reservoir Simulation Library in Python`

![](/images/logo.png)
<!--- 
Petroleum Reservoir Simulation using Scientific Computing and Machine Learning With Python developed by [Zakariya ABUGRIN](https://github.com/zakgrin).
--->

![](/images/five_spot_single_phase.gif)
*Example: Pressure Distribution of Single Phase Flow in Five Spot Wells Patterns*

**Table of Content:**

- [Introduction](#introduction)
- [Installation](#installation)
- [Import Convention](#import-convention)
- [Version](#version)
- [Capabilities](#capabilities)
- [Vision](#vision)
- [Motivation](#motivation)
- [Sponsorship](#sponsorship)
- [How Can I Support This Project?](#how-can-i-support-this-project)
- [Convention and Rules](#convention-and-rules)
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

- Install Python:

    For this project, `python==3.7.9` is recommended. To download the executable installer files directly:

  - Windows: [Windows x86-64 executable installer](https://www.python.org/ftp/python/3.7.9/python-3.7.9-amd64.exe)
  - macOS: [macOS 64-bit installer](https://www.python.org/ftp/python/3.7.9/python-3.7.9-macosx10.9.pkg)
  - Linux/UNIX: [Gzipped source tarball](https://www.python.org/ftp/python/3.7.9/Python-3.7.9.tgz).

  More download options can be found at the official website [python.org](https://www.python.org/downloads/release/python-379/).

- Load from GitHub (for developers):

    You can load project files directly from the official GitHub repository. Skip this step if you want to install this library using `pip` command directly.

  - Download [reservoirflow](https://github.com/zakgrin/reservoirflow) repository from GitHub:

    ```bash
    git clone https://github.com/zakgrin/reservoirflow.git
    ```

  - Navigate to the project folder:

    ```bash
    cd reservoirflow
    ```

  - Setup a Python virtual environment:

      ```bash
      python -m venv .venv
      ```

  - Activate your environment:

    - Windows:

        ```bash
        source .venv/Scripts/Activate
        ```

    - Linux or Mac:

        ```bash
        source .venv/bin/activate
        ```

  - Update `pip`:

      ```bash
      pip install --upgrade pip
      ```

  - Setup the library:

    ```bash
    pip install -e .
    ```

    or

    ```bash
    python setup.py install
    ```

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

## Capabilities

A long-term plan with high ambitions is already in place for this library. As of now, this tool is still very small relative to the entire scope of work (i.e. only about 5%). Supported capabilities for far are shown below:

|**Feature**|**Type**|**Support**|**Starting From**|**Comment**|
|-----------|--------|-----------|-----------------|-----------|
|**Model Type**|Black Oil|YES|`v0.1.0`|Only single phase on Cartesian grid.|
||Compositional|No|`-`||
|**Flow Dimensions**|1D|Yes|`v0.1.0`|Flow dimension is defined based on grid shape.|
||2D|Yes|`v0.1.0`|Flow dimension is defined based on grid shape.|
||3D|Yes|`v0.1.0`|Flow dimension is defined based on grid shape.|
|**Grid Types**|Regular Cartesian|Yes|`v0.1.0`||
||Single Well Radial|No|`-`||
||Cartesian|No|`-`||
|**Fluid Phases**| Single Phase|Yes|`v0.1.0`||
||Two Phases|No|`-`||
||Three Phases|No|`-`||
|**Fluid Compressibility**| Incompressible |Yes|`v0.1.0`||
|| Slightly Compressible|Yes|`v0.1.0`||
|| Compressible|No|`-`||
|**Rock Compressibility**| Incompressible |Yes|`v0.1.0`||
|| Slightly Compressible|Yes|`v0.1.0`||
|**Experiments**|Core-Flooding|No|`-`||
|| Slim-Tube|No|`-`||
|**Well Types**|Single Cell|Yes|`v0.1.0`||
||Multiple Cells|No|`-`||
||Horizontal|No|`-`||
||Slanted|No|`-`||
|**History Matching**|Conventional|No| `-`||
||Machine Learning|No|`-`||
|**Production Optimization**|Reinforcement Learning|No|`-`||
|**Solutions**|Numerical [FDM]|Yes|`v0.1.0`|Requires an `Initializer` and a `Solver`|
||Analytical [1D]|Yes|`v0.1.0`|Available only for `1D Single Phase` problems.|
||Neurical [1D]|Yes|`v0.1.0`|Neural networks based solution (e.g., PINNs).|
|**Numerical Initializers**|Vectorized|Yes|`v0.1.0`|Used for efficient computing.|
||Symbolic|Yes|`v0.1.0`|Used for debugging.|
|**Numerical Solvers**|Iterative Solvers|Yes|`v0.1.0`|Requires sparse matrices.|
||Neurical Solvers*|No|`-`|Solving a system of linear equations using neural networks.|
|**Efficient Computing**|Sparse Matrices|Yes|`v0.1.0`|Required for iterative solvers.|
||Threading|Yes|`v0.1.0`|Used for numerical symbolic initializers.|
||Concurrent|Yes|`v0.1.0`|Used for numerical symbolic initializers.|
||GPU Computing|No|`-`||
||TPU Computing|No|`-`||
||Quantum Computing|No|`-`||
|

(*) Innovative solution introduced within this work.

## Vision

This library brings reservoir simulation to the Python ecosystem including data analytics and machine learning tools which is in itself a huge benefit. Having this option allows for better integration and support with intelligent fields (e.g. software containers based on this library can be built and deployed to improve production and reservoir management such as real-time well production optimization). In contrast to commercial software where reservoir simulation studies are relatively isolated, important actions can be immediately predicted and made available for the field hardware to execute.

## Motivation

Accessible and reproducible research is very scarce in the oil & gas industry. This unfortunate reality led many reservoir engineers to be limited and dependent on commercial tools that are closed black-boxes which led to severe limitations when it comes to skills and innovation. For example, most reservoir engineers today can not easily express their scientific knowledge by coding and normally are not trained to do so. This perhaps is the main reason why the oil & gas industry is still behind in the AI revolution.

The author believes that this sad reality has to be changed as soon as possible. Reservoir engineers face difficult challenges in getting careers in or outside the oil & gas industry mainly due to the lack of soft skills. Otherwise, a reservoir engineer, who has a great mix of science (i.e. mathematics, statistics, scientific computing, programming, chemistry, physics, fluid mechanics and thermodynamics, geology, etc.), should have no reason to stay unemployed.

This sad reality can be changed if reservoir engineers are trained to express their scientific and engineering knowledge with coding. In addition, this will bring more excitement and job satisfaction. Those who learn to do reservoir simulation in Python should be familiar with the scientific computing and machine learning stack in Python. Therefore, it will be easy for them to change careers and work as data scientists, data analysts, or ML engineers in other industries.

On a larger scale, an open-source reservoir simulation library within the Python ecosystem will be very important to students, universities, and researchers. Therefore, the growth of this tool can only be taken as a positive growth for a new and healthy oil & gas community that we try to create. However, this brings us to the next topic of [Sponsorship](#sponsorship).

Unfortunately, this work is not sponsored. As a result, the current license does not allow for commercial applications. Companies who are interested to use and support this tool for commercial applications must get in touch to get a proper license.

## Sponsorship

Any open-source tool requires resources that cost both time and money. Although much of open-source work is driven by the good intention of having accessible and reproducible research, open-source has also proved to be a very successful business model especially when the tool grows to a production level. Compare for example Python with other commercial tools (e.g. Matlab). Python is open-source but also much more successful!

To reach a stable production level, we will need enough investment to develop, upgrade, and maintain this tool. Here comes the role of oil & gas companies (i.e. especially operation companies) who may benefit a lot from this tool to sponsor this project. This will help us to dedicate more resources and to add new features (i.e. let it be Quantum Computing!) or fix some bugs.

If you are an individual who may benefit from this work, consider supporting this project by [Sponsorship on Patreon](https://www.patreon.com/zakgrin) or by any means. Even by just starting to use this tool and reporting some errors or insights. If you can't then that is fine. Our main duty is to always keep this tool open-source.

## How Can I Support This Project?

There are two ways to support this project. One way is by [Sponsorship](#sponsorship). The other way is to offer some technical support as explained below.

**Normal Users:**

As a normal user, you can support this project in many ways:

- Give a star in GitHub for [`reservoirflow`](https://github.com/zakgrin/reservoirflow) (you need a GitHub Account).
- Start using `reservoirflow` for your projects or thesis.
- Compare `reservoirflow` with other commercial tools and openly share your feedback.
- Share your progress and experience on [LinkedIn](https://www.linkedin.com/feed/) or any other platform. Use `#reservoirflow` hashtag, so the community can interact with your work. Mention [Zakariya Abugrin](https://www.linkedin.com/in/zakariya-abugrin/), so the author can give you quick feedback.
- Report some issues you face under [Issues](https://github.com/zakgrin/reservoirflow/issues).

**Developers:**

If you are a developer or want to be one, you need additionally to:

- [Fork](https://docs.github.com/en/get-started/quickstart/fork-a-repo) this repository in your machine and start learning about the source-code.
- Create your feature branch:

    ```Bash
    git checkout -b my_new_feature
    ```

- Commit your changes:

    ```Bash
    git commit -am "Add my new feature"
    ```

- Push to the branch:

    ```Bash
    git push origin my_new_feature
    ```

- [Contribute](https://docs.github.com/en/get-started/quickstart/contributing-to-projects) by creating a [Pull Request](https://github.com/zakgrin/reservoirflow/pulls) to add features or solve bugs under [Issues](https://github.com/zakgrin/reservoirflow/issues). Please keep in mind that you need to stick with the project [Convention and Rules](#convention-and-rules).

**Companies:**

If you are a company that may use this tool for commercial application, then the best way is to [Sponsorship](#sponsorship) this project. Additionally, you can do the following:

- Use this tool for your internal projects.
- In case new features are needed, then open a new request under [issues](https://github.com/zakgrin/reservoirflow/issues). In contrast, we also expect some [Sponsorship](#sponsorship), so we have the required resources to develop this feature specifically for you.
- Encourage your RE teams to participate in this project as developers and not only as normal users. If you do not have any, hire one who can do so.
- Use our future offers that we will announce later once the community grows enough (i.e. training, consultation, special use cases, etc.).

## Convention and Rules

An excellent library needs two things: clean code and clean documentation. The author tried to follow most of the Python enhancement proposals mentioned in [PEP0](https://peps.python.org/pep-0000/). Mainly:

- [PEP20](https://peps.python.org/pep-0020/): guiding principles for Python’s design.
- [PEP8](https://peps.python.org/pep-0008/): a clean code that can be easily understood by others.
- [PEP257](https://peps.python.org/pep-0257/): a clean documentation.
- [PEP484](https://peps.python.org/pep-0484/): a standard syntax for type annotations.

Future developers and engineers who intend to contribute to this library should always keep in mind these standards to keep the source code clean and easy to follow. This should allow developing complex reservoir simulation models that are both innovative and efficient.

## License

This work is licensed under a [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License](https://creativecommons.org/licenses/by-nc-sa/4.0/). Detailed license can be found [here](https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode).

The current license does not allow commercial use. The author would like to change the license to [BSD 3-Clause License](https://opensource.org/licenses/BSD-3-Clause) which allows for commercial use. However, this depends on the level of sponsorship and support that is offered for this project especially from oil and gas operation companies.
