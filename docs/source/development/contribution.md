# Contributing to ReservoirFlow

## Code of Conduct

An excellent library needs two things: a clean code with a clean documentation. The author tried to follow most of the Python enhancement proposals mentioned in [PEP0](https://peps.python.org/pep-0000/). Mainly:

- [PEP20](https://peps.python.org/pep-0020/): guiding principles for Python’s design.
- [PEP8](https://peps.python.org/pep-0008/): a clean code that can be easily understood by others.
- [PEP257](https://peps.python.org/pep-0257/): a clean documentation.
- [PEP484](https://peps.python.org/pep-0484/): a standard syntax for type annotations.

Future developers and engineers who intend to contribute to this library should always keep in mind these standards to keep the source code clean and easy to follow. This should allow developing complex reservoir simulation models that are both innovative and efficient.

## Load Repository

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

## Contribute

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

- [Contribute](https://docs.github.com/en/get-started/quickstart/contributing-to-projects) by creating a [Pull Request](https://github.com/zakgrin/reservoirflow/pulls) to add features or solve bugs under [Issues](https://github.com/zakgrin/reservoirflow/issues). Please keep in mind that you need to stick with the project [Code of Conduct](#code-of-conduct).
