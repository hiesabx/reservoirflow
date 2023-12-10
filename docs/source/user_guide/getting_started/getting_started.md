# Getting Started

## Installation

- Install `reservoirflow` directly from [PyPi](https://pypi.org/):

    ```bash
    pip install reservoirflow
    ```

    You need to install python. In addition, you may need to create a project directory with a dedicated virtual environment. If you are new to the python ecosystem, follow the steps below. Note that the steps below require command-line tool. You can execute these commands in IDE, `cmd`, or `bash` terminals.

- Install Python:

    For this project, `python>=3.8` is recommended. More download options can be found at the official website [python.org](https://www.python.org/downloads/release/python-379/). Make sure that you add Python to the `Path` variable of your operating system.

- Create a new folder:

    ```bash
    mkdir myrf
    ```

- Navigate to the project folder:

    ```bash
    cd myrf
    ```

- Setup a Python virtual environment:

    ```bash
    python -m venv .venv
    ```

- Activate your environment:

  - Windows:

      ```bash
      .venv/Scripts/Activate.bat
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
    pip install reservoirflow
    ```

    This library will also install `numpy`, `pandas`, and `pyvista`.



```{include} /_static/comments_section.md
```
