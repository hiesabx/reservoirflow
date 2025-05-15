# Getting Started 🐤

## Installation

- Install Python:

    For this project, You need to install python. We recommend `python==3.11`. More download options can be found at the official website [python.org](https://www.python.org/downloads/release/python-379/). Make sure that you add Python to the `Path` variable of your operating system.

- Install `reservoirflow` directly from [PyPi](https://pypi.org/project/reservoirflow/):

    ```console
    $ pip install reservoirflow
    ```

    ```{important}
    You may need to create a project directory with a dedicated virtual environment. If you are new to the python ecosystem, follow the steps below. Note that the steps below require command-line tool. You can execute these commands in IDE, `cmd`, or `bash` terminals..
    ```

## Installation in Virtual Environment

- Create a new folder:

    ```console
    $ mkdir myrf
    ```

- Navigate to the project folder:

    ```console
    $ cd myrf
    ```

- Setup a Python virtual environment:

    ```console
    $ python -m venv .venv
    ```

- Activate your environment:

  - Windows:

      ```console
      $ .venv/Scripts/Activate.bat
      ```

  - Linux or Mac:

      ```console
      $ source .venv/Scripts/Activate
      ```

- Update `pip`:

    ```console
    $ pip install --upgrade pip
    ```

    or

    ```console
    $ python -m pip install --upgrade pip
    ```

- Install `reservoirflow` in your virtual environment:

    ```console
    $ pip install reservoirflow
    ```

    This library will also install `numpy`, `pandas`, and `pyvista`.


```{include} /_static/comments_section.md
```
