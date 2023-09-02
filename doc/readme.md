# Doc ReadMe


## Environment

You can create a separate environment for documentation using the following commands:

```bash
python -m venv .denv
source .denv/Scripts/Activate
pip install --upgrade pip
pip install -r requirements.txt
```

If you are on windows, you need to open __cmd__ and run the following command:

```bash
make.bat
```

Windows:

```bash
sphinx-quickstart
sphinx-apidoc -o . ..\reservoirflow -f
sphinx-apidoc -o doc reservoirflow/ -f --ext-autodoc
```
