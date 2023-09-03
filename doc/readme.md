# Doc ReadMe

## Environment

You can create a separate environment for documentation using the following commands:

```bash
mrdir doc
cd doc
python -m venv .denv
source .denv/Scripts/Activate
pip install --upgrade pip
pip install -r requirements.txt
sphinx-quickstart
```

## Generate rst and html files

### bash

´´´bash
sphinx-apidoc -o modules/ ../reservoirflow/ -f -M -T
sphinx-build -b html ./ _build/html/
```

```
sphinx-apidoc -o doc/modules reservoirflow/ -f -M -T
sphinx-build -b html doc/ doc/_build/html/
```

### cmd 

´´´bash
sphinx-apidoc -o modules\ ..\reservoirflow\ -f -M -T
make.bat clean html | make.bat html #(69 warnings) or
sphinx-build -b html .\ _build\html # (69 warnings)

sphinx-apidoc -o doc\modules\ reservoirflow\ -f -M -T
doc\make.bat clean html | doc\make.bat html 
sphinx-build -b html doc\ doc\_build\html\ # (6 warnings)
```

## Run Documentation

at doc location run livereload by `python run_livereload.py`