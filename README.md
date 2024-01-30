# Glenn_Bfx_repo README
## Tools


## Setup
- Clone repository
  - git clone git@github.com:glenncarson/Glenn-Bfx_Repo
- Download python 3.9.18 (Alternatively, use pyenv)
  - `wget https://www.python.org/ftp/python/3.9.18/Python-3.9.18.tgz`
- Extract python
  - `tar -xvf Python-3.9.18.tgz`
  - `/configure `
  - `make `
  - `make test `
  -` sudo make install`
- find path to python
  - `which python3`
- install poetry
  - `curl -sSL https://install.python-poetry.org | <path-to-python3> - --version 1.7.1`
- navigate to Glenn-Bfx-Repo
  - `cd <Glenn-Bfx-Repo>`
- Create poetry virtual env
  - `poetry shell`
  - Note: use this command to launch the ve in the future
    - exit with: `exit` command
- Install dependencies
  - `poetry install`
- Setup python interpreter for PyCharm
  - Find virtualenv path
    - `which python`
  - navigate to venv dir
    - e.g., `cd /Users/glenn.carson/Coding_Workspace/PycharmProjects/Glenn_Bfx_Repo/.venv`
  - create an excecutable to launch the interpreter
    - `nano python_interp.sh`
    - copy paste these contents:
    - `#!/bin/bash`
    - `source <venv dir>/bin/activate`
      - e.g., `<venv dir>/bin/python "$@"`
    - make the file executable: 
      - `chmod +x python_interp.sh`
    - test the file:
      - `./python_interp.sh`
    - 