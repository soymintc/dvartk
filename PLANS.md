# Development plans :construction_worker:

## Guidelines
[ ] essentially make a SNV/SV comparison toolkit (with extension potential for Indels)
[ ] make empty functions and test loading
[ ] down-to-up development: endpoint functions first, then their inputs

## Details - SNV
[o] get maf as input
[-] get vcf as input
[x] get type column name (e.g. “TYPE”) as a parameter
[x] get variant type as input
[ ] allow plug-and-play filtering options
[ ] have functions that can be reused for MMCTM count file data
[ ] on a separate branch, make modules that make MMCTM input files

## Details - SV
[ ] get maf as input
[ ] get vcf as input
[ ] get type column name (e.g. “TYPE”) as a parameter
[ ] get variant type as input
[ ] allow plug-and-play filtering options
[ ] have functions that can be reused for MMCTM count file data
[ ] on a separate branch, make modules that make MMCTM input files

# Record of how this package was made
```
curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python -
source $HOME/.poetry/env
mkdir dvartk; cd dvartk
poetry new .
poetry install
poetry add -D pre-commit  # for PEP formatting
# remove README.rst and add README.md
poetry shell
# edit, git add
pre-commit run all-files
# git commit
poetry config repositories.testpypi https://test.pypi.org/legacy/
poetry config http-basic.testpypi __token__ $PYPITOKEN # in $HOME/.pypirc
poetry build
poetry publish -r testpypi
poetry add -D scriv[toml]  # for changelogs
vi pyproject.toml # add: [tool.scriv]\nformat = "md"
# would like to do pytest but bloody libffi version in CentOS is 6 not 7
scriv create
# then edit the changelog.d/*.md file
scriv collect  # creates CHANGELOG
poetry config pypi-token.pypi $PYPITOKEN

for package in $MY_PACKAGES; do poetry add $package; done
poetry install
poetry publish --build
```
