SHELL = /bin/zsh

.PHONY: help grid comparison build upload test_upload

help:
	@echo "Use this make file to either run the model grid or model comparison."

grid:
	@python "scripts/terminal scripts/run_model_grid.py"

comparison:
	@python "scripts/terminal scripts/run_comparison.py"

build:
	@rm -rf dist
	@pip install --upgrade build
	@python -m build

# Ensure you update first the version number in `setup.py`, 
# `pyproject.py`, and `docs/_source/conf.py`
upload: build
	@twine upload --repository pypi dist/*
	@rm -rf dist

test_upload: build
	@twine upload --repository testpypi dist/*
	@rm -rf dist

