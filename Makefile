SHELL = /bin/zsh

.PHONY: help grid comparison build upload

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

upload: build
	@twine upload --repository pypi dist/*
	@rm -rf dist
