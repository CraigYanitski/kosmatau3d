SHELL = /bin/zsh

.PHONY: help grid comparison

help:
	@echo "Use this make file to either run the model grid or model comparison."

grid:
	@python "scripts/terminal scripts/run_model_grid.py"

comparison:
	@python "scripts/terminal scripts/run_comparison.py"
