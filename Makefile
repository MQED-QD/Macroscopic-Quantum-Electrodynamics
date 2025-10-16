# This is a special target that lists commands that are not files.
# It's good practice to include all your command targets here.
.PHONY: docs docs-clean

##@ Documentation

docs: ## Generate the HTML documentation from the docs/ directory
	@echo ">>> Generating Sphinx HTML documentation..."
	$(MAKE) -C docs html

docs-clean: ## Remove the built documentation directory
	@echo ">>> Cleaning documentation build directory..."
	$(MAKE) -C docs clean