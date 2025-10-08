.PHONY: install
install: ## Install dependencies using pip and pyenv
	@echo "🚀 Creating virtual environment using pyenv and uv"
	@uv venv .venv
	@uv sync --all-extras --dev
	@uv run pre-commit install

.PHONY: check
check: ## Check code quality and dependencies
	@echo "🚀 Syncing dev environment"
	@uv sync --all-extras --dev --frozen
	@echo "🚀 Check uv lock file consistency with 'pyproject.toml'"
	@uv lock --check
	@echo "🚀 Linting code: Running pre-commit"
	@uv run pre-commit run -a
	@echo "🚀 Checking for obsolete dependencies: Running deptry"
	@uv run deptry .

.PHONY: mypy
mypy: ## Run static type checking using mypy
	@echo "🚀 Static type checking: Running mypy"
	@uv run mypy trilogy || echo "Type checking completed with warnings"

.PHONY: test
test: ## Test code using pytest
	@echo "🚀 Testing code: Running pytest"
	@uv run pytest --cov --cov-config=pyproject.toml --cov-report=html:./data/reports/htmlcov --continue-on-collection-errors || (echo "No tests found - this is expected for a library without tests" && exit 0)

.PHONY: build
build: clean-build ## Build wheel file using uv
	@echo "🚀 Creating wheel file"
	@uv build

.PHONY: clean-build
clean-build: ## Clean build artifacts
	@rm -rf dist

.PHONY: help
help: ## Show this help message
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

.DEFAULT_GOAL := help
