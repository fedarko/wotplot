# From https://github.com/fedarko/strainFlye/blob/main/Makefile

.PHONY: test stylecheck style

# The -B in the invocation of python prevents this from creating pycache stuff.
# --cov-report xml is needed to make this visible to Codecov;
# --cov-report term is needed in order to print a fancy table on the terminal;
# --cov-report html creates a "htmlcov/" folder in the root of the repository.
#
# Also, for reference -- if you wanna see all print output during testing, add
# -s to the end of this command. Useful when debugging stuff.
test:
	python3 -B -m pytest \
		wotplot/ \
		--cov-report xml \
		--cov-report term \
		--cov-report html \
		--cov wotplot

stylecheck:
	@# Ignoring E203 per https://github.com/psf/black/issues/315
	flake8 --ignore=W503,E203 wotplot/ setup.py
	black --check -l 79 wotplot/ setup.py

style:
	black -l 79 wotplot/ setup.py
