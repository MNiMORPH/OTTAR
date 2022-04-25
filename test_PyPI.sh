#! /bin/sh
twine upload --repository-url https://test.pypi.org/legacy/ dist/*
firefox https://test.pypi.org/project/OTTAR/
