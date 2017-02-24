#!/usr/bin/env python

from distutils.core import setup

setup(
    name='grond',
    description='What do you want to bust today?!',
    version='0.1',
    author='Sebastian Heimann',
    author_email='sebastian.heimann@gfz-potsdam.de',
    packages=['grond'],
    package_dir={'grond': 'src'},
    scripts=['apps/grond'],
    package_data={'grond': []},
    data_files=[('/etc/bash_completion.d', ['extras/grond'])])
