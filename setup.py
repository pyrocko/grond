#!/usr/bin/env python

from distutils.core import setup

setup(
    name='grond',
    description='What do you want to bust today?!',
    version='0.1',
    author='Sebastian Heimann',
    author_email='sebastian.heimann@gfz-potsdam.de',
    packages=['grond', 'grond.baraddur'],
    scripts=['apps/grond'],
    package_dir={'grond': 'src'},
    package_data={'grond': [],
                  'grond': ['baraddur/templates/*.html',
                            'baraddur/res/*']},
    data_files=[('/etc/bash_completion.d', ['extras/grond'])],
    )
