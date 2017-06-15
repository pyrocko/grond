#!/usr/bin/env python
import os
from distutils.core import setup


def grond_completion():
    if os.access('/etc/bash_completion.d/', os.W_OK):
        return [('/etc/bash_completion.d', ['extras/grond'])]
    return []


setup(
    name='grond',
    description='What do you want to bust today?!',
    version='0.1',
    author='Sebastian Heimann',
    author_email='sebastian.heimann@gfz-potsdam.de',
    packages=['grond', 'grond.baraddur', 'grond.problems', 'grond.solver'],
    scripts=['apps/grond'],
    package_dir={'grond': 'src'},
    package_data={'grond': ['baraddur/templates/*.html',
                            'baraddur/res/*']},
    data_files=[] + grond_completion(),
    )
