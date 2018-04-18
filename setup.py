#!/usr/bin/env python
import os
from setuptools import setup


def grond_completion():
    if os.access('/etc/bash_completion.d/', os.W_OK):
        return [('/etc/bash_completion.d', ['extras/grond'])]
    return []


setup(
    name='grond',
    description='What do you want to bust today?!',
    version='0.2',
    author='Sebastian Heimann',
    author_email='sebastian.heimann@gfz-potsdam.de',
    packages=[
        'grond',
        'grond.apps',
        'grond.targets',
        'grond.targets.waveform',
        'grond.targets.satellite',
        'grond.targets.gnss_campaign',
        'grond.problems',
        'grond.problems.cmt',
        'grond.problems.double_dc',
        'grond.problems.rectangular',
        'grond.optimizers',
        'grond.optimizers.highscore',
        'grond.analysers',
        'grond.listeners',
        'grond.report',
    ],
    python_requires='>=3.5',
    entry_points={
        'console_scripts': [
            'grond = grond.apps.grond:main',
        ]
    },
    package_dir={'grond': 'src'},
    package_data={'grond': ['baraddur/templates/*.html',
                            'baraddur/res/*',
                            'report/app/css/*.css',
                            'report/app/js/*.js',
                            'report/app/*.html']},
    data_files=[] + grond_completion(),
    )
