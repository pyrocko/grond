#!/usr/bin/env python
from setuptools import setup


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
        'grond.targets.waveform_phase_ratio',
        'grond.targets.waveform_oac',
        'grond.targets.satellite',
        'grond.targets.gnss_campaign',
        'grond.problems',
        'grond.problems.cmt',
        'grond.problems.double_dc',
        'grond.problems.rectangular',
        'grond.optimisers',
        'grond.optimisers.highscore',
        'grond.analysers',
        'grond.analysers.noise_analyser',
        'grond.analysers.target_balancing',
        'grond.report',
        'grond.plot',
    ],
    python_requires='>=3.5',
    entry_points={
        'console_scripts': [
            'grond = grond.apps.grond:main',
        ]
    },
    package_dir={'grond': 'src'},

    package_data={
        'grond': ['report/app/*.html',
                  'report/app/favicon.png',
                  'report/app/css/*.css',
                  'report/app/css/*.map',
                  'report/app/js/*.js',
                  'report/app/js/*.map']},

    data_files=[],
    )
