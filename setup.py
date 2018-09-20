#!/usr/bin/env python

from setuptools import setup
from setuptools.command.install import install


class CustomInstallCommand(install):
    def run(self):
        install.run(self)


setup(
    cmdclass={
        'install': CustomInstallCommand,
    },

    name='grond',

    description='A probabilistic earthquake source inversion framework. '
                'Designed and crafted in Mordor.',

    version='1.0',

    author='The Grond Developers',

    author_email='info@pyrocko.org',

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
        'grond.problems.multirectangular',
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
        'grond': [
            'report/app/*.html',
            'report/app/favicon.png',
            'report/app/templates/*.html',
            'report/app/css/*.css',
            'report/app/js/*.js']},

    data_files=[],

    license='GPLv3',

    classifiers=[
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: Implementation :: CPython',
        'Operating System :: POSIX',
        'Operating System :: MacOS',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Visualization',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Software Development :: Libraries :: Application Frameworks',
        ],

    keywords=[
        'seismology, waveform analysis, earthquake modelling, geophysics,'
        ' geophysical inversion'],
    )
