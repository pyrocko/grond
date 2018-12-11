from pyrocko.guts import Object, String, Dict, Bool

from grond.setup_info import git_sha1, local_modifications, \
    version, long_version, installed_date

guts_prefix = 'grond'


class VersionInfo(Object):
    grond_version = String.T(yamlstyle="'")
    grond_long_version = String.T(yamlstyle="'")
    git_sha1 = String.T(optional=True, yamlstyle="'")
    local_modifications = Bool.T(optional=True)
    installed_date = String.T(optional=True, yamlstyle="'")
    dependencies = Dict.T(String.T(), String.T(yamlstyle="'"))


def version_info():

    deps = {}

    try:
        import pyrocko
        deps['pyrocko'] = pyrocko.long_version
    except ImportError:
        pass

    try:
        import numpy
        deps['numpy'] = numpy.__version__
    except ImportError:
        pass

    try:
        import scipy
        deps['scipy'] = scipy.__version__
    except ImportError:
        pass

    try:
        import matplotlib
        deps['matplotlib'] = matplotlib.__version__
    except ImportError:
        pass

    try:
        from pyrocko.gui.qt_compat import Qt
        deps['PyQt'] = Qt.PYQT_VERSION_STR
        deps['Qt'] = Qt.QT_VERSION_STR
    except ImportError:
        pass

    import sys
    deps['python'] = '%s.%s.%s' % sys.version_info[:3]

    vi = VersionInfo(
        grond_version=version,
        grond_long_version=long_version,
        git_sha1=git_sha1,
        local_modifications=local_modifications,
        installed_date=installed_date,
        dependencies=deps)

    return vi


__all__ = ['VersionInfo']
