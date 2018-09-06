from pyrocko.guts import Object, String, Dict

from grond.version import __version__

guts_prefix = 'grond'


class VersionInfo(Object):
    grond_version = String.T(yamlstyle="'")
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
        grond_version=__version__,
        dependencies=deps)

    return vi


__all__ = ['VersionInfo']
