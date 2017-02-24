import os.path as op
from string import Template
from pyrocko.guts import Object, String


def xjoin(basepath, path):
    if path is None and basepath is not None:
        return basepath
    elif op.isabs(path) or basepath is None:
        return path
    else:
        return op.join(basepath, path)


def xrelpath(path, start):
    if op.isabs(path):
        return path
    else:
        return op.relpath(path, start)


class Forbidden(Exception):
    pass


class Path(String):
    pass


class GrondError(Exception):
    pass


def expand_template(template, d):
    try:
        return Template(template).substitute(d)
    except KeyError as e:
        raise GrondError(
            'invalid placeholder "%s" in template: "%s"' % (str(e), template))
    except ValueError:
        raise GrondError(
            'malformed placeholder in template: "%s"' % template)


class HasPaths(Object):
    path_prefix = Path.T(optional=True)

    def __init__(self, *args, **kwargs):
        Object.__init__(self, *args, **kwargs)
        self._basepath = None
        self._parent_path_prefix = None

    def set_basepath(self, basepath, parent_path_prefix=None):
        self._basepath = basepath
        self._parent_path_prefix = parent_path_prefix
        for (prop, val) in self.T.ipropvals(self):
            if isinstance(val, HasPaths):
                val.set_basepath(
                    basepath, self.path_prefix or self._parent_path_prefix)

    def get_basepath(self):
        assert self._basepath is not None
        return self._basepath

    def change_basepath(self, new_basepath, parent_path_prefix=None):
        assert self._basepath is not None

        self._parent_path_prefix = parent_path_prefix
        if self.path_prefix or not self._parent_path_prefix:

            self.path_prefix = op.normpath(xjoin(xrelpath(
                self._basepath, new_basepath), self.path_prefix))

        for val in self.T.ivals(self):
            if isinstance(val, HasPaths):
                val.change_basepath(
                    new_basepath, self.path_prefix or self._parent_path_prefix)

        self._basepath = new_basepath

    def expand_path(self, path, extra=None):
        assert self._basepath is not None

        if extra is None:
            def extra(path):
                return path

        path_prefix = self.path_prefix or self._parent_path_prefix

        if path is None:
            return None
        elif isinstance(path, basestring):
            return extra(
                op.normpath(xjoin(self._basepath, xjoin(path_prefix, path))))
        else:
            return [
                extra(
                    op.normpath(xjoin(self._basepath, xjoin(path_prefix, p))))
                for p in path]
