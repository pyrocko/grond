import logging
import os.path as op
from string import Template
from pyrocko.guts import Object, String, Float, Unicode

try:
    newstr = unicode
except NameError:
    newstr = str


logger = logging.getLogger('grond.meta')


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


class ADict(dict):
    def __getattr__(self, k):
        return self[k]

    def __setattr__(self, k, v):
        self[k] = v


class Parameter(Object):
    name__ = String.T()
    unit = Unicode.T(optional=True)
    scale_factor = Float.T(default=1., optional=True)
    scale_unit = Unicode.T(optional=True)
    label = Unicode.T(optional=True)

    def __init__(self, *args, **kwargs):
        if len(args) >= 1:
            kwargs['name'] = args[0]
        if len(args) >= 2:
            kwargs['unit'] = newstr(args[1])

        self.groups = [None]
        self._name = None

        Object.__init__(self, **kwargs)

    def get_label(self, with_unit=True):
        lbl = [self.label or self.name]
        if with_unit:
            unit = self.get_unit_label()
            if unit:
                lbl.append('[%s]' % unit)

        return ' '.join(lbl)

    def set_groups(self, groups):
        if not isinstance(groups, list):
            raise AttributeError('Groups must be a list of strings.')
        self.groups = groups

    def _get_name(self):
        if None not in self.groups:
            return '%s:%s' % (':'.join(self.groups), self._name)
        return self._name

    def _set_name(self, value):
        self._name = value

    name = property(_get_name, _set_name)

    @property
    def name_nogroups(self):
        return self._name

    def get_value_label(self, value, format='%(value)g%(unit)s'):
        value = self.scaled(value)
        unit = self.get_unit_suffix()
        return format % dict(value=value, unit=unit)

    def get_unit_label(self):
        if self.scale_unit is not None:
            return self.scale_unit
        elif self.unit:
            return self.unit
        else:
            return None

    def get_unit_suffix(self):
        unit = self.get_unit_label()
        if not unit:
            return ''
        else:
            return ' %s' % unit

    def scaled(self, x):
        if isinstance(x, tuple):
            return tuple(v/self.scale_factor for v in x)
        if isinstance(x, list):
            return list(v/self.scale_factor for v in x)
        else:
            return x/self.scale_factor

    def inv_scaled(self, x):
        if isinstance(x, tuple):
            return tuple(v*self.scale_factor for v in x)
        if isinstance(x, list):
            return list(v*self.scale_factor for v in x)
        else:
            return x*self.scale_factor


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
        elif isinstance(path, str):
            return extra(
                op.normpath(xjoin(self._basepath, xjoin(path_prefix, path))))
        else:
            return [
                extra(
                    op.normpath(xjoin(self._basepath, xjoin(path_prefix, p))))
                for p in path]
