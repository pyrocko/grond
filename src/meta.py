import logging
import math
import re
import numpy as num
import os.path as op
from string import Template
from pyrocko.guts import Object, String, Float, Unicode, StringPattern, Bool
from pyrocko import util

guts_prefix = 'grond'


try:
    newstr = unicode
except NameError:
    newstr = str


logger = logging.getLogger('grond.meta')
km = 1e3

classes_with_have_get_plot_classes = []


def has_get_plot_classes(cls):
    classes_with_have_get_plot_classes.append(cls)
    return cls


class StringID(StringPattern):
    pattern = r'^[A-Za-z][A-Za-z0-9._-]{0,64}$'


StringID.regex = re.compile(StringID.pattern)


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


def ordersort(x):
    isort = num.argsort(x)
    iorder = num.empty(isort.size, dtype=num.int)
    iorder[isort] = num.arange(isort.size)
    return iorder


def nextpow2(i):
    return 2**int(math.ceil(math.log(i) / math.log(2.)))


def gather(seq, key, sort=None, filter=None):
    d = {}
    for x in seq:
        if filter is not None and not filter(x):
            continue

        k = key(x)
        if k not in d:
            d[k] = []

        d[k].append(x)

    if sort is not None:
        for v in d.values():
            v.sort(key=sort)

    return d


def str_dist(dist):
    if dist < 10.0:
        return '%g m' % dist
    elif 10. <= dist < 1. * km:
        return '%.0f m' % dist
    elif 1. * km <= dist < 10. * km:
        return '%.1f km' % (dist / km)
    else:
        return '%.0f km' % (dist / km)


def str_duration(t):
    s = ''
    if t < 0.:
        s = '-'

    t = abs(t)

    if t < 10.0:
        return s + '%.2g s' % t
    elif 10.0 <= t < 3600.:
        return s + util.time_to_str(t, format='%M:%S min')
    elif 3600. <= t < 24 * 3600.:
        return s + util.time_to_str(t, format='%H:%M h')
    else:
        return s + '%.1f d' % (t / (24. * 3600.))


try:
    nanmedian = num.nanmedian
except AttributeError:
    def nanmedian(a, axis=None):
        if axis is None:
            return num.median(a[num.isfinite(a)])
        else:
            shape_out = list(a.shape)
            shape_out.pop(axis)
            out = num.empty(shape_out, dtype=a.dtype)
            out[...] = num.nan
            for iout in num.ndindex(tuple(shape_out)):
                iin = list(iout)
                iin[axis:axis] = [slice(0, a.shape[axis])]
                b = a[tuple(iin)]
                out[iout] = num.median(b[num.isfinite(b)])

            return out


class Forbidden(Exception):
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
    optional = Bool.T(default=True, optional=True)

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
            return '%s.%s' % ('.'.join(self.groups), self._name)
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


class Path(String):
    pass


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

    def rel_path(self, path):
        return xrelpath(path, self.get_basepath())


def nslc_to_pattern(s):
    toks = s.split('.')
    if len(toks) == 1:
        return '*.%s.*.*' % s
    elif len(toks) == 2:
        return '%s.*.*' % s
    elif len(toks) == 3:
        return '%s.*' % s
    elif len(toks) == 4:
        return s
    else:
        raise GrondError('Invalid net.sta.loc.cha pattern: %s' % s)


def nslcs_to_patterns(seq):
    return [nslc_to_pattern(s) for s in seq]


class SelectionError(GrondError):
    pass


# --select="magnitude_min:5 tag_contains:a,b "

g_conditions = {}

selected_operators_1_1 = {
    'min': lambda data, key, value: data[key] >= value,
    'max': lambda data, key, value: data[key] <= value,
    'is': lambda data, key, value: data[key] == value}

selected_operators_1_n = {
    'in': lambda data, key, values:
        data[key] in values}

selected_operators_n_1 = {}

selected_operators_n_n = {
    'contains': lambda data, key, values:
        any(d in values for d in data[key])}


selected_operators = set()

for s in (selected_operators_1_1, selected_operators_1_n,
          selected_operators_n_1, selected_operators_n_n):
    for k in s:
        selected_operators.add(k)


def _parse_selected_expression(expression):
    for condition in expression.split():
        condition = condition.strip()
        if condition not in g_conditions:
            try:
                argument, value = condition.split(':', 1)
            except ValueError:
                raise SelectionError(
                    'Invalid condition in selection expression: '
                    '"%s", must be "ARGUMENT:VALUE"' % condition)

            argument = argument.strip()

            try:
                key, operator = argument.rsplit('_', 1)
            except ValueError:
                raise SelectionError(
                    'Invalid argument in selection expression: '
                    '"%s", must be "KEY_OPERATOR"' % argument)

            if operator not in selected_operators:
                raise SelectionError(
                    'Invalid operator in selection expression: '
                    '"%s", available: %s' % (
                        operator, ', '.join('"%s"' % s for s in sorted(
                            list(selected_operators)))))

            g_conditions[condition] = key, operator, value

        yield g_conditions[condition]


def selected(expression, data, types):
    results = []
    for (key, operator, value) in _parse_selected_expression(expression):
        if key not in data:
            raise SelectionError(
                'Invalid key in selection expression: '
                '"%s", available:\n  %s' % (
                    key, '\n  '.join(sorted(data.keys()))))

        typ = types[key]
        if not isinstance(typ, tuple):
            if operator in selected_operators_1_1:
                results.append(
                    selected_operators_1_1[operator](data, key, typ(value)))
            elif operator in selected_operators_1_n:
                values = list(typ(v) for v in value.split(','))
                results.append(
                    selected_operators_1_n[operator](data, key, values))
            else:
                raise SelectionError(
                    'Cannot use operator "%s" with argument "%s".' % (
                        operator,
                        key))

        else:
            if operator in selected_operators_n_1:
                results.append(
                    selected_operators_n_1[operator](data, key, typ[1](value)))
            elif operator in selected_operators_n_n:
                values = typ[0](typ[1](v) for v in value.split(','))
                results.append(
                    selected_operators_n_n[operator](data, key, values))
            else:
                raise SelectionError(
                    'Cannot use operator "%s" with argument "%s".' % (
                        operator,
                        key))

    return all(results)


__all__ = '''
    Forbidden
    GrondError
    Path
    HasPaths
    Parameter
    StringID
'''.split()
