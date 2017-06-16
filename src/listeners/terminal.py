from pyrocko import util
import progressbar as pbar
from .base import Listener


class color:
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'


class TerminalListener(Listener):

    col_width = 15
    row_name = color.BOLD + '{:<{col_param_width}s}' + color.END
    parameter_fmt = '{:{col_width}s}'

    def __init__(self):
        self.current_state = None
        self.pbars = {}

    def progress_start(self, name, niter):
        pbar = util.progressbar('analysing problem', niter)
        self.pbars[name] = pbar
        pbar.start()

    def progress_update(self, name, iiter):
        self.pbars[name].update(iiter)

    def progress_finish(self, name):
        self.pbars[name].finish()
        del self.pbars[name]

    def state(self, state):
        lines = []
        self.current_state = state

        def l(t):
            lines.append(t)

        out_ln = self.row_name +\
            ''.join([self.parameter_fmt] * len(state.parameter_sets))
        col_param_width = max([len(p) for p in state.parameter_names]) + 2

        l('Problem name: {s.problem_name}'
          '\t({s.runtime:s} - remaining {s.runtime_remaining})'
            .format(s=state))
        l('Iteration {s.iiter} / {s.niter}\t\t({s.iter_per_second:.1f} iter/s)'
          .format(s=state))

        l(out_ln.format(
            *['Parameter'] + state.parameter_sets.keys(),
            col_param_width=col_param_width,
            col_width=self.col_width,
            type='s'))

        def fmt(s):
            return util.gform(s, significant_digits=(self.col_width-1-6)/2)

        for ip, parameter_name in enumerate(state.parameter_names):
            l(out_ln.format(
                parameter_name,
                *[fmt(v[ip]) for v in state.parameter_sets.values()],
                col_param_width=col_param_width,
                col_width=self.col_width))

        l(state.extra_text.format(
            col_param_width=col_param_width,
            col_width=self.col_width,))

        lines[0:0] = ['\033[2J']
        l('')
        print '\n'.join(lines)
