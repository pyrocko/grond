import curses


class State(object):
    iiter = 0
    niter = 0
    iter_sec = 0.
    problem_name = ''
    parameter_names = []
    column_names = []
    values = []
    text = ''


class _CursesPad(object):
    def __init__(self, pad):
        self.pad = pad
        self.rows, self.cols = self.pad.getyx()

    def resize_pad(self):
        return
        self.pad.resize(self.rows, self.cols)


class CursesListener(object):

    class ParameterTable(_CursesPad):

        value_fmt = '{0:8.4g}'
        column_padding = 2

        def update(self, state):
            pad = self.pad
            pad.clear()
            if not state:
                return

            parameter_names = ['Parameters'] + state.parameters
            col = 0

            for icol in range(len(state.values) + 1):
                row = 0

                if icol == 0:
                    col_width = max([len(p) for p in parameter_names])
                    for name in parameter_names:
                        pad.addstr(
                            row, col,
                            '{:<{width}}'.format(
                                name, width=col_width),
                            curses.A_BOLD)
                        row += 1

                else:
                    igroup = icol - 1
                    col_heading = state.column_names[igroup]
                    col_width = max(
                        len(col_heading),
                        len(self.value_fmt.format(0.)) + self.column_padding)

                    pad.addstr(row, col,
                               '{:>{width}}'.format(
                                    col_heading, width=col_width),
                               curses.A_BOLD)

                    for iv, v in enumerate(state.values[igroup]):
                        row += 1
                        vstr = ' ' * self.column_padding +\
                            self.value_fmt.format(v)
                        pad.addstr(row, col, vstr)

                col += col_width
            self.rows = row
            self.resize_pad()
            pad.noutrefresh()

    class Footer(_CursesPad):

        def update(self, state):
            pad = self.pad
            pad.clear()
            if not state:
                return

            pad.addstr(0, 0, 'Performance:')
            pad.addstr(0, 14, '%.1f iter/s' % state.iter_sec)
            pad.addstr(1, 0, state.text)
            self.rows = 3
            self.resize_pad()
            pad.noutrefresh()

    class Header(_CursesPad):

        def update(self, state):
            pad = self.pad
            pad.clear()
            if not state:
                return

            pad.addstr(0, 0, 'Problem Name:')
            pad.addstr(0, 14, state.problem_name,
                       curses.A_BOLD)
            pad.addstr(1, 0, 'Iteration:')
            pad.addstr(1, 14, '%d / %d' % (state.iiter, state.niter),
                       curses.A_BOLD)
            self.rows = 3
            self.resize_pad()
            pad.noutrefresh()

    def __init__(self):
        self.scr = None
        self.state = None
        curses.wrapper(self.set_screen)

        self.header_pad = self.Header(self.scr.subpad(3, 100, 0, 0))
        self.parameter_pad = self.ParameterTable(self.scr.subpad(3, 0))
        self.footer_pad = self.Footer(self.scr.subpad(3, 100, 5, 0))

    def set_screen(self, scr):
        self.scr = scr

    def state(self, state):
        self.state = state
        self.parameter_pad.update(state)
        self.header_pad.update(state)
        self.footer_pad.update(state)

        self.footer_pad.pad.mvwin(
            self.parameter_pad.rows + self.parameter_pad.pad.getparyx()[0] + 2,
            0)
        self.scr.refresh()
