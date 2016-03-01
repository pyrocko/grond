import unittest
from pyrocko import util
from grond import HasPaths, Path


class PathTestCase(unittest.TestCase):

    def test_pathstuff(self):

        class B(HasPaths):
            p1 = Path.T()
            p2 = Path.T()

        class A(HasPaths):
            p1 = Path.T()
            p2 = Path.T()
            b = B.T()

        for path_prefix_a in (
                None, 'relative_prefix_a', '/absolute_prefix_a'):
            for path_prefix_b in (
                    None, 'relative_prefix_b', '/absolute_prefix_b'):
                a = A(
                    path_prefix=path_prefix_a,
                    p1='abc/x.txt',
                    p2='/absolute/x.txt',
                    b=B(
                        path_prefix=path_prefix_b,
                        p1='abc/y.txt',
                        p2='/absolute/y.txt'))

                a.set_basepath('rundir')

                t1 = a.expand_path(a.p1)
                t2 = a.expand_path(a.p2)
                t3 = a.b.expand_path(a.b.p1)
                t4 = a.b.expand_path(a.b.p2)

                a.change_basepath('resultdir')

                assert t1 == a.expand_path(a.p1)
                assert t2 == a.expand_path(a.p2)
                assert t3 == a.b.expand_path(a.b.p1)
                assert t4 == a.b.expand_path(a.b.p2)


if __name__ == '__main__':
    util.setup_logging('test_path', 'warning')
    unittest.main()
