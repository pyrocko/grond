import time
import tempfile
from grond.problems.base import RandomStateManager


def test_get_state():
    rstate = RandomStateManager()
    rs = rstate.get_rstate('test')
    rs1 = rstate.get_rstate('test', seed=123123123)

    assert rs is rs1
    assert rstate.nstates == 1

    rs2 = rstate.get_rstate('test2')
    assert rstate.nstates == 2


def test_save_load_state():
    rstate = RandomStateManager()
    rstate.get_rstate('test')
    rs = rstate.get_rstate('test2', seed=4414412)
    rs.uniform()

    with tempfile.NamedTemporaryFile() as f:
        rstate.save_state(f.name)
        unif = rs.uniform()

        rstate_load = RandomStateManager()
        rstate_load.load_state(f.name)

        rs_loaded = rstate_load.get_rstate('test2')
        assert unif == rs_loaded.uniform()
        assert rstate_load.nstates == 2


def test_performance():
    rstate_perf = RandomStateManager()
    for r in range(1000):
        rstate_perf.get_rstate('test-%d' % r)

    with tempfile.NamedTemporaryFile() as f:
        t0 = time.time()
        rstate_perf.save_state(f.name)
        print('saved %d states in %f s'
              % (rstate_perf.nstates, time.time()-t0))

        rstate_load = RandomStateManager()
        t0 = time.time()
        rstate_load.load_state(f.name)
        print('loaded %d states in %f s'
              % (rstate_load.nstates, time.time()-t0))


if __name__ == '__main__':
    test_get_state()
    test_save_load_state()
    test_performance()
