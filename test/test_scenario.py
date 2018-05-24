from grond import scenario as grond_scenario
import tempfile

km = 1e3

STORE_STATIC = 'crust2_ib_static'
STORE_WAVEFORMS = 'crust2_ib'


def test_scenario():

    observations = [
        grond_scenario.WaveformObservation(
            nstations=20,
            store_id=STORE_WAVEFORMS),
        grond_scenario.InSARObservation(
            store_id=STORE_STATIC),
        grond_scenario.GNSSCampaignObservation(
            nstations=20,
            store_id=STORE_STATIC)
        ]

    problems = [
        grond_scenario.DCSourceProblem(
            nevents=1),
        grond_scenario.RectangularSourceProblem()
    ]

    for prob in problems:
        with tempfile.TemporaryDirectory(prefix='grond') as project_dir:

            scenario = grond_scenario.GrondScenario(
                project_dir,
                center_lat=41., center_lon=33.3,
                radius=200.*km)
            for obs in observations:
                scenario.add_observation(obs)

            scenario.set_problem(prob)

            scenario.build(force=True, interactive=False)
