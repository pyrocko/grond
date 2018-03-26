
class PhaseRatioTargetGroup(TargetGroup):
    distance_min = Float.T(optional=True)
    distance_max = Float.T(optional=True)
    distance_3d_min = Float.T(optional=True)
    distance_3d_max = Float.T(optional=True)
    depth_min = Float.T(optional=True)
    depth_max = Float.T(optional=True)
    limit = Int.T(optional=True)
    measure_a = FeatureMeasure.T()
    measure_b = FeatureMeasure.T()

    def get_targets(self, ds, event, default_path):
        logger.debug('Selecting waveform targets...')
        origin = event
        targets = []

        for st in ds.get_stations():
            for cha in self.channels:
                if ds.is_blacklisted((st.nsl() + (cha,))):
                    continue

                target = WaveformMisfitTarget(
                    quantity='displacement',
                    codes=st.nsl() + (cha,),
                    lat=st.lat,
                    lon=st.lon,
                    depth=st.depth,
                    interpolation=self.interpolation,
                    store_id=self.store_id,
                    misfit_config=self.misfit_config,
                    manual_weight=self.weight,
                    normalisation_family=self.normalisation_family,
                    path=self.path or default_path)

                if self.distance_min is not None and \
                   target.distance_to(origin) < self.distance_min:
                    continue

                if self.distance_max is not None and \
                   target.distance_to(origin) > self.distance_max:
                    continue

                if self.distance_3d_min is not None and \
                   target.distance_3d_to(origin) < self.distance_3d_min:
                    continue

                if self.distance_3d_max is not None and \
                   target.distance_3d_to(origin) > self.distance_3d_max:
                    continue

                if self.depth_min is not None and \
                   target.depth < self.depth_min:
                    continue

                if self.depth_max is not None and \
                   target.depth > self.depth_max:
                    continue

                azi, _ = target.azibazi_to(origin)
                if cha == 'R':
                    target.azimuth = azi - 180.
                    target.dip = 0.
                elif cha == 'T':
                    target.azimuth = azi - 90.
                    target.dip = 0.
                elif cha == 'Z':
                    target.azimuth = 0.
                    target.dip = -90.

                target.set_dataset(ds)
                targets.append(target)

        if self.limit:
            return weed(origin, targets, self.limit)[0]
        else:
            return targets



class PhaseRatioResult(MisfitResult):
    pass


class PhaseRatioTarget(MisfitTarget):

    def prepare_modelling(self):
        pass
        # provide gf.Targets for modelling

    def finalize_modelling(self, results):
        # calc ratio
        result = PhaseRatioResult()
        return result


__all__ = '''
    PhaseRatioTargetGroup
    PhaseRatioTarget
    PhaseRatioResult
'''.split()
