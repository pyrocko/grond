import logging
import numpy as num
from matplotlib import cm, gridspec

from grond.plot.config import PlotConfig
from grond.plot.collection import PlotItem

from matplotlib import pyplot as plt
from matplotlib import patches
from pyrocko.guts import Tuple, Float, String, Int, Bool
import utm

logger = logging.getLogger('grond.targets.satellite.plot')

km = 1e3
d2r = num.pi/180.
guts_prefix = 'grond'


def scale_axes(axis, scale, offset=0.):
    from matplotlib.ticker import ScalarFormatter

    class FormatScaled(ScalarFormatter):

        @staticmethod
        def __call__(value, pos):
            return '{:,.1f}'.format((offset + value) * scale).replace(',', ' ')

    axis.set_major_formatter(FormatScaled())


class SatelliteTargetDisplacement(PlotConfig):
    ''' Maps showing surface displacements from satellite and modelled data '''

    name = 'satellite'
    dpi = Int.T(
        default=250)
    size_cm = Tuple.T(
        2, Float.T(),
        default=(22., 12.))
    colormap = String.T(
        default='RdBu',
        help='Colormap for the surface displacements')
    relative_coordinates = Bool.T(
        default=False,
        help='Show relative coordinates, initial location centered at 0N, 0E')

    def make(self, environ):
        cm = environ.get_plot_collection_manager()
        history = environ.get_history()
        optimiser = environ.get_optimiser()
        ds = environ.get_dataset()

        environ.setup_modelling()

        cm.create_group_mpl(
            self,
            self.draw_static_fits(ds, history, optimiser),
            title=u'InSAR Displacements',
            section='fits',
            feather_icon='navigation',
            description=u' Maps showing subsampled surface displacements as'
                        u' observed, modelled and the residual (observed minus'
                        u' modelled).\n The displacement values predicted by'
                        u' the orbit-ambiguity ramps are added to the modelled'
                        u' displacements (middle panels). The color shows the'
                        u' LOS displacement values associated with, and the'
                        u' extent of, every quadtree box. The light grey dots'
                        u' show the focal point of pixels combined in the'
                        u' quadtree box. This point corresponds to the'
                        u' position of the modelled data point.\n The large dark'
                        u' grey dot shows the reference source position. The'
                        u' grey filled box shows the surface projection of the'
                        u' modelled source, with the thick-lined edge marking'
                        u' the upper fault edge. '
                        u' Complete data extent is shown.')

    def draw_static_fits(self, ds, history, optimiser, closeup=False):
        from pyrocko.orthodrome import latlon_to_ne_numpy
        problem = history.problem

        sat_targets = problem.satellite_targets
        for target in sat_targets:
            target.set_dataset(ds)

        gms = problem.combine_misfits(
            history.misfits,
            extra_correlated_weights=optimiser.get_correlated_weights(problem))
        isort = num.argsort(gms)
        gms = gms[isort]
        models = history.models[isort, :]
        xbest = models[0, :]
        # nsources = problem.nsources #help
        nsources = 2
        if nsources is not None:
            sources = []
            for i in range(nsources):
                source_i = problem.get_source(xbest, i)
                sources.append(source_i)
            source = sources[0]
        else:
            source = problem.get_source(xbest)
        results = problem.evaluate(xbest, targets=sat_targets)

        def initAxes(ax, scene, title, last_axes=False):
            ax.set_title(title)
            ax.tick_params(length=2)

            if scene.frame.isMeter():
                ax.set_xlabel('Easting [km]')
                scale_x = {'scale': 1./km}
                scale_y = {'scale': 1./km}
                if not self.relative_coordinates:
                    utm_E, utm_N, utm_zone, utm_zone_letter =\
                        utm.from_latlon(source.lat,
                                        source.lon)
                    scale_x['offset'] = utm_E
                    scale_y['offset'] = utm_N

                    if last_axes:
                        ax.text(0.975, 0.025,
                                'UTM Zone %d%s' % (utm_zone, utm_zone_letter),
                                va='bottom', ha='right',
                                fontsize=8, alpha=.7,
                                transform=ax.transAxes)
                ax.set_aspect('equal')

            elif scene.frame.isDegree():
                ax.set_xlabel('Lon [°]')
                scale_x = {'scale': 1.}
                scale_y = {'scale': 1.}
                if not self.relative_coordinates:
                    scale_x['offset'] = source.lat
                    scale_y['offset'] = source.lon
                ax.set_aspect(1./num.cos(source.lat*d2r))

            scale_axes(ax.get_xaxis(), **scale_x)
            scale_axes(ax.get_yaxis(), **scale_y)

        def drawSource(ax, scene):
            if scene.frame.isMeter():
                if nsources is not None:
                    fns = []
                    fes = []
                    for subsource in sources:
                        fn_sub, fe_sub = subsource.outline(cs='xy').T
                        fes.append(fe_sub)
                        fns.append(fn_sub)
                else:
                    fn, fe = source.outline(cs='xy').T
                    fn -= fn.mean()
                    fe -= fe.mean()

            elif scene.frame.isDegree():
                if nsources is not None:
                    fns = []
                    fes = []
                    for subsource in sources:
                        fn_sub, fe_sub = subsource.outline(cs='latlon').T
                        fes.append(fe_sub)
                        fns.append(fn_sub)

                else:
                    fn, fe = source.outline(cs='latlon').T
                    fn -= source.effective_lat
                    fe -= source.effective_lon

            # source is centered
            ax.scatter(0., 0., color='black', s=3, alpha=.5, marker='o')# shows references
            if nsources is not None:
                for fe, fn in zip(fes, fns):
                    ax.fill(fe, fn,
                            edgecolor=(0., 0., 0.),
                            facecolor=(.5, .5, .5), alpha=0.5)
                    ax.plot(fe[0:2], fn[0:2], 'k', linewidth=1.3)

            else:

                ax.fill(fe, fn,
                        edgecolor=(0., 0., 0.),
                        facecolor=(.1, .2, .4), alpha=0.5)
                ax.plot(fe[0:2], fn[0:2], 'k', linewidth=1.3)

        def mapDisplacementGrid(displacements, scene):
            arr = num.full_like(scene.displacement, fill_value=num.nan)
            qt = scene.quadtree

            for syn_v, l in zip(displacements, qt.leaves):
                arr[l._slice_rows, l._slice_cols] = syn_v

            arr[scene.displacement_mask] = num.nan
            return arr

        def drawLeaves(ax, scene, offset_e=0., offset_n=0.):
            rects = scene.quadtree.getMPLRectangles()
            for r in rects:
                r.set_edgecolor((.4, .4, .4))
                r.set_linewidth(.5)
                r.set_facecolor('none')
                r.set_x(r.get_x() - offset_e)
                r.set_y(r.get_y() - offset_n)
            map(ax.add_artist, rects)

            ax.scatter(scene.quadtree.leaf_coordinates[:, 0] - offset_e,
                       scene.quadtree.leaf_coordinates[:, 1] - offset_n,
                       s=.25, c='black', alpha=.1)

        def addArrow(ax, scene):
            phi = num.nanmean(scene.phi)
            los_dx = num.cos(phi + num.pi) * .0625
            los_dy = num.sin(phi + num.pi) * .0625

            az_dx = num.cos(phi - num.pi/2) * .125
            az_dy = num.sin(phi - num.pi/2) * .125

            anchor_x = .9 if los_dx < 0 else .1
            anchor_y = .85 if los_dx < 0 else .975

            az_arrow = patches.FancyArrow(
                x=anchor_x-az_dx, y=anchor_y-az_dy,
                dx=az_dx, dy=az_dy,
                head_width=.025,
                alpha=.5, fc='k',
                head_starts_at_zero=False,
                length_includes_head=True,
                transform=ax.transAxes)

            los_arrow = patches.FancyArrow(
                x=anchor_x-az_dx/2, y=anchor_y-az_dy/2,
                dx=los_dx, dy=los_dy,
                head_width=.02,
                alpha=.5, fc='k',
                head_starts_at_zero=False,
                length_includes_head=True,
                transform=ax.transAxes)

            ax.add_artist(az_arrow)
            ax.add_artist(los_arrow)

        urE, urN, llE, llN = (0., 0., 0., 0.)
        for target in sat_targets:

            if target.scene.frame.isMeter():
                off_n, off_e = map(float, latlon_to_ne_numpy(
                    target.scene.frame.llLat, target.scene.frame.llLon,
                    source.lat, source.lon))
            if target.scene.frame.isDegree():
                off_n = source.lat - target.scene.frame.llLat
                off_e = source.lon - target.scene.frame.llLon

            turE, turN, tllE, tllN = zip(
                *[(l.gridE.max()-off_e,
                   l.gridN.max()-off_n,
                   l.gridE.min()-off_e,
                   l.gridN.min()-off_n)
                  for l in target.scene.quadtree.leaves])

            turE, turN = map(max, (turE, turN))
            tllE, tllN = map(min, (tllE, tllN))
            urE, urN = map(max, ((turE, urE), (urN, turN)))
            llE, llN = map(min, ((tllE, llE), (llN, tllN)))

        def generate_plot(sat_target, result, ifig):

            scene = sat_target.scene

            fig = plt.figure()
            fig.set_size_inches(*self.size_inch)
            gs = gridspec.GridSpec(
                2, 3,
                wspace=.05, hspace=.2,
                left=.1, right=.975, top=.95,
                height_ratios=[12, 1])

            item = PlotItem(
                name='fig_%i' % ifig,
                attributes={'targets': [sat_target.path]},
                title=u'Satellite Surface Displacements - %s'
                      % scene.meta.scene_title,
                description=u'''
Surface displacements derived from satellite data.
(Left) the input data, (center) the modelled
data and (right) the model residual.
'''.format(meta=scene.meta))

            stat_obs = result.statics_obs
            stat_syn = result.statics_syn['displacement.los']
            res = stat_obs - stat_syn

            if scene.frame.isMeter():
                offset_n, offset_e = map(float, latlon_to_ne_numpy(
                    scene.frame.llLat, scene.frame.llLon,
                    source.lat, source.lon))
            elif scene.frame.isDegree():
                offset_n = source.lat - scene.frame.llLat
                offset_e = source.lon - scene.frame.llLon

            im_extent = (scene.frame.E.min() - offset_e,
                         scene.frame.E.max() - offset_e,
                         scene.frame.N.min() - offset_n,
                         scene.frame.N.max() - offset_n)

            abs_displ = num.abs([stat_obs.min(), stat_obs.max(),
                                 stat_syn.min(), stat_syn.max(),
                                 res.min(), res.max()]).max()

            cmw = cm.ScalarMappable(cmap=self.colormap)
            cmw.set_clim(vmin=-abs_displ, vmax=abs_displ)
            cmw.set_array(stat_obs)

            axes = [fig.add_subplot(gs[0, 0]),
                    fig.add_subplot(gs[0, 1]),
                    fig.add_subplot(gs[0, 2])]

            ax = axes[0]
            ax.imshow(mapDisplacementGrid(stat_obs, scene),
                      extent=im_extent, cmap=self.colormap,
                      vmin=-abs_displ, vmax=abs_displ,
                      origin='lower')
            drawLeaves(ax, scene, offset_e, offset_n)
            drawSource(ax, scene)
            addArrow(ax, scene)
            initAxes(ax, scene, 'Observed')

            ax.text(.025, .025, 'Scene ID: %s' % scene.meta.scene_id,
                    fontsize=8, alpha=.7,
                    va='bottom', transform=ax.transAxes)
            if scene.frame.isDegree():
                ax.set_ylabel('Lat [°]')
            elif scene.frame.isMeter():
                ax.set_ylabel('Northing [km]')

            ax = axes[1]
            ax.imshow(mapDisplacementGrid(stat_syn, scene),
                      extent=im_extent, cmap=self.colormap,
                      vmin=-abs_displ, vmax=abs_displ,
                      origin='lower')
            drawLeaves(ax, scene, offset_e, offset_n)
            drawSource(ax, scene)
            addArrow(ax, scene)
            initAxes(ax, scene, 'Model')
            ax.get_yaxis().set_visible(False)

            ax = axes[2]
            ax.imshow(mapDisplacementGrid(res, scene),
                      extent=im_extent, cmap=self.colormap,
                      vmin=-abs_displ, vmax=abs_displ,
                      origin='lower')
            drawLeaves(ax, scene, offset_e, offset_n)
            drawSource(ax, scene)
            addArrow(ax, scene)
            initAxes(ax, scene, 'Residual', last_axes=True)
            ax.get_yaxis().set_visible(False)

            for ax in axes:
                ax.set_xlim(llE, urE)
                ax.set_ylim(llN, urN)

            if closeup:
                if scene.frame.isMeter():
                    if nsources is not None:
                        fns = []
                        fes = []
                        for subsource in sources:
                            fn_sub, fe_sub = subsource.outline(cs='xy').T
                            fes.append(fe_sub)
                            fns.append(fn_sub)
                        fnv = list(fns[0])
                        fnv.extend(list(fns[1]))
                        fnv = num.array(fnv)
                        fev = list(fes[0])
                        fev.extend(list(fes[1]))
                        fev = num.array(fev)
                        off_n = (fns[0][0] + fns[1][1]) / 2
                        off_e = (fes[0][0] + fes[1][1]) / 2
                    else:
                        fn, fe = source.outline(cs='xy').T
                        if fn.size > 1:
                            off_n = (fn[0] + fn[1]) / 2
                            off_e = (fe[0] + fe[1]) / 2
                        else:
                            off_n = fn[0]
                            off_e = fe[0]
                        fnv = fn
                        fev = fe

                elif scene.frame.isDegree():
                    if nsources is not None:
                        fns = []
                        fes = []
                        for subsource in sources:
                            fn_sub, fe_sub = subsource.outline(cs='latlon').T
                            fn_sub -= subsource.effective_lat
                            fe_sub -= subsource.effective_lon
                            fes.append(fe_sub)
                            fns.append(fn_sub)
                        fnv = list(fns[0])
                        fnv.extend(list(fns[1]))
                        fnv = num.array(fnv)
                        fev = list(fes[0])
                        fev.extend(list(fes[1]))
                        fev = num.array(fev)
                        off_n = (fns[0][0] + fns[1][1]) / 2
                        off_e = (fes[0][0] + fes[1][1]) / 2
                    else:
                        fn, fe = source.outline(cs='latlon').T
                        fn -= source.effective_lat
                        fe -= source.effective_lon
                        if fn.size > 1:
                            off_n = (fn[0] + fn[1]) / 2
                            off_e = (fe[0] + fe[1]) / 2
                        else:
                            off_n = fn[0]
                            off_e = fe[0]
                        fnv = fn
                        fev = fe

                fault_size = 2*num.sqrt(max(abs(fnv-off_n))**2
                                        + max(abs(fev-off_e))**2)
                fault_size *= self.map_scale
                if fault_size == 0.0:
                    if scene.frame.isMeter():
                        fault_size = 1000.0
                    elif scene.frame.isDegree():
                        fault_size = 1.0

                for ax in axes:
                    ax.set_xlim(-fault_size/2 + off_e, fault_size/2 + off_e)
                    ax.set_ylim(-fault_size/2 + off_n, fault_size/2 + off_n)

            cax = fig.add_subplot(gs[1, :])
            cbar = fig.colorbar(cmw, cax=cax, orientation='horizontal',
                                use_gridspec=True)

            cbar.set_label('LOS Displacement [m]')

            return (item, fig)

        for ifig, (sat_target, result) in enumerate(zip(sat_targets, results)):
            yield generate_plot(sat_target, result, ifig)


class SatelliteTargetDisplacementCloseup(SatelliteTargetDisplacement):
    ''' Close-up of satellite surface displacements and modelled data. '''
    name = 'satellite_closeup'

    map_scale = Float.T(
        default=2.,
        help='Scale the map surroundings, larger value zooms out.')

    def make(self, environ):
        cm = environ.get_plot_collection_manager()
        history = environ.get_history()
        optimiser = environ.get_optimiser()
        ds = environ.get_dataset()

        environ.setup_modelling()

        cm.create_group_mpl(
            self,
            self.draw_static_fits(ds, history, optimiser, closeup=True),
            title=u'InSAR Displacements (Closeup)',
            section='fits',
            feather_icon='zoom-in',
            description=u'''
Maps showing subsampled surface displacements as observed, modelled and the
residual (observed minus modelled).

The displacement values predicted by the orbit-ambiguity ramps are added to the
modelled displacements (middle panels). The color shows the LOS displacement
values associated with, and the extent of, every quadtree box. The light grey
dots show the focal point of pixels combined in the quadtree box. This point
corresponds to the position of the modelled data point.

The large dark grey dot shows the reference source position. The grey filled
box shows the surface projection of the modelled source, with the thick-lined
edge marking the upper fault edge. Map is focused around the fault's extent.
''')

class SatelliteCustom(PlotConfig):
    ''' Maps showing surface displacements from satellite and modelled data '''

    name = 'satellite_custom'
    dpi = Int.T(
        default=250)
    size_cm = Tuple.T(
        2, Float.T(),
        default=(22., 10.))
    colormap = String.T(
        default='RdBu',
        help='Colormap for the surface displacements')
    relative_coordinates = Bool.T(
        default=False,
        help='Show relative coordinates, initial location centered at 0N, 0E')
    color_parameter = String.T(default='misfit')
    #dem_path = String.T()

    def make(self, environ):
        cm = environ.get_plot_collection_manager()
        #history = environ.get_history()
        history = environ.get_history(subset='harvest')

        optimiser = environ.get_optimiser()
        ds = environ.get_dataset()

        environ.setup_modelling()

        cm.create_group_mpl(
            self,
            self.draw_static_fits(ds, history, optimiser),
            title=u'Satellite Surface Displacements',
            section='fits',
            feather_icon='navigation',
            description=u'Maps showing surface displacements'
                        u' from satellite and modelled data. '
                        u'The ensemble\'s best fault model is outlined,'
                        u' the gray dot indicates the initial location'
                        u' before the optimisation. '
                        u'Complete data extent is shown.')

    def draw_static_fits(self, ds, history, optimiser, closeup=False):
        from pyrocko.orthodrome import latlon_to_ne_numpy
        problem = history.problem
        color_parameter = self.color_parameter
        sat_targets = problem.satellite_targets
        for target in sat_targets:
            target.set_dataset(ds)

        # models = history.models
        #
        # if optimiser.nbootstrap == 0:
        #     gms = problem.combine_misfits(history.misfits)
        # else:
        #     gms = history.bootstrap_misfits[:, 0]
        gms = problem.combine_misfits(
            history.misfits,
            extra_correlated_weights=optimiser.get_correlated_weights(problem))
        isort = num.argsort(gms)
        gms = gms[isort]
        models = history.models[isort, :]
        xbest = models[0, :]

        #isort = num.argsort(gms)[::-1]
        #gms = gms[isort]
        print("isort", len(isort))
        #models = models[isort, :]

        #xbest = models[0, :]
        nmodels = models.shape[0]
        print("nmodels", nmodels)
        # ifs überflüssig?
        if color_parameter == 'misfit':
            iorder = num.arange(nmodels)
            icolor = iorder

        if color_parameter in problem.parameter_names:
            ind = problem.name_to_index(color_parameter)
            icolor = problem.extract(models, ind)

        #source = problem.get_source(xbest)
        nsources = 2
        sources = []
        for i in range(nsources):
                source_i = problem.get_source(xbest, i)
                sources.append(source_i)
        source = sources[0]
        results = problem.evaluate(xbest, targets=sat_targets)

        # src_cnt = []
        # for sourc in sources:
        #     src_nn, src_ee, src_zz = sourc.outline(cs='xyz').T
        #     src_cnt.append([num.mean(src_ee), num.mean(src_nn), num.mean(src_zz)])

        # centroids relative to reference location
        centroids_xyz = []
        # relative position of the nucleation point (off the centroid)
        nucs_xy_xyz = []
        outlines_e = []
        outlines_n = []
        for mods in models:
            for i in range(nsources):
                srcx = problem.get_source(mods,i)

                nuc_x_xyz = [num.sin(num.deg2rad(srcx.strike)) \
                              * 0.5*srcx.length*srcx.nucleation_x,
                             num.cos(num.deg2rad(srcx.strike)) \
                              * 0.5*srcx.length*srcx.nucleation_x,
                             0]
                nuc_y_xyz = [num.cos(num.deg2rad(srcx.strike)) \
                             * num.cos(num.deg2rad(srcx.dip)) \
                             * 0.5*srcx.width*srcx.nucleation_y,
                             -num.sin(num.deg2rad(srcx.strike)) \
                             * num.cos(num.deg2rad(srcx.dip)) \
                             * 0.5*srcx.width*srcx.nucleation_y,
                             num.sin(num.deg2rad(srcx.dip)) \
                             * 0.5*srcx.width*srcx.nucleation_y]



                xy = num.array([nuc_x_xyz[0] + nuc_y_xyz[0],
                                nuc_x_xyz[1] + nuc_y_xyz[1],
                                nuc_x_xyz[2] + nuc_y_xyz[2]])
                centroid =  num.array([srcx.east_shift \
                             + num.cos(num.deg2rad(srcx.strike)) \
                             * num.cos(num.deg2rad(srcx.dip)) * 0.5*srcx.width,
                             srcx.north_shift - num.sin(num.deg2rad(srcx.strike)) \
                             * num.cos(num.deg2rad(srcx.dip))* 0.5*srcx.width,
                             srcx.depth + num.sin(num.deg2rad(srcx.dip)) \
                             * 0.5*srcx.width])

                fn, fe = srcx.outline(cs='xy').T



                if num.shape(nucs_xy_xyz)[0]==0:
                    # relative position of the nucleation point (off the centroid)
                    nucs_xy_xyz =  xy
                    # centroids relative to reference location
                    centroids_xyz = centroid
                    outlines_e = fe
                    print("fe",fe)
                    outlines_n = fn
                else:
                    nucs_xy_xyz = num.vstack((nucs_xy_xyz, xy))
                    centroids_xyz = num.vstack((centroids_xyz, centroid))
                    outlines_e = num.vstack((outlines_e, fe))
                    outlines_n = num.vstack((outlines_n, fn))

        for ipar in range(problem.ncombined):
            checkpar = problem.combined[ipar]
            if checkpar.name=='nucleation_x1':
                xpar1 = ipar
                nuc_x_pars1 = problem.combined[ipar]
            elif checkpar.name=='nucleation_y1':
                ypar1 = ipar
                nuc_y_pars1 = problem.combined[ipar]
            elif checkpar.name=='depth1':
                zpar1 = ipar
                src_z_pars1 = problem.combined[ipar]
            elif checkpar.name=='east_shift1':
                epar1 = ipar
                src_e_pars1 = problem.combined[ipar]
            elif checkpar.name=='north_shift1':
                npar1 = ipar
                src_n_pars1 = problem.combined[ipar]
            elif checkpar.name=='nucleation_x2':
                xpar2 = ipar
                nuc_x_pars2 = problem.combined[ipar]
            elif checkpar.name=='nucleation_y2':
                ypar2 = ipar
                nuc_y_pars2 = problem.combined[ipar]
            elif checkpar.name=='depth2':
                zpar2 = ipar
                src_z_pars2 = problem.combined[ipar]
            elif checkpar.name=='east_shift2':
                epar2 = ipar
                src_e_pars2 = problem.combined[ipar]
            elif checkpar.name=='north_shift2':
                npar2 = ipar
                src_n_pars2 = problem.combined[ipar]

        fx1 = problem.extract(models, xpar1)
        fy1 = problem.extract(models, ypar1)

        fx2 = problem.extract(models, xpar2)
        fy2 = problem.extract(models, ypar2)

        print("len(nuc_x_pars1.scaled ",len(nuc_x_pars1.scaled(fx1)))
        def initAxes(ax, scene, title, last_axes=False):
            ax.set_title(title)
            ax.tick_params(length=2)

            if scene.frame.isMeter():
                ax.set_xlabel('Easting [km]')
                scale_x = {'scale': 1./km}
                scale_y = {'scale': 1./km}
                if not self.relative_coordinates:
                    import utm
                    utm_E, utm_N, utm_zone, utm_zone_letter =\
                        utm.from_latlon(source.lat, source.lon)
                    scale_x['offset'] = utm_E
                    scale_y['offset'] = utm_N

                    if last_axes:
                        ax.text(0.975, 0.025,
                                'UTM Zone %d%s' % (utm_zone, utm_zone_letter),
                                va='bottom', ha='right',
                                fontsize=8, alpha=.7,
                                transform=ax.transAxes)

            elif scene.frame.isDegree():
                ax.set_xlabel('Lon [°]')
                scale_x = {'scale': 1.}
                scale_y = {'scale': 1.}
                if not self.relative_coordinates:
                    scale_x['offset'] = source.lon
                    scale_y['offset'] = source.lat

            scale_axes(ax.get_xaxis(), **scale_x)
            scale_axes(ax.get_yaxis(), **scale_y)
            ax.set_aspect('equal')

        def drawSource(ax, scene):
            if scene.frame.isMeter():
                fns = []
                fes = []
                for subsource in sources:
                    fn_sub, fe_sub = subsource.outline(cs='xy').T
                    fes.append(fe_sub)
                    fns.append(fn_sub)
                #fn, fe = source.outline(cs='xy').T
            elif scene.frame.isDegree():
                fns = []
                fes = []
                for subsource in sources:
                    fn_sub, fe_sub = subsource.outline(cs='latlon').T
                    fes.append(fe_sub)
                    fns.append(fn_sub)
                # fn, fe = source.outline(cs='latlon').T
                # fn -= source.lat
                # fe -= source.lon
            off_nuc_x = []
            off_nuc_y = []
            for subsource in sources:
                off_nuc_x.append([num.sin(num.deg2rad(subsource.strike)) \
                             *0.5*subsource.length*subsource.nucleation_x,
                             num.cos(num.deg2rad(subsource.strike)) \
                             *0.5*subsource.length*subsource.nucleation_x])
                off_nuc_y.append([num.cos(num.deg2rad(subsource.strike)) \
                             * num.sin(num.deg2rad(subsource.dip)) \
                             *0.5*subsource.width*subsource.nucleation_y,
                             num.cos(num.deg2rad(subsource.strike)) \
                             * num.sin(num.deg2rad(subsource.dip)) \
                             *0.5*subsource.width*subsource.nucleation_y])
            # nucs_x_map = [num.sin(num.deg2rad(source.strike)) \
            #               * 0.5*source.length*nuc_x_pars.scaled(fx),
            #              num.cos(num.deg2rad(source.strike)) \
            #               * 0.5*source.length*nuc_x_pars.scaled(fx)]
            # nucs_y_map = [num.cos(num.deg2rad(source.strike)) \
            #              * num.cos(num.deg2rad(source.dip)) \
            #              * 0.5*source.width*nuc_y_pars.scaled(fy),
            #              num.sin(num.deg2rad(source.strike)) \
            #              * num.cos(num.deg2rad(source.dip)) \
            #              * 0.5*source.width*nuc_y_pars.scaled(fy)]

            #num.shape(centroids_xyz)
            # source is centered
            ax.scatter(0., 0., color='black', s=3, alpha=.5, marker='o')
            ###ax.scatter(nuc_x_pars, nuc_y_pars,  color='black', s=3, alpha=.5, marker='o')
            for fe,fn,off_nuc_xx,off_nuc_yy in zip(fes, fns, off_nuc_x, off_nuc_y):
                ax.fill(fe, fn,
                        edgecolor=(0., 0., 0.),
                        facecolor=(.5, .5, .5), alpha=0.3)
                ax.scatter(fe[:-1].mean()+off_nuc_xx[0]+off_nuc_yy[0],
                           fn[:-1].mean()+off_nuc_xx[1]+off_nuc_yy[1], color = 'g',
                           s=17, marker='o')
            #ax.scatter(centroids_xyz[:,0],
            #           centroids_xyz[:,1], color = 'm',
            #           s=2, marker='*')
            #ax.scatter(nucs_xy_xyz[:,0] + centroids_xyz[:,0] ,
            #           nucs_xy_xyz[:,1] + centroids_xyz[:,1] , color = 'k',
            #           s=0.3, marker='o')
                ax.scatter(fe[:-1].mean(),fn[:-1].mean(), s=5, marker='x')

                ax.plot(fe[0:2], fn[0:2], 'k', linewidth=1.3)

        def drawSource_xz(ax, scene):
            if scene.frame.isMeter():
                #fn, fe = source.outline(cs='xy').T
                fns = []
                fes = []
                for subsource in sources:
                    fn_sub, fe_sub = subsource.outline(cs='xy').T
                    fes.append(fe_sub)
                    fns.append(fn_sub)
            elif scene.frame.isDegree():
                # fn, fe = source.outline(cs='latlon').T
                # fn -= source.lat
                # fe -= source.lon
                fns = []
                fes = []
                for subsource in sources:
                    fn_sub, fe_sub = subsource.outline(cs='latlon').T
                    fes.append(fe_sub)
                    fns.append(fn_sub)

            off_nuc_x = []
            off_nuc_y = []
            for subsource in sources:
                off_nuc_x.append([num.sin(num.deg2rad(subsource.strike))*\
                                0.5*subsource.length*subsource.nucleation_x,
                             num.cos(num.deg2rad(subsource.strike))*\
                                0.5*subsource.length*subsource.nucleation_x])
                off_nuc_y.append([num.cos(num.deg2rad(subsource.strike)) \
                             * num.sin(num.deg2rad(subsource.dip)) \
                             *0.5*subsource.width*subsource.nucleation_y,
                             num.cos(num.deg2rad(subsource.strike)) \
                             * num.sin(num.deg2rad(subsource.dip)) \
                             *0.5*subsource.width*subsource.nucleation_y])
                # nucs_x_map = [num.sin(num.deg2rad(source.strike)) \
                #               * 0.5*source.length*nuc_x_pars.scaled(fx),
                #              num.cos(num.deg2rad(source.strike)) \
                #               * 0.5*source.length*nuc_x_pars.scaled(fx)]
                # nucs_y_map = [num.cos(num.deg2rad(source.strike)) \
                #              * num.cos(num.deg2rad(source.dip)) \
                #              * 0.5*source.width*nuc_y_pars.scaled(fy),
                #              -num.sin(num.deg2rad(source.strike)) \
                #              * num.cos(num.deg2rad(source.dip)) \
                #              * 0.5*source.width*nuc_y_pars.scaled(fy)]

            #num.shape(centroids_xyz)
            # source is centered
            ax.scatter(0., 0., color='black', s=3, alpha=.5, marker='o')
            ###ax.scatter(nuc_x_pars, nuc_y_pars,  color='black', s=3, alpha=.5, marker='o')
            for fe,fn,off_nuc_xx,off_nuc_yy in zip(fes, fns, off_nuc_x, off_nuc_y):
                ax.fill(fe, fn,
                        edgecolor=(0., 0., 0.),
                        facecolor=(.5, .5, .5), alpha=0.3)
                ax.scatter(fe[:-1].mean()+off_nuc_xx[0]+off_nuc_yy[0],
                           fn[:-1].mean()+off_nuc_xx[1]+off_nuc_yy[1], color = 'g',
                           s=17, marker='o')
                # ax.scatter(centroids_xyz[:,0],
                #            centroids_xyz[:,1], color = 'm',
                #            s=2, marker='*')
                # ax.scatter(nucs_xy_xyz[:,0] + centroids_xyz[:,0] ,
                #            nucs_xy_xyz[:,1] + centroids_xyz[:,1] , color = 'b',
                #            s=0.3, marker='o')
                ax.scatter(fe[:-1].mean(),fn[:-1].mean(), s=5, marker='x')

                ax.plot(fe[0:2], fn[0:2], 'k', linewidth=1.3)

        def mapDisplacementGrid(displacements, scene):
            arr = num.full_like(scene.displacement, fill_value=num.nan)
            qt = scene.quadtree

            for syn_v, l in zip(displacements, qt.leaves):
                arr[l._slice_rows, l._slice_cols] = syn_v

            arr[scene.displacement_mask] = num.nan
            return arr

        def drawLeaves(ax, scene, offset_e=0, offset_n=0):
            rects = scene.quadtree.getMPLRectangles()
            for r in rects:
                r.set_edgecolor((.4, .4, .4))
                r.set_linewidth(.5)
                r.set_facecolor('none')
                r.set_x(r.get_x() - offset_e)
                r.set_y(r.get_y() - offset_n)
            map(ax.add_artist, rects)

            ax.scatter(scene.quadtree.leaf_coordinates[:, 0] - offset_e,
                       scene.quadtree.leaf_coordinates[:, 1] - offset_n,
                       s=.25, c='black', alpha=.1)

        def addArrow(ax, scene):
            phi = num.nanmean(scene.phi)
            los_dx = num.cos(phi + num.pi) * .0625
            los_dy = num.sin(phi + num.pi) * .0625

            az_dx = num.cos(phi - num.pi/2) * .125
            az_dy = num.sin(phi - num.pi/2) * .125

            anchor_x = .9 if los_dx < 0 else .1
            anchor_y = .85 if los_dx < 0 else .975

            az_arrow = patches.FancyArrow(
                x=anchor_x-az_dx, y=anchor_y-az_dy,
                dx=az_dx, dy=az_dy,
                head_width=.025,
                alpha=.5, fc='k',
                head_starts_at_zero=False,
                length_includes_head=True,
                transform=ax.transAxes)

            los_arrow = patches.FancyArrow(
                x=anchor_x-az_dx/2, y=anchor_y-az_dy/2,
                dx=los_dx, dy=los_dy,
                head_width=.02,
                alpha=.5, fc='k',
                head_starts_at_zero=False,
                length_includes_head=True,
                transform=ax.transAxes)

            ax.add_artist(az_arrow)
            ax.add_artist(los_arrow)

        urE, urN, llE, llN = (0., 0., 0., 0.)
        for target in sat_targets:

            if target.scene.frame.isMeter():
                off_n, off_e = map(float, latlon_to_ne_numpy(
                    target.scene.frame.llLat, target.scene.frame.llLon,
                    source.lat, source.lon))
            if target.scene.frame.isDegree():
                off_n = source.lat - target.scene.frame.llLat
                off_e = source.lon - target.scene.frame.llLon

            turE, turN, tllE, tllN = zip(
                *[(l.gridE.max()-off_e,
                   l.gridN.max()-off_n,
                   l.gridE.min()-off_e,
                   l.gridN.min()-off_n)
                  for l in target.scene.quadtree.leaves])

            turE, turN = map(max, (turE, turN))
            tllE, tllN = map(min, (tllE, tllN))
            urE, urN = map(max, ((turE, urE), (urN, turN)))
            llE, llN = map(min, ((tllE, llE), (llN, tllN)))

        def generate_plot(sat_target, result, ifig):

            scene = sat_target.scene

            fig = plt.figure()
            fig.set_size_inches(*self.size_inch)
            gs = gridspec.GridSpec(
                3, 3,
                wspace=.05, hspace=.2,
                left=.1, right=.975, top=.95,
                width_ratios=[2,2,2],
                height_ratios=[12, 8, 1])

            item = PlotItem(
                name='fig_%i' % ifig,
                attributes={'targets': [sat_target.path]},
                title=u'Satellite Surface Displacements - %s'
                      % scene.meta.scene_title,
                description=u'''Surface displacements derived from
satellite data, Scene {meta.scene_title} (id: {meta.scene_id}).
 (Left) the input data, (center) the
modelled data and (right) the model residual.'''.format(meta=scene.meta))

            stat_obs = result.statics_obs
            stat_syn = result.statics_syn['displacement.los']
            res = stat_obs - stat_syn

            if scene.frame.isMeter():
                offset_n, offset_e = map(float, latlon_to_ne_numpy(
                    scene.frame.llLat, scene.frame.llLon,
                    source.lat, source.lon))
            elif scene.frame.isDegree():
                offset_n = source.lat - scene.frame.llLat
                offset_e = source.lon - scene.frame.llLon

            im_extent = (scene.frame.E.min() - offset_e,
                         scene.frame.E.max() - offset_e,
                         scene.frame.N.min() - offset_n,
                         scene.frame.N.max() - offset_n)

            im_extent_n_cut = (0., 30.,
                         scene.frame.N.min() - offset_n,
                         scene.frame.N.max() - offset_n)

            im_extent_e_cut = (scene.frame.E.min() - offset_e,
                         scene.frame.E.max() - offset_e,
                         0.,
                         30.)

            abs_displ = num.abs([stat_obs.min(), stat_obs.max(),
                                 stat_syn.min(), stat_syn.max(),
                                 res.min(), res.max()]).max()

            cmw = cm.ScalarMappable(cmap=self.colormap)
            minmax_col=0.7
            cmw.set_clim(vmin=-minmax_col, vmax=minmax_col)
            cmw.set_array(stat_obs)

            axes = [fig.add_subplot(gs[0, 0]),
                    fig.add_subplot(gs[0, 1]),
                    fig.add_subplot(gs[0, 2])]

            ax = axes[0]
            ax.imshow(mapDisplacementGrid(stat_obs, scene),
                      extent=im_extent, cmap=self.colormap,
                      vmin=-abs_displ, vmax=abs_displ,
                      origin='lower')
            drawLeaves(ax, scene, offset_e, offset_n)
            drawSource(ax, scene)
            addArrow(ax, scene)
            initAxes(ax, scene, 'Observed')

            ax.text(.025, .025, 'Scene ID: %s' % scene.meta.scene_id,
                    fontsize=8, alpha=.7,
                    va='bottom', transform=ax.transAxes)
            if scene.frame.isDegree():
                ax.set_ylabel('Lat [°]')
            elif scene.frame.isMeter():
                ax.set_ylabel('Northing [km]')

            #ax = axes[1]
            #ax.fill(src_zz, src_nn, color=(0.5, 0.5, 0.5), alpha=0.5)
            #ax.plot(src_zz[:2], src_nn[:2], linewidth=2., color='black', alpha=0.5)

            #ax.plot(src_cnt[2], src_cnt[1],'ob')
            #ax.set_xlim(0., 30.)
            #ax.invert_xaxis()
            #ax.set_ylim(llN, urN)
            #ax.set_aspect(1/km)

            ax = axes[1]
            ax.imshow(mapDisplacementGrid(stat_syn, scene),
                      extent=im_extent, cmap=self.colormap,
                      vmin=-abs_displ, vmax=abs_displ,
                      origin='lower')
            drawLeaves(ax, scene, offset_e, offset_n)
            drawSource(ax, scene)
            addArrow(ax, scene)
            initAxes(ax, scene, 'Model')
            ax.get_yaxis().set_visible(False)

            ax = axes[2]
            ax.imshow(mapDisplacementGrid(res, scene),
                      extent=im_extent, cmap=self.colormap,
                      vmin=-abs_displ, vmax=abs_displ,
                      origin='lower')
            drawLeaves(ax, scene, offset_e, offset_n)
            drawSource(ax, scene)

            print("shape: outlines_e",num.shape(outlines_e))
            for i in num.arange(outlines_e.shape[0]):
                ax.fill(outlines_e[i, :], outlines_n[ i,:],
                        color='grey', linewidth = 0., alpha = 0.1)
            addArrow(ax, scene)
            initAxes(ax, scene, 'Residual', last_axes=True)
            ax.get_yaxis().set_visible(False)

            for ax in axes:
                ax.set_xlim(llE, urE)
                ax.set_ylim(llN, urN)

            if closeup:
                print("no closeup available for this plot configuration")
            #     closeup_zoom_factor = 3
            #     if scene.frame.isMeter():
            #         fn, fe = source.outline(cs='xy').T
            #     elif scene.frame.isDegree():
            #         fn, fe = source.outline(cs='latlon').T
            #         fn -= source.lat
            #         fe -= source.lon
            #
            #     off_n = (fn[0] + fn[1]) / closeup_zoom_factor
            #     off_e = (fe[0] + fe[1]) / closeup_zoom_factor
            #
            #     fault_size = 2*num.sqrt(max(abs(fn-off_n))**2
            #                             + max(abs(fe-off_e))**2)
            #     fault_size *= 1.25  # add 25%
            #
            #     for ax in axes:
            #         ax.set_xlim(-fault_size/closeup_zoom_factor + off_e, fault_size/closeup_zoom_factor + off_e)
            #         ax.set_ylim(-fault_size/closeup_zoom_factor + off_n, fault_size/closeup_zoom_factor + off_n)

            # ab hier noch für 2 Sources bauen (hier wird standardmäßig 1. Source gezeigt)
            xaxes = [fig.add_subplot(gs[1, 0]),
                     fig.add_subplot(gs[1, 2])]
            # strike-plane view of fault
            # best model
            ax = xaxes[0]

            ax.fill(xbest[0]+[-xbest[3]/2, xbest[3]/2, xbest[3]/2, -xbest[3]/2],
                    [xbest[2], xbest[2] , xbest[2] + xbest[4], xbest[2] + xbest[4]], color=(0.5, 0.5, 0.5), alpha=0.5)
            ax.plot(xbest[9], xbest[2]+xbest[4]/2 + xbest[10], linewidth=2., color='black', alpha=0.5)

            # bootstrap models (from harvest)
            # fault outline
            ax.plot(models[:,0]+[-models[:,3]/2, models[:,3]/2,
                      models[:,3]/2, -models[:,3]/2,
                      -models[:,3]/2],
                    [models[:,2], models[:,2],
                                   models[:,2] + models[:, 4],
                                   models[:,2] + models[:, 4],
                                  models[:,2]],
                    linewidth=.7, color='grey', alpha=0.2)
            # nucleation point
            ax.plot(models[:,0]+ models[:,3]/2*models[:,9],
                    models[:,2] + models[:, 4]/2 + models[:, 4]/2* models[:, 10], '.r', color='b', alpha=0.5)
            #print(src_zz,src_ee, llE, urE, src_cnt)
            #ax.invert_yaxis()
            #ax.set_xlim(-xbest[3]/2-xbest[3]*0.1, xbest[3]/2+xbest[3]*0.1)
            #ax.set_ylim(0, xbest[2]+xbest[4]+xbest[4]*0.4)
            #ax.plot(src_cnt[0]/km, src_cnt[2]/km,'ob')
            #ax.set_aspect(km)
            ax.invert_yaxis()

            # dip-plane view of fault plane
            ax = xaxes[1]
            #ax.fill([-xbest[3]/2, xbest[3]/2, xbest[3]/2, -xbest[3]/2],
                    #[xbest[2], xbest[2] , xbest[2] + xbest[4], xbest[2] + xbest[4]], color=(0.5, 0.5, 0.5), alpha=0.5)
            #ax.plot(xbest[9], xbest[2]+xbest[4]/2 + xbest[10], linewidth=2., color='black', alpha=0.5)

            # in-strike nukleation
            # x
            nuk_x_strplane = \
            num.cos(num.deg2rad(models[:, 7]))* models[:, 4]/2 +\
            num.cos(num.deg2rad(models[:, 7]))* models[:, 4]/2 * models[:,10]

            nuk_y_strplane = models[:,2] +\
            num.sin(num.deg2rad(models[:, 7]))* models[:, 4]/2+\
            num.sin(num.deg2rad(models[:, 7]))* models[:, 4]/2 * models[:,10]

            f_line_x = num.array([models[:, 1] ,
             models[:, 1]+ num.cos(num.deg2rad(models[:, 7]))* models[:, 4]])
            f_line_y = num.array([models[:, 2], models[:, 2] + \
                num.sin(num.deg2rad(models[:, 7]))* models[:, 4]])

            print("shape: f_line_x", f_line_x.shape, "shape: f_line_y", f_line_y.shape)
            #ax.plot(nuk_x_strplane,  nuk_y_strplane,
            #        linewidth=.7, color='k', alpha=0.5)

            ax.plot(f_line_x, f_line_y, color='grey', alpha=0.5)
            #print(src_zz,src_ee, llE, urE, src_cnt)
            #ax.invert_yaxis()
            #ax.set_xlim(-xbest[3]/2-xbest[3]*0.1, xbest[3]/2+xbest[3]*0.1)
            #ax.set_ylim(0, xbest[2]+xbest[4]+xbest[4]*0.4)
            #ax.plot(src_cnt[0]/km, src_cnt[2]/km,'ob')
            ax.set_aspect('equal')
            ax.invert_yaxis()
            # bis hier

            cax = fig.add_subplot(gs[2, :])
            cbar = fig.colorbar(cmw, cax=cax, orientation='horizontal',
                                use_gridspec=True)
                                #aspect=1000)
            cbar.set_label('LOS Displacement [m]')

            return (item, fig)

        for ifig, (sat_target, result) in enumerate(zip(sat_targets, results)):
            yield generate_plot(sat_target, result, ifig)

def get_plot_classes():
    return [SatelliteTargetDisplacement, SatelliteTargetDisplacementCloseup,
            SatelliteCustom]
