import logging
import numpy as num
from matplotlib import cm, gridspec

from grond.plot.config import PlotConfig
from grond.plot.collection import PlotItem

from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator, StrMethodFormatter
from matplotlib.ticker import FormatStrFormatter
from matplotlib import patches
from pyrocko.guts import Tuple, Float, String, Int, Bool, StringChoice
import matplotlib.colors as mcolors


logger = logging.getLogger('grond.targets.satellite.plot')

km = 1e3
d2r = num.pi/180.
guts_prefix = 'grond'

def PlotAnyColor(col_data, shad_data, cmap_in, shad_lim, tick_step = 500,
                 contrast = 1., mask = 0, data_limits = (-0.5, 0.5)):
    '''Plot color on shaded relief.'''
    
    from scipy.ndimage import convolve as im_conv
    cmap = plt.get_cmap(cmap_in)
    def norm_0_1(data, data_limits = (0, 0)):
        if num.sum(num.abs(data_limits))==0:
            data_out = data - num.nanmin(num.nanmin(data))
            data_out /= float(num.nanmax(num.nanmax(data_out)))
        else:
            data_out = num.where(data<data_limits[0], data_limits[0], data)
            data_out = num.where(data>data_limits[1], data_limits[1], data_out)

            data_out = data_out - data_limits[0]
            data_out /= (data_limits[1]-data_limits[0])
            
        return data_out

    # 1) Creating the shading
    # Light source from somewhere above - psychologically the best choice
    # from upper left
    ramp = num.array([[1, 0], [0, -1.]]) * contrast

    # convolution of two 2-dimensional arrays    
    shad = im_conv(shad_data*1e-3, ramp.T)
    shad = -1. * shad
   
    # if there are strong artifical edges in the data, shades get
    # dominated by them. Cutting off the largest and smallest 2% of
    # shades helps
    percentile2 = num.quantile(shad, 0.02)
    percentile98 = num.quantile(shad, 0.98)
    
    shad = num.where( shad > percentile98, percentile98, shad)
    shad = num.where( shad < percentile2, percentile2, shad)

    #normalize to range [0 1]

    shad = norm_0_1(shad)
    shad = num.where(mask == 1., shad, num.nan) 
    
    # apply reduced range to shad_lim to balance gray color
    shad *= (shad_lim[1] - shad_lim[0])
    shad += shad_lim[0]
    
    # 2) Create colormap from data and colormap
    # normalize
    norm_col_data =  norm_0_1(col_data, data_limits)

    # create ticks for plotting - real values for the labels
    # and their position in normed data for the ticks
    
    # range of heights in topo map
    if  num.sum(num.abs(data_limits))==0:
        range_col_data = num.nanmax(col_data) - num.nanmin(col_data)
        n_ticks_col = num.floor(range_col_data / tick_step) + 1
        ticks_col = num.arange(n_ticks_col) * tick_step + num.nanmin(col_data)
        if num.nanmin(col_data)<0:
            min_tick = num.nonzero((ticks_col>=0) & (ticks_col<tick_step))
            ticks_col_sh = ticks_col[1:] - ticks_col[min_tick[0][0]]
            ticks_col_shift = num.hstack((ticks_col_sh, ticks_col_sh[-1] + tick_step))
            ticks_col = ticks_col_shift
            
        ticks_col2 = num.hstack((ticks_col, num.nanmax(col_data)))
        ticks_col = num.hstack((num.nanmin(col_data), ticks_col2))
    else:
        range_col_data = data_limits[1] - data_limits[0]
        
        n_ticks_col = num.floor(range_col_data / tick_step) + 1
        ticks_col = num.arange(n_ticks_col) * tick_step + data_limits[0]
        if data_limits[0]<0:
            min_tick = num.nonzero((ticks_col>=0) & (ticks_col<tick_step))
            ticks_col_sh = ticks_col[1:] - ticks_col[min_tick[0][0]]
            ticks_col_shift = num.hstack((ticks_col_sh, ticks_col_sh[-1] + tick_step))
            ticks_col = ticks_col_shift
            
        ticks_col2 = num.hstack((ticks_col, data_limits[1]))
        ticks_col = num.hstack((data_limits[0], ticks_col2))

    norm_pos_ticks = norm_0_1(ticks_col)
    norm_pos_ticks = norm_pos_ticks[1:-1]
    ticks_col = ticks_col[1:-1]
    
    # 3) Combine color and shading 
    rgb_map = cmap(norm_col_data)
    for chan in num.arange(3):
        rgb_map[:,:,chan]= num.where(num.isnan(norm_col_data),1., rgb_map[:,:,chan])

    for chan in num.arange(3):
        rgb_map[:,:,chan] = num.multiply(rgb_map[:,:,chan], shad)

        
    return rgb_map, norm_pos_ticks, ticks_col

def scale_axes(axis, scale, offset=0., suffix=''):
    from matplotlib.ticker import ScalarFormatter

    class FormatScaled(ScalarFormatter):

        @staticmethod
        def __call__(value, pos):
            return '{:,.1f}{:}'.format((offset + value) * scale, suffix)\
                .replace(',', ' ')

    axis.set_major_formatter(FormatScaled())

class SatelliteTargetDisplacement(PlotConfig):
    ''' Maps showing surface displacements from satellite and modelled data '''

    name = 'satellite'
    font_size = Float.T(
        default = 10.)
        
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
    fit = StringChoice.T(
        default='best', choices=['best', 'mean'],
        help='Show the \'best\' or \'mean\' fits and source model from the'
         ' ensamble.')
    
    shown_disp_unit = StringChoice.T(
        default='m', choices=['m', 'mm', 'cm'],
        help='Show results in \'m\' , \'cm\' or \'mm\' ')
    color_tick_step = Float.T(
        default = 10.,
        help='Separation of tick at colorbar in the unit that is to be shown')
    
    show_topo = Bool.T(
        default=False,
        help='show topography')
    
    show_leaf_centres = Bool.T(
        default=False,
        help='show the center points of Quadtree leaves')
    
    set_map_limits = Tuple.T(
        4, Float.T(),
        default=(0., 0., 0., 0.))
    
    common_color_scale = Bool.T(
        default=False,
        help='Results shown with common color scale for all satellite data sets (based on the data)')
    
    source_outline_color = String.T(
        default='grey',
        help='Choose color of source outline from named matplotlib Colors')
    
    

    def make(self, environ):
        cm = environ.get_plot_collection_manager()
        history = environ.get_history(subset='harvest')
        optimiser = environ.get_optimiser()
        ds = environ.get_dataset()

        environ.setup_modelling()

        cm.create_group_mpl(
            self,
            self.draw_static_fits(ds, history, optimiser),
            title=u'InSAR Displacements',
            section='fits',
            feather_icon='navigation',
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
edge marking the upper fault edge. Complete data extent is shown.
''')

    def draw_static_fits(self, ds, history, optimiser, closeup=False):
        from pyrocko.orthodrome import latlon_to_ne_numpy
        problem = history.problem

        sat_targets = problem.satellite_targets
        for target in sat_targets:
            target.set_dataset(ds)

        if self.fit == 'best':
            source = history.get_best_source()
            model = history.get_best_model()
        elif self.fit == 'mean':
            source = history.get_mean_source()
            model = history.get_mean_model()

        results = problem.evaluate(model, targets=sat_targets)

        def initAxes(ax, scene, title, last_axes=False):
            ax.set_title(title, fontsize = self.font_size)
            ax.tick_params(length=2)
            
            if scene.frame.isMeter():
                ax.set_xlabel('Easting [km]', fontsize = self.font_size)
                scale_x = dict(scale=1./km)
                scale_y = dict(scale=1./km)
                if not self.relative_coordinates:
                    import utm
                    utm_E, utm_N, utm_zone, utm_zone_letter =\
                        utm.from_latlon(source.effective_lat,
                                        source.effective_lon)
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
                ax.set_aspect(1./num.cos(source.effective_lat*d2r))
                #ax.apply_aspect()

                scale_x = dict(scale=1.00, suffix='°')
                scale_y = dict(scale=1., suffix='°')
                #if not self.relative_coordinates:
                #    scale_x['offset'] = source.effective_lon
                #    scale_y['offset'] = source.effective_lat

            # reducing the ticks to 4 or 5 
            # caused problems with the position of yticks, which
            # where offset (looks like the stretching of the yaxis
            # does not go well together with the relabeling, but can't
            # be sure.
            
            #nticks_lon = 4 if abs(scene.frame.llLon) >= 100 else 5
            #ax.xaxis.set_major_locator(MaxNLocator(nticks_lon))
            #ax.yaxis.set_major_locator(MaxNLocator(5))
            
            #scale_axes(ax.get_xaxis(), **scale_x)
            #scale_axes(ax.get_yaxis(), **scale_y)
            

            

        def drawSource(ax, scene):
            if scene.frame.isMeter():
                fn, fe = source.outline(cs='xy').T
                fn -= fn.mean()
                fe -= fe.mean()
            elif scene.frame.isDegree():
                fn, fe = source.outline(cs='latlon').T
                if self.relative_coordinates:
                    fn -= source.effective_lat
                    fe -= source.effective_lon
            

            # source is centered
            if self.relative_coordinates:
                ax.scatter(0., 0., color='black', s=3, alpha=.5, marker='o')
            else:
                ax.scatter(source.effective_lon, source.effective_lat, color='black', s=3, alpha=.5, marker='o')
            
            if self.source_outline_color in mcolors.cnames:
                ax.fill(fe, fn,
                    facecolor=self.source_outline_color, alpha=0.7)
            else:
                #Todo: Logger message "choosen source color not a named matplotlib color
                # defaulting to grey 
                ax.fill(fe, fn,
                    facecolor='gray', alpha=0.7)
                
            ax.plot(fe[0:2], fn[0:2], 'k', linewidth=1.)

        def mapDisplacementGrid(displacements, scene, show_topo = False, 
                                abs_displ = 0., tick_step = 0.01):
            radar_map = num.full_like(scene.displacement, fill_value=num.nan)
            qt = scene.quadtree
            for syn_v, l in zip(displacements, qt.leaves):
                radar_map[l._slice_rows, l._slice_cols] = syn_v

            radar_map[scene.displacement_mask] = num.nan
                
            if show_topo:
                    
                elevation = scene.get_elevation()
                shad_lim = [0.4, .99]
                mask = num.where(elevation == 0., 0., 1. )
              
                    
                radar_map, norm_pos_ticks, ticks_col = PlotAnyColor(radar_map, 
                                                                    elevation,
                                                                    self.colormap, 
                                                                    shad_lim, 
                                                                    tick_step = tick_step,
                                                                    contrast = 1, 
                                                                    mask = mask,
                                                                    data_limits = (-abs_displ, abs_displ))

                
            
            return radar_map, norm_pos_ticks, ticks_col

        def drawLeaves(ax, scene, offset_e=0., offset_n=0.):
            rects = scene.quadtree.getMPLRectangles()
            for r in rects:
                r.set_edgecolor((.4, .4, .4))
                r.set_linewidth(.5)
                r.set_facecolor('none')
                r.set_x(r.get_x() - offset_e)
                r.set_y(r.get_y() - offset_n)
            map(ax.add_artist, rects)

            if self.show_leaf_centres:
                #ax.scatter(scene.quadtree.leaf_coordinates[:, 0] - offset_e,
                        #scene.quadtree.leaf_coordinates[:, 1] - offset_n,
                        #s=.25, c='black', alpha=.1)
                ax.scatter(scene.quadtree.leaf_coordinates[:, 0] ,
                        scene.quadtree.leaf_coordinates[:, 1] ,
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
                alpha=.9, fc='k',
                head_starts_at_zero=False,
                length_includes_head=True,
                transform=ax.transAxes)

            los_arrow = patches.FancyArrow(
                x=anchor_x-az_dx/2, y=anchor_y-az_dy/2,
                dx=los_dx, dy=los_dy,
                head_width=.02,
                alpha=.9, fc='k',
                head_starts_at_zero=False,
                length_includes_head=True,
                transform=ax.transAxes)

            ax.add_artist(az_arrow)
            ax.add_artist(los_arrow)

        urE, urN, llE, llN = (0., 0., 0., 0.)
        
        targets_displ_max = 0.
        targets_displ_min = 0.
        for target in sat_targets:
            if num.sum(self.set_map_limits)==0:
                        #ToDo: catch lower limit above upper limit
                        #Todo: catch limits entirely outside of map extent 
                if target.scene.frame.isMeter():
                    off_n, off_e = map(float, latlon_to_ne_numpy(
                        target.scene.frame.llLat, target.scene.frame.llLon,
                        source.effective_lat, source.effective_lon))
                if target.scene.frame.isDegree():
                    off_n = source.effective_lat - target.scene.frame.llLat
                    off_e = source.effective_lon - target.scene.frame.llLon

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
                    
            elif num.sum(self.set_map_limits)!=0:
                if target.scene.frame.isMeter():
                    off_n, off_e = map(float, latlon_to_ne_numpy(
                        self.set_map_limits[2], self.set_map_limits[0],
                        source.effective_lat, source.effective_lon))
                if target.scene.frame.isDegree():
                    off_n = source.effective_lat - self.set_map_limits[2]
                    off_e = source.effective_lon - self.set_map_limits[0]
                    
                if self.relative_coordinates:    
                    urE = self.set_map_limits[1] - source.effective_lon
                    urN = self.set_map_limits[3] - source.effective_lat
                    llE = -off_e
                    llN = -off_n
                else:
                    urE = self.set_map_limits[1]
                    urN = self.set_map_limits[3]
                    llE = self.set_map_limits[0]
                    llN = self.set_map_limits[2]
                    
          
            if targets_displ_min>num.nanmin(target.scene.displacement):
                targets_displ_min = num.nanmin(target.scene.displacement)
                
            if targets_displ_max<num.nanmax(target.scene.displacement):
                targets_displ_max = num.nanmax(target.scene.displacement)  
                
            if self.shown_disp_unit == 'cm':
                tick_step = self.color_tick_step / 100.
            if self.shown_disp_unit == 'mm':
                tick_step = self.color_tick_step / 1000.
 
        def generate_plot(sat_target, result, ifig):

            scene = sat_target.scene
            parameters = {'xtick.labelsize': self.font_size,
                          'ytick.labelsize': self.font_size}
            

            plt.rcParams.update(parameters)
            fig = plt.figure()
            fig.set_size_inches(*self.size_inch)
            gs = gridspec.GridSpec(
                2, 3,
                wspace=.15, hspace=.2,
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
                    source.effective_lat, source.effective_lon))
            elif scene.frame.isDegree():
                offset_n = source.effective_lat - scene.frame.llLat
                offset_e = source.effective_lon - scene.frame.llLon
            if self.relative_coordinates:
                im_extent = (scene.frame.E.min() - offset_e,
                            scene.frame.E.max() - offset_e,
                            scene.frame.N.min() - offset_n,
                            scene.frame.N.max() - offset_n)
            else:
                im_extent = (scene.frame.E.min() + scene.frame.llLon,
                            scene.frame.E.max() + scene.frame.llLon,
                            scene.frame.N.min() + scene.frame.llLat,
                            scene.frame.N.max() + scene.frame.llLat)
                
            if not self.common_color_scale:
                abs_displ = num.abs([num.nanmin(stat_obs), num.nanmax(stat_obs),
                                    num.nanmin(stat_syn), num.nanmax(stat_syn),
                                    num.nanmin(res), num.nanmax(res)]).max()
            else:
                abs_displ = num.max([num.abs(targets_displ_min), 
                                    num.abs(targets_displ_max)])
                
            if self.shown_disp_unit == 'm':
                shown_abs_displ = abs_displ
            elif self.shown_disp_unit=='mm':
                shown_abs_displ = abs_displ * 1e3
            elif self.shown_disp_unit=='cm':
                shown_abs_displ = abs_displ * 1e2

            cmw = cm.ScalarMappable(cmap=self.colormap)
            cmw.set_clim(vmin=-shown_abs_displ, vmax=shown_abs_displ)
            cmw.set_array(stat_obs)
            print(im_extent)
            axes = [fig.add_subplot(gs[0, 0]),
                    fig.add_subplot(gs[0, 1]),
                    fig.add_subplot(gs[0, 2])]

            ax = axes[0]

            disp_map, norm_pos_ticks, ticks_col = mapDisplacementGrid(stat_obs, scene, 
                                          self.show_topo, 
                                          abs_displ = abs_displ,
                                          tick_step = tick_step)
            ax.imshow(disp_map,
                      extent=im_extent, cmap=self.colormap,
                      vmin=-abs_displ, vmax=abs_displ,
                      origin='lower')

            drawLeaves(ax, scene, offset_e, offset_n)
            drawSource(ax, scene)
            addArrow(ax, scene)
            
            
            ax.set_xlim(llE, urE)
            ax.set_ylim(llN, urN)
            initAxes(ax, scene, 'Observed')
            
            # Problems:
            # maybe too many ticks, because the reduction causes trouble
            # the second digit after the comma is needed for small areas
            # but doesn;t show :/
            # these I tried - but they cause the relative coordinates to
            # show up 
            #ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            #plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:.2f}'))

        
            ax.text(.025, .025, 'Scene ID: %s' % scene.meta.scene_id,
                    fontsize=8, alpha=.7,
                    va='bottom', transform=ax.transAxes)
            if scene.frame.isMeter():
                ax.set_ylabel('Northing [km]', fontsize = self.font_size)

            ax = axes[1]

            disp_map, norm_pos_ticks, ticks_col = mapDisplacementGrid(stat_syn, scene, 
                                          self.show_topo, 
                                          abs_displ = abs_displ,
                                          tick_step = tick_step)
            ax.imshow(disp_map,  
                      extent=im_extent, 
                      cmap=self.colormap,
                      vmin=-abs_displ, vmax=abs_displ,
                      origin='lower')
            drawLeaves(ax, scene, offset_e, offset_n)
            drawSource(ax, scene)
            addArrow(ax, scene)
            ax.set_xlim(llE, urE)
            ax.set_ylim(llN, urN)
            initAxes(ax, scene, 'Model')
            ax.get_yaxis().set_visible(False)

            ax = axes[2]
            disp_map, norm_pos_ticks, ticks_col = mapDisplacementGrid(res, scene, 
                                          self.show_topo, 
                                          abs_displ = abs_displ,
                                          tick_step = tick_step)
            ax.imshow(disp_map,   
                      extent=im_extent, cmap=self.colormap,
                      vmin=-abs_displ, vmax=abs_displ,
                      origin='lower')
            drawLeaves(ax, scene, offset_e, offset_n)


            drawSource(ax, scene)
            addArrow(ax, scene)
            ax.get_yaxis().set_visible(False)
            ax.set_xlim(llE, urE)
            ax.set_ylim(llN, urN)
            initAxes(ax, scene, 'Residual', last_axes=True)
           

            if closeup:
                if scene.frame.isMeter():
                    fn, fe = source.outline(cs='xy').T
                elif scene.frame.isDegree():
                    fn, fe = source.outline(cs='latlon').T
                    fn -= source.effective_lat
                    fe -= source.effective_lon

                if fn.size > 1:
                    off_n = (fn[0] + fn[1]) / 2
                    off_e = (fe[0] + fe[1]) / 2
                else:
                    off_n = fn[0]
                    off_e = fe[0]

                fault_size = 2*num.sqrt(max(abs(fn-off_n))**2
                                        + max(abs(fe-off_e))**2)
                fault_size *= self.map_scale
                if fault_size == 0.0:
                    extent = (scene.frame.N[-1] + scene.frame.E[-1]) / 2
                    fault_size = extent * .25

                for ax in axes:
                    ax.set_xlim(-fault_size/2 + off_e, fault_size/2 + off_e)
                    ax.set_ylim(-fault_size/2 + off_n, fault_size/2 + off_n)

           
            cax = fig.add_subplot(gs[1, 1])
            
            if self.shown_disp_unit == 'mm':
                ticks_col *= 1e3 
            cbar = fig.colorbar(cmw, cax=cax, ticks=ticks_col, orientation='horizontal',
                                use_gridspec=True)
            cbar.ax.set_xticklabels(ticks_col)
            
            if self.shown_disp_unit == 'm':
                cbar.set_label('LOS Displacement [m]', fontsize = self.font_size)
                plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}')) # 2 decimal places

            if self.shown_disp_unit == 'cm':
                cbar.set_label('LOS Displacement [cm]', fontsize = self.font_size)
                plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}')) 

            if self.shown_disp_unit == 'mm':
                cbar.set_label('LOS Displacement [mm]' , fontsize = self.font_size)
                plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}')) # 2 decimal places

                
            return (item, fig)

        for ifig, (sat_target, result) in enumerate(zip(sat_targets, results)):
            yield generate_plot(sat_target, result, ifig)

class SatelliteTargetDisplacementCloseup(SatelliteTargetDisplacement):
    ''' Close-up of satellite surface displacements and modelled data. '''
    name = 'satellite_closeup'
    map_scale = Float.T(
        default=2.,
        help='Scale the map surroundings, larger value zooms out.')
    
    font_size = Float.T(
        default = 10.)
        
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
    fit = StringChoice.T(
        default='best', choices=['best', 'mean'],
        help='Show the \'best\' or \'mean\' fits and source model from the'
         ' ensamble.')
    
    shown_disp_unit = StringChoice.T(
        default='m', choices=['m', 'mm', 'cm'],
        help='Show results in \'m\' , \'cm\' or \'mm\' ')
    color_tick_step = Float.T(
        default = 10.,
        help='Separation of tick at colorbar in the unit that is to be shown')
    
    show_topo = Bool.T(
        default=False,
        help='show topography')
    
    show_leaf_centres = Bool.T(
        default=False,
        help='show the center points of Quadtree leaves')
    
    set_map_limits = Tuple.T(
        4, Float.T(),
        default=(0., 0., 0., 0.))
    
    common_color_scale = Bool.T(
        default=False,
        help='Results shown with common color scale for all satellite data sets (based on the data)')
    

    def make(self, environ):
        cm = environ.get_plot_collection_manager()
        history = environ.get_history(subset='harvest')
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
        
def get_plot_classes():
    return [SatelliteTargetDisplacement, SatelliteTargetDisplacementCloseup]
