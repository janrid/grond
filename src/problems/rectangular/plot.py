from ..plot import SourceLocationMap
from pyrocko import orthodrome as od
import numpy as num
import logging
from pyrocko.plot import automap
from grond.plot.collection import PlotItem

km = 1e3
logger = logging.getLogger('grond.problem.rectangular.plot')

class SourceLocationMap_Rectangular(SourceLocationMap):
    ''' Map showing source outlines of the ensemble.'''

    name = 'source_location'

    def make(self, environ):
        cm = environ.get_plot_collection_manager()
        history = environ.get_history()
        optimiser = environ.get_optimiser()
        ds = environ.get_dataset()

        environ.setup_modelling()

        cm.create_group_automap(
            self,
            self.draw_map(environ),
            title=u'Source Outline Map',
            section='fits',
            feather_icon='map',
            description=u'''
Map showing the best source outlines of all bootstrap chains on a geographical
 map.

Colours give the relative misfits of the source models. Black colours give low
 misfit, white colours high misfit.
''')

    def draw_map(self, environ):
        item = PlotItem(name=self.name,
                title=u'Source location map',
                description=u'''
    Map showing the best source outlines (rectangular source).
    ''')

        history = environ.get_history(subset='harvest')
        optimiser = environ.get_optimiser()
        ds = environ.get_dataset()

        environ.setup_modelling()
        problem = history.problem

        gms = problem.combine_misfits(history.misfits,
            extra_correlated_weights=optimiser.get_correlated_weights(problem))
        isort = num.argsort(gms)
        gms = gms[isort]
        models = history.models[isort, :]
        xbest = models[0, :]

        nmodels = models.shape[0]

        source = problem.get_source(xbest)

        outlines_e0 = []
        outlines_n0 = []
        centroids_xy = []

        for mods in models:
            srcx = problem.get_source(mods)
            fe, fn = srcx.outline(cs='lonlat').T
            # centroids relative to reference location
            centroid =  num.array([srcx.east_shift \
                     + num.cos(num.deg2rad(srcx.strike)) \
                     * num.cos(num.deg2rad(srcx.dip)) * 0.5*srcx.width,
                     srcx.north_shift - num.sin(num.deg2rad(srcx.strike)) \
                     * num.cos(num.deg2rad(srcx.dip))* 0.5*srcx.width])

            if num.shape(centroids_xy)[0]==0:
                centroids_xy = centroid
                outlines_e0 = fe
                outlines_n0 = fn
            else:
                centroids_xy = num.vstack((centroids_xy, centroid))
                outlines_e0 = num.vstack((outlines_e0, fe))
                outlines_n0 = num.vstack((outlines_n0, fn))

        ref_pos = num.array([source.lat, source.lon])
        centroids_x = num.array([centroids_xy[i][0] for i in range(len(centroids_xy))])
        centroids_y = num.array([centroids_xy[i][1] for i in range(len(centroids_xy))])
        centroids_latlon = od.ne_to_latlon(ref_pos[0], ref_pos[1], centroids_y,
            centroids_x) # get geographical coordinates of centroids
        ref_loc = ref_pos.reshape([2,1])
        locations = num.hstack((centroids_latlon, ref_loc))
        lat, lon = od.geographic_midpoint(locations[0], locations[1])

        if self.radius is None:
                radius = od.distance_accurate50m_numpy(
                            lat[num.newaxis], lon[num.newaxis],
                            locations[0], locations[1]).max()
                radius *= 1.5 #1.1

        if radius < 20.*km:
                logger.warn(
                    'Radius of map defaulting to 20 km')
                radius = 20*km

        m = automap.Map(
                width=self.size_cm[0],
                height=self.size_cm[1],
                lat=lat,
                lon=lon,
                radius=radius,
                show_topo=self.show_topo,
                show_grid=self.show_grid,
                show_rivers=self.show_rivers,
                color_wet=(216, 242, 254),
                color_dry=(238, 236, 230))

        # add outlines/centroids to map
        color, alpha = self.get_shading(nmodels)

        for index in num.arange(nmodels):
            m.gmt.psxy(
                    in_rows=num.array([outlines_e0[index],
                        outlines_n0[index]]).T,
                    L='',
                    G=color[index],
                    t=alpha[index],
                    *m.jxyr)

        yield (item, m)
