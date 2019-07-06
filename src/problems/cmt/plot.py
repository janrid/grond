from ..plot import SourceLocationMap
from pyrocko import gmtpy, orthodrome as od
import numpy as num
import logging
from pyrocko.plot import automap
from grond.plot.collection import PlotItem

km = 1e3
logger = logging.getLogger('grond.problem.cmt.plot')

class SourceLocationMap_CMT(SourceLocationMap):
    ''' Map showing source positions of the ensemble.'''

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
            title=u'Source Position Map',
            section='fits',
            feather_icon='map',
            description=u'''
Map showing the best source centroid positions of all bootstrap chains on
 a geographical map.

Colours give the relative misfits of the source models. Black colours give low
 misfit, white colours high misfit.
''')

    def draw_map(self, environ):
        item = PlotItem(name=self.name,
                title=u'Source location map',
                description=u'''
    Map showing the best source centroid positions.
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

        centroids_xy = []

        for mods in models:
            srcx = problem.get_source(mods)
            centroid = num.array([srcx.east_shift, srcx.north_shift])

            if num.shape(centroids_xy)[0]==0:
                centroids_xy = centroid
            else:
                centroids_xy = num.vstack((centroids_xy, centroid))

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

        factor_symbl_size = 2.0
        mag=6.4
        ev_symb = 'c'+str((mag*factor_symbl_size)*8 / gmtpy.cm)+'p'
        for index in num.arange(nmodels):
            m.gmt.psxy(
                in_rows=[[centroids_latlon[1][index],
                    centroids_latlon[0][index]]],
                S=ev_symb,
                G=color[index],
                t=alpha[index],
                #W='1p,darkred',
                *m.jxyr)

        yield (item, m)
