import numpy as num
import logging
#import sys
from pyrocko import gf, util
from pyrocko.guts import String, Float, Dict, Int
from optparse import OptionParser

from grond.meta import expand_template, Parameter

from ..base import Problem, ProblemConfig

guts_prefix = 'grond'
logger = logging.getLogger('grond.problems.multirectangular.problem')
km = 1e3
as_km = dict(scale_factor=km, scale_unit='km')


#(options, args) = parser.parse_args(sys.argv[1:])

class MultiRectangularProblemConfig(ProblemConfig):

    ranges = Dict.T(String.T(), gf.Range.T())
    decimation_factor = Int.T(default=1)
    distance_min = Float.T(default=0.)
    nsources = Int.T(default=1)

    def get_problem(self, event, target_groups, targets):
        base_source = gf.RectangularSource.from_pyrocko_event(
            event,
            anchor='top',
            decimation_factor=self.decimation_factor)

        subs = dict(
            event_name=event.name,
            event_time=util.time_to_str(event.time))

        problem = MultiRectangularProblem(
            name=expand_template(self.name_template, subs),
            base_source=base_source,
            distance_min=self.distance_min, 
            target_groups=target_groups,
            targets=targets,
            ranges=self.ranges,
            norm_exponent=self.norm_exponent)

        return problem


class MultiRectangularProblem(Problem):
    nsources = 2
    #for i in range(0, 100):
        #if "--nsources="+str(i) in sys.argv:
            #nsources = int(i)
    #if nsources is None:
        #print('input --nsources= to go command missing')

    problem_parameters = []
    problem_waveform_parameters = []

    for i in range(1,nsources+1):
        problem_parameters.append(Parameter('north_shift%s' % i,
                                            'm',
                                            label='Northing %s' %i,
                                            **as_km))
        problem_parameters.append(Parameter('east_shift%s' % i,
                                            'm',
                                            label='Easting %s' %i,
                                            **as_km))
        problem_parameters.append(Parameter('depth%s' % i, 'm',
                                            label='Depth %s' %i,
                                            **as_km))
        problem_parameters.append(Parameter('length%s' % i,
                                            'm',
                                            label='Length %s' %i,
                                            **as_km))
        problem_parameters.append(Parameter('width%s' % i, 'm',
                                            label='Width %s' %i,
                                            **as_km))
        problem_parameters.append(Parameter('dip%s' % i, 'deg',
                                            label='Dip %s' %i))
        problem_parameters.append(Parameter('strike%s' % i,
                                            'deg',
                                            label='Strike %s' %i))
        problem_parameters.append(Parameter('rake%s' % i, 'deg',
                                            label='Rake %s' %i))
        problem_parameters.append(Parameter('slip%s' % i, 'm',
                                            label='Slip %s' %i))

        problem_parameters.append(Parameter('nucleation_x%s' % i,
                                            'offset',
                                            label='Nucleation X %s' %i))
        problem_parameters.append(Parameter('nucleation_y%s' % i,
                                            'offset',
                                            label='Nucleation Y %s' %i))
        problem_parameters.append(Parameter('time%s' % i, 's',
                                             label='Time %s' %i))

    dependants = []
    distance_min = Float.T(default=0.0)

    def pack(self, source):
        arr = self.get_parameter_array(source)
        for ip, p in enumerate(self.parameters):
            if p.name == 'time':
                arr[ip] -= self.base_source.time
        return arr

    def get_source(self, x, i):
        d = self.get_parameter_dict(x[0+12*i:12+i*12], nsources=2)#x, nsources=self.nsources) # not looking nice but needed for correct usage of params #from branch multisource_new
        i = i+1 # for realistic numbers
        p = {}
        for k in self.base_source.keys():
            if k in d:
                p[k] = float(
                    self.ranges[k+str(i)].make_relative(self.base_source[k],
                                                        d[k]))
        source = self.base_source.clone(**p)

        return source

    def random_uniform(self, xbounds):
        x = num.zeros(self.nparameters)
        for i in range(self.nparameters):
            x[i] = num.random.uniform(xbounds[i, 0], xbounds[i, 1])

        return x

    def preconstrain(self, x):
        # source = self.get_source(x)
        # if any(self.distance_min > source.distance_to(t)
        #        for t in self.targets):
            # raise Forbidden()
        return x

    @classmethod
    def get_plot_classes(cls):
        plots = super(MultiRectangularProblem, cls).get_plot_classes()
        return plots


__all__ = '''
    MultiRectangularProblem
    MultiRectangularProblemConfig
'''.split()
