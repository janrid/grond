#!/usr/bin/env python

import sys
import os.path as op
import logging
from optparse import OptionParser

from pyrocko import util, marker

import grond

logger = logging.getLogger('grond.main')
km = 1e3


def d2u(d):
    if isinstance(d, dict):
        return dict((k.replace('-', '_'), v) for (k, v) in d.items())
    else:
        return d.replace('-', '_')


subcommand_descriptions = {
    'scenario': 'create an scenario project',
    'init': 'initialise new project structure or print configuration',
    'events': 'print available event names for given configuration',
    'check': 'check data and configuration',
    'go': 'run Grond optimisation',
    'forward': 'run forward modelling',
    'harvest': 'manually run harvesting',
    'plot': 'plot optimisation result',
    'movie': 'visualize optimiser evolution',
    'export': 'export results',
    'report': 'create result report',
    'diff': 'compare two configs or other normalized Grond YAML files',
    'qc-polarization': 'check sensor orientations with polarization analysis',
    'upgrade-config': 'upgrade config file to the latest version of Grond',
}

subcommand_usages = {
    'scenario': 'scenario [options] <project_dir>',
    'init': 'init [options] <project_dir>',
    'events': 'events <configfile>',
    'check': 'check <configfile> <eventnames> ... [options]',
    'go': 'go <configfile> <eventnames> ... [options]',
    'forward': (
        'forward <rundir> [options]',
        'forward <configfile> <eventnames> ... [options]'),
    'harvest': 'harvest <rundir> [options]',
    'plot': 'plot <plotnames> <rundir> [options]',
    'movie': 'movie <rundir> <xpar> <ypar> <filetemplate> [options]',
    'export': 'export (best|mean|ensemble|stats) <rundirs> ... [options]',
    'report': (
        'report <rundir> ... [options]',
        'report <configfile> <eventnames> ...'),
    'diff': 'diff <left_path> <right_path>',
    'qc-polarization': 'qc-polarization <configfile> <eventname> '
                       '<target_group_path> [options]',
    'upgrade-config': 'upgrade-config <configfile>',
}

subcommands = subcommand_descriptions.keys()

program_name = 'grond'

usage_tdata = d2u(subcommand_descriptions)
usage_tdata['program_name'] = program_name

usage = '''%(program_name)s <subcommand> [options] [--] <arguments> ...

Subcommands:

    scenario        %(scenario)s
    init            %(init)s
    events          %(events)s
    check           %(check)s
    go              %(go)s
    forward         %(forward)s
    harvest         %(harvest)s
    plot            %(plot)s
    movie           %(movie)s
    export          %(export)s
    report          %(report)s
    diff            %(diff)s
    qc-polarization %(qc_polarization)s
    upgrade-config  %(upgrade_config)s

To get further help and a list of available options for any subcommand run:

    %(program_name)s <subcommand> --help

''' % usage_tdata


def main(args=None):
    if not args:
        args = sys.argv

    args = list(sys.argv)
    if len(args) < 2:
        sys.exit('Usage: %s' % usage)

    args.pop(0)
    command = args.pop(0)

    if command in subcommands:
        globals()['command_' + d2u(command)](args)

    elif command in ('--help', '-h', 'help'):
        if command == 'help' and args:
            acommand = args[0]
            if acommand in subcommands:
                globals()['command_' + acommand](['--help'])

        sys.exit('Usage: %s' % usage)

    else:
        die('no such subcommand: %s' % command)


def add_common_options(parser):
    parser.add_option(
        '--loglevel',
        action='store',
        dest='loglevel',
        type='choice',
        choices=('critical', 'error', 'warning', 'info', 'debug'),
        default='info',
        help='set logger level to '
             '"critical", "error", "warning", "info", or "debug". '
             'Default is "%default".')


def process_common_options(options):
    util.setup_logging(program_name, options.loglevel)


def cl_parse(command, args, setup=None, details=None):
    usage = subcommand_usages[command]
    descr = subcommand_descriptions[command]

    if isinstance(usage, str):
        usage = [usage]

    susage = '%s %s' % (program_name, usage[0])
    for s in usage[1:]:
        susage += '\n%s%s %s' % (' '*7, program_name, s)

    description = descr[0].upper() + descr[1:] + '.'

    if details:
        description = description + '\n\n%s' % details

    parser = OptionParser(usage=susage, description=description)

    if setup:
        setup(parser)

    add_common_options(parser)
    (options, args) = parser.parse_args(args)
    process_common_options(options)
    return parser, options, args


def die(message, err=''):
    if err:
        sys.exit('%s failed: %s \n %s' % (program_name, message, err))
    else:
        sys.exit('%s failed: %s' % (program_name, message))


def help_and_die(parser, message):
    parser.print_help(sys.stderr)
    sys.stderr.write('\n')
    die(message)


def command_scenario(args):

    STORE_STATIC = 'crust2_ib_static'
    STORE_WAVEFORMS = 'crust2_ib'

    def setup(parser):
        parser.add_option(
            '--waveforms', dest='waveforms', action='store_true',
            help='add waveform configuration. '
                 '(default)')
        parser.add_option(
            '--insar', dest='insar', action='store_true',
            help='add InSAR displacement scenes using kite containers. '
                 '(see https://pyrocko.org)')
        parser.add_option(
            '--gnss', dest='gnss', action='store_true',
            help='add GNSS campaign data using kite containers. '
                 '(see https://pyrocko.org)')
        parser.add_option(
            '--nstations', dest='nstations', type=int, default=20,
            help='number of seismic stations to create (default: %default)')
        parser.add_option(
            '--gnss_nstations', dest='gnss_nstations', type=int, default=20,
            help='number of GNSS campaign stations to create'
                 ' (default: %default)')
        parser.add_option(
            '--nevents', dest='nevents', type=int, default=1,
            help='number of events to create (default: %default)')
        parser.add_option(
            '--lat', dest='lat', type=float, default=41.0,
            help='center latitude of the scenario (default: %default)')
        parser.add_option(
            '--lon', dest='lon', type=float, default=33.3,
            help='center latitude of the scenario (default: %default)')
        parser.add_option(
            '--radius', dest='radius', type=float, default=200.,
            help='radius of the the scenario in [km] (default: %default)')
        parser.add_option(
            '--source', dest='source', type=str, default='dc',
            help='source to generate \'dc\' (double couple)'
                 ' or\'rectangular\' (rectangular finite fault)'
                 ' (default: \'%default\')')
        parser.add_option(
            '--gf_waveforms', dest='store_waveforms', type=str,
            default=STORE_WAVEFORMS,
            help='Green\'s function store for waveform modelling, '
                 '(default: %default)')
        parser.add_option(
            '--gf_static', dest='store_statics', type=str,
            default=STORE_STATIC,
            help='Green\'s function store for static modelling, '
                 '(default: %default)')
        parser.add_option(
            '--force', dest='force', action='store_true',
            help='overwrite existing project folder.')

    parser, options, args = cl_parse('scenario', args, setup)

    if len(args) == 1:
        project_dir = op.join(op.curdir, args[0])
    else:
        parser.print_help()
        sys.exit(1)

    from grond import scenario as grond_scenario
    scenario = grond_scenario.GrondScenario(
        project_dir,
        center_lat=options.lat, center_lon=options.lon,
        radius=options.radius*km)

    if not options.waveforms and not options.insar and not options.gnss:
        options.waveforms = True

    if options.waveforms:
        obs = grond_scenario.WaveformObservation(
            nstations=options.nstations,
            store_id=options.store_waveforms)
        scenario.add_observation(obs)

    if options.insar:
        obs = grond_scenario.InSARObservation(
            store_id=options.store_statics)
        scenario.add_observation(obs)

    if options.gnss:
        obs = grond_scenario.GNSSCampaignObservation(
            nstations=options.gnss_nstations,
            store_id=options.store_statics)
        scenario.add_observation(obs)

    if options.source == 'dc':
        problem = grond_scenario.DCSourceProblem(
            nevents=options.nevents)
    elif options.source == 'rectangular':
        problem = grond_scenario.RectangularSourceProblem(
            nevents=options.nevents)
    scenario.set_problem(problem)

    scenario.build(force=options.force, interactive=True)


def command_init(args):

    from . import cmd_init as init

    def setup(parser):
        parser.add_option(
            '--waveforms', dest='waveforms', action='store_true',
            help='add waveform configuration '
                 '(default)')
        parser.add_option(
            '--insar', dest='insar', action='store_true',
            help='add InSAR displacement scenes using kite containers'
                 '(https://pyrocko.org')
        parser.add_option(
            '--gnss', dest='gnss', action='store_true',
            help='add GNSS campaign configuration')
        parser.add_option(
            '--full', dest='full', action='store_true',
            help='create a full configuration, from targets above')
        parser.add_option(
            '--problem-cmt', dest='cmt', action='store_true',
            help='add a CMT source problem')
        parser.add_option(
            '--problem-rectangular', dest='rectangular',
            action='store_true',
            help='add a finite rectangular source problem')
        parser.add_option(
            '--force', dest='force', action='store_true',
            help='overwrite existing project folder')

    parser, options, args = cl_parse('init', args, setup)

    project_dir = None
    if len(args) == 1:
        project_dir = op.join(op.curdir, args[0])

    if not options.insar and not options.waveforms:
        options.waveforms = True

    project = init.GrondProject()

    if options.waveforms:
        project.add_waveforms()
        project.set_cmt_source()
    if options.insar:
        project.add_insar()
        project.set_rectangular_source()
    if options.gnss:
        project.add_gnss()
        project.set_rectangular_source()

    if options.full:
        project = init.GrondProject()

        project.add_waveforms()
        project.add_insar()
        project.add_gnss()
        project.set_rectangular_source()

    if options.cmt:
        project.set_cmt_source()
    if options.rectangular:
        project.set_rectangular_source()

    if len(args) == 1:
        project_dir = op.join(op.curdir, args[0])
        project.build(project_dir, options.force)
    else:
        sys.stdout.write(project.dump())


def command_events(args):
    def setup(parser):
        pass

    parser, options, args = cl_parse('events', args, setup)
    if len(args) != 1:
        help_and_die(parser, 'missing arguments')

    config_path = args[0]
    config = grond.read_config(config_path)

    for event_name in grond.get_event_names(config):
        print(event_name)


def command_check(args):
    def setup(parser):
        parser.add_option(
            '--target-ids', dest='target_string_ids', metavar='TARGET_IDS',
            help='process only selected targets. TARGET_IDS is a '
                 'comma-separated list of target IDs. Target IDs have the '
                 'form SUPERGROUP.GROUP.NETWORK.STATION.LOCATION.CHANNEL.')

        parser.add_option(
            '--waveforms', dest='show_waveforms', action='store_true',
            help='show raw, restituted, projected, and processed waveforms')

        parser.add_option(
            '--nrandom', dest='n_random_synthetics', metavar='N', type='int',
            default=10,
            help='set number of random synthetics to forward model (default: '
                 '10). If set to zero, create synthetics for the reference '
                 'solution.')

    parser, options, args = cl_parse('check', args, setup)
    if len(args) < 2:
        help_and_die(parser, 'missing arguments')

    config_path = args[0]
    event_names = args[1:]

    config = grond.read_config(config_path)

    target_string_ids = None
    if options.target_string_ids:
        target_string_ids = options.target_string_ids.split(',')

    grond.check(
        config,
        event_names=event_names,
        target_string_ids=target_string_ids,
        show_waveforms=options.show_waveforms,
        n_random_synthetics=options.n_random_synthetics)


def command_go(args):
    def setup(parser):
        parser.add_option(
            '--force', dest='force', action='store_true',
            help='overwrite existing run directory')
        parser.add_option(
            '--preserve', dest='preserve', action='store_true',
            help='preserve old rundir')
        parser.add_option(
            '--status', dest='status', default='state',
            type='choice',
            choices=['state', 'quiet'],
            help='status output selection (choices: state, quiet, default: '
                 'state)')
        parser.add_option(
            '--parallel', dest='nparallel', type='int', default=1,
            help='set number of events to process in parallel')

    parser, options, args = cl_parse('go', args, setup)

    if len(args) == 1:
        config_path = args[0]
        config = grond.read_config(config_path)

        event_names = config.get_event_names()
        if len(event_names) == 1:
            args.append(event_names[0])
        else:
            help_and_die(
                parser,
                'missing <eventnames>, candidates are:\n\n%s' % '\n'.join(
                    event_names))

    if len(args) < 2:
        help_and_die(parser, 'missing arguments')

    config_path = args[0]
    event_names = args[1:]

    config = grond.read_config(config_path)

    grond.go(
        config,
        event_names=event_names,
        force=options.force,
        preserve=options.preserve,
        status=options.status,
        nparallel=options.nparallel)


def command_forward(args):
    def setup(parser):
        pass

    parser, options, args = cl_parse('forward', args, setup)
    if len(args) < 1:
        help_and_die(parser, 'missing arguments')

    event_names = args[1:]

    if not event_names:
        help_and_die(parser, 'no event names given')

    run_path = args[0]
    grond.forward(
        run_path,
        event_names=event_names)


def command_harvest(args):
    def setup(parser):
        parser.add_option(
            '--force', dest='force', action='store_true',
            help='overwrite existing harvest directory')
        parser.add_option(
            '--neach', dest='neach', type=int, default=10,
            help='take NEACH best samples from each chain (default: 10)')
        parser.add_option(
            '--weed', dest='weed', type=int, default=0,
            help='weed out bootstrap samples with bad global performance. '
                 '0: no weeding (default), '
                 '1: only bootstrap chains where all NEACH best samples '
                 'global misfit is less than the global average misfit of all '
                 'NEACH best in all chains plus one standard deviation are '
                 'included in the harvest ensemble, '
                 '2: same as 1 but additionally individual samples are '
                 'removed if their global misfit is greater than the global '
                 'average misfit of all NEACH best in all chains, '
                 '3: harvesting is done on the global chain only, bootstrap '
                 'chains are excluded')

    parser, options, args = cl_parse('harvest', args, setup)
    if len(args) != 1:
        help_and_die(parser, 'no rundir')

    run_path, = args
    grond.harvest(
        run_path,
        force=options.force,
        nbest=options.neach,
        weed=options.weed)


def command_plot(args):

    import matplotlib
    matplotlib.use('Agg')

    from grond.environment import Environment

    def setup(parser):
        pass

    details = '''Available <plotnames> can be listed with

    `grond plot list ( <rundir> | <configfile> <eventname> )`
    `grond plot <plotname> ( <rundir> | <configfile> <eventname> )`
    `grond plot config ( <rundir> | <configfile> <eventname> )`
'''

    parser, options, args = cl_parse('plot', args, setup, details)

    if len(args) not in (2, 3):
        help_and_die(parser, 'two or three arguments required')

    env = Environment(*args[1:])
    from grond import plot
    if args[0] == 'list':
        print('Usage `grond plot <plot_name>'
              ' ( <rundir> | <configfile> <eventname> )`')
        print('Available plots for this rundir:')
        plot_names, plot_doc = zip(*[(pc.name, pc.__doc__)
                                     for pc in env.get_plot_classes()])
        plot_descs = [doc.split('\n')[0].strip() for doc in plot_doc]
        left_spaces = max([len(pn) for pn in plot_names])

        print()
        for name, desc in zip(plot_names, plot_descs):
            print('{name:<{ls}} - {desc}'.format(
                ls=left_spaces, name=name, desc=desc))

    elif args[0] == 'config':
        plot_config_collection = plot.get_plot_config_collection(env)
        print(plot_config_collection)

    elif args[0] == 'all':
        plots = plot.get_plot_names(env)
        plot.make_plots(env, plots)

    elif op.exists(args[0]):
        plots = plot.PlotConfigCollection.load(args[0])
        plot.make_plots(env, plots)

    else:
        plot_names = [name.strip() for name in args[0].split(',')]
        plot.make_plots(env, plot_names=plot_names)


def command_movie(args):

    import matplotlib
    matplotlib.use('Agg')

    def setup(parser):
        pass

    parser, options, args = cl_parse('movie', args, setup)

    if len(args) != 4:
        help_and_die(parser, 'four arguments required')

    run_path, xpar_name, ypar_name, movie_filename_template = args

    from grond import plot

    movie_filename = movie_filename_template % {
        'xpar': xpar_name,
        'ypar': ypar_name}

    try:
        plot.make_movie(run_path, xpar_name, ypar_name, movie_filename)

    except grond.GrondError as e:
        die(str(e))


def command_export(args):

    def setup(parser):
        parser.add_option(
            '--type', dest='type', metavar='TYPE',
            choices=('event', 'event-yaml', 'source', 'vector'),
            help='select type of objects to be exported. Choices: '
                 '"event" (default), "event-yaml", "source", "vector".')

        parser.add_option(
            '--parameters', dest='parameters', metavar='PLIST',
            help='select parameters to be exported. PLIST is a '
                 'comma-separated list where each entry has the form '
                 '"<parameter>[.<measure>]". Available measures: "best", '
                 '"mean", "std", "minimum", "percentile16", "median", '
                 '"percentile84", "maximum".')

        parser.add_option(
            '--output', dest='filename', metavar='FILE',
            help='write output to FILE')

    parser, options, args = cl_parse('export', args, setup)
    if len(args) < 2:
        help_and_die(parser, 'arguments required')

    what = args[0]

    dirnames = args[1:]

    what_choices = ('best', 'mean', 'ensemble', 'stats')

    if what not in what_choices:
        help_and_die(
            parser,
            'invalid choice: %s (choose from %s)' % (
                repr(what), ', '.join(repr(x) for x in what_choices)))

    if options.parameters:
        pnames = options.parameters.split(',')
    else:
        pnames = None

    try:
        grond.export(
            what,
            dirnames,
            filename=options.filename,
            type=options.type,
            pnames=pnames)

    except grond.GrondError as e:
        die(str(e))


def command_report(args):

    import matplotlib
    matplotlib.use('Agg')

    from grond.environment import Environment
    from grond.report import report, report_index

    def setup(parser):
        parser.add_option(
            '--index-only', dest='index_only', action='store_true',
            help='create index only')
        parser.add_option(
            '--open', dest='open', action='store_true',
            help='open webpage')
        parser.add_option(
            '--config', dest='config',
            help='configuration file to use')
        parser.add_option(
            '--update-without-plotting',
            dest='update_without_plotting',
            action='store_true',
            help='quick-and-dirty update parameter files without plotting')

    parser, options, args = cl_parse('report', args, setup)

    if options.config:
        conf = report.read_config(options.config)
    else:
        conf = None

    if options.index_only:
        report_index(conf)
        sys.exit(0)

    if len(args) < 1:
        help_and_die(parser, 'arguments required')

    if all(op.isdir(rundir) for rundir in args):
        rundirs = args
        try:
            for rundir in rundirs:
                env = Environment(rundir)
                report(
                    env, conf,
                    update_without_plotting=options.update_without_plotting)

        except grond.GrondError as e:
            die(str(e))

    else:
        if len(args) < 2:
            help_and_die(parser, 'arguments required')

        config_path = args[0]
        event_names = args[1:]

        try:
            env = Environment(config_path)
            for event_name in event_names:
                env.set_event_name(event_name)
                report(
                    env, conf,
                    update_without_plotting=options.update_without_plotting)

        except grond.GrondError as e:
            die(str(e))

    if options.open:
        import webbrowser
        webbrowser.open(op.join(conf.reports_base_path, 'index.html'))


def command_qc_polarization(args):

    def setup(parser):
        parser.add_option(
            '--time-factor-pre', dest='time_factor_pre', type=float,
            metavar='NUMBER',
            default=0.5,
            help='set duration to extract before synthetic P phase arrival, '
                 'relative to 1/fmin. fmin is taken from the selected target '
                 'group in the config file (default=%default)')
        parser.add_option(
            '--time-factor-post', dest='time_factor_post', type=float,
            metavar='NUMBER',
            default=0.5,
            help='set duration to extract after synthetic P phase arrival, '
                 'relative to 1/fmin. fmin is taken from the selected target '
                 'group in the config file (default=%default)')
        parser.add_option(
            '--distance-min', dest='distance_min', type=float,
            metavar='NUMBER',
            help='minimum event-station distance [m]')
        parser.add_option(
            '--distance-max', dest='distance_max', type=float,
            metavar='NUMBER',
            help='maximum event-station distance [m]')
        parser.add_option(
            '--depth-min', dest='depth_min', type=float,
            metavar='NUMBER',
            help='minimum station depth [m]')
        parser.add_option(
            '--depth-max', dest='depth_max', type=float,
            metavar='NUMBER',
            help='maximum station depth [m]')
        parser.add_option(
            '--picks', dest='picks_filename',
            metavar='FILENAME',
            help='add file with P picks in Snuffler marker format')
        parser.add_option(
            '--save', dest='output_filename',
            metavar='FILENAME.FORMAT',
            help='save output to file FILENAME.FORMAT')
        parser.add_option(
            '--dpi', dest='output_dpi', type=float, default=120.,
            metavar='NUMBER',
            help='DPI setting for raster formats (default=120)')

    parser, options, args = cl_parse('qc-polarization', args, setup)
    if len(args) != 3:
        help_and_die(parser, 'missing arguments')

    if options.output_filename:
        import matplotlib
        matplotlib.use('Agg')

    import grond.qc

    config_path, event_name, target_group_path = args

    config = grond.read_config(config_path)

    ds = config.get_dataset(event_name)

    engine = config.engine_config.get_engine()

    nsl_to_time = None
    if options.picks_filename:
        markers = marker.load_markers(options.picks_filename)
        marker.associate_phases_to_events(markers)

        nsl_to_time = {}
        for m in markers:
            if isinstance(m, marker.PhaseMarker):
                ev = m.get_event()
                if ev is not None and ev.name == event_name:
                    nsl_to_time[m.one_nslc()[:3]] = m.tmin

        if not nsl_to_time:
            help_and_die(
                parser,
                'no markers associated with event "%s" found in file "%s"' % (
                    event_name, options.picks_filename))

    target_group_paths_avail = []
    for target_group in config.target_groups:
        name = target_group.path
        if name == target_group_path:
            imc = target_group.misfit_config
            fmin = imc.fmin
            fmax = imc.fmax
            ffactor = imc.ffactor

            store = engine.get_store(target_group.store_id)
            timing = '{cake:P|cake:p|cake:P\\|cake:p\\}'

            grond.qc.polarization(
                ds, store, timing, fmin=fmin, fmax=fmax, ffactor=ffactor,
                time_factor_pre=options.time_factor_pre,
                time_factor_post=options.time_factor_post,
                distance_min=options.distance_min,
                distance_max=options.distance_max,
                depth_min=options.depth_min,
                depth_max=options.depth_max,
                nsl_to_time=nsl_to_time,
                output_filename=options.output_filename,
                output_dpi=options.output_dpi)

            return

        target_group_paths_avail.append(name)

        die('no target group with path "%s" found. Available: %s' % (
            target_group_path, ', '.join(target_group_paths_avail)))


def command_upgrade_config(args):
    def setup(parser):
        parser.add_option(
            '--diff', dest='diff', action='store_true',
            help='create diff between normalized old and new versions')

    parser, options, args = cl_parse('upgrade-config', args, setup)
    if len(args) != 1:
        help_and_die(parser, 'missing argument <configfile>')

    from grond import upgrade
    upgrade.upgrade_config_file(args[0], diff=options.diff)


def command_diff(args):
    def setup(parser):
        pass

    parser, options, args = cl_parse('diff', args, setup)
    if len(args) != 2:
        help_and_die(parser, 'requires exactly two arguments')

    from grond.config import diff_configs
    diff_configs(*args)


if __name__ == '__main__':
    main()
