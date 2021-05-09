'''
Marc Coll 05/2021

Add support for passing command line arguments. In order to use this feature you
just need to add the following two lines to your script:

    from mudpy import parsecmdargs
    parsecmdargs.parse_cmdline_arguments(globals())

They must be placed right after all the global configuration variables (ncpus,
model_name, fault_name, NFFT...) have been declared.
'''

import argparse
import sys
import re
from obspy.core import UTCDateTime


def parse_cmdline_arguments(gobal_vars):
    if (len(sys.argv) <= 1):
        return

    '''
    Add here the names of the variables you want to be able to modify from the command line,
    and then modify the parse function accordingly
    '''
    params = [
        'init',
        'make_ruptures',
        'make_GFs',
        'make_synthetics',
        'make_waveforms',
        'load_distances',
        'G_from_file',
        'ncpus',
        'use_gpu',
        'model_name',
        'fault_name',
        'slab_name',
        'mesh_name',
        'distances_name',
        'UTM_zone',
        'scaling_law',
        'dynamic_GFlist',
        'dist_threshold',
        'Nrealizations',
        'max_slip',
        'hurst',
        'Ldip',
        'Lstrike',
        'lognormal',
        'slip_standard_deviation',
        'num_modes',
        'rake',
        'force_magnitude',
        'force_area',
        'no_random',
        'time_epi',
        'hypocenter',
        'force_hypocenter',
        'mean_slip',
        'center_subfault',
        'use_hypo_fraction',
        'source_time_function',
        'rise_time_depths',
        'shear_wave_fraction',
        'GF_list',
        'G_name',
        'NFFT',
        'dt',
        'dk',
        'pmin',
        'pmax',
        'kmax',
        'custom_stf'
    ]
    dict_params = {}
    for p in params:
        #Get values from existing variables using their names
        dict_params[p] = gobal_vars[p]
    ok,actions=parse(dict_params)
    if (ok == False):
        exit(1)
    if (len(actions) > 0):
        print("Actions: {}\nConfiguration:".format(actions))
    max_length = max(len(p) for p in params)
    for p in params:
        padding = ''.ljust(max_length - len(p) + 1)
        print("\t{}{}= {}".format(p, padding, dict_params[p]))
        #Set values to existing variables using their names
        gobal_vars[p] = dict_params[p]
    print("\n")


def parse(params):
    prog = sys.argv[0]
    description = "MudPy's command line arguments"
    epilog = '''Examples:
    {} init
    {} make_ruptures -load_distances=T
    {} make_waveforms -ncpus=16 -use_gpu=F -nfft=128 -source_time_function=dreger
    {} make_ruptures make_g_files make_waveforms -load_distances=F -ncpus=16 -use_gpu=T
    '''.format(prog, prog, prog, prog)

    bool_choices=['T', 'F', 't', 'f', '0', '1']
    ldip_lstrike_choices=['auto', 'MB2002', 'MH2019']

    # Parse command line parameters
    parser = argparse.ArgumentParser(description=description, epilog=epilog,
                                     formatter_class=argparse.RawTextHelpFormatter)

    # What do you want to do?
    parser.add_argument('actions', nargs='*',
                        metavar='[{init,make_ruptures,make_g_files,make_waveforms}]',
                        choices=['init', 'make_ruptures', 'make_g_files', 'make_waveforms'],
                        help='What do you want to do?\n'
                             ' * init: initialize project folders\n'
                             ' * make_ruptures: generate rupture models\n'
                             ' * make_g_files: prepare waveforms and synthetics\n'
                             ' * make_waveforms: synthesize the waveforms')

    parser.add_argument('-load_distances', default=str(params['load_distances']), metavar='[TF]',
                        choices=bool_choices,
                        help='load the distances instead of calculating them')
    parser.add_argument('-g_from_file', default=str(params['G_from_file']), metavar='[TF]',
                        choices=bool_choices,
                        help='read the LARGE matrix (must run make_waveforms first)')

    # Run-time parameters
    parser.add_argument('-ncpus', type=int, default=params['ncpus'], metavar='N',
                        help='number of CPUs to use')
    parser.add_argument('-use_gpu', default=str(params['use_gpu']), metavar='[TF]',
                        choices=bool_choices,
                        help='use the GPU, if available')
    parser.add_argument('-model_name', default=params['model_name'], metavar='NAME',
                        help='velocity model')
    parser.add_argument('-fault_name', default=params['fault_name'], metavar='NAME',
                        help='fault geometry')
    parser.add_argument('-slab_name', default=params['slab_name'], metavar='NAME',
                        help='slab 1.0 Ascii file (only used for 3D fault)')
    parser.add_argument('-mesh_name', default=params['mesh_name'], metavar='NAME',
                        help='GMSH output file (only used for 3D fault)')
    parser.add_argument('-distances_name', default=params['distances_name'], metavar='NAME',
                        help='name of distance matrix')
    parser.add_argument('-utm_zone', default=params['UTM_zone'], metavar='ZONE',
                        help='as defined in the Universal Transverse Mercator coordinate system')
    parser.add_argument('-scaling_law', default=params['scaling_law'], metavar='[TSN]',
                        choices=['T', 'S', 'N'],
                        help='T for thrust, S for strike-slip, N for normal')

    # Dynamic GFlist options
    parser.add_argument('-dynamic_gflist', default=str(params['dynamic_GFlist']), metavar='[TF]',
                        choices=bool_choices,
                        help='use dynamic GFlist')
    parser.add_argument('-dist_threshold', type=float, default=params['dist_threshold'], metavar='N',
                        help='(degree) station to the closest subfault must be closer to this distance')

    # Slip parameters
    parser.add_argument('-nrealizations', type=int, default=params['Nrealizations'], metavar='N',
                        help='fake ruptures to generate per magnitude bin. Let Nrealizations %% ncpus=0')
    parser.add_argument('-max_slip', type=int, default=params['max_slip'], metavar='N',
                        help='maximum slip (m) allowed in the model')

    # Correlation function parameters
    parser.add_argument('-hurst', type=float, default=params['hurst'], metavar='N',
                        help='0.4~0.7 is reasonable')
    parser.add_argument('-ldip', default=params['Ldip'],
                        choices=ldip_lstrike_choices,
                        help='Correlation length scaling')
    parser.add_argument('-lstrike', default=params['Lstrike'],
                        choices=ldip_lstrike_choices,
                        help='Lstrike value (correlation function parameter)')
    parser.add_argument('-lognormal', default=str(params['lognormal']), metavar='[TF]',
                        choices=bool_choices,
                        help='')
    parser.add_argument('-slip_standard_deviation', type=float, default=params['slip_standard_deviation'], metavar='N',
                        help='')
    parser.add_argument('-num_modes', type=int, default=params['num_modes'], metavar='N',
                        help='modes in K-L expantion (max#= number of subfaults)')
    parser.add_argument('-rake', type=float, default=params['rake'], metavar='N',
                        help='')

    # Rupture parameters
    parser.add_argument('-force_magnitude', default=str(params['force_magnitude']), metavar='[TF]',
                        choices=bool_choices,
                        help='make the magnitudes EXACTLY the value in target_Mw')
    parser.add_argument('-force_area', default=str(params['force_area']), metavar='[TF]',
                        choices=bool_choices,
                        help='forces using the entire fault area defined by the .fault file')
    parser.add_argument('-no_random', default=str(params['no_random']), metavar='[TF]',
                        choices=bool_choices,
                        help='if true uses median length/width if false draws from prob. distribution')
    parser.add_argument('-time_epi', default=params['time_epi'].isoformat(), metavar='TIME',
                        help='defines the hypocentral time')
    parser.add_argument('-hypocenter', default=','.join(str(s) for s in params['hypocenter']), metavar='N,N,N',
                        help='defines the specific hypocenter location if force_hypocenter=True')
    parser.add_argument('-force_hypocenter', default=str(params['force_hypocenter']), metavar='[TF]',
                        choices=bool_choices,
                        help='forces hypocenter to occur at specified location as opposed to random')
    parser.add_argument('-mean_slip', default=params['mean_slip'], metavar='NAME',
                        help='provide path to file name of .rupt to be used as mean slip pattern')
    parser.add_argument('-center_subfault', type=int, default=params['center_subfault'], metavar='N',
                        help='use this subfault as center for defining rupt area')
    parser.add_argument('-use_hypo_fraction', default=str(params['use_hypo_fraction']), metavar='[TF]',
                        choices=bool_choices,
                        help='if true use hypocenter PDF positions from Melgar & Hayes 2019, if false then selects at random')

    # Kinematic parameters
    parser.add_argument('-source_time_function', default=params['source_time_function'],
                        choices=['triangle', 'cosine', 'dreger'],
                        help='calculation method used in source time function')
    parser.add_argument('-rise_time_depths', default=','.join(str(s) for s in params['rise_time_depths']), metavar='N,N',
                        help='transition depths for rise time scaling')
    parser.add_argument('-shear_wave_fraction', type=float, default=params['shear_wave_fraction'], metavar='N',
                        help='fraction of shear wave speed to use as mean rupture velocity')

    # Station information (only used when syntehsizing waveforms)
    parser.add_argument('-gf_list', default=params['GF_list'], metavar='NAME',
                        help='name of the GF list file')
    parser.add_argument('-g_name', default=params['G_name'], metavar='NAME',
                        help='name of the GF folder')

    # Displacement and velocity waveform parameters
    parser.add_argument('-nfft', type=int, default=params['NFFT'], metavar='N',
                        help='NFFT (displacement and velocity waveform parameter)')
    parser.add_argument('-dt', type=float, default=params['dt'], metavar='N',
                        help='dt (displacement and velocity waveform parameter)')

    # fk-parameters
    parser.add_argument('-dk', type=float, default=params['dk'], metavar='N',
                        help='dk (fk parameter)')
    parser.add_argument('-pmin', type=int, default=params['pmin'], metavar='N',
                        help='pmin (fk parameter)')
    parser.add_argument('-pmax', type=int, default=params['pmax'], metavar='N',
                        help='pmax (fk parameter)')
    parser.add_argument('-kmax', type=int, default=params['kmax'], metavar='N',
                        help='kmax (fk parameter)')
    parser.add_argument('-custom_stf', default=params['custom_stf'], metavar='STF',
                        help='custom_stf (fk parameter)')


    # Default initialization of actions
    params['init'] = 0
    params['make_ruptures'] = 0
    params['make_GFs'] = 0
    params['make_synthetics'] = 0
    params['make_waveforms'] = 0


    # Get parameter values
    args = parser.parse_args()
    if ('init' in args.actions):
        params['init'] = 1
    if ('make_ruptures' in args.actions):
        params['make_ruptures'] = 1
        params['load_distances'] = get_intbool(args.load_distances)
    if ('make_g_files' in args.actions):
        params['make_GFs'] = 1
        params['make_synthetics'] = 1
    if ('make_waveforms' in args.actions):
        params['make_waveforms'] = 1
        params['G_from_file'] = get_intbool(args.g_from_file)

    params['ncpus'] = args.ncpus
    params['use_gpu'] = get_bool(args.use_gpu)
    params['model_name'] = args.model_name
    params['fault_name'] = args.fault_name
    params['slab_name'] = args.slab_name
    params['mesh_name'] = args.mesh_name
    params['distances_name'] = args.distances_name
    params['UTM_zone'] = args.utm_zone
    params['scaling_law'] = args.scaling_law
    params['dynamic_GFlist'] = get_bool(args.dynamic_gflist)
    params['dist_threshold'] = args.dist_threshold
    params['Nrealizations'] = args.nrealizations
    params['max_slip'] = args.max_slip
    params['hurst'] = args.hurst
    params['Ldip'] = args.ldip
    params['Lstrike'] = args.lstrike
    params['lognormal'] = get_bool(args.lognormal)
    params['slip_standard_deviation'] = args.slip_standard_deviation
    params['num_modes'] = args.num_modes
    params['rake'] = args.rake
    params['force_magnitude'] = args.force_magnitude
    params['force_area'] = args.force_area
    params['no_random'] = args.no_random
    params['time_epi'] = UTCDateTime(args.time_epi)
    params['hypocenter'] = get_list(args.hypocenter, float)
    params['force_hypocenter'] = args.force_hypocenter
    params['mean_slip'] = args.mean_slip
    params['center_subfault'] = args.center_subfault
    params['use_hypo_fraction'] = args.use_hypo_fraction
    params['source_time_function'] = args.source_time_function
    params['rise_time_depths'] = get_list(args.rise_time_depths, int)
    params['shear_wave_fraction'] = args.shear_wave_fraction
    params['GF_list'] = args.gf_list
    params['G_name'] = args.g_name
    params['NFFT'] = args.nfft
    params['dt'] = args.dt
    params['dk'] = args.dk
    params['pmin'] = args.pmin
    params['pmax'] = args.pmax
    params['kmax'] = args.kmax
    params['custom_stf'] = args.custom_stf

    # Check if some of the parameters make sense
    if (args.ncpus < 1):
        print('Invalid number of CPUs')
        return False,None

    regex_timeepi = re.compile('^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}(.[0-9]+Z?)?$')
    if (not regex_timeepi.match(args.time_epi)):
        print('Invalid hypocentral time')
        return False,None

    if (len(params['hypocenter']) != 3):
        print('Invalid number of hypocenter coordinates')
        return False,None

    if (len(params['rise_time_depths']) != 2):
        print('Invalid number of rise time depths')
        return False,None

    regex_utmzone = re.compile('^(A|B|Y|Z|([1-5][0-9]|60|0[1-9])[C-X])$')
    if (not regex_utmzone.match(args.utm_zone)):
        print('Invalid UTM zone')
        return False,None

    return True,args.actions


def get_list(value, dtype):
    return list(map(dtype, list(value.split(','))))


def get_bool(value):
    v = value.upper()
    return (v == 'T' or value == '1')


def get_intbool(value):
    if get_bool(value):
        return 1
    else:
        return 0
