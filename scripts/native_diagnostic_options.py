"""Optional native diagnostics shared by discovery and targeted batch CLIs."""


def add_arguments(parser):
    parser.add_argument('--performance-profile', choices=('on', 'off'), default='off',
                        help='opt-in native input/producer timing diagnostics (requires native 0.8.5+)')
    parser.add_argument('--quality-diagnostics', choices=('off', 'phred33'), default='off',
                        help='explicit FASTQ quality encoding declaration; observation only, never a filter')
    parser.add_argument('--quality-max-record-bases', type=int,
                        help='retain complete quality for reads up to this length (1..1000000); requires phred33')


def native_diagnostic_options(args):
    quality = getattr(args, 'quality_diagnostics', 'off')
    profile = getattr(args, 'performance_profile', 'off')
    limit = getattr(args, 'quality_max_record_bases', None)
    if quality not in ('off', 'phred33') or profile not in ('on', 'off'):
        raise ValueError('Invalid native diagnostic option')
    if limit is not None and (quality != 'phred33' or type(limit) is not int or not 1 <= limit <= 1_000_000):
        raise ValueError('quality-max-record-bases requires phred33 and must be 1..1000000')
    result = []
    if profile == 'on':
        result += ['--performance-profile', 'on']
    if quality == 'phred33':
        result += ['--quality-diagnostics', 'phred33']
    if limit is not None:
        result += ['--quality-max-record-bases', str(limit)]
    return result


def validate_diagnostic_output(result, args):
    if getattr(args, 'performance_profile', 'off') == 'on' and not isinstance(result.get('performance_profile'), dict):
        raise ValueError('Requested native performance profile is missing')
    if getattr(args, 'quality_diagnostics', 'off') == 'phred33':
        from quality_observation import describe_quality
        context = describe_quality(result, None)
        if context['diagnostics'] is None:
            raise ValueError('Requested native quality diagnostics are missing')
        expected = getattr(args, 'quality_max_record_bases', None) or 100_000
        if context['diagnostics'].get('max_record_bases') != expected:
            raise ValueError('Native quality retention limit differs from requested configuration')
