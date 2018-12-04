#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# Copyright © 2018 rgregoir <rgregoir@laurier>
#
# Distributed under terms of the MIT license.

import math
import os
import os.path
import pyBigWig
from sty import fg, bg, ef, rs
import sys
from subprocess import call

dirname = os.path.dirname(os.path.abspath(__file__))

inputFiles = [
    os.path.join(dirname, './in/25574.Blueprint.ERS487305.RNA-Seq.signal_reverse.bigWig'),
    os.path.join(dirname, './in/25575.Blueprint.ERS487305.RNA-Seq.signal_forward.bigWig'),
    os.path.join(dirname, './in/25582.Blueprint.ERS487306.RNA-Seq.signal_reverse.bigWig'),
]

outputFile = os.path.join(dirname, 'output.bw')
deviationFile = os.path.join(dirname, 'standard-deviation.bw')

outputChrom = 'chr1'
outputPosition = [0, 500000]
outputRange = [0, 1000]
deviationDefault = 0.0

position = '%s:%d-%d' % (outputChrom, outputPosition[0], outputPosition[1])
executable = os.path.normpath(os.path.join(dirname, '../bigWigMergePlus'))


def main():
    tests = [
        test_deviation_with_default,
        test_deviation_without_default,
    ]

    print('')
    for test in tests:
        print_title('Running', test.__name__)
        test()


#########
# Tests #
#########

def test_deviation_without_default():
    command = [
        executable,
        '-compress',
        '-range=%d-%d' % (outputRange[0], outputRange[1]),
        '-position=%s' % position,
        '-deviation=%s' % deviationFile,
    ] + inputFiles + [
        outputFile
    ]
    print('\nCommand: \x1b[1m%s\x1b[0m\n' % ' '.join(command))
    call(command)

    tracks = [pyBigWig.open(f) for f in inputFiles]

    dev = pyBigWig.open(deviationFile)
    out = pyBigWig.open(outputFile)

    summaries = [t.header() for t in tracks]

    print('')
    print('Checking standard-deviation on position chr1:%d-%d' % tuple(outputPosition))
    print('outputRange = [%d, %d]' % tuple(outputRange))
    print('%d tracks' % len(tracks))
    #  print('t1: minVal = %d, maxVal = %d' % (s1['minVal'], s1['maxVal']))
    #  print('t2: minVal = %d, maxVal = %d' % (s2['minVal'], s2['maxVal']))
    #  print('t3: minVal = %d, maxVal = %d' % (s3['minVal'], s3['maxVal']))
    print('')

    for i in range(outputPosition[1]):
        all_values = [get_value_transformed(tracks[n], summaries[n], i) for n in range(len(tracks))]
        vdev = get_value_at(dev, i)
        vout = get_value_at(out, i)

        if math.isnan(vdev):
            vdev = 0
        if math.isnan(vout):
            vout = 0

        values = [val for val in all_values if not math.isnan(val)]

        if len(values) <= 1:
            if vdev == 0:
                continue
            print('Error: vdev should be 0')
            print({'i': i, 'vdev': vdev, 'all_values': all_values})
            sys.exit(1)

        n = len(values)
        mean = sum(values) / n
        stdDeviation = math.sqrt(
            sum(map(lambda x: math.pow(x - mean, 2), values)) / (n - 1)
        )

        merged = sum([v / len(inputFiles) for v in values])

        if not isclose(merged, vout, 1e-02):
            print_error_message('merge value mismatch')
            print({'i': i, 'all_values': all_values})
            print('')
            print('actual:   %.16f' % vout)
            print('expected: %.16f' % merged)
            print('')
            sys.exit(1)

        #  print({'i': i, 'v1': v1, 'v2': v2, 'v3': v3, 'vdev': vdev, 'stdDeviation': stdDeviation})
        #  print('%d: %f %f' % (i, round_to_n(vdev, 3), round_to_n(stdDeviation, 3)))
        if not isclose(vdev, stdDeviation, 1e-01):
            print_error_message('standard-deviation value mismatch')
            print({'i': i, 'vdev': vdev, 'stdDeviation': stdDeviation})
            print({'all_values': all_values})
            print({'raw_values': raw_values})
            print({'____values': values})

            print('')
            print('actual:   %.4f' % round_to_n(vdev, 4))
            print('expected: %.4f' % round_to_n(stdDeviation, 4))
            print('')
            sys.exit(1)

def test_deviation_with_default():
    command = [
        executable,
        '-compress',
        '-range=%d-%d' % (outputRange[0], outputRange[1]),
        '-position=%s' % position,
        '-deviation=%s' % deviationFile,
        '-deviationDefault=%s' % deviationDefault,
    ] + inputFiles + [
        outputFile
    ]
    print('\nCommand: \x1b[1m%s\x1b[0m\n' % ' '.join(command))
    call(command)

    tracks = [pyBigWig.open(f) for f in inputFiles]

    dev = pyBigWig.open(deviationFile)
    out = pyBigWig.open(outputFile)

    summaries = [t.header() for t in tracks]

    print('')
    print('Checking standard-deviation on position chr1:%d-%d' % tuple(outputPosition))
    print('outputRange = [%d, %d]' % tuple(outputRange))
    print('%d tracks' % len(tracks))
    #  print('t1: minVal = %d, maxVal = %d' % (s1['minVal'], s1['maxVal']))
    #  print('t2: minVal = %d, maxVal = %d' % (s2['minVal'], s2['maxVal']))
    #  print('t3: minVal = %d, maxVal = %d' % (s3['minVal'], s3['maxVal']))
    print('')

    for i in range(outputPosition[1]):
        raw_values = [get_value_at(tracks[n], i) for n in range(len(tracks))]
        all_values = [get_value_transformed(tracks[n], summaries[n], i) for n in range(len(tracks))]
        vdev = get_value_at(dev, i)
        vout = get_value_at(out, i)

        if math.isnan(vdev):
            vdev = 0
        if math.isnan(vout):
            vout = 0

        values = [val if not math.isnan(val) else deviationDefault for val in all_values]

        if len(values) <= 1:
            if vdev == 0:
                continue
            print('Error: vdev should be 0')
            print({'i': i, 'vdev': vdev, 'all_values': all_values})
            sys.exit(1)

        n = len(values)
        mean = sum(values) / n
        stdDeviation = math.sqrt(
            sum(map(lambda x: math.pow(x - mean, 2), values)) / (n - 1)
        )

        merged = sum([v / len(inputFiles) for v in values])

        if not isclose(merged, vout, 1e-02):
            print_error_message('merge value mismatch')
            print({'i': i})
            print({'raw_values': raw_values})
            print({'all_values': all_values})
            print({'____values': values})
            print('')
            print('actual:   %.16f' % vout)
            print('expected: %.16f' % merged)
            print('')
            sys.exit(1)

        #  print({'i': i, 'v1': v1, 'v2': v2, 'v3': v3, 'vdev': vdev, 'stdDeviation': stdDeviation})
        #  print('%d: %f %f' % (i, round_to_n(vdev, 3), round_to_n(stdDeviation, 3)))
        if not isclose(vdev, stdDeviation, 1e-01):
            print_error_message('standard-deviation value mismatch')
            print({'i': i, 'vdev': vdev, 'stdDeviation': stdDeviation})
            print({'all_values': all_values})
            print({'raw_values': raw_values})
            print({'____values': values})

            print('')
            print('actual:   %.4f' % round_to_n(vdev, 4))
            print('expected: %.4f' % round_to_n(stdDeviation, 4))
            print('')
            sys.exit(1)


###########
# Helpers #
###########

def print_title(title, text=''):
    print(fg.black + bg.white + ef.bold + ' ' + title + ' ' + rs.all + ef.bold + ' ' + text + rs.all)

def print_error_message(message):
    print(fg.white + bg.red + ef.bold + ' Error ' + rs.all + ef.bold + ' ' + message + rs.all)

def get_value_at(track, position):
    val = track.values(outputChrom, position, position + 1)[0]
    return val

def get_value_transformed(track, summary, position):
    val = get_value_at(track, position)
    if math.isnan(val):
        return val
    val = \
        ((((val - summary['minVal']) / (summary['maxVal'] - summary['minVal'])) * (outputRange[1] - outputRange[0])) + outputRange[0])
    return val

def round_to_n(x, n):
    if x == 0:
        return x
    try:
        return round(x, -int(math.floor(math.log10(x))) + (n - 1))
    except Exception as e:
        print(e)
        print(x)
        print(n)
        sys.exit(1)

def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


if __name__ == '__main__':
    main()
