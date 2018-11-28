#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2018 rgregoir <rgregoir@laurier>
#
# Distributed under terms of the MIT license.

import math
import os
import os.path
import pyBigWig
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

outputPosition = [0, 500000]
outputRange = [0, 1000]
position = 'chr1:%d-%d' % (outputPosition[0], outputPosition[1])
executable = os.path.normpath(os.path.join(dirname, '../bigWigMergePlus'))


def main():
    command = [
        executable,
        '-compress',
        '-range=%d-%d' % (outputRange[0], outputRange[1]),
        '-position=%s' % position,
        '-deviation=%s' % deviationFile,
        inputFiles[0],
        inputFiles[1],
        inputFiles[2],
        outputFile
    ]
    print '\nRunning: \x1b[1m%s\x1b[0m\n' % ' '.join(command)
    call(command)

    t1 = pyBigWig.open(inputFiles[0])
    t2 = pyBigWig.open(inputFiles[1])
    t3 = pyBigWig.open(inputFiles[2])

    dev = pyBigWig.open(deviationFile)

    s1 = t1.header()
    s2 = t2.header()
    s3 = t3.header()

    print ''
    print 'Checking standard-deviation on position chr1:%d-%d' % tuple(outputPosition)
    print 'outputRange = [%d, %d]' % tuple(outputRange)
    print 't1: minVal = %d, maxVal = %d' % (s1['minVal'], s1['maxVal'])
    print 't2: minVal = %d, maxVal = %d' % (s2['minVal'], s2['maxVal'])
    print 't3: minVal = %d, maxVal = %d' % (s3['minVal'], s3['maxVal'])
    print ''

    for i in range(outputPosition[1]):
        v1 = get_value_transformed(t1, s1, i)
        v2 = get_value_transformed(t2, s2, i)
        v3 = get_value_transformed(t3, s3, i)
        vdev = get_value_at(dev, i)

        if math.isnan(vdev):
            vdev = 0

        values = []

        if not math.isnan(v1):
            values.append(v1)
        if not math.isnan(v2):
            values.append(v2)
        if not math.isnan(v3):
            values.append(v3)

        if len(values) <= 1:
            if vdev == 0:
                continue
            print 'Error: vdev should be 0'
            print {'i': i, 'v1': v1, 'v2': v2, 'v3': v3, 'vdev': vdev}
            sys.exit(1)

        n = len(values)
        mean = sum(values) / n
        stdDeviation = math.sqrt(
            sum(map(lambda x: math.pow(x - mean, 2), values)) / (n - 1)
        )

        #  print {'i': i, 'v1': v1, 'v2': v2, 'v3': v3, 'vdev': vdev, 'stdDeviation': stdDeviation}
        #  print '%d: %f %f' % (i, round_to_n(vdev, 3), round_to_n(stdDeviation, 3))
        if not isclose(vdev, stdDeviation):
            print 'Error: value mismatch'
            print {'i': i, 'v1': v1, 'v2': v2, 'v3': v3, 'vdev': vdev, 'stdDeviation': stdDeviation}
            print ''
            print 'Rounded:'
            print {'vdev': round_to_n(vdev, 4), 'stdDeviation': round_to_n(stdDeviation, 4)}
            print get_value_at(t1, i), get_value_transformed(t1, s1, i)
            print get_value_at(t2, i), get_value_transformed(t2, s2, i)
            print get_value_at(t3, i), get_value_transformed(t3, s3, i)
            sys.exit(1)


def get_value_at(track, position):
    val = track.values('chr1', position, position + 1)[0]
    return val

def get_value_transformed(track, summary, position):
    val = get_value_at(track, position)
    if math.isnan(val):
        return val
    val = \
        ((((val - summary['minVal']) / (summary['maxVal'] - summary['minVal'])) * (outputRange[1] - outputRange[0])) + outputRange[0]) \
        / len(inputFiles)
    return val

def round_to_n(x, n):
    if x == 0:
        return x
    try:
        return round(x, -int(math.floor(math.log10(x))) + (n - 1))
    except Exception as e:
        print e
        print x
        print n
        sys.exit(1)

def isclose(a, b, rel_tol=1e-02, abs_tol=0.0):
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


if __name__ == '__main__':
    main()
