#!/usr/bin/env python
import sys
import csv

import framereader

def frame_callback(frame, pvars):
    # only care about active mass, delete other keys
    frame['x_t'] = [p['x_t'] for p in frame['particles'] if (p['active'] != 0 and p['y'] > 0.95)][-1]
    frame['particles'] = None
    return frame

if (len(sys.argv) <= 1):
    print 'usage: FRAME_DATA [outfile]'
    print '       outfile consists of [frameID, frametime, framemass] rows in csv format.'
    exit(127)

infile = sys.argv[1]
if (len(sys.argv) > 2):
    outfile = sys.argv[2]
else:
    outfile = ".".join(infile.split('.')[:-1]) + '_mass.csv'

print 'Using input file:', infile
print 'Creating output file:', outfile

i = 0
with open(infile, 'r') as f_in:
    with open(outfile, 'w') as f_csv:
        wcsv = csv.writer(f_csv)
        row = ['Frame', 't', 'x_t']
        wcsv.writerow(row)
        while True:
            frame = framereader.read_frame_cb(f_in,
                        lambda x: frame_callback(x, 'x_t'),
                        None)
            if (frame is None):
                break
            else:
                print "Read Frame:", i
                row = [i, frame['time'], frame['x_t']]
                wcsv.writerow(row)
                i = i + 1

print "Bye."

