#!/usr/bin/env python
FRAME = 0
PARTICLE = 1

# Reads a frame of data output from mpm_2d
# f_in - already opened file to read from
# f_cb - callback to perform when a frame is read. Takes the frame dictionary
#           and is presumed to return a (possibly modified) frame dictionary.
# p_cb - callback to perform when a particle is read. Takes the
#           particle dictionary and is presumed to return a (possibly modified)
#           particle dictionary.
def read_frame_cb(f_in, f_cb, p_cb):
    frame = None
    state = FRAME
    i = 0
    num_particles = 0
    for line in f_in:
        if (state == FRAME):
            tok = line.split(' ')
            num_particles = int(tok[2])
            state = PARTICLE
            frame = {'time':float(tok[1]), 'particles':[0]*num_particles}
            i = 0
        elif (state == PARTICLE):
            tok = line.split(' ')
            tokf = map(float, tok)
            particle = {'m':tokf[0], 'v':tokf[1], 'x':tokf[2], 'y':tokf[3],
                        'x_t':tokf[4], 'y_t':tokf[5], 'sxx':tokf[6],
                        'sxy':tokf[7], 'syy':tokf[8], 'ux':tokf[9],
                        'uy':tokf[10], 'gammap':tokf[11]}
            if (p_cb is None):
                frame['particles'][i] = particle
            else:
                frame['particles'][i] = p_cb(particle)
            i = i + 1
            if (i >= num_particles):
                if (f_cb is not None):
                    frame = f_cb(frame)
                break
        else:
            print 'Entered unknown state.'
            return None
    return frame


def read_frame(f_in):
    return read_frame_cb(f_in, None, None)
