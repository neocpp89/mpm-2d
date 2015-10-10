#!/bin/bash
parallel --slf machines.slf --progress "cd $(pwd); ./frame_to_csv.py {}/frame_particle_data.txt x,y,x_t,y_t,m 0:60 {}/fpd" ::: drop/*.output
parallel --slf machines.slf --progress "cd $(pwd); rm {}/com; for x in {0..59}; do ./com_finder_simple.py {}/fpd_\$x.csv >> {}/com; done;" ::: drop/*.output
parallel "paste -d ' ' drop/{}_*.output/com > {}_com" ::: `seq 0 0.5 3`
