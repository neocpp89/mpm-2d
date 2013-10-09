/**
    \file main.c
    \author Sachith Dunatunga
    \date 02.06.12

    mpm_2d -- An implementation of the Material Point Method in 2D.
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <getopt.h>
#include <pthread.h>
#include <confuse.h>

#include "material.h"
#include "particle.h"
#include "point.h"
#include "process.h"
#include "reader.h"
#include "writer.h"
#include "visualization.h"

#define dispg(x) printf(#x " = %g\n", x)
#define dispd(x) printf(#x " = %d\n", x)

#define SAMPLE_HZ 60.0f

static struct s_state {
    char *outputdir;

    char *h5file;
    char *fpfile;
    char *fefile;

    char *statefile;

    double tmax;
    int restart;        /* is this a restart? */
} g_state;

static const char* g_optstring = "o:r:";
char *outputdir = "./";
const char *h5file = "frame_data.h5";
const char *fpfile = "frame_particle_data.txt";
const char *fefile = "frame_element_data.txt";
const char *statefile = "state.txt";
const char *cfgfile = "simulation.cfg";

volatile int want_sigterm = 0;

/* loading.c (user defined) */
void initial_loads(job_t *job);
void time_varying_loads(job_t *job);

/*----------------------------------------------------------------------------*/
void signal_callback_handler(int signum)
{
   want_sigterm = 1;
   return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void usage(char *program_name)
{
    printf("%s: [OPTIONS] grid_file particle_file [t_max]\n", program_name);
    printf("\tdefault t_max is 1.0\n");
    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
int main(int argc, char **argv)
{
    /* libconfuse config parsing variables */
    cfg_opt_t timestep_opts[] =
    {
        CFG_FLOAT("dt-max", 1e-2, CFGF_NONE),
        CFG_FLOAT("dt-min", 1e-10, CFGF_NONE),
        CFG_FLOAT("dt", 1e-6, CFGF_NONE),
        CFG_INT("automatic-dt", 1, CFGF_NONE),
        CFG_INT("allow-dt-increase", 0, CFGF_NONE),
        CFG_INT("stable-dt-threshold", 4, CFGF_NONE),
        CFG_END()
    };
    cfg_opt_t implicit_opts[] =
    {
        CFG_FLOAT("displacement-norm-ratio", 1e-2, CFGF_NONE),
        CFG_FLOAT("residual-norm-ratio", 1e-1, CFGF_NONE),
        CFG_FLOAT("converged-displacement-norm", 1e-12, CFGF_NONE),
        CFG_INT("unstable-iteration-count", 10, CFGF_NONE),
        CFG_END()
    };
    cfg_opt_t opts[] =
    {
        CFG_SEC("timestep", timestep_opts, CFGF_NONE),
        CFG_SEC("implicit", implicit_opts, CFGF_NONE),
        CFG_END()
    };
    cfg_t *cfg;
    cfg_t *cfg_timestep;
    cfg_t *cfg_implicit;

    const int num_threads = 1;

    double h = 0.1;
    int N = (1+ceil(1/h));
    int i;
    int j;
    int len;
    double t_stop = 1;

    grid_t g;
    particle_t *pdata = NULL;
    int plen;

    job_t *job = NULL;
    FILE *frame_fd;
    FILE *state_fd;
    FILE *felement_fd;

    pthread_t *threads;
    threadtask_t *tasks;

    h5writer_t *h5;

    int opt;
    int leftover_argc;
    char **leftover_argv;

    /* allows us to trap interrupt and dump data before exit. */
    signal(SIGINT, signal_callback_handler);

    /* set default command line state */
    g_state.outputdir = outputdir;
    g_state.h5file = NULL;
    g_state.fpfile = NULL;
    g_state.fefile = NULL;
    g_state.restart = 0;
    g_state.tmax = 0;

    /* parse command line options */
    opt = getopt(argc, argv, g_optstring);

    while (opt != -1) {
        switch (opt) {
            case 'r':
                g_state.restart = 1;
                state_fd = fopen(optarg, "r");
                job = read_state(state_fd);
                fclose(state_fd);
                state_fd = NULL;
                break;
            case 'o':
                g_state.outputdir = optarg;
                break;
            default:
                break;
        }
        opt = getopt(argc, argv, g_optstring);
    }


    /* Add an additional byte in case we need to add a '/'. */
    len = 2 + strlen(g_state.outputdir);
    printf("Using output directory \"%s\".\n", g_state.outputdir);
    g_state.h5file = (char *)malloc(len + strlen(h5file));
    g_state.fpfile = (char *)malloc(len + strlen(fpfile));
    g_state.fefile = (char *)malloc(len + strlen(fefile));
    g_state.statefile = (char *)malloc(len + strlen(statefile));
    if (g_state.outputdir[len-2] == '/') {
        sprintf(g_state.h5file, "%s%s", g_state.outputdir, h5file);
        sprintf(g_state.fpfile, "%s%s", g_state.outputdir, fpfile);
        sprintf(g_state.fefile, "%s%s", g_state.outputdir, fefile);
        sprintf(g_state.statefile, "%s%s", g_state.outputdir, statefile);
    } else {
        sprintf(g_state.h5file, "%s/%s", g_state.outputdir, h5file);
        sprintf(g_state.fpfile, "%s/%s", g_state.outputdir, fpfile);
        sprintf(g_state.fefile, "%s/%s", g_state.outputdir, fefile);
        sprintf(g_state.statefile, "%s/%s", g_state.outputdir, statefile);
    }

    leftover_argv = argv + optind;
    leftover_argc = argc - optind;

    if (g_state.restart) {
        goto job_start;
    }

    if (leftover_argc < 2) {
        usage(argv[0]);
        exit(0);
    }

    read_grid_params(&g, leftover_argv[0]);
    printf("Finished reading grid file \"%s\".\n", leftover_argv[0]);
    read_particles(&pdata, &plen, leftover_argv[1]);
    printf("Finished reading particle file \"%s\".\n", leftover_argv[1]);

    if (leftover_argc >= 3) {
        sscanf(leftover_argv[2], "%lf", &t_stop);
    } else {
        t_stop = 1;
    }

    N = g.N;
    h = g.len / (N-1);

    printf("Grid parameters: N = %d, h = %g.\n", N, h);

    job = mpm_init(N, h, pdata, plen, t_stop);

    /* parse config file */
    cfg = cfg_init(opts, CFGF_NONE);

    if (cfg_parse(cfg, cfgfile) == CFG_PARSE_ERROR) {
        fprintf(stderr, "Fatal -- cannot parse configuration file '%s'.\n",
            cfgfile);
        exit(127);
    }

    /* section for timestep */
    cfg_timestep = cfg_getsec(cfg, "timestep");

    job->timestep.dt_max = cfg_getfloat(cfg_timestep, "dt-max");
    job->timestep.dt_min = cfg_getfloat(cfg_timestep, "dt-min");
    job->timestep.dt = cfg_getfloat(cfg_timestep, "dt");

    job->timestep.automatic_dt =
        cfg_getint(cfg_timestep, "automatic-dt");

    if (job->timestep.automatic_dt != 0) {
        job->timestep.dt = job->dt;
    }

    job->timestep.allow_dt_increase =
        cfg_getint(cfg_timestep, "allow-dt-increase");
    job->timestep.stable_dt_threshold =
        cfg_getint(cfg_timestep, "stable-dt-threshold");

    fprintf(stderr, "Timestep options set:\n");
    fprintf(stderr, "dt_max: %e\n", job->timestep.dt_max);
    fprintf(stderr, "dt_min: %e\n", job->timestep.dt_min);
    fprintf(stderr, "dt: %e\n", job->timestep.dt);
    fprintf(stderr, "automatic_dt: %d\n", job->timestep.automatic_dt);
    fprintf(stderr, "allow_dt_increase: %d\n", job->timestep.allow_dt_increase);
    fprintf(stderr, "stable_dt_threshold: %d\n", job->timestep.stable_dt_threshold);

    /* section for implicit solver */
    cfg_implicit = cfg_getsec(cfg, "implicit");

    job->implicit.du_norm_ratio =
        cfg_getfloat(cfg_implicit, "displacement-norm-ratio");
    job->implicit.q_norm_ratio =
        cfg_getfloat(cfg_implicit, "residual-norm-ratio");
    job->implicit.du_norm_converged =
        cfg_getfloat(cfg_implicit, "converged-displacement-norm");
    job->implicit.unstable_iteration_count =
        cfg_getint(cfg_implicit, "unstable-iteration-count");

    fprintf(stderr, "Implicit options set:\n");
    fprintf(stderr, "du_norm_ratio: %e\n", job->implicit.du_norm_ratio);
    fprintf(stderr, "q_norm_ratio: %e\n", job->implicit.q_norm_ratio);
    fprintf(stderr, "du_norm_converged: %e\n", job->implicit.du_norm_converged);


job_start:
    dispg(job->t_stop);
    dispg(job->dt);
    dispg(job->h);
    dispd(job->num_particles);
    dispd(job->num_nodes);
    dispd(job->num_elements);

    threads = (pthread_t *)malloc(sizeof(pthread_t) * num_threads);
    tasks = (threadtask_t *)malloc(sizeof(threadtask_t) * num_threads);
    for (i = 0; i < num_threads; i++) {
        tasks[i].id = i;
        tasks[i].num_threads = num_threads;
        tasks[i].job = job;
    }
    printf("Using %d %s.\n", num_threads, (num_threads > 1)?"threads":"thread");

    initial_loads(job);

    j = floor(job->t * SAMPLE_HZ);
    dispg(job->t);
    dispd(j);

    frame_fd = fopen(g_state.fpfile, "w");
    felement_fd = fopen(g_state.fefile, "w");
    state_fd = fopen(g_state.statefile, "w");
    h5 = h5_init(g_state.h5file, job);

    while (job->t < job->t_stop && !want_sigterm) {
        mpm_step(job);

        if (job->t >= (j / SAMPLE_HZ)) {
            write_frame(frame_fd, j, job->t, job);
/*            h5_write_frame(h5, j, job);*/
            write_element_frame(felement_fd, j, job->t, job);

            j++;
            printf("\rt = %05.3fs\t[%3d%%]\t", job->t, (int)((100 * job->t) / (job->t_stop)));
            fflush(stdout);
        }

        time_varying_loads(job);
    }

    /* dump state to file */
    write_state(state_fd, job);
/*    h5_write_state("state.h5", job);*/

    fclose(frame_fd);
    fclose(felement_fd);
    fclose(state_fd);

    printf("\n");
    printf("Freeing allocated memory.\n");

    cfg_free(cfg);
    mpm_cleanup(job);
    h5_cleanup(h5);

    free(pdata);
    pdata = NULL;

    /* kill all threads */
    pthread_exit(NULL);
    return 0;
}
/*----------------------------------------------------------------------------*/

