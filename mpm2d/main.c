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

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <dlfcn.h>

#include "exitcodes.h"

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

static struct state_s {
    char *outputdir;

    char *gridfile;
    char *particlefile;

    double tmax;
    int restart;        /* is this a restart? */
} g_state;

static const char* g_optstring = "o:r:p:g:";
char *outputdir = "./";
const char *cfgfile = "simulation.cfg";

volatile int want_sigterm = 0;

/* function pointer for type of mpm step (implicit or explicit). */
void (*mpm_step)(job_t *);

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
int validate_path(cfg_t *cfg, cfg_opt_t *opt)
{
    char *str = cfg_opt_getnstr(opt, cfg_opt_size(opt) - 1);
    struct stat s;

    /* check if the directory supplied exists and is actually a directory. */
    if(stat(str, &s) != 0 || !S_ISDIR(s.st_mode))
    {
        cfg_error(cfg, "Unable to open '%s' for option '%s.%s'.",
            str, cfg->name, opt->name);
        return -1;
    }

    return 0;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
int validate_readonly_file(cfg_t *cfg, cfg_opt_t *opt)
{
    char *str = cfg_opt_getnstr(opt, cfg_opt_size(opt) - 1);
    struct stat s;

    /* check if the file supplied exists and is a regular file. */
    if(stat(str, &s) != 0 || !S_ISREG(s.st_mode))
    {
        cfg_error(cfg, "Unable to open '%s' for option '%s.%s'.",
            str, cfg->name, opt->name);
        return -1;
    }

    return 0;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
int set_solver_type(cfg_t *cfg, cfg_opt_t *opt, const char *value, void *result)
{
    if (strcmp(value, "implicit") == 0) {
        *(enum solver_e *)result = IMPLICIT_SOLVER;
    } else if (strcmp(value, "explicit") == 0) {
        *(enum solver_e *)result = EXPLICIT_SOLVER;
    } else {
        cfg_error(cfg, "Invalid value for option '%s': %s", opt->name, value);
        return -1;
    }
    return 0;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void usage(char *program_name)
{
    printf("%s: [OPTIONS] grid_file particle_file [t_max]\n", program_name);
    printf("\tOPTIONS are any of:\n");
    printf("\t\t-o DIR, specify output directory. Overrides config file value.\n");
    printf("\t\t-r STATE, restart analysis using STATE as input.\n");
    printf("\t\t-p PFILE, particle file to use. Overrides config file value.\n");
    printf("\t\t-g GFILE, grid file to use. Overrides config file value.\n");
    printf("\n\tdefault t_max is 1.0\n");
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
    cfg_opt_t solver_opts[] =
    {
        CFG_INT_CB("solver-type", IMPLICIT_SOLVER, CFGF_NONE, &set_solver_type),
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
    cfg_opt_t explicit_opts[] =
    {
        CFG_END()
    };
    cfg_opt_t output_opts[] =
    {
        CFG_STR("job-name", "Default", CFGF_NONE),
        CFG_STR("job-description", "None", CFGF_NONE),
        CFG_INT("job-id", 0, CFGF_NONE),
        CFG_STR("directory", "jobs/", CFGF_NONE),
        CFG_STR("user", "unknown", CFGF_NONE),
        CFG_INT("override-directory-with-id", 0, CFGF_NONE),
        CFG_INT("prepend-date", 1, CFGF_NONE),
        CFG_STR("particle-file", "frame_particle_data.txt", CFGF_NONE),
        CFG_INT("enable-particle-output", 1, CFGF_NONE),
        CFG_STR("element-file", "frame_element_data.txt", CFGF_NONE),
        CFG_INT("enable-element-output", 0, CFGF_NONE),
        CFG_STR("state-file", "state.txt", CFGF_NONE),
        CFG_INT("trap-terminate-interrupt", 1, CFGF_NONE),
        CFG_INT("save-state-on-terminate", 1, CFGF_NONE),
        CFG_STR("log-file", "job.log", CFGF_NONE),
        CFG_INT("log-level", 1, CFGF_NONE),
        CFG_END()
    };
    cfg_opt_t input_opts[] =
    {
        CFG_STR("initial-particle-file", "particles.txt", CFGF_NONE),
        CFG_STR("grid-file", "grid.txt", CFGF_NONE),
        CFG_END()
    };
    cfg_opt_t material_opts[] =
    {
        CFG_STR("material-file", "builtin.so", CFGF_NONE),
        CFG_INT("use-builtin", 1, CFGF_NONE),
        CFG_END()
    };
    cfg_opt_t opts[] =
    {
        CFG_SEC("timestep", timestep_opts, CFGF_NONE),
        CFG_SEC("solver", solver_opts, CFGF_NONE),
        CFG_SEC("implicit", implicit_opts, CFGF_NONE),
        CFG_SEC("explicit", explicit_opts, CFGF_NONE),
        CFG_SEC("output", output_opts, CFGF_NONE),
        CFG_SEC("input", input_opts, CFGF_NONE),
        CFG_SEC("material", material_opts, CFGF_NONE),
        CFG_END()
    };
    cfg_t *cfg;
    cfg_t *cfg_solver;
    cfg_t *cfg_timestep;
    cfg_t *cfg_implicit;
    cfg_t *cfg_output;
    cfg_t *cfg_input;
    cfg_t *cfg_material;

    const char *solver_names[] = {
        "Implicit",
        "Explicit",
        "N/A"
    };

    const int num_threads = 1;
    int command_line_outputdir = 0;
    int command_line_gridfile = 0;
    int command_line_particlefile = 0;
    char *s;
    char *s_dlerror;

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
    FILE *state_fd;

    pthread_t *threads;
    threadtask_t *tasks;

    void *material_so_handle;

    int opt;
    int leftover_argc;
    char **leftover_argv;


    /* parse config file */
    cfg = cfg_init(opts, CFGF_NONE);
    cfg_set_validate_func(cfg, "output|directory", validate_path);

    if (cfg_parse(cfg, cfgfile) == CFG_PARSE_ERROR) {
        fprintf(stderr, "Fatal -- cannot parse configuration file '%s'.\n",
            cfgfile);
        exit(EXIT_ERROR_CFG_PARSE);
    }

    /* section for input */
    cfg_input = cfg_getsec(cfg, "input");
    g_state.particlefile = cfg_getstr(cfg_input, "initial-particle-file");
    g_state.gridfile = cfg_getstr(cfg_input, "grid-file");

    /* allows us to trap interrupt and dump data before exit. */
    signal(SIGINT, signal_callback_handler);

    /* set default command line state */
    g_state.outputdir = outputdir;
/*    g_state.h5file = NULL;*/
/*    g_state.fpfile = NULL;*/
/*    g_state.fefile = NULL;*/
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
                command_line_outputdir = 1;
                break;
            case 'p':
                g_state.particlefile = optarg;
                command_line_particlefile = 1;
                break;
            case 'g':
                g_state.gridfile = optarg;
                command_line_gridfile = 1;
                break;
            default:
                break;
        }
        opt = getopt(argc, argv, g_optstring);
    }


    /* Add an additional byte if we need to add a '/'. */
    if (command_line_outputdir != 0) {
        len = strlen(g_state.outputdir);
        printf("Using output directory \"%s\".\n", g_state.outputdir);
        if (g_state.outputdir[len-1] != '/') {
            s = g_state.outputdir;
            len += 2;
            g_state.outputdir = (char *)malloc(len);
            snprintf(g_state.outputdir, len, "%s/", s);
        }
    }

    leftover_argv = argv + optind;
    leftover_argc = argc - optind;

    if (g_state.restart) {
        goto job_start;
    }

/*    if (leftover_argc < 2) {*/
/*        usage(argv[0]);*/
/*        exit(0);*/
/*    }*/

    read_grid_params(&g, g_state.gridfile);
    printf("Finished reading grid file \"%s\".\n", g_state.gridfile);
    read_particles(&pdata, &plen, g_state.particlefile);
    printf("Finished reading particle file \"%s\".\n", g_state.particlefile);

    if (leftover_argc >= 1) {
        sscanf(leftover_argv[0], "%lf", &t_stop);
    } else {
        t_stop = 1;
    }

    N = g.N;
    h = g.len / (N-1);

    printf("Grid parameters: N = %d, h = %g.\n", N, h);

    job = mpm_init(N, h, pdata, plen, t_stop);

    /* section for material options */
    cfg_material = cfg_getsec(cfg, "material");
    fprintf(stderr, "old function pointer %p\n", job->material.calculate_stress); 
    job->material.use_builtin = cfg_getint(cfg_material, "use-builtin");
    if (job->material.use_builtin == 0) {
        job->material.material_filename =
            cfg_getstr(cfg_material, "material-file");
        material_so_handle =
            dlopen(job->material.material_filename, RTLD_LAZY);
        if (material_so_handle == NULL) {
            fprintf(stderr, "Can't dlopen() material file '%s': %s.\n",
                job->material.material_filename, dlerror());
            exit(EXIT_ERROR_MATERIAL_FILE);
        }
        *(void **)(&(job->material.material_init)) =
            dlsym(material_so_handle, "material_init");
        if ((s_dlerror = dlerror()) != NULL) {
            fprintf(stderr, "Error loading symbol 'material_init': %s.\n",
                s_dlerror);
            exit(EXIT_ERROR_MATERIAL_FILE);
        }
        *(void **)(&(job->material.calculate_stress)) =
            dlsym(material_so_handle, "calculate_stress");
        if ((s_dlerror = dlerror()) != NULL) {
            fprintf(stderr, "Error loading symbol 'calculate_stress': %s.\n",
                s_dlerror);
            exit(EXIT_ERROR_MATERIAL_FILE);
        }
    }
    fprintf(stderr, "new function pointer %p\n", job->material.calculate_stress); 
    fprintf(stderr, "\nMaterial options set:\n");
    fprintf(stderr, "material_filename: %s\n", job->material.material_filename);
    fprintf(stderr, "use_builtin: %d\n", job->material.use_builtin);
    

    /* section for solver options */
    cfg_solver = cfg_getsec(cfg, "solver");
    job->solver = cfg_getint(cfg_solver, "solver-type");
    fprintf(stderr, "\nSolver options set:\n");
    fprintf(stderr, "solver: %d (%s)\n",
        job->solver, solver_names[(int)job->solver]);
    if (job->solver == IMPLICIT_SOLVER) {
        mpm_step = implicit_mpm_step;
    } else if (job->solver == EXPLICIT_SOLVER) {
        mpm_step = explicit_mpm_step;
    } else {
        fprintf(stderr, "Unknown solver type.\n");
        exit(-1);
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

    fprintf(stderr, "\nTimestep options set:\n");
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

    fprintf(stderr, "\nImplicit options set:\n");
    fprintf(stderr, "du_norm_ratio: %e\n", job->implicit.du_norm_ratio);
    fprintf(stderr, "q_norm_ratio: %e\n", job->implicit.q_norm_ratio);
    fprintf(stderr, "du_norm_converged: %e\n", job->implicit.du_norm_converged);
    fprintf(stderr, "unstable_iteration_count: %d\n", job->implicit.unstable_iteration_count);

    /* section for output */
    cfg_output = cfg_getsec(cfg, "output");

    /*
        We can specify the output directory on the commandline, use that to
        override the configuration file.
    */
    if (command_line_outputdir) {
        job->output.directory = g_state.outputdir;
    } else {
        job->output.directory = cfg_getstr(cfg_output, "directory");
    }
    job->output.user = cfg_getstr(cfg_output, "user");
    job->output.particle_filename = cfg_getstr(cfg_output, "particle-file");
    job->output.element_filename = cfg_getstr(cfg_output, "element-file");
    job->output.state_filename = cfg_getstr(cfg_output, "state-file");
    job->output.log_filename = cfg_getstr(cfg_output, "log-file");

    job->output.modified_directory = 0;
    len = strlen(job->output.directory);
    if (job->output.directory[len - 1] != '/') {
        job->output.directory = (char *)malloc(len + 2);
        strcpy(job->output.directory, cfg_getstr(cfg_output, "directory"));
        job->output.directory[len] = '/';
        job->output.directory[len + 1] = 0;
        job->output.modified_directory = 1;
    }

    len = strlen(job->output.directory) + strlen(job->output.particle_filename) + 1;
    job->output.particle_filename_fullpath = (char *)malloc(len);
    snprintf(job->output.particle_filename_fullpath, len, "%s%s",
        job->output.directory, job->output.particle_filename);
    len = strlen(job->output.directory) + strlen(job->output.element_filename) + 1;
    job->output.element_filename_fullpath = (char *)malloc(len);
    snprintf(job->output.element_filename_fullpath, len, "%s%s",
        job->output.directory, job->output.element_filename);
    len = strlen(job->output.directory) + strlen(job->output.state_filename) + 1;
    job->output.state_filename_fullpath = (char *)malloc(len);
    snprintf(job->output.state_filename_fullpath, len, "%s%s",
        job->output.directory, job->output.state_filename);
    len = strlen(job->output.directory) + strlen(job->output.log_filename) + 1;
    job->output.log_filename_fullpath = (char *)malloc(len);
    snprintf(job->output.log_filename_fullpath, len, "%s%s",
        job->output.directory, job->output.log_filename);

    job->output.particle_fd = fopen(job->output.particle_filename_fullpath, "w");
    if (job->output.particle_fd == NULL) {
      printf( "Failed to open particle_fd file: %s\n", job->output.particle_filename_fullpath );
      exit(-1);
    }
    job->output.element_fd = fopen(job->output.element_filename_fullpath, "w");
    if (job->output.particle_fd == NULL) {
      printf( "Failed to open element_fd file: %s\n", job->output.element_filename_fullpath );
      exit(-1);
    }
    job->output.state_fd = fopen(job->output.state_filename_fullpath, "w");
    if (job->output.particle_fd == NULL) {
      printf( "Failed to open state_fd file: %s\n", job->output.state_filename_fullpath );
      exit(-1);
    }
    job->output.log_fd = fopen(job->output.log_filename_fullpath, "w");
    if (job->output.particle_fd == NULL) {
      printf( "Failed to open log_fd file: %s\n", job->output.log_filename_fullpath );
      exit(-1);
    }

    fprintf(stderr, "\nOutput options set:\n");
    fprintf(stderr, "output_directory: %s\n", job->output.directory);
    fprintf(stderr, "user: %s\n", job->output.user);
    fprintf(stderr, "particle_filename: %s\n", job->output.particle_filename);
    fprintf(stderr, "element_filename: %s\n", job->output.element_filename);
    fprintf(stderr, "state_filename: %s\n", job->output.state_filename);

    fprintf(stderr, "particle_filename_fullpath: %s\n",
        job->output.particle_filename_fullpath);
    fprintf(stderr, "element_filename_fullpath: %s\n",
        job->output.element_filename_fullpath);
    fprintf(stderr, "state_filename_fullpath: %s\n",
        job->output.state_filename_fullpath);
    fprintf(stderr, "log_filename_fullpath: %s\n",
        job->output.log_filename_fullpath);

/*    exit(0);*/

    job->dt = job->timestep.dt;
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
    (*(job->material.material_init))(job);

    j = floor(job->t * SAMPLE_HZ);
    dispg(job->t);
    dispd(j);

/*    frame_fd = fopen(g_state.fpfile, "w");*/
/*    felement_fd = fopen(g_state.fefile, "w");*/
/*    state_fd = fopen(g_state.statefile, "w");*/
/*    h5 = h5_init(g_state.h5file, job);*/

    while (job->t < job->t_stop && !want_sigterm) {
        mpm_step(job);

        if (job->t >= (j / SAMPLE_HZ)) {
            write_frame(job->output.particle_fd, j, job->t, job);
/*            h5_write_frame(h5, j, job);*/
            write_element_frame(job->output.element_fd, j, job->t, job);

            j++;
            printf("\rt = %05.3fs\t[%3d%%]\t", job->t, (int)((100 * job->t) / (job->t_stop)));
            fflush(stdout);
        }

        time_varying_loads(job);
    }

    /* dump state to file */
    write_state(job->output.state_fd, job);
/*    h5_write_state("state.h5", job);*/

    printf("Closing files.\n");
    fclose(job->output.particle_fd);
    fclose(job->output.element_fd);
    fclose(job->output.state_fd);

/*    fclose(frame_fd);*/
/*    fclose(felement_fd);*/
/*    fclose(state_fd);*/

    printf("\n");
    printf("Freeing allocated memory.\n");

    cfg_free(cfg);
    if (job->output.modified_directory) {
        free(job->output.directory);
    }
    free(job->output.particle_filename_fullpath);
    free(job->output.element_filename_fullpath);
    free(job->output.state_filename_fullpath);
    mpm_cleanup(job);
/*    h5_cleanup(h5);*/

    free(pdata);
    pdata = NULL;

    /* kill all threads */
    pthread_exit(NULL);
    return 0;
}
/*----------------------------------------------------------------------------*/

