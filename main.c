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

#include <time.h>

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

#define MIN(x,y) (((x) > (y))?(y):(x))
#define MAX(x,y) (((x) < (y))?(y):(x))

#define DEFAULT_SAMPLE_HZ 60.0f

/* A bit ugly, but useful for unwinding error handling. */
#define JUMP_IF_NULL(val, gotolabel, msg) \
    do { if (val == NULL) { fprintf(stderr, "%s", msg); goto gotolabel; } } while(0)
/* If 'cond' is true, print msg to stderr and jump to gotolabel. */
#define JUMP_IF(cond, gotolabel, msg) \
    do { if (cond) { fprintf(stderr, "%s", msg); goto gotolabel; } } while(0)

static struct state_s {
    char *outputdir;

    char *gridfile;
    char *particlefile;
    char *materialso;

    double tmax;
    int restart;        /* is this a restart? */
} g_state;

static const char* g_optstring = "hc:o:r:p:g:u:t:";
char *outputdir = "./";
const char *default_cfgfile = "simulation.cfg";

volatile int want_sigterm = 0;

/* function pointer for type of mpm step (implicit or explicit). */
void (*mpm_step)(job_t *);

/* loading.c (user defined) */
void initial_loads(job_t *job);
void time_varying_loads(job_t *job);

/* threaded helper function */
void *mpm_run_until(void *_task);

/*----------------------------------------------------------------------------*/
void signal_callback_handler(int signum)
{
   want_sigterm = 1;
   fprintf(stderr, "\nRecieved sigterm, terminating.\n");
   exit(0);
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
    } else if (strcmp(value, "explicit-usf") == 0) {
        *(enum solver_e *)result = EXPLICIT_SOLVER_USF;
    } else if (strcmp(value, "explicit-usl") == 0) {
        *(enum solver_e *)result = EXPLICIT_SOLVER_USL;
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
    printf("%s: [OPTIONS] [t_max]\n", program_name);
    printf("\tOPTIONS are any of:\n");
    printf("\t\t-o DIR, specify output directory. Overrides config file value.\n");
    printf("\t\t-r STATE, restart analysis using STATE as input.\n");
    printf("\t\t-p PFILE, particle file to use. Overrides config file value.\n");
    printf("\t\t-g GFILE, grid file to use. Overrides config file value.\n");
    printf("\t\t-u MATERIAL, material shared object to use. Overrides config file value.\n");
    printf("\t\t-c CFGFILE, configuration file to use. Default is '%s'.\n", default_cfgfile);
    printf("\t\t-t THREADS, number of threads to use. Default is 1.\n");
    printf("\t\t-h This help message.");
    printf("\n\tdefault t_max is 1.0\n");
    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
int main(int argc, char **argv)
{
    /* Configuration file used by libconfuse. */
    char *cfgfile = default_cfgfile;
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
        CFG_FLOAT("sample-rate", DEFAULT_SAMPLE_HZ, CFGF_NONE),
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
        CFG_FLOAT_LIST("properties", "{}", CFGF_NONE),
        CFG_INT_LIST("integer-properties", "{}", CFGF_NONE),
        CFG_END()
    };
    cfg_opt_t boundary_opts[] =
    {
        CFG_STR("boundary-conditions-file", "builtin.so", CFGF_NONE),
        CFG_INT("use-builtin", 1, CFGF_NONE),
        CFG_FLOAT_LIST("properties", "{}", CFGF_NONE),
        CFG_INT_LIST("integer-properties", "{}", CFGF_NONE),
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
        CFG_SEC("boundary-conditions", boundary_opts, CFGF_NONE),
        CFG_FUNC("include", cfg_include),
        CFG_END()
    };
    cfg_t *cfg;
    cfg_t *cfg_solver;
    cfg_t *cfg_timestep;
    cfg_t *cfg_implicit;
    cfg_t *cfg_output;
    cfg_t *cfg_input;
    cfg_t *cfg_material;
    cfg_t *cfg_boundary;

    const char *solver_names[] = {
        "Implicit",
        "Explicit - Update Stress First",
        "Explicit - Update Stress Last",
        "N/A"
    };

    int num_threads = 1;
    int command_line_outputdir = 0;
    int command_line_gridfile = 0;
    int command_line_particlefile = 0;
    int command_line_materialso = 0;
    char *s;
    char *s_dlerror;

    double h = 0.1;
    int N = (1+ceil(1/h));
    int i;
    int j;
    int len;
    double t_stop = 1;

    struct timespec wallstart, wallstop;
    long ns;
    long stepcount;

    grid_t g;
    particle_t *pdata = NULL;
    int plen;

    job_t *job = NULL;
    FILE *state_fd;

    void *material_so_handle;

    int opt;
    int leftover_argc;
    char **leftover_argv;

    threadtask_t *tasks;
    pthread_t *threads;

    size_t psplit, nsplit, esplit;

    /* set default command line state */
    g_state.outputdir = outputdir;
    g_state.restart = 0;
    g_state.tmax = 0;

    /* parse command line options */
    opt = getopt(argc, argv, g_optstring);

    while (opt != -1) {
        switch (opt) {
            case 'c':
                cfgfile = optarg;
                /* Make sure file exists and is readable. */
                if (access(cfgfile, R_OK | F_OK) != 0) {
                    JUMP_IF_NULL(NULL, _cfgfile_error, "Can't read config file.\n");
                }
                break;
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
            case 'u':
                g_state.materialso = optarg;
                command_line_materialso = 1;
                break;
            case 't':
                num_threads = atoi(optarg);
                if (num_threads <= 0) {
                    num_threads = 1;
                }
                break;
            case 'h':
            default:
                usage(argv[0]);
                goto _commandline_error;
        }
        opt = getopt(argc, argv, g_optstring);
    }

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
    if (command_line_particlefile != 1) {
        g_state.particlefile = cfg_getstr(cfg_input, "initial-particle-file");
    }
    if (command_line_gridfile != 1) {
        g_state.gridfile = cfg_getstr(cfg_input, "grid-file");
    }

    /* allows us to trap interrupt and dump data before exit. */
    signal(SIGINT, signal_callback_handler);

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

    JUMP_IF(read_grid_params(&g, g_state.gridfile) != 0,
        _read_grid_file_error, "Error reading grid file.\n");
    printf("Finished reading grid file \"%s\".\n", g_state.gridfile);
    JUMP_IF(read_particles(&pdata, &plen, g_state.particlefile) != 0,
        _read_particle_file_error, "Error reading particle file.\n");
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
    job->material.calculate_stress_threaded = NULL;
    if (command_line_materialso != 0) {
        job->material.use_builtin = 0;
        job->material.material_filename = g_state.materialso;
    } else {
        job->material.use_builtin = cfg_getint(cfg_material, "use-builtin");
        job->material.material_filename = cfg_getstr(cfg_material, "material-file");
    }
    if (job->material.use_builtin == 0) {
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
        *(void **)(&(job->material.calculate_stress_threaded)) =
            dlsym(material_so_handle, "calculate_stress_threaded");
        if ((s_dlerror = dlerror()) != NULL) {
            fprintf(stderr, "Error loading symbol 'calculate_stress_threaded': %s.\n",
                s_dlerror);
            exit(EXIT_ERROR_MATERIAL_FILE);
        }
    }

    job->material.num_fp64_props = cfg_size(cfg_material, "properties");
    job->material.num_int_props = cfg_size(cfg_material, "integer-properties");
    job->material.fp64_props = (double *)malloc(sizeof(double) * job->material.num_fp64_props);
    job->material.int_props = (int *)malloc(sizeof(int) * job->material.num_int_props);
    for (i = 0; i < job->material.num_fp64_props; i++) {
        job->material.fp64_props[i] = cfg_getnfloat(cfg_material, "properties", i);
    }
    for (i = 0; i < job->material.num_int_props; i++) {
        job->material.int_props[i] = cfg_getnint(cfg_material, "integer-properties", i);
    }

    fprintf(stderr, "\nMaterial options set:\n");
    fprintf(stderr, "material_filename: %s\n", job->material.material_filename);
    fprintf(stderr, "use_builtin: %d\n", job->material.use_builtin);
    fprintf(stderr, "num_fp64_props: %d\n", job->material.num_fp64_props);
    fprintf(stderr, "num_int_props: %d\n", job->material.num_int_props);

    fprintf(stderr, "fp64_props: { ");
    for (i = 0; i < job->material.num_fp64_props; i++) {
        fprintf(stderr, "%.3g ", job->material.fp64_props[i]);
    }
    fprintf(stderr, "}\n");
    fprintf(stderr, "int_props: { ");
    for (i = 0; i < job->material.num_int_props; i++) {
        fprintf(stderr, "%d ", job->material.int_props[i]);
    }
    fprintf(stderr, "}\n");

    /* section for boundary condition options */
    /* XXX: fix implementation! need to do dlopen/dlsym etc */
    cfg_boundary = cfg_getsec(cfg, "boundary-conditions");

    job->boundary.num_fp64_props = cfg_size(cfg_boundary, "properties");
    job->boundary.num_int_props = cfg_size(cfg_boundary, "integer-properties");
    job->boundary.fp64_props = (double *)malloc(sizeof(double) * job->boundary.num_fp64_props);
    job->boundary.int_props = (int *)malloc(sizeof(int) * job->boundary.num_int_props);
    for (i = 0; i < job->boundary.num_fp64_props; i++) {
        job->boundary.fp64_props[i] = cfg_getnfloat(cfg_boundary, "properties", i);
    }
    for (i = 0; i < job->boundary.num_int_props; i++) {
        job->boundary.int_props[i] = cfg_getnint(cfg_boundary, "integer-properties", i);
    }

    fprintf(stderr, "\nBoundary options set:\n");
    fprintf(stderr, "num_fp64_props: %d\n", job->boundary.num_fp64_props);
    fprintf(stderr, "num_int_props: %d\n", job->boundary.num_int_props);

    fprintf(stderr, "fp64_props: { ");
    for (i = 0; i < job->boundary.num_fp64_props; i++) {
        fprintf(stderr, "%.3g ", job->boundary.fp64_props[i]);
    }
    fprintf(stderr, "}\n");
    fprintf(stderr, "int_props: { ");
    for (i = 0; i < job->boundary.num_int_props; i++) {
        fprintf(stderr, "%d ", job->boundary.int_props[i]);
    }
    fprintf(stderr, "}\n");

    /* section for solver options */
    cfg_solver = cfg_getsec(cfg, "solver");
    job->solver = cfg_getint(cfg_solver, "solver-type");
    fprintf(stderr, "\nSolver options set:\n");
    fprintf(stderr, "solver: %d (%s)\n",
        job->solver, solver_names[(int)job->solver]);
    if (job->solver == IMPLICIT_SOLVER) {
        mpm_step = implicit_mpm_step;
    } else if (job->solver == EXPLICIT_SOLVER_USF) {
        mpm_step = explicit_mpm_step_usf;
    } else if (job->solver == EXPLICIT_SOLVER_USL) {
        mpm_step = explicit_mpm_step_usl;
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

    /*
        Modify the output directory to add a trailing slash if it doesn't
        exist.
    */
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
        JUMP_IF_NULL(job->output.particle_fd, _particle_fd_error,
            "Can't open particle file for output.\n");
    job->output.element_fd = fopen(job->output.element_filename_fullpath, "w");
        JUMP_IF_NULL(job->output.particle_fd, _element_fd_error,
            "Can't open element file for output.\n");
    job->output.state_fd = fopen(job->output.state_filename_fullpath, "w");
        JUMP_IF_NULL(job->output.particle_fd, _state_fd_error,
            "Can't open state file for output.\n");
    job->output.log_fd = fopen(job->output.log_filename_fullpath, "w");
        JUMP_IF_NULL(job->output.particle_fd, _log_fd_error,
            "Can't open log file for output.\n");

    /* sampling rate */
    job->output.sample_rate_hz = cfg_getfloat(cfg_output, "sample-rate");
    if (job->output.sample_rate_hz < 0) {
        job->output.sample_rate_hz = DEFAULT_SAMPLE_HZ;
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

    fprintf(stderr, "sample_rate_hz: %5.4f\n", job->output.sample_rate_hz);
/*    exit(0);*/

    job->dt = job->timestep.dt;
job_start:
    dispg(job->t_stop);
    dispg(job->dt);
    dispg(job->h);
    dispd(job->num_particles);
    dispd(job->num_nodes);
    dispd(job->num_elements);

    job->threads = (pthread_t *)malloc(sizeof(pthread_t) * num_threads);
    tasks = (threadtask_t *)malloc(sizeof(threadtask_t) * num_threads);
    job->step_barrier = (pthread_barrier_t *)malloc(sizeof(pthread_barrier_t));
    job->serialize_barrier = (pthread_barrier_t *)malloc(sizeof(pthread_barrier_t));
    psplit = (job->num_particles / num_threads) + ((job->num_particles % num_threads != 0)?(1):(0));
    nsplit = (job->num_nodes / num_threads) + ((job->num_nodes % num_threads != 0)?(1):(0));
    esplit = (job->num_elements / num_threads) + ((job->num_elements % num_threads != 0)?(1):(0));
    for (i = 0; i < num_threads; i++) {
        tasks[i].id = i;
        tasks[i].num_threads = num_threads;
        tasks[i].job = job;
        /* split over particles */
        tasks[i].offset = i * psplit;
        tasks[i].blocksize = MIN(psplit, job->num_particles - tasks[i].offset);

        /* split over nodes */
        tasks[i].n_offset = i * nsplit;
        tasks[i].n_blocksize = MIN(nsplit, job->num_nodes - tasks[i].n_offset);
        
        /* split over elements */
        tasks[i].e_offset = i * esplit;
        tasks[i].e_blocksize = MIN(esplit, job->num_elements - tasks[i].e_offset);

        tasks[i].gidx_max = job->num_particles;
        tasks[i].gidx_min = 0;
        if (tasks[i].blocksize < 256 && num_threads > 1 && i != (num_threads -1)) {
            fprintf(stderr, "WARNING: blocksize for task is < 256! This may result in lower performance compared to a single threaded instance.\n");
        }
    }
    job->num_threads = num_threads;
    threads = (pthread_t *)malloc(sizeof(pthread_t) * job->num_threads);
    printf("Using %d %s.\n", num_threads, (num_threads > 1)?"threads":"thread");
    if(pthread_barrier_init(job->step_barrier, NULL, num_threads))
    {
        fprintf(stderr, "Could not create pthread step_barrier!\n");
        exit(EXIT_ERROR_THREADING);
    }
    if(pthread_barrier_init(job->serialize_barrier, NULL, num_threads))
    {
        fprintf(stderr, "Could not create pthread serialize_barrier!\n");
        exit(EXIT_ERROR_THREADING);
    }

    /* create element color lists on first step. */
    job->update_elementlists = (int *)malloc(sizeof(int) * job->num_threads);
    for (i = 0; i < job->num_threads; i++) {
        job->update_elementlists[i] = 1;
    }

    /* Create structure to sort particle ids for easier parallel mapping. */
/*    job->particle_by_element_color_offsets = (size_t *)malloc(sizeof(size_t) * job->num_colors * job->num_threads);*/
/*    job->particle_by_element_color_lengths = (size_t *)malloc(sizeof(size_t) * job->num_colors * job->num_threads);*/
/*    job->particle_by_element_color_list = (size_t *)malloc(sizeof(size_t) * job->num_particles);*/

    job->particle_by_element_color_lengths = (size_t *)malloc(sizeof(size_t) * job->num_colors * job->num_threads);
    job->particle_by_element_color_lists = (size_t **)malloc(sizeof(size_t *) * job->num_colors * job->num_threads);
    for (i = 0; i < (job->num_colors * job->num_threads); i++) {
        job->particle_by_element_color_lists[i] = (size_t *)malloc(sizeof(size_t) * job->num_particles);
    }

    /* Actually find filled elements for creating parallel list. */
    for (i = 0; i < job->num_elements; i++) {
        job->elements[i].filled = 0;
        job->elements[i].n = 0;
        job->elements[i].m = 0;
    }
    find_filled_elements(job);

    initial_loads(job);
    (*(job->material.material_init))(job);

    j = floor(job->t * job->output.sample_rate_hz);
    job->frame = j;
    dispg(job->t);
    dispd(j);

    if (validate_bc_properties(job) == 0) {
        fprintf(stderr, "Exiting due to error with BC properties.\n");
        exit(EXIT_ERROR_BC_PROPS);
    }

    job->frame = floor(job->t * job->output.sample_rate_hz);

    fprintf(stderr, "Starting timer...\n");
    clock_gettime(CLOCK_REALTIME, &wallstart);
    clock_gettime(CLOCK_REALTIME, &(job->tic));
    job->stepcount = 0;

    for (i = 0; i < (job->num_threads - 1); i++) {
        pthread_create(&(threads[i]), NULL, &mpm_run_until, &(tasks[i]));
    }
    mpm_run_until(&(tasks[job->num_threads-1]));

    for (i = 0; i < (job->num_threads - 1); i++) {
        pthread_join(threads[i], NULL);
    }

    fprintf(stderr, "\nStopping timer...\n");
    clock_gettime(CLOCK_REALTIME, &wallstop);
    ns = 1E9 * (wallstop.tv_sec - wallstart.tv_sec) + (wallstop.tv_nsec - wallstart.tv_nsec);
    printf("Elapsed Time: %.3fs\n", ns / 1E9);

    /* dump state to file */
    write_state(job->output.state_fd, job);
/*    h5_write_state("state.h5", job);*/

    printf("Closing files.\n");
    fclose(job->output.log_fd);
_log_fd_error:
    fclose(job->output.state_fd);
_state_fd_error:
    fclose(job->output.element_fd);
_element_fd_error:
    fclose(job->output.particle_fd);
_particle_fd_error:

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

_read_particle_file_error:
_read_grid_file_error:
_commandline_error:
_cfgfile_error:
    fprintf(stderr, "Exiting.\n");

    /* kill all threads */
    pthread_exit(NULL);
    return 0;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void *mpm_run_until(void *_task)
{
    threadtask_t *task = (threadtask_t *)_task;
    job_t *job = task->job;
    long ns = 0;
    int rc = 0;

    fprintf(stderr, "Starting thread: id=%zu, offset=%zu, blocksize=%zu, noff=%zu, nblk=%zu, eoff=%zu, eblk=%zu.\n",
        task->id, task->offset, task->blocksize,
        task->n_offset, task->n_blocksize,
        task->e_offset, task->e_blocksize);

    while (job->t < job->t_stop && !want_sigterm) {
        explicit_mpm_step_usl_threaded(task);

        /* have one thread write out the file */
        rc = pthread_barrier_wait(job->serialize_barrier);
        if (rc == PTHREAD_BARRIER_SERIAL_THREAD) {
            job->stepcount++;

            if (job->t >= (job->frame / job->output.sample_rate_hz)) {
/*                v2_write_frame(job->output.directory, NULL, job, v2_write_particle, NULL);*/
                write_frame(job->output.particle_fd, job->frame, job->t, job);
                write_element_frame(job->output.element_fd, job->frame, job->t, job);

                job->frame++;
                clock_gettime(CLOCK_REALTIME, &(job->toc));
                ns = 1E9 * (job->toc.tv_sec - job->tic.tv_sec) + (job->toc.tv_nsec - job->tic.tv_nsec);
                printf("\rt = %05.3fs\t[%3d%%]\tframe took %05.3fs (%d steps)",
                    job->t,
                    (int)((100 * job->t) / (job->t_stop)),
                    ns / 1e9,
                    job->stepcount
                );
                job->stepcount = 0;
                memcpy(&(job->tic), &(job->toc), sizeof(struct timespec));
                fflush(stdout);
                fflush(job->output.log_fd);
            }

            time_varying_loads(job);
        }

        /* XXX: mpm_step must wait for threads BEFORE starting! */
    }

    return NULL;
}
/*----------------------------------------------------------------------------*/

