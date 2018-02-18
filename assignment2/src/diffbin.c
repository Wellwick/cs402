#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <errno.h>
#include <math.h>

static void print_usage(void);
static void print_version(void);
static void print_help(void);

static char *progname;

#define PACKAGE "diffbin"
#define VERSION "1.0"

/* Command line options */
static struct option long_opts[] = {
    { "help",    0, NULL, 'h' },
    { "version", 0, NULL, 'V' },
    { "epsilon", 1, NULL, 'e' },
    { "mode",    1, NULL, 'm' },
    { 0,           0, 0,   0  }
};

#define GETOPTS "e:hm:V"

#define MODE_DIFF 0
#define MODE_OUTPUT_U 1
#define MODE_OUTPUT_V 2 
#define MODE_OUTPUT_P 3 
#define MODE_OUTPUT_FLAGS 4 

int main(int argc, char **argv)
{
    FILE *f1, *f2;
    int imax, jmax, i, j;

    float *u1, *u2, *v1, *v2, *p1, *p2;
    char *flags1, *flags2;
    float epsilon = 1e-7;
    int mode = MODE_DIFF;
    int show_help = 0, show_usage = 0, show_version = 0;
    progname = argv[0];

    int optc;
    while ((optc = getopt_long(argc, argv, GETOPTS, long_opts, NULL)) != -1) {
        switch (optc) {
            case 'h':
                show_help = 1;
                break;
            case 'V':
                show_version = 1;
                break;
            case 'e':
                epsilon = atof(optarg);
                break;
            case 'm':
                if (strcasecmp(optarg, "diff") == 0) {
                    mode = MODE_DIFF;
                } else if (strcasecmp(optarg, "plot-u") == 0) {
                    mode = MODE_OUTPUT_U;
                } else if (strcasecmp(optarg, "plot-v") == 0) {
                    mode = MODE_OUTPUT_V;
                } else if (strcasecmp(optarg, "plot-p") == 0) {
                    mode = MODE_OUTPUT_P;
                } else if (strcasecmp(optarg, "plot-flags") == 0) {
                    mode = MODE_OUTPUT_FLAGS;
                } else {
                    fprintf(stderr, "%s: Invalid mode '%s'\n", progname, optarg);
                    show_usage = 1;
                }
                break;
            default:
                show_usage = 1;
        }
    }

    if (show_version) {
        print_version();
        if (!show_help) {
            return 0;
        }
    }

    if (show_help) {
        print_help();
        return 0;
    }
    
    if (show_usage || optind != (argc - 2)) {
        print_usage();
        return 1;
    }


    if ((f1 = fopen(argv[optind], "rb"))  == NULL) {
        fprintf(stderr, "Could not open '%s': %s\n", argv[optind],
            strerror(errno));
        return 1;
    }
    if ((f2 = fopen(argv[optind+1], "rb"))  == NULL) {
        fprintf(stderr, "Could not open '%s': %s\n", argv[optind+1],
            strerror(errno));
        return 1;
    }

    fread(&imax, sizeof(int), 1, f1);
    fread(&jmax, sizeof(int), 1, f1);
    fread(&i, sizeof(int), 1, f2);
    fread(&j, sizeof(int), 1, f2);
    if (i != imax || j != jmax) {
        printf("Number of cells differ! (%dx%d vs %dx%d)\n", imax, jmax, i, j);
        return 1;
    }

    float xlength1, ylength1, xlength2, ylength2;
    fread(&xlength1, sizeof(float), 1, f1);
    fread(&ylength1, sizeof(float), 1, f1);
    fread(&xlength2, sizeof(float), 1, f2);
    fread(&ylength2, sizeof(float), 1, f2);
    if (xlength1 != xlength2 || ylength1 != ylength2) {
        printf("Image domain dimensions differ! (%gx%g vs %gx%g)\n",
            xlength1, ylength1, xlength2, ylength2);
        return 1;
    }

    u1 = malloc(sizeof(float) * (jmax + 2));
    u2 = malloc(sizeof(float) * (jmax + 2));
    v1 = malloc(sizeof(float) * (jmax + 2));
    v2 = malloc(sizeof(float) * (jmax + 2));
    p1 = malloc(sizeof(float) * (jmax + 2));
    p2 = malloc(sizeof(float) * (jmax + 2));
    flags1 = malloc(jmax + 2);
    flags2 = malloc(jmax + 2);
    if (!u1 || !u2 || !v1 || !v2 || !p1 || !p2 || !flags1 || !flags2) {
        fprintf(stderr, "Couldn't allocate enough memory.\n");
        return 1;
    }

    int diff_found = 0;
    for (i = 0; i < imax + 2 && !diff_found; i++) {
        fread(u1, sizeof(float), jmax + 2, f1);
        fread(v1, sizeof(float), jmax + 2, f1);
        fread(p1, sizeof(float), jmax + 2, f1);
        fread(flags1, 1, jmax + 2, f1);
        fread(u2, sizeof(float), jmax + 2, f2);
        fread(v2, sizeof(float), jmax + 2, f2);
        fread(p2, sizeof(float), jmax + 2, f2);
        fread(flags2, 1, jmax + 2, f2);
        for (j = 0; j < jmax + 2 && !diff_found; j++) {
            float du, dv, dp;
            int dflags;
            du = u1[j] - u2[j];
            dv = v1[j] - v2[j];
            dp = p1[j] - p2[j];
            dflags = flags1[j] - flags2[j];
            switch (mode) {
                case MODE_DIFF:
                    if (fabs(du) > epsilon || fabs(dv) > epsilon ||
                        fabs(dp) > epsilon || fabs(dflags) > epsilon) {
                        diff_found = 1; 
                    }
                    break;
                case MODE_OUTPUT_U:
                    printf("%g%c", du, (j==jmax+1)?'\n':' ');
                    break;
                case MODE_OUTPUT_V:
                    printf("%g%c", dv, (j==jmax+1)?'\n':' ');
                    break;
                case MODE_OUTPUT_P:
                    printf("%g%c", dp, (j==jmax+1)?'\n':' ');
                    break;
                case MODE_OUTPUT_FLAGS:
                    printf("%d%c", dflags, (j==jmax+1)?'\n':' ');
                    break;
            }
        }
    }
    if (diff_found) {
        printf("Files differ.\n");
        return 1;
    }
    if (mode == MODE_DIFF) {
        printf("Files identical.\n");
    }
    return 0;
}

static void print_usage(void)
{
    fprintf(stderr, "Try '%s --help' for more information.\n", progname);
}

static void print_version(void)
{
    fprintf(stderr, "%s %s\n", PACKAGE, VERSION);
}

static void print_help(void)
{
    fprintf(stderr, "%s. A utility to compare karman state files.\n\n",
        PACKAGE);
    fprintf(stderr, "Usage %s [OPTIONS] FILE1 FILE2\n\n", progname);
    fprintf(stderr, "  -h, --help            Print a summary of the options\n");
    fprintf(stderr, "  -V, --version         Print the version number\n");
    fprintf(stderr, "  -e, --epsilon=EPSILON Set epsilon: the maximum allowed difference\n");
    fprintf(stderr, "  -m, --mode=MODE       Set the mode, may be one of 'diff', 'plot-u',\n");
    fprintf(stderr, "                        'plot-v', 'plot-p', or 'plot-flags'. The plot\n");
    fprintf(stderr, "                        modes produce output ready to be used by the\n");
    fprintf(stderr, "                        'splot matrix' command in gnuplot.\n");
}

