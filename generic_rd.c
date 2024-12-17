#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <gsl/gsl_interp.h>
#include <time.h>
#include <string.h>
#include <float.h>
#include <unistd.h>
#include <pthread.h>

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERRCODE 2
#define ERR(e)                                             \
    {                                                      \
        printf("Application Error: %s\n", nc_strerror(e)); \
        exit(ERRCODE);                                     \
    }

#define TYPE_VAR_NC_BYTE signed char
#define TYPE_VAR_NC_CHAR char
#define TYPE_VAR_NC_SHORT short
#define TYPE_VAR_NC_INT int
#define TYPE_VAR_NC_FLOAT float
#define TYPE_VAR_NC_DOUBLE double
#define TYPE_VAR_NC_UBYTE unsigned char
#define TYPE_VAR_NC_USHORT unsigned short
#define TYPE_VAR_NC_UINT unsigned int
#define TYPE_VAR_NC_INT64 long long
#define TYPE_VAR_NC_UINT64 unsigned long long
#define TYPE_VAR_NC_STRING char *

#define ALLOCATE_MEMORY(TYPE, LENGTH)                    \
    case TYPE:                                           \
        data = malloc(LENGTH * sizeof(TYPE_VAR_##TYPE)); \
        break;

#define PRINT_VALUES(TYPE, FMT)                                       \
    case TYPE:                                                        \
        for (int j = 0; j < file->var_len; j++)                       \
        {                                                             \
            printf(FMT "\n", ((TYPE_VAR_##TYPE *)data[i].values)[j]); \
        };                                                            \
        break;

#define PRINT_LN_VALUES(TYPE, FMT)                                \
    case TYPE:                                                    \
    {                                                             \
        printf(FMT "\t", ((TYPE_VAR_##TYPE *)data[j].values)[i]); \
    };                                                            \
    break;

#define W_VALUE 9969209968386869046778552952102584320.0
#define R_VALUE NAN

#define IF_CHANGE(TYPE)                                                         \
    case TYPE:                                                                  \
    {                                                                           \
        if (isnan((float)((TYPE_VAR_##TYPE *)data[i].values)[j]))               \
        {                                                                       \
            count++;                                                            \
            continue;                                                           \
        }                                                                       \
        if (((TYPE_VAR_##TYPE *)data[i].values)[j] == (TYPE_VAR_##TYPE)W_VALUE) \
        {                                                                       \
            ((TYPE_VAR_##TYPE *)data[i].values)[j] = (TYPE_VAR_##TYPE)R_VALUE;  \
            count++;                                                            \
        }                                                                       \
    };                                                                          \
    break;

#define K 20                     // Interval size
#define WINDOW (2 * K) + 1       // Window size
#define WINDOW_MINUS (2 * K) - 1 // Maximum size of interpolation range

#define IF_ISNAN(TYPE)                                                    \
    case TYPE:                                                            \
    {                                                                     \
        if (isnan((float)((TYPE_VAR_##TYPE *)data[i].values)[j]))         \
        {                                                                 \
            if (!d_time[0] && j != 0)                                     \
            {                                                             \
                d_time[0] = ((TYPE_VAR_##TYPE *)data[0].values)[j - 1];   \
                d_interp[0] = ((TYPE_VAR_##TYPE *)data[i].values)[j - 1]; \
                d_index[0] = j;                                           \
                count++;                                                  \
                continue;                                                 \
            }                                                             \
            count++;                                                      \
            continue;                                                     \
        }                                                                 \
        if (d_time[0] && !d_time[1])                                      \
        {                                                                 \
            if (count > WINDOW_MINUS)                                     \
            {                                                             \
                count = 0;                                                \
                d_time[0] = 0;                                            \
                d_interp[0] = 0;                                          \
                d_index[0] = 0;                                           \
                continue;                                                 \
            };                                                            \
            d_time[1] = ((TYPE_VAR_##TYPE *)data[0].values)[j];           \
            d_interp[1] = ((TYPE_VAR_##TYPE *)data[i].values)[j];         \
            d_index[1] = j;                                               \
        }                                                                 \
    };                                                                    \
    break;

#define INTERP_DATA(TYPE)                                                   \
    case TYPE:                                                              \
    {                                                                       \
        gsl_interp_init(interp, d_time, d_interp, 2);                       \
        for (int l = d_index[0]; l < d_index[1]; l++)                       \
        {                                                                   \
            if (isnan((float)((TYPE_VAR_##TYPE *)data[i].values)[l]))       \
            {                                                               \
                ((TYPE_VAR_##TYPE *)data[i].values)[l] =                    \
                    gsl_interp_eval(interp,                                 \
                                    d_time,                                 \
                                    d_interp,                               \
                                    ((TYPE_VAR_##TYPE *)data[0].values)[l], \
                                    acc);                                   \
            }                                                               \
        }                                                                   \
        d_changes += count;                                                 \
        count = 0;                                                          \
        d_time[0] = 0;                                                      \
        d_time[1] = 0;                                                      \
        d_interp[0] = 0;                                                    \
        d_interp[1] = 0;                                                    \
        d_index[0] = 0;                                                     \
        d_index[1] = 0;                                                     \
    };                                                                      \
    break;

// DATE_ISO - format "YYYY-MM-DDTHH:MM:00"
#define PREDICTION_INIT "2019-01-01T00:00:00"
#define PREDICTION_END "2019-12-31T23:54:00" //"2019-12-31T23:54:00"
#define TRAINING_INIT "2018-01-01T00:00:00"
#define TRAINING_END "2018-12-31T23:54:00"

#define MONACHE(TYPE)                                                     \
    case TYPE:                                                            \
    {                                                                     \
        for (int x = 0; x < WINDOW; x++)                                  \
        {                                                                 \
            double diff =                                                 \
                ((TYPE_VAR_##TYPE *)data[i].values)[(forecast - K) + x] - \
                ((TYPE_VAR_##TYPE *)data[i].values)[(analog - K) + x];    \
            sum += diff * diff;                                           \
        }                                                                 \
    };                                                                    \
    break;

#define N_NEIGHBORS 25
#define NUM_THREADS 4

typedef struct Kdtree
{
    unsigned int window_id;
    struct Kdtree *left;
    struct Kdtree *right;
} Kdtree;

// Definindo a estrutura para armazenar um ponto com sua distância
typedef struct
{
    unsigned int window_id;
    double distance;
} NearPoint;

typedef struct
{
    char field[NC_MAX_NAME + 1];
    int id;
    int nvar_available;
    int nvar_interp;
    float nvar_unavailable;
    nc_type type;
    void *values;
} VarData;

typedef struct
{
    int ncId;
    int nvars;
    size_t var_len;
    VarData *data;
} BaseApp;

typedef struct
{
    int p_init;
    int p_end;
    int t_init;
    int t_end;
    int n_pthreads;
    int t;
    int i;
    VarData *predictor_s;
    Kdtree *root;
} Params;

void *create_base_reconstruction(int, char *[]);
void read_header_data(BaseApp *, char *);
void read_body_data(BaseApp *, char *);
void *allocate_memory(char, size_t);
void deallocate_memory(BaseApp *, int);
void read_data(int, VarData *);
void print_io_netcdf_test(BaseApp *);
void print_available_data   (BaseApp *, int);
void print_data             (BaseApp *, int);
void count_data_unavailable (BaseApp *, int);
void interp_data            (BaseApp *, int);
time_t convert_iso_to_time(char *);
int binary_search           (BaseApp *, int);
double monacheMetric(VarData *, int, int, int, int);
void process_data_kdtree(BaseApp *, int, int, int);
void process_data_brute_force(BaseApp *, int, Params);
Kdtree *create_kdt_node(float);
Kdtree *insert_kdt_node(Kdtree *, float *, int, int);
void deallocate_kdtree(Kdtree *);
double euclidean_distance(Kdtree *, float *, int, int);
int compare_near_point(const void *, const void *);
void search_nearest_points(float *, Kdtree *, int, int, NearPoint *, int *);
double kmeans_euclidean_distance(float *, float *, int);
void kmeans_sort_centroids(float **, float *, int, int, NearPoint *);
void kmeans_find_closest_points(float **, float *, int *, int, NearPoint *, int, int);
void kmeans(float *, int, int, float **, int *, int);
void process_data_kmeans(BaseApp *, int, Params);
void *thread_brute_force_function(void *);
void process_data_brute_force_parallel(BaseApp *, int, Params);
void *thread_kdtree_function(void *);
void process_data_kdtree_parallel(BaseApp *, int, Params);

int main(int argc, char *argv[])
{
    time_t start, end;
    start = time(NULL);
    // double elapsed_time;

    BaseApp *data = create_base_reconstruction(argc, &argv[0]);
    Params params;

    params.p_init = binary_search(data, convert_iso_to_time(PREDICTION_INIT));
    params.p_end = binary_search(data, convert_iso_to_time(PREDICTION_END));
    params.t_init = binary_search(data, convert_iso_to_time(TRAINING_INIT));
    params.t_end = binary_search(data, convert_iso_to_time(TRAINING_END));
    params.n_pthreads = NUM_THREADS;

    // printf("prediction_init: %i\n", params.p_init);
    // printf("prediction_end: %i\n", params.p_end);
    // printf("training_init: %i\n", params.t_init);
    // printf("training_end: %i\n", params.t_end);
    // printf("%i,%i,%i,%i,%i,", params.n_pthreads, params.t_init, params.t_end, params.p_init, params.p_end);
    // print_io_netcdf_test(data);
    count_data_unavailable(data, argc);
    // print_available_data(data, argc);
    interp_data(data, argc);
    count_data_unavailable(data, argc);
    // print_available_data(data, argc);
    // process_data_kdtree(data, argc, params.p_init, params.p_end);
    process_data_kdtree_parallel(data, argc, params);
    // process_data_brute_force(data, argc, params);
    // process_data_brute_force_parallel(data, argc, params);
    // process_data_kmeans(data, argc, params);

    deallocate_memory(data, argc);
    end = time(NULL);
    // printf("Process %d: thread main: program completed: exiting\n", getpid());

    // printf("Tempo de execução: %.2f segundos\n", difftime(end, start));
    printf("%.0f", difftime(end, start));
    pthread_exit(NULL);
    // return 0;
}

void *create_base_reconstruction(int argc, char *argv[])
{
    // printf("Files Number: %i\n", (argc - 1));
    printf("%i,", (argc - 1));
    for (int i = 1; i < argc; i++)
    {
        // printf("File %i: %s\n", i, argv[i]);
    }

    BaseApp *files = (BaseApp *)malloc((argc - 1) * sizeof(BaseApp));
    if (!files)
        exit(EXIT_FAILURE);
    for (int i = 0; i < (argc - 1); i++)
    {
        read_header_data(&files[i], argv[i + 1]);
        read_body_data(&files[i], argv[i + 1]);
    }

    return files;
}

void read_header_data(BaseApp *file, char *argv)
{
    int retVal;

    /* Open the file. NC_NOWRITE tells netCDF we want read-only access
     * to the file.*/
    if ((retVal = nc_open(argv, NC_NOWRITE, &file->ncId)))
        ERR(retVal);

    /* Get base dimension for variables */
    if ((retVal = nc_inq_dimlen(file->ncId, 0, &file->var_len)))
        ERR(retVal)

    /* Get the number of variables in the file */
    if ((retVal = nc_inq_nvars(file->ncId, &file->nvars)))
        ERR(retVal);
}

void read_body_data(BaseApp *file, char *argv)
{
    int retVal;

    /* Allocate memory for variables names */
    VarData *varData = (VarData *)malloc(file->nvars * sizeof(VarData));

    /* List all variables */
    for (int i = 0; i < file->nvars; i++)
    {
        /* Get var name */
        if ((retVal = nc_inq_varname(file->ncId, i, varData[i].field)))
            ERR(retVal);
        /* Get id type var */
        if ((retVal = nc_inq_varid(file->ncId, varData[i].field, &varData[i].id)))
            ERR(retVal);
        /* Get type var */
        if ((retVal = nc_inq_vartype(file->ncId, i, &varData[i].type)))
            ERR(retVal);

        /* Allocate memory with specify dimension */
        varData[i].values = allocate_memory(varData[i].type, file->var_len);

        read_data(file->ncId, &varData[i]);
        // printf("vn: %s\t vt: %i\t vl: %li\t id: %i\n", varData[i].field, varData[i].type, file->var_len, varData[i].id);
    }
    file->data = varData;
    varData = NULL;
}

void *allocate_memory(char type, size_t length)
{
    void *data = NULL;
    switch (type)
    {
        ALLOCATE_MEMORY(NC_BYTE, length);
        ALLOCATE_MEMORY(NC_CHAR, length);
        ALLOCATE_MEMORY(NC_SHORT, length);
        ALLOCATE_MEMORY(NC_INT, length);
        ALLOCATE_MEMORY(NC_FLOAT, length);
        ALLOCATE_MEMORY(NC_DOUBLE, length);
        ALLOCATE_MEMORY(NC_UBYTE, length);
        ALLOCATE_MEMORY(NC_USHORT, length);
        ALLOCATE_MEMORY(NC_UINT, length);
        ALLOCATE_MEMORY(NC_INT64, length);
        ALLOCATE_MEMORY(NC_UINT64, length);
        ALLOCATE_MEMORY(NC_STRING, length);
    default:
        fprintf(stderr, "Error: Unknown data type.\n");
        exit(EXIT_FAILURE);
    }
    if (data == NULL)
    {
        fprintf(stderr, "Error: Failed to allocate memory.\n");
        exit(EXIT_FAILURE);
    }
    return data;
}

void deallocate_memory(BaseApp *file, int argc)
{
    int retVal;

    VarData *varData = NULL;
    for (int i = 0; i < (argc - 1); i++)
    {
        varData = file[i].data;
        for (int j = 0; j < file[i].nvars; j++)
        {
            free(varData[j].values);
            varData[j].values = NULL;
        }
        free(varData);
        varData = NULL;
        file[i].data = NULL;

        /* Close the file, freeing all resources. */
        if ((retVal = nc_close(file[i].ncId)))
            ERR(retVal);
    }
    free(file);
}

void read_data(int ncId, VarData *data)
{
    switch (data->type)
    {
    case NC_BYTE:
    case NC_CHAR:
        nc_get_var_text(ncId, data->id, data->values);
        break;
    case NC_SHORT:
        nc_get_var_short(ncId, data->id, data->values);
        break;
    case NC_INT:
        nc_get_var_int(ncId, data->id, data->values);
        break;
    case NC_FLOAT:
        nc_get_var_float(ncId, data->id, data->values);
        break;
    case NC_DOUBLE:
        nc_get_var_double(ncId, data->id, data->values);
        break;
    case NC_UBYTE:
        nc_get_var_uchar(ncId, data->id, data->values);
        break;
    case NC_USHORT:
        nc_get_var_ushort(ncId, data->id, data->values);
        break;
    case NC_UINT:
        nc_get_var_uint(ncId, data->id, data->values);
        break;
    case NC_INT64:
        nc_get_var_longlong(ncId, data->id, data->values);
        break;
    case NC_UINT64:
        nc_get_var_ulonglong(ncId, data->id, data->values);
        break;
    case NC_STRING:
        nc_get_var_string(ncId, data->id, data->values);
        break;
    default:
        fprintf(stderr, "Unknown variable type.\n");
        exit(1);
    }
}

void print_io_netcdf_test(BaseApp *file)
{
    if (file == NULL)
    {
        printf("No data to print for this file.\n");
        return;
    }

    VarData *data = file->data;

    for (int i = 0; i < file->nvars; i++)
    {
        if (data == NULL || data[i].values == NULL)
        {
            printf("No data to print.\n");
            return;
        }

        switch (data[i].type)
        {
            PRINT_VALUES(NC_BYTE, "%hhd");
            PRINT_VALUES(NC_CHAR, "%c");
            PRINT_VALUES(NC_SHORT, "%hd");
            PRINT_VALUES(NC_INT, "%i");
            PRINT_VALUES(NC_FLOAT, "%.1f");
            PRINT_VALUES(NC_DOUBLE, "%lf");
            PRINT_VALUES(NC_UBYTE, "%hhu");
            PRINT_VALUES(NC_USHORT, "%hu");
            PRINT_VALUES(NC_UINT, "%u");
            PRINT_VALUES(NC_INT64, "%lld");
            PRINT_VALUES(NC_UINT64, "%llu");
            PRINT_VALUES(NC_STRING, "%s");
        }
    }
}

void print_available_data(BaseApp *file, int argc)
{
    for (int f = 0; f < (argc - 1); f++)
    {
        if (file == NULL || file[f].data == NULL)
        {
            printf("No data to print for this file.\n");
            return;
        }

        printf("available data in file %i:\n", f + 1);

        for (int t = 0; t < file[f].nvars; t++)
        {
            t == 0 ? printf("\t\t")
                   : printf("%.2f%%\t", file[f].data[t].nvar_unavailable);
        }

        printf("\n");
    }
}

void print_data(BaseApp *file, int argc)
{
    VarData *data = NULL;
    for (int f = 0; f < (argc - 1); f++)
    {
        if (file == NULL || file[f].data == NULL)
        {
            printf("No data to print for this file.\n");
            return;
        }

        data = file[f].data;

        for (int i = 0; i < file[f].var_len; i++)
        // for (int i = 0; i < 1; i++)
        {

            if (i == 0)
            {
                printf("file data %i:\n", f + 1);

                for (int t = 0; t < file[f].nvars; t++)
                {
                    t == 0 ? printf("%s\t\t", file[f].data[t].field)
                           : printf("%s\t", file[f].data[t].field);
                }
                printf("\n");
            }

            for (int j = 0; j < file[f].nvars; j++)
            {
                switch (data[j].type)
                {
                    PRINT_LN_VALUES(NC_BYTE, "%hhd");
                    PRINT_LN_VALUES(NC_CHAR, "%c");
                    PRINT_LN_VALUES(NC_SHORT, "%hd");
                    PRINT_LN_VALUES(NC_INT, "%i");
                    PRINT_LN_VALUES(NC_FLOAT, "%.1f");
                    PRINT_LN_VALUES(NC_DOUBLE, "%lf");
                    PRINT_LN_VALUES(NC_UBYTE, "%hhu");
                    PRINT_LN_VALUES(NC_USHORT, "%hu");
                    PRINT_LN_VALUES(NC_UINT, "%u");
                    PRINT_LN_VALUES(NC_INT64, "%lld");
                    PRINT_LN_VALUES(NC_UINT64, "%llu");
                    PRINT_LN_VALUES(NC_STRING, "%s");
                }
            }
            printf("\n");
        }
    }
}

void count_data_unavailable(BaseApp *file, int argc)
{
    VarData *data = NULL;

    for (int f = 0; f < (argc - 1); f++)
    {
        if (file == NULL || file[f].data == NULL)
        {
            printf("No data to print for this file.\n");
            return;
        }

        data = file[f].data;

        for (int i = 0; i < file[f].nvars; i++)
        {
            int count = 0;

            if (data == NULL || data[i].values == NULL)
            {
                printf("No data %s to print.\n", data[i].field);
                return;
            }

            for (int j = 0; j < file[f].var_len; j++)
            // for (int j = 0; j < 1; j++)
            {

                switch (data[i].type)
                {
                    IF_CHANGE(NC_BYTE);
                    IF_CHANGE(NC_CHAR);
                    IF_CHANGE(NC_SHORT);
                    IF_CHANGE(NC_INT);
                    IF_CHANGE(NC_FLOAT);
                    IF_CHANGE(NC_DOUBLE);
                    IF_CHANGE(NC_UBYTE);
                    IF_CHANGE(NC_USHORT);
                    IF_CHANGE(NC_UINT);
                    IF_CHANGE(NC_INT64);
                    IF_CHANGE(NC_UINT64);
                }
            }
            
            if (!data[i].nvar_available)
            {
                data[i].nvar_available = file[f].var_len - count;
            }
            else
                data[i].nvar_interp = file[f].var_len - count;

            data[i].nvar_unavailable = 100 - (((float)count / (float)file[f].var_len) * 100);
        }
    }
}

void interp_data(BaseApp *file, int argc)
{
    VarData *data = NULL;
    int d_changes = 0;
    double d_time[2] = {0, 0};
    double d_interp[2] = {0, 0};
    int d_index[2] = {0, 0};

    // Initialize the interpolation objects
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_interp *interp = gsl_interp_alloc(gsl_interp_linear, 2);

    // TYPE_VAR_NC_FLOAT v_interp[5] = {0};
    for (int f = 0; f < (argc - 1); f++)
    {
        if (file == NULL || file[f].data == NULL)
        {
            printf("Err: interp_data: No data to print for this file.\n");
            return;
        }

        data = file[f].data;

        // printf("\t\t");
        for (int i = 1; i < file[f].nvars; i++)
        {
            int count = 0;
            if ((int)file[f].data[i].nvar_unavailable != 0 && (int)file[f].data[i].nvar_unavailable != 100)
            {
                for (int j = 0; j < file[f].var_len; j++)
                {
                    switch (data[i].type)
                    {
                        IF_ISNAN(NC_BYTE);
                        IF_ISNAN(NC_CHAR);
                        IF_ISNAN(NC_SHORT);
                        IF_ISNAN(NC_INT);
                        IF_ISNAN(NC_FLOAT);
                        IF_ISNAN(NC_DOUBLE);
                        IF_ISNAN(NC_UBYTE);
                        IF_ISNAN(NC_USHORT);
                        IF_ISNAN(NC_UINT);
                        IF_ISNAN(NC_INT64);
                        IF_ISNAN(NC_UINT64);
                    }
                    if (d_time[1])
                    {
                        switch (data[i].type)
                        {
                            INTERP_DATA(NC_BYTE);
                            INTERP_DATA(NC_CHAR);
                            INTERP_DATA(NC_SHORT);
                            INTERP_DATA(NC_INT);
                            INTERP_DATA(NC_FLOAT);
                            INTERP_DATA(NC_DOUBLE);
                            INTERP_DATA(NC_UBYTE);
                            INTERP_DATA(NC_USHORT);
                            INTERP_DATA(NC_UINT);
                            INTERP_DATA(NC_INT64);
                            INTERP_DATA(NC_UINT64);
                        }
                    }
                };
                // i == 0 ? printf("\t\t") : printf("%i\t", d_changes);
                d_changes = 0;
                continue;
            }
            // printf("0\t");
        }
        // printf("\n");
    }
    // free resorcer
    gsl_interp_free(interp);
    gsl_interp_accel_free(acc);
}

time_t convert_iso_to_time(char *rawtime)
{
    struct tm time_info;
    memset(&time_info, 0, sizeof(struct tm));

    // Convert an ISO 8601 string to struct tm
    if (sscanf(rawtime, "%4d-%2d-%2dT%2d:%2d:%2d",
               &time_info.tm_year, &time_info.tm_mon, &time_info.tm_mday,
               &time_info.tm_hour, &time_info.tm_min, &time_info.tm_sec) != 6)
    {
        printf("Error converting ISO 8601 string.\n");
        return -1; // Conversion error
    }

    // Adjust the values ​​of the tm struct for the mktime function
    time_info.tm_year -= 1900; // Year since 1900
    time_info.tm_mon -= 1;     // Months 0 to 11

    return mktime(&time_info) / 60; // Convert struct tm to time_t
}

int binary_search(BaseApp *file, int target)
{
    int *values = file[0].data[0].values;

    int init = 0, end = file[0].var_len - 1;

    while (init <= end)
    {
        int middle = init + (end - init) / 2;

        if (values[middle] == target)
            return middle;

        if (values[middle] < target)
            init = middle + 1;

        else
            end = middle - 1;
    }

    return -1;
}

double monacheMetric(VarData *data, int forecast, int analog, int t_init, int i)
{
    double sum = 0.0;
    int flag = 0;

    switch (data[i].type)
    {
        MONACHE(NC_BYTE);
        MONACHE(NC_CHAR);
        MONACHE(NC_SHORT);
        MONACHE(NC_INT);
        MONACHE(NC_FLOAT);
        MONACHE(NC_DOUBLE);
        MONACHE(NC_UBYTE);
        MONACHE(NC_USHORT);
        MONACHE(NC_UINT);
        MONACHE(NC_INT64);
        MONACHE(NC_UINT64);
    }

    return sqrt(sum);
}

// Function to create a new node in the Kdtree
Kdtree *create_kdt_node(float window_id)
{
    Kdtree *node = (Kdtree *)malloc(sizeof(Kdtree));
    node->window_id = window_id;
    node->left = NULL;
    node->right = NULL;
    return node;
}

// Function to insert a new point in the Kdtree
Kdtree *insert_kdt_node(Kdtree *root, float *values, int window_id, int depth)
{
    if (root == NULL)
        return create_kdt_node(window_id);

    int axis = depth % WINDOW; // Set axis based on depth

    // Compare if window axis is greater than stored value.
    if (values[(window_id - K) + axis] <
        values[(root->window_id - K) + axis])
        root->left = insert_kdt_node(root->left, values, window_id, depth + 1);
    else
        root->right = insert_kdt_node(root->right, values, window_id, depth + 1);

    return root;
}

// free memory kdtree
void deallocate_kdtree(Kdtree *root)
{
    if (root == NULL)
        return;

    deallocate_kdtree(root->left);
    deallocate_kdtree(root->right);
    free(root);
}

// Function to calculate euclidean distance (monache)
double euclidean_distance(Kdtree *root, float *values, int window_id, int target_id)
{
    double sum = 0.0;
    for (int i = 0; i < WINDOW; i++)
    {
        double diff =
            (float)values[(target_id - K) + i] -
            (float)values[(root->window_id - K) + i];
        sum += diff * diff;
    }

    return sqrt(sum);
}

// Function to compare distances for heap
int compare_near_point(const void *x, const void *y)
{
    NearPoint *point1 = (NearPoint *)x;
    NearPoint *point2 = (NearPoint *)y;

    if (point1->distance < point2->distance)
        return 1;
    if (point1->distance > point2->distance)
        return -1;
    return 0;
}

// Function to search the n nearest points
void search_nearest_points(float *values, Kdtree *root, int target_id, int depth, NearPoint *nearest, int *found)
{
    if (root == NULL)
        return;

    int axis = depth % WINDOW; // Set axis based on depth

    double distance = euclidean_distance(root, values, root->window_id, target_id);

    if (!isnan(distance))
    {
        if (*found < N_NEIGHBORS)
        {
            nearest[*found].window_id = root->window_id;
            nearest[*found].distance = distance;
            (*found)++;
            if (*found == N_NEIGHBORS)
                qsort(nearest, N_NEIGHBORS, sizeof(NearPoint), compare_near_point);
        }
        else if (distance < nearest[0].distance)
        {
            nearest[0].window_id = root->window_id;
            nearest[0].distance = distance;
            qsort(nearest, N_NEIGHBORS, sizeof(NearPoint), compare_near_point);
        }
    }

    Kdtree *next = NULL;
    Kdtree *other = NULL;
    // Compare if window axis is greater than stored value.
    if (values[(target_id - K) + axis] < values[(root->window_id - K) + axis])
    {
        next = root->left;
        other = root->right;
    }
    else
    {
        next = root->right;
        other = root->left;
    }

    // search_nearest_points(values, root->left, target_id, depth + 1, nearest, found);
    // search_nearest_points(values, root->right, target_id, depth + 1, nearest, found);
    search_nearest_points(values, next, target_id, depth + 1, nearest, found);

    if (*found < N_NEIGHBORS || fabs(values[(target_id - K) + axis] -
                                     values[(root->window_id - K) + axis]) < nearest[0].distance)
        search_nearest_points(values, other, target_id, depth + 1, nearest, found);
    // printf("distance: %.2lf\n", distance);
}

void process_data_kdtree(BaseApp *file, int argc, int p_init, int p_end)
{
    float *predictor_s = NULL;
    float *predicted_s = NULL;

    if (file == NULL || file[0].data == NULL)
    {
        printf("Err: process_data: No data to print for this file.\n");
        return;
    }

    for (int i = 1; i - 1 < file[0].nvars - 13; i++)
    {
        for (int f = 0; f < (argc - 1); f++)
        {
            if (f == 0)
            {
                if (predicted_s == NULL)
                {
                    predicted_s = file[f].data[i].values;
                    continue;
                }
                continue;
            }

            predictor_s = (float *)file[f].data[i].values;

            // printf("%i ", (int)file[f].data[i].nvar_unavailable);
            if ( //(int)file[f].data[i].nvar_unavailable != 0
                (int)file[f].data[i].nvar_unavailable >= 85 &&
                (int)file[f].data[i].nvar_unavailable != 100)
            {
                Kdtree *root = NULL;
                for (int analog = K; analog < p_init; analog++)
                {
                    root = insert_kdt_node(root, predictor_s, analog, 0);
                }

                for (int forecast = p_init; forecast <= p_end; forecast++)
                {
                    int found = 0;
                    NearPoint *nearest = (NearPoint *)malloc(N_NEIGHBORS * sizeof(NearPoint));
                    search_nearest_points(predictor_s, root, forecast, 0, nearest, &found);

                    // printf("------------------------------------------------------------------\n");
                    // printf("Posição: %i\n", forecast);

                    // for (int i = 0; i < found; i++)
                    // {
                    //     printf("| P%i D: %.2f ", i, nearest[i].distance);
                    // }
                    // printf("\n");

                    free(nearest);
                }

                deallocate_kdtree(root);
            }
            predicted_s = NULL;
            predictor_s = NULL;
        }
    }
}

void process_data_brute_force(BaseApp *file, int argc, Params params)
{
    VarData *predictor_s = NULL;
    VarData *predicted_s = NULL;

    if (file == NULL || file[0].data == NULL)
    {
        printf("Err: process_data: No data to print for this file.\n");
        return;
    }

    for (int i = 1; i - 1 < file[0].nvars - 13; i++)
    {
        for (int f = 0; f < (argc - 1); f++)
        {
            if (f == 0)
            {
                predicted_s = file[f].data;
                continue;
            }

            predictor_s = file[f].data;

            // printf("%i ", (int)file[f].data[i].nvar_unavailable);
            if ((int)file[f].data[i].nvar_unavailable >= 85 &&
                (int)file[f].data[i].nvar_unavailable != 100)
            {
                for (int forecast = params.p_init; forecast <= params.p_end; forecast++)
                {
                    int found = 0;
                    NearPoint *nearest = (NearPoint *)malloc(N_NEIGHBORS * sizeof(NearPoint));
                    for (int analog = params.t_init + K; analog < params.t_end; analog++)
                    {
                        double distance = monacheMetric(predictor_s,
                                                        forecast,
                                                        analog,
                                                        params.t_end,
                                                        i);

                        if (found < N_NEIGHBORS)
                        {
                            nearest[found].window_id = analog;
                            nearest[found].distance = distance;
                            found++;
                            if (found == N_NEIGHBORS)
                                qsort(nearest, N_NEIGHBORS, sizeof(NearPoint), compare_near_point);
                        }
                        else if (distance < nearest[0].distance)
                        {
                            nearest[0].window_id = analog;
                            nearest[0].distance = distance;
                            qsort(nearest, N_NEIGHBORS, sizeof(NearPoint), compare_near_point);
                        }
                    }
                    printf("------------------------------------------------------------------\n");
                    printf("Posição: %i\n", forecast);

                    for (int i = 0; i < found; i++)
                    {
                        printf("| P%i D: %.2f ", i, nearest[i].distance);
                    }
                    printf("\n");

                    free(nearest);
                }
            }
        }
    }
}

// Função para calcular a distância euclidiana entre dois vetores de tamanho `X`, ignorando valores NaN
double kmeans_euclidean_distance(float *a, float *b, int X)
{
    float sum = 0.0;
    int valid_values = 0;

    for (int i = 0; i < X; i++)
    {
        if (!isnan(a[i]) && !isnan(b[i]))
        { // Verifica se ambos os valores não são NaN
            sum += pow(a[i] - b[i], 2);
            valid_values++;
        }
    }

    // Se nenhum valor válido for encontrado, retorna uma distância "infinita"
    if (valid_values == 0)
    {
        return INFINITY;
    }

    return sqrt(sum);
}

// Função para ordenar centróides com base na distância ao ponto P
void kmeans_sort_centroids(float **centroids, float *values, int P, int num_clusters, NearPoint *distances)
{

    for (int i = 0; i < num_clusters; i++)
    {
        distances[i].window_id = i;
        distances[i].distance = kmeans_euclidean_distance(centroids[i], &values[i - (WINDOW - 1) / 2], WINDOW);
    }

    // Ordenação simples usando Bubble Sort
    for (int i = 0; i < num_clusters - 1; i++)
    {
        for (int j = 0; j < num_clusters - i - 1; j++)
        {
            if (distances[j].distance > distances[j + 1].distance)
            {
                NearPoint temp = distances[j];
                distances[j] = distances[j + 1];
                distances[j + 1] = temp;
            }
        }
    }
}

// Função para encontrar os n pontos mais próximos de todos centróides, ignorando NaN
void kmeans_find_closest_points(float **centroids, float *data, int *labels, int n, NearPoint *centroid_distances, int n_best, int num_clusters)
{
    // printf("\nOs %i pontos mais próximos (ignorando NaN):\n", n_best);

    // Calcular as distâncias dos pontos e armazenar
    NearPoint *points_distances = malloc(n * sizeof(NearPoint));
    int total_points = 0;

    for (int cluster_rank = 0; cluster_rank < num_clusters && total_points < n_best; cluster_rank++)
    {
        int closest_centroid_index = centroid_distances[cluster_rank].window_id;
        float *closest_centroid = centroids[closest_centroid_index];

        int count = 0;
        for (int i = 0; i < n; i++)
        {
            // Ignorar pontos com NaN e que não pertençam ao cluster atual
            if (labels[i] == closest_centroid_index)
            {
                int has_nan = 0;

                // Verifica se a janela centrada no ponto contém algum NaN
                for (int j = i - (WINDOW - 1) / 2; j <= i + (WINDOW - 1) / 2; j++)
                {
                    if (isnan(data[j]))
                    {
                        has_nan = 1;
                        break;
                    }
                }

                // Se a janela não tiver NaN, calcula a distância
                if (!has_nan)
                {
                    points_distances[total_points + count].window_id = i;
                    points_distances[total_points + count].distance = kmeans_euclidean_distance(&data[i - (WINDOW - 1) / 2], closest_centroid, WINDOW);
                    count++;
                }
            }
        }

        // Ordena os pontos deste cluster pelo mais próximo ao centróide
        for (int i = total_points; i < total_points + count - 1; i++)
        {
            for (int j = total_points; j < total_points + count - 1 - (i - total_points); j++)
            {
                if (points_distances[j].distance > points_distances[j + 1].distance)
                {
                    NearPoint temp = points_distances[j];
                    points_distances[j] = points_distances[j + 1];
                    points_distances[j + 1] = temp;
                }
            }
        }

        // Adiciona os pontos até atingir n_best ou o máximo de pontos disponíveis neste cluster
        total_points += count;
    }

    // printf("count: %i\n", total_points);
    // Exibe os n melhores pontos mais próximos
    for (int i = 0; i < n_best && i < total_points; i++)
    {
        printf("| P%i D: %.2f ", i, points_distances[i].distance);
        // printf("Ponto %i (Distância: %.2f): ", points_distances[i].window_id, points_distances[i].distance);
        // for (int j = 0; j < WINDOW; j++)
        // {
        //     printf("%.2f ", data[points_distances[i].window_id - (WINDOW - 1) / 2 + j]);
        // }
    }
    printf("\n");

    free(points_distances);
}

// Função para aplicar o K-means sem armazenar janelas, ignorando NaN
void kmeans(float *data, int n, int num_clusters, float **centroids, int *labels, int max_iters)
{
    int changed;
    float min_dist, dist;

    // Inicializa a janela deslizante com o primeiro ponto
    float window[WINDOW];
    int has_nan = 0;

    // Inicializa a janela centrada no ponto K
    for (int j = 0; j <= 2 * K; j++)
    {
        if (isnan(data[j]))
        {
            has_nan = 1;
            break;
        }
        window[j] = data[j]; // Inicializa a janela
    }

    for (int iter = 0; iter < max_iters; iter++)
    {
        changed = 0;

        // Atribui pontos aos clusters
        for (int i = K; i < n - K; i++)
        {
            // Atualiza a janela deslizante apenas para o novo ponto
            if (i > K)
            {
                // Remova o valor à esquerda e adicione o valor à direita
                if (!isnan(data[i - K - 1]) && !isnan(data[i + K]))
                {
                    // Shift para a esquerda
                    for (int j = 0; j < WINDOW - 1; j++)
                    {
                        window[j] = window[j + 1];
                    }
                    // Adiciona o novo valor à direita
                    window[WINDOW - 1] = data[i + K];
                }
                else
                {
                    // Se houver NaN, recria a janela completamente
                    has_nan = 0;
                    for (int j = i - K; j <= i + K; j++)
                    {
                        if (isnan(data[j]))
                        {
                            has_nan = 1;
                            break;
                        }
                        window[j - (i - K)] = data[j]; // Recria a janela
                    }
                    if (has_nan)
                    {
                        continue; // Ignora esta janela se houver NaN
                    }
                }
            }

            // Encontra o centróide mais próximo
            min_dist = kmeans_euclidean_distance(window, centroids[0], WINDOW);
            labels[i] = 0;
            for (int c = 1; c < num_clusters; c++)
            {
                dist = kmeans_euclidean_distance(window, centroids[c], WINDOW);
                if (dist < min_dist)
                {
                    min_dist = dist;
                    labels[i] = c;
                }
            }
        }

        // Atualiza os centróides
        for (int c = 0; c < num_clusters; c++)
        {
            int count = 0;
            float new_centroid[WINDOW] = {0};

            // Calcula a média para atualizar o centróide
            for (int i = K; i < n - K; i++)
            {
                if (labels[i] == c)
                {
                    int idx = 0;
                    int has_nan = 0;
                    for (int j = i - K; j <= i + K; j++)
                    {
                        if (isnan(data[j]))
                        {
                            has_nan = 1;
                            break;
                        }
                        new_centroid[idx++] += data[j];
                    }
                    if (has_nan)
                        continue;
                    count++;
                }
            }

            // Atualiza o centróide se houver pontos válidos
            if (count > 0)
            {
                for (int j = 0; j < WINDOW; j++)
                {
                    new_centroid[j] /= count;
                    if (fabs(new_centroid[j] - centroids[c][j]) > 1e-4)
                    {
                        changed = 1;
                    }
                    centroids[c][j] = new_centroid[j];
                }
            }
        }

        // Se nenhum centróide mudar, parar
        if (!changed)
            break;
    }
}

void process_data_kmeans(BaseApp *file, int argc, Params params)
{
    // clock_t start = clock();
    // clock_t end;
    // double elapsed_time;

    float *predictor_s = NULL;
    float *predicted_s = NULL;

    if (file == NULL || file[0].data == NULL)
    {
        printf("Err: process_data: No data to print for this file.\n");
        return;
    }

    for (int i = 1; i - 1 < file[0].nvars - 13; i++)
    {
        for (int f = 0; f < (argc - 1); f++)
        {
            int count = 0;
            if (f == 0)
            {
                if (predicted_s == NULL)
                {
                    predicted_s = file[f].data[i].values;
                    printf("line 1255 - f: %i\n", f);
                    continue;
                }
                printf("line 1258 - f: %i", f);
                continue;
            }

            predictor_s = (float *)file[f].data[i].values;

            // printf("%i ", (int)file[f].data[i].nvar_unavailable);
            if ((int)file[f].data[i].nvar_unavailable >= 85 &&
                (int)file[f].data[i].nvar_unavailable != 100)
            {
                int n = params.t_init; // - params.t_end;

                // int num_clusters = sqrt(file[f].data[i].nvar_available); // Número de clusters
                int num_clusters = sqrt((params.t_end - params.t_init)); // Número de clusters
                int n_best = N_NEIGHBORS;                                // Número de melhores pontos a serem retornados

                // Inicializando os centróides com os primeiros k pontos
                float **centroids = (float **)malloc(num_clusters * sizeof(float *));
                for (int i = 0; i < num_clusters; i++)
                {
                    centroids[i] = (float *)malloc((WINDOW) * sizeof(float));
                    for (int j = 0; j < WINDOW; j++)
                    {
                        centroids[i][j] = predictor_s[i + j];
                    }
                }

                // Rótulos dos clusters
                int *labels = (int *)malloc(n * sizeof(int));

                // clock_t start = clock();
                // clock_t end;
                // double elapsed_time;

                // Aplicando o K-means, ignorando NaN
                kmeans(predictor_s, n, num_clusters, centroids, labels, 100);

                // end = clock();
                // elapsed_time = (double)(end - start) / CLOCKS_PER_SEC;
                // printf("Tempo de execução: %.2f segundos\n", elapsed_time);

                // Distâncias dos centróides
                NearPoint *centroid_distances = (NearPoint *)malloc(num_clusters * sizeof(NearPoint));

                for (int forecast = params.p_init; forecast <= params.p_end /*params.p_init + 1*/; forecast++)
                {
                    // start = clock();
                    printf("Posição: %i\n", forecast);
                    // Ordenando os centróides com base na distância ao ponto P
                    kmeans_sort_centroids(centroids, predictor_s, forecast, num_clusters, centroid_distances);

                    // Encontrando os n melhores pontos mais próximos dos centróides, ignorando NaN
                    kmeans_find_closest_points(centroids, predictor_s, labels, n, centroid_distances, n_best, num_clusters);
                    // end = clock();
                    // elapsed_time = (double)(end - start) / CLOCKS_PER_SEC;
                    // printf("Tempo de execução: %.2f segundos\n", elapsed_time);
                }
                // Liberando memória
                for (int i = 0; i < num_clusters; i++)
                {
                    free(centroids[i]);
                }
                free(centroids);
                free(labels);
                free(centroid_distances);
            }
        }
    }
}

void *thread_brute_force_function(void *arg)
{
    Params *params = (Params *)arg;
    // printf("p_init = %i p_end = %i t_init = %i t_end = %i n_pthreads = %i t = %i i = %i\n",
    //        params->p_init,
    //        params->p_end,
    //        params->t_init,
    //        params->t_end,
    //        params->n_pthreads,
    //        params->t,
    //        params->i);

    for (int forecast = params->p_init + params->t; forecast <= params->p_end; forecast += params->n_pthreads)
    {
        int found = 0;
        NearPoint *nearest = (NearPoint *)malloc(N_NEIGHBORS * sizeof(NearPoint));
        for (int analog = params->t_init + K; analog < params->t_end; analog++)
        {
            double distance = monacheMetric(params->predictor_s,
                                            forecast,
                                            analog,
                                            params->t_end,
                                            params->i);

            if (found < N_NEIGHBORS)
            {
                nearest[found].window_id = analog;
                nearest[found].distance = distance;
                found++;
                if (found == N_NEIGHBORS)
                    qsort(nearest, N_NEIGHBORS, sizeof(NearPoint), compare_near_point);
            }
            else if (distance < nearest[0].distance)
            {
                nearest[0].window_id = analog;
                nearest[0].distance = distance;
                qsort(nearest, N_NEIGHBORS, sizeof(NearPoint), compare_near_point);
            }
        }
        // printf("------------------------------------------------------------------\n");
        // printf("Posição: %i\n", forecast);

        // for (int i = 0; i < found; i++)
        // {
        //     printf("| P%i D: %.2f ", i, nearest[i].distance);
        // }
        // printf("\n");

        free(nearest);
    }
    // printf("Process %d: thread %ld: Hello World !\n", getpid(), params->t);
    pthread_exit(NULL);
}

void process_data_brute_force_parallel(BaseApp *file, int argc, Params args)
{
    VarData *predictor_s = NULL;
    VarData *predicted_s = NULL;
    Params params[args.n_pthreads];
    time_t start, end;
    int sum_tsearche = 0;

    if (file == NULL || file[0].data == NULL)
    {
        printf("Err: process_data: No data to print for this file.\n");
        return;
    }

    for (int i = 1; i - 1 < file[0].nvars - 13; i++)
    {
        for (int f = 0; f < (argc - 1); f++)
        {
            if (f == 0)
            {
                predicted_s = file[f].data;
                continue;
            }

            predictor_s = file[f].data;

            // printf("%i ", (int)file[f].data[i].nvar_unavailable);
            if ((int)file[f].data[i].nvar_unavailable >= 85 &&
                (int)file[f].data[i].nvar_unavailable != 100)
            {
                pthread_t threads[NUM_THREADS];

                for (int t = 0; t < NUM_THREADS; t++)
                {
                    params[t].p_init = args.p_init;
                    params[t].p_end = args.p_end;
                    params[t].t_init = args.t_init;
                    params[t].t_end = args.t_end;
                    params[t].n_pthreads = NUM_THREADS;
                    params[t].t = t;
                    params[t].i = i;
                    params[t].predictor_s = predictor_s;
                }
                start = time(NULL);
                for (int t = 0; t < NUM_THREADS; t++)
                {
                    // printf("Process %d: thread main: creating thread %ld\n", getpid(), t);
                    pthread_create(&threads[t], NULL, thread_brute_force_function, (void *)&params[t]);
                }

                // Aguarda todas as threads terminarem
                for (int t = 0; t < NUM_THREADS; t++)
                {
                    pthread_join(threads[t], NULL);
                }
                end = time(NULL);
                printf("%.0f,", difftime(end, start));
                // printf("Process %d: thread main: program completed: exiting\n", getpid());
            }
        }
    }
}

void *thread_kdtree_function(void *arg)
{
    Params *params = (Params *)arg;
    // printf("p_init = %i p_end = %i t_init = %i t_end = %i n_pthreads = %i t = %i i = %i\n",
    //        params->p_init,
    //        params->p_end,
    //        params->t_init,
    //        params->t_end,
    //        params->n_pthreads,
    //        params->t,
    //        params->i);

    for (int forecast = params->p_init + params->t; forecast <= params->p_end; forecast += params->n_pthreads)
    {
        int found = 0;
        NearPoint *nearest = (NearPoint *)malloc(N_NEIGHBORS * sizeof(NearPoint));
        search_nearest_points((float *)params->predictor_s, (Kdtree *)params->root, forecast, 0, nearest, &found);

        // printf("------------------------------------------------------------------\n");
        // printf("Posição: %i\n", forecast);

        // for (int i = 0; i < found; i++)
        // {
        //     printf("| P%i D: %.2f ", i, nearest[i].distance);
        // }
        // printf("\n");

        free(nearest);
    }
    pthread_exit(NULL);
}

void process_data_kdtree_parallel(BaseApp *file, int argc, Params args)
{
    time_t start1, end1, start2, end2;
    float sum_tsearch = 0, sum_tprocess = 0;
    float *predictor_s = NULL;
    float *predicted_s = NULL;
    Params params[args.n_pthreads];

    if (file == NULL || file[0].data == NULL)
    {
        printf("Err: process_data: No data to print for this file.\n");
        return;
    }

    for (int i = 1; i - 1 < file[0].nvars - 13; i++)
    {
        for (int f = 0; f < (argc - 1); f++)
        {
            if (f == 0)
            {
                if (predicted_s == NULL)
                {
                    predicted_s = file[f].data[i].values;
                    continue;
                }
                continue;
            }

            predictor_s = (float *)file[f].data[i].values;

            // printf("%i ", (int)file[f].data[i].nvar_unavailable);
            if ( //(int)file[f].data[i].nvar_unavailable != 0
                (int)file[f].data[i].nvar_unavailable >= 85 &&
                (int)file[f].data[i].nvar_unavailable != 100)
            {
                pthread_t threads[NUM_THREADS];
                Kdtree *root = NULL;
                int count = 0;

                // printf("start = %i,", start1);
                start1 = time(NULL);
                for (int analog = args.t_init + K; analog < args.t_end; analog++)
                {
                    root = insert_kdt_node(root, predictor_s, analog, 0);
                    count++;
                }
                end1 = time(NULL);
                sum_tprocess += difftime(end1, start1);
                // printf("end = %i,", end1);
                // printf("count = %i\n", count);

                for (int t = 0; t < NUM_THREADS; t++)
                {
                    params[t].p_init = args.p_init;
                    params[t].p_end = args.p_end;
                    params[t].t_init = args.t_init;
                    params[t].t_end = args.t_end;
                    params[t].n_pthreads = NUM_THREADS;
                    params[t].t = t;
                    params[t].i = i;
                    params[t].predictor_s = (VarData *)predictor_s;
                    params[t].root = root;
                }

                start2 = time(NULL);
                for (int t = 0; t < NUM_THREADS; t++)
                {
                    // printf("Process %d: thread main: creating thread %ld\n", getpid(), t);
                    pthread_create(&threads[t], NULL, thread_kdtree_function, (void *)&params[t]);
                }

                // Aguarda todas as threads terminarem
                for (int t = 0; t < NUM_THREADS; t++)
                {
                    pthread_join(threads[t], NULL);
                }
                // printf("Process %d: thread main: program completed: exiting\n", getpid());

                end2 = time(NULL);
                deallocate_kdtree(root);
            }

            sum_tsearch += difftime(end2, start2);
            predicted_s = NULL;
            predictor_s = NULL;
        }
    }

    printf("%.0f,%.0f,", sum_tprocess, sum_tsearch);
    // %.0f,
}
