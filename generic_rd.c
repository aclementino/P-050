#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <gsl/gsl_interp.h>
#include <time.h>
#include <string.h>

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

#define DATE_ISO_INIT "2017-01-01T00:00:00" // format "YYYY-MM-DDTHH:MM:00"
#define DATE_ISO_END "2017-12-31T00:00:00"  // format "YYYY-MM-DDTHH:MM:00"

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

void *create_base_reconstruction(int, char *[]);
void read_header_data(BaseApp *, char *);
void read_body_data(BaseApp *, char *);
void *allocate_memory(char, size_t);
void deallocate_memory(BaseApp *, int, char *[]);
void read_data(int, VarData *);
void print_io_netcdf_test(BaseApp *);
void print_available_data(BaseApp *, int);
void print_data(BaseApp *, int);
void count_data_unavailable(BaseApp *, int);
void interp_data(BaseApp *, int);
time_t convert_iso_to_time(char *);
int binary_search(BaseApp *, int);
double monacheMetric(VarData *, int, int, int, int);
void process_data_kdtree(BaseApp *, int, int, int);
void process_data_brute_force(BaseApp *, int, int, int);
Kdtree *create_kdt_node(float);
Kdtree *insert_kdt_node(Kdtree *, float *, int, int);
void deallocate_kdtree(Kdtree *);
double euclidean_distance(Kdtree *, float *, int, int);
int compare_near_point(const void *, const void *);
void search_nearest_points(float *, Kdtree *, int, int, NearPoint *, int *);

int main(int argc, char *argv[])
{
    int time_init = 0,
        time_end = 0;
    // time_init = convert_iso_to_time(DATE_ISO_INIT);
    // time_end = convert_iso_to_time(DATE_ISO_END);
    BaseApp *data = create_base_reconstruction(argc, &argv[0]);
    time_init = binary_search(data, convert_iso_to_time(DATE_ISO_INIT));
    time_end = binary_search(data, convert_iso_to_time(DATE_ISO_END));

    printf("time_init: %i\n", time_init);
    printf("time_end: %i\n", time_end);
    // print_io_netcdf_test(data);
    count_data_unavailable(data, argc);
    // print_available_data(data, argc);
    interp_data(data, argc);
    count_data_unavailable(data, argc);
    // process_data_kdtree(data, argc, time_init, time_end);
    // process_data_brute_force(data, argc, time_init, time_end);
    // print_available_data(data, argc);
    // print_data(data, argc);

    deallocate_memory(data, argc, &argv[0]);
    data = NULL;
    return 0;
}

void *create_base_reconstruction(int argc, char *argv[])
{
    printf("Files Number: %i\n", (argc - 1));

    for (int i = 1; i < argc; i++)
    {
        printf("File %i: %s\n", i, argv[i]);
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

void deallocate_memory(BaseApp *file, int argc, char *argv[])
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
        printf("Unknown variable type\n");
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

// Function ro search the n nearest points
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

void process_data_kdtree(BaseApp *file, int argc, int t_init, int t_end)
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
            int count = 0;
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
                for (int analog = K; analog < t_init; analog++)
                {
                    root = insert_kdt_node(root, predictor_s, analog, 0);
                }

                for (int forecast = t_init; forecast <= t_end /*t_init + 2000*/; forecast++)
                {
                    int found = 0;
                    NearPoint *nearest = (NearPoint *)malloc(N_NEIGHBORS * sizeof(NearPoint));
                    search_nearest_points(predictor_s, root, forecast, 0, nearest, &found);

                    printf("------------------------------------------------------------------\n");
                    printf("Posição: %i\n", forecast);

                    for (int i = 0; i < found; i++)
                    {
                        // printf("Ponto (");
                        // for (int j = (nearest[i].window_id - K); j < WINDOW; j++)
                        // {
                        //     printf("%.2f, ", (predictor_s)[j]);
                        //     if (j % 10 == 0)
                        //         printf("\n");
                        // }
                        // printf(") com distância %.2f\n", nearest[i].distance);
                        printf("| P%i W_ID: %i D: %.2f ", i, nearest[i].window_id, nearest[i].distance);
                    }
                    printf("\n");

                    // for (int i = 0; i < N_NEIGHBORS; i++)
                    // {
                    //     nearest[i].window_id = 0;
                    //     nearest[i].distance = 0.0;
                    // }
                    free(nearest);
                }

                deallocate_kdtree(root);
            }
            predicted_s = NULL;
            predictor_s = NULL;
        }
    }
}

void process_data_brute_force(BaseApp *file, int argc, int t_init, int t_end)
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
            int count = 0;
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
                for (int forecast = t_init; forecast <= t_end /*t_init + 2000*/; forecast++)
                {
                    int found = 0;
                    NearPoint *nearest = (NearPoint *)malloc(N_NEIGHBORS * sizeof(NearPoint));
                    for (int analog = K; analog < t_init; analog++)
                    {
                        double distance = monacheMetric(predictor_s,
                                                        forecast,
                                                        analog,
                                                        t_init,
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
                        printf("| P%i W_ID: %i D: %.2f ", i, nearest[i].window_id, nearest[i].distance);
                    }
                    printf("\n");

                    free(nearest);
                }
            }
        }
    }
}