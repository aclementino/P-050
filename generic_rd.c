#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <gsl/gsl_interp.h>

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

#define K 20                       // Interval size
#define WIN_SIZE (2 * K) + 1       // Window size
#define WIN_SIZE_MINUS (2 * K) - 1 // Maximum size of interpolation range

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
            if (count > WIN_SIZE_MINUS)                                   \
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

int main(int argc, char *argv[])
{
    BaseApp *data = create_base_reconstruction(argc, &argv[0]);

    // print_io_netcdf_test(data);
    count_data_unavailable(data, argc);
    print_available_data(data, argc);
    interp_data(data, argc);
    count_data_unavailable(data, argc);
    print_available_data(data, argc);
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