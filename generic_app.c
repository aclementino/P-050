#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>

#define VALUE_ERR 9969209968386869046778552952102584320.0
#define VALUE_NAN NAN

typedef struct dimensions
{
    char name[NC_MAX_NAME + 1];
    size_t len;
} Dimension;

typedef struct attributes
{
    char name[NC_MAX_NAME + 1];
} Attribute;

typedef struct variables
{
    char name[NC_MAX_NAME + 1];
    nc_type type;
    int id;
    int natts;
    float nvar_unavailable;
    Attribute *atts;
    void *data;
} Variable;

typedef struct netcdf
{
    char file_name[NC_MAX_NAME + 1];
    int ncid_in;
    int ndims;
    int nvars;
    int natts;
    Dimension *dim;
    Attribute *g_atts; // global_attributes
    Variable *var;
} Netcdf;

void handle_error(int);
void *create_struct();
void *allocate_memory();
void desallocate_memory();
void read_file(Netcdf *, char *);
void read_dimensions(Netcdf *);
void read_variables(Netcdf *);
void read_body();

int main()
{
    printf("Loading...\n");
    return 0;
}

/*
 * Function to handle NetCDF errors, by printing an error message and exiting
 * with a non-zero status.
 */
void handle_error(int status)
{
    if (status != NC_NOERR)
    {
        fprintf(stderr, "%s\n", nc_strerror(status));
        exit(1);
    }
}

void *create_struct() {}
void *allocate_memory() {}
void desallocate_memory() {}
void read_file(Netcdf *file, char *argv)
{
    /*
     * Open the input NetCDF file.
     * NC_NOWRITE tells netCDF we want read-only access to the file.
     */
    handle_error(nc_open(argv, NC_NOWRITE, &file->ncid_in));

    /*
     * Get information from the input file
     */
    handle_error(nc_inq(file->ncid_in,
                        &file->ndims,
                        &file->nvars,
                        &file->natts,
                        NULL));

    /* Get base dimension for variables */
    read_dimensions(file);

    /* Get the number of variables in the file */
    // TODO
}
void read_dimensions(Netcdf *file)
{
    Dimension *dim = file->dim;
    /* Get dimensions */
    for (int i = 0; i < file->ndims; i++)
    {
        handle_error(nc_inq_dim(file->ncid_in, i, dim[i].name, &dim[i].len));
    }
}
void read_variable(Netcdf *file)
{
    Variable *var = file->var;
    // Get variables and yours attributes
    for (int i = 0; i < file->nvars; i++)
    {
        handle_error(nc_inq_var(file->ncid_in,
                                i,
                                var[i].name,
                                &var[i].type,
                                NULL,
                                NULL,
                                &var[i].natts));

        Attribute *atts = var[i].atts;
        // Get variable's attributes
        for (int j = 0; j < var[i].natts; j++)
        {
            char att_name[NC_MAX_NAME + 1];
            handle_error(nc_inq_attname(file->ncid_in, i, j, atts[j].name));
        }
    }
}
void read_body() {}