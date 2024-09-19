#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>

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

#define VALUE_ERR 9969209968386869046778552952102584320.0
#define VALUE_NAN NAN

#define K 20                       // Interval size
#define WIN_SIZE (2 * K) + 1       // Window size
#define WIN_SIZE_MINUS (2 * K) - 1 // Maximum size of interpolation range

typedef struct
{
    char name[NC_MAX_NAME + 1];
    size_t len;
} Dimension;

typedef struct
{
    char name[NC_MAX_NAME + 1];
} Attribute;

typedef struct
{
    char name[NC_MAX_NAME + 1];
    nc_type type;
    int id;
    int natts;
    float nvar_unavailable;
    Attribute *atts;
    void *data;
} Variable;

typedef struct
{
    char file_name[NC_MAX_NAME + 1];
    int ncid_in;
    int ndims;
    int nvars;
    int natts;
    Dimension *dim;
    Attribute *g_atts; // global_attributes
    Variable *var;
} NetCDF;

typedef struct
{
    unsigned int window_id;
    double distance;
} NearPoint;

typedef struct Kdtree
{
    unsigned int window_id;
    struct Kdtree *left;
    struct Kdtree *right;
} Kdtree;

void handle_error(int);
void *create_struct();
void *allocate_memory();
void desallocate_memory();
void read_file(NetCDF *, char *);
void read_dimensions(NetCDF *);
void read_variables(NetCDF *);
void read_body();
void write_file();

Kdtree *create_kdt_node(float);                         // Function to create a new node in the Kdtree
Kdtree *insert_kdt_node(Kdtree *, float *, int, int);   // Function to insert a new point in the Kdtree
void deallocate_kdtree(Kdtree *);                       // Free memory kdtree
double euclidean_distance(Kdtree *, float *, int, int); // Function to calculate euclidean distance (monache)
int compare_near_point(const void *, const void *);     // Function to compare distances for heap

int main(int argc, char *argv[])
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
void read_file(NetCDF *file, char *argv)
{
    /*
     * Open the input NetCDF file.
     * NC_NOWRITE tells NetCDF we want read-only access to the file.
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
void read_dimensions(NetCDF *file)
{
    Dimension *dim = file->dim;
    /* Get dimensions */
    for (int i = 0; i < file->ndims; i++)
    {
        handle_error(nc_inq_dim(file->ncid_in, i, dim[i].name, &dim[i].len));
    }
}
void read_variable(NetCDF *file)
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
void write_file() {}

Kdtree *create_kdt_node(float window_id)
{
    Kdtree *node = (Kdtree *)malloc(sizeof(Kdtree));
    node->window_id = window_id;
    node->left = NULL;
    node->right = NULL;
    return node;
}
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
void deallocate_kdtree(Kdtree *root)
{
    if (root == NULL)
        return;

    deallocate_kdtree(root->left);
    deallocate_kdtree(root->right);
    free(root);
}
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

    search_nearest_points(values, next, target_id, depth + 1, nearest, found);

    if (*found < N_NEIGHBORS || fabs(values[(target_id - K) + axis] -
                                     values[(root->window_id - K) + axis]) < nearest[0].distance)
        search_nearest_points(values, other, target_id, depth + 1, nearest, found);
}
