#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>

#define VALUE_ERR 9969209968386869046778552952102584320.0
#define VALUE_NAN NAN

typedef struct variables
{
    char name[NC_MAX_NAME + 1];
    int id;
    float nvar_unavailable;
    nc_type type;
    void *data;
} variable;

typedef struct global_attributes
{
    // id ;
    // summary ;
    // title ;
    // time_coverage_start ;
    // time_coverage_final ;
    // station ;
    // longitude ;
    // latitude ;
} attribute;

typedef struct netcdf
{
    int ncid;
    int nvars;
    size_t dimension;
    attribute global_attributes;
    variable *variables;
} netcdf;

void *create_struct() {};
void *allocate_memory() {};
void desallocate_memory() {};
void read_file() {};
void read_header() {};
void read_body() {};

int main()
{
    printf("Loading...\n");
    return 0;
}

void *create_struct() {}
void *allocate_memory() {}
void desallocate_memory() {}
void read_file() {}
void read_header() {}
void read_body() {}