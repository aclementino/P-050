#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>

#define FILE_NAME_IN "wdsv2_2011_2019.nc"  // Nome do arquivo de entrada
#define FILE_NAME_OUT "output.nc" // Nome do arquivo de saída

// Função para lidar com erros do NetCDF
void handle_error(int status) {
    if (status != NC_NOERR) {
        fprintf(stderr, "%s\n", nc_strerror(status));
        exit(1);
    }
}

int main() {
    int ncid_in, ncid_out;
    int ndims, nvars, natts, unlimdimid;

    // Abrindo o arquivo NetCDF de entrada
    handle_error(nc_open(FILE_NAME_IN, NC_NOWRITE, &ncid_in));

    // Obtendo informações do arquivo de entrada
    handle_error(nc_inq(ncid_in, &ndims, &nvars, &natts, &unlimdimid));

    // Criando o arquivo NetCDF de saída
    handle_error(nc_create(FILE_NAME_OUT, NC_CLOBBER, &ncid_out));

    // Copiando as dimensões
    int dimids[NC_MAX_DIMS];
    for (int i = 0; i < ndims; i++) {
        char dim_name[NC_MAX_NAME + 1];
        size_t dimlen;
        handle_error(nc_inq_dim(ncid_in, i, dim_name, &dimlen));
        handle_error(nc_def_dim(ncid_out, dim_name, dimlen, &dimids[i]));
    }

    // Copiando as variáveis e seus atributos
    for (int i = 0; i < nvars; i++) {
        char var_name[NC_MAX_NAME + 1];
        nc_type var_type;
        int var_ndims, var_natts, var_dimids[NC_MAX_VAR_DIMS];

        handle_error(nc_inq_var(ncid_in, i, var_name, &var_type, &var_ndims, var_dimids, &var_natts));
        handle_error(nc_def_var(ncid_out, var_name, var_type, var_ndims, var_dimids, NULL));

        // Copiando os atributos da variável
        for (int j = 0; j < var_natts; j++) {
            char att_name[NC_MAX_NAME + 1];
            handle_error(nc_inq_attname(ncid_in, i, j, att_name));
            handle_error(nc_copy_att(ncid_in, i, att_name, ncid_out, i));
        }
    }

    // Copiando os atributos globais
    for (int i = 0; i < natts; i++) {
        char att_name[NC_MAX_NAME + 1];
        handle_error(nc_inq_attname(ncid_in, NC_GLOBAL, i, att_name));
        handle_error(nc_copy_att(ncid_in, NC_GLOBAL, att_name, ncid_out, NC_GLOBAL));
    }

    // Finalizando a definição do arquivo de saída
    handle_error(nc_enddef(ncid_out));

    // Copiando os dados das variáveis
    for (int i = 0; i < nvars; i++) {
        char var_name[NC_MAX_NAME + 1];
        nc_type var_type;
        int var_ndims;
        int var_dimids[NC_MAX_VAR_DIMS];
        size_t var_dimlen[NC_MAX_VAR_DIMS];

        handle_error(nc_inq_var(ncid_in, i, var_name, &var_type, &var_ndims, var_dimids, NULL));

        // Obtendo o tamanho das dimensões de cada variável
        for (int j = 0; j < var_ndims; j++) {
            handle_error(nc_inq_dimlen(ncid_in, var_dimids[j], &var_dimlen[j]));
        }

        // Alocando memória para os dados e copiando
        switch (var_type) {
            case NC_BYTE: {
                signed char *data = (signed char *)malloc(sizeof(signed char) * var_dimlen[0]);
                handle_error(nc_get_var(ncid_in, i, data));
                handle_error(nc_put_var(ncid_out, i, data));
                free(data);
                break;
            }
            case NC_CHAR: {
                char *data = (char *)malloc(sizeof(char) * var_dimlen[0]);
                handle_error(nc_get_var(ncid_in, i, data));
                handle_error(nc_put_var(ncid_out, i, data));
                free(data);
                break;
            }
            case NC_SHORT: {
                short *data = (short *)malloc(sizeof(short) * var_dimlen[0]);
                handle_error(nc_get_var(ncid_in, i, data));
                handle_error(nc_put_var(ncid_out, i, data));
                free(data);
                break;
            }
            case NC_INT: {
                int *data = (int *)malloc(sizeof(int) * var_dimlen[0]);
                handle_error(nc_get_var(ncid_in, i, data));
                handle_error(nc_put_var(ncid_out, i, data));
                free(data);
                break;
            }
            case NC_FLOAT: {
                float *data = (float *)malloc(sizeof(float) * var_dimlen[0]);
                handle_error(nc_get_var(ncid_in, i, data));
                handle_error(nc_put_var(ncid_out, i, data));
                free(data);
                break;
            }
            case NC_DOUBLE: {
                double *data = (double *)malloc(sizeof(double) * var_dimlen[0]);
                handle_error(nc_get_var(ncid_in, i, data));
                handle_error(nc_put_var(ncid_out, i, data));
                free(data);
                break;
            }
            default:
                fprintf(stderr, "Tipo de variável não suportado.\n");
                exit(1);
        }
    }

    // Fechando os arquivos NetCDF
    handle_error(nc_close(ncid_in));
    handle_error(nc_close(ncid_out));

    printf("Arquivo NetCDF copiado com sucesso para %s!\n", FILE_NAME_OUT);

    return 0;
}