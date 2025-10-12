//
// Created by 고기호 on 25. 10. 10.
//

#ifndef KOH_H
#define KOH_H

#endif //KOH_H
#include <stdbool.h>

struct MPI_info {
    int nx, ny, px, py, local_xsize, local_ysize, coord_px, coord_py;
};

struct RES {
    double local_sum;
    double local_max;
};

void print_info(struct MPI_info info);

void check_bound(int rank, struct MPI_info info);

void exchange_bound(int rank, struct MPI_info info, double **T_old);

void exchange_res(int rank, struct MPI_info info, double **T_new);

void make_message(int i, struct MPI_info info, double **T_old, double *message);

void read_message(int i, struct MPI_info info, double **T_old, double *message);

void insulate(struct MPI_info info, double **T_old);//아 단열 까먹었네 뮤솜ㄴ아ㅣ러매ㅑㅔㄷ저래ㅑ머ㅐㅑㅔㄷ조ㅜ래ㅑㅔ뮤

void gather_total_heat(struct MPI_info info, double **T_new, int rank);

int get_dest(int i, int rank, struct MPI_info info);

int get_origin(int i, int rank, struct MPI_info info);
