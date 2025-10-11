#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "koh.h"


//array extended version

int main(int argc, char *argv[]) {
    int local_xsize, local_ysize;
    int global_x, global_y;
    int steps = 50;
    int rank, size;
    double **T_old;
    double **T_new;
    static double alpha;
    static double g_sum = 0;
    static double g_max = 0;

    struct MPI_info info;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 6) {
        if (rank == 0)
            printf("사용법: mpirun -np 프로세스수 ./main NX NY PX PY ALPHA STEPS\n");
        return 1;
    }

    info.nx = atoi(argv[1]);
    info.ny = atoi(argv[2]);
    info.px = atoi(argv[3]);
    info.py = atoi(argv[4]);
    alpha = atof(argv[5]);
    steps = atoi(argv[6]);

    if (info.px * info.py != size && info.nx % info.px != 0 || info.ny % info.py != 0) {
        printf("프로세스 수가 일치하지 않습니다. 프로그램 종료합니다\n");
        return 0;
    }

    info.pex = info.nx / (info.nx / info.px); //pe num in x
    info.pey = info.ny / (info.ny / info.py); //pe num in y


    // 아니지 여기다가 if (마지막 줄)이면 배열 따로 설정하게 각 열, 행의 마지막 프로세스
    // 밑에 식 맞음 ㅇㅇ 일단은 이거 뺴고 진행
    if ((rank + 1) % info.pex == 0 && (rank + 1) > info.pex * (info.pey - 1)) {
    }

    local_xsize = info.nx / info.pex;
    local_ysize = info.ny / info.pey;

    global_x = local_xsize * (rank % info.px);
    global_y = local_ysize * (rank / info.px);

    // printf("rank %d global_x %d global_y %d\n", rank,global_x,global_y);


    T_old = (double **) malloc((local_xsize + 2) * sizeof(double *)); //양 옆의 다른 PE의 정보 저장하기 위해
    T_new = (double **) malloc((local_xsize + 2) * sizeof(double *));
    for (int i = 0; i < local_xsize + 2; i++) {
        T_old[i] = (double *) malloc((local_ysize + 2) * sizeof(double));
        T_new[i] = (double *) malloc((local_ysize + 2) * sizeof(double));
        for (int j = 0; j < local_ysize + 2; j++) {
            T_old[i][j] = 0.0;
            T_new[i][j] = 0.0;
        }
    }

    //발원 초기 위치 설정
    if (rank == 0) {
        //step 0
        T_new[1][1] = 100.0;
        T_old[1][1] = 100.0;
        // print_info(info);
    }


    int arr_x[4] = {-1, 1, 0, 0};
    int arr_y[4] = {0, 0, -1, 1};

    check_bound(rank, info); //본인 랭크의 끝인지 확인

    for (int i = 1; i < steps; i++) {
        //먼저 교환을 하고
        exchange_bound(rank, global_x, global_y, info, T_old);

        for (int x = 1; x < local_xsize + 1; x++) {
            for (int y = 1; y < local_ysize + 1; y++) {
                double arr[4] = {0}; //동서남북
                for (int t = 0; t < 4; t++) {
                    int local_tx = x + arr_x[t];
                    int local_ty = y + arr_y[t];
                    // //범위 벗어남 여기서 안함
                    // if (global_tx<0||global_ty<0||global_tx>=info.nx||global_ty>=info.ny) {
                    //     arr[t]=0;
                    //     printf("out of range\n");
                    //     continue;
                    // }
                    // //다른 프로세스로 부터 받아와야함
                    // if (local_tx<0||local_ty<0||local_tx>=local_xsize||local_ty>=local_ysize) {)

                    arr[t] = T_old[local_tx][local_ty];
                    // printf("local_tx %d local_ty %d global_tx %d global_ty %d arr[t] %f\n",local_tx,local_ty,global_tx,global_ty,arr[t]);
                }
                double sum = arr[0] + arr[1] + arr[2] + arr[3];
                T_new[x][y] = T_old[x][y] + (sum - 4 * T_old[x][y]) * alpha;
                // printf("x %d y %d old %f new %f sum %f alpha %f",x,y,T_old[x][y],T_new[x][y],sum,alpha);
            }
        }

        //T_old에 있는 복사본 초기화
        for (int j = 0; j < local_xsize + 2; j++) {
            T_old[j][0] = 0.0;
            T_old[j][local_ysize + 2 - 1] = 0.0;
        }

        for (int j = 0; j < local_ysize + 2; j++) {
            T_old[0][j] = 0.0;
            T_old[local_xsize + 2 - 1][j] = 0.0;
        }

        double **tmp = T_old;
        T_old = T_new;
        T_new = tmp;
    }

    //완
    for (int i = 0; i < local_xsize; i++) {
        free(T_old[i]);
        free(T_new[i]);
    }
    free(T_old);
    free(T_new);

    MPI_Finalize();

    return 0;
}
