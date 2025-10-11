#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>



struct MPI_info{
    int nx, ny, px, py;
};

struct Message {
    int x, y, temp;
};

void print_info(struct MPI_info info);

int main(int argc, char *argv[]) {
    int local_xsize, local_ysize;
    int global_x, global_y;
    int steps = 50;
    int rank, size;
    double**T_old;
    double**T_new;
    static double alpha;
    static double g_sum=0;
    static double g_max=0;

    struct MPI_info info;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 6) {
        if (rank == 0)
            printf("사용법: mpirun -np [프로세스수] ./main NX NY PX PY ALPHA STEPS\n");
        return 1;
    }

    info.nx = atoi(argv[1]);
    info.ny = atoi(argv[2]);
    info.px = atoi(argv[3]);
    info.py = atoi(argv[4]);
    alpha = atof(argv[5]);
    steps = atoi(argv[6]);

    if (info.px*info.py!=size) {
        printf("프로세스 수가 일치하지 않습니다. 프로그램 종료합니다\n");
        return 0;
    }

    local_xsize=info.nx/info.px;
    local_ysize=info.ny/info.py;

    T_old=(double**)malloc(local_xsize*sizeof(double*));
    T_new=(double**)malloc(local_xsize*sizeof(double*));
    for (int i = 0; i < local_xsize; i++) {
        T_old[i] = (double*)malloc(local_ysize*sizeof(double));
        T_new[i] = (double*)malloc(local_ysize*sizeof(double));
    }
    for (int i = 0; i < local_xsize; ++i) {
        for (int j = 0; j < local_ysize; ++j) {
            T_old[i][j] = 0.0;
            T_new[i][j] = 0.0;
        }
    }
    if (rank==0) {
        //step 0
        T_new[0][0]=100.0;
        T_old[0][0]=100.0;
        // print_info(info);
    }

    global_x=local_xsize*(rank%info.px);
    global_y=local_ysize*(rank/info.px);

    printf("rank %d global_x %d global_y %d\n", rank,global_x,global_y);

    int arr_x[4] =  {-1,1,0,0};
    int arr_y[4] =  {0,0,-1,1};

    for (int i = 1; i < steps; i++) {
        //먼저 교환을 하고 rank
        if (rank==0) {
            for (int x=0; x < local_xsize; x++) {
                for (int y=0; y < local_ysize; y++) {
                    double arr[4] = {0};//동서남북
                    for (int t=0; t < 4; t++) {
                        int local_tx = x+arr_x[t];
                        int global_tx=global_x+local_tx;
                        int local_ty = y+arr_y[t];
                        int global_ty=global_y+local_ty;
                        printf("local_tx %d local_ty %d global_tx %d global_ty %d\n",local_tx,local_ty,global_tx,global_ty);
                        //범위 벗어남
                        if (global_tx<0||global_ty<0||global_tx>=info.nx||global_ty>=info.ny) {
                            arr[t]=0;
                            continue;
                        }

                        //다른 프로세스로 부터 받아와야함
                        if (local_tx<0||local_ty<0||local_tx>=local_xsize||local_ty>=local_ysize) {

                        }else {//내 안에서 해결
                            arr[t]=T_old[local_tx][local_ty];
                        }

                    }
                    double sum = arr[0]+arr[1]+arr[2]+arr[3];
                    T_new[x][y]=T_old[x][y]+(sum-4*T_old[x][y])*alpha;
                    printf("x %d y %d new %f",x,y,T_new[x][y]);
                }
            }
        }else {

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

void print_info(struct MPI_info info) {
    int x_size = info.nx/info.px;
    int y_size = info.ny/info.py;
    for (int y=0; y<info.ny; y++) {
        if (y%y_size==0&&y!=0) {
            for (int i=0;i<info.nx+(info.nx/x_size);i++) {
                printf("-");
            }
            printf("\n");
        }

        for (int x=0; x<info.nx; x++) {
            if (x%x_size==0&&x!=0) {
                printf("|");
            }
            printf("0");
        }
        printf("\n");
    }
}