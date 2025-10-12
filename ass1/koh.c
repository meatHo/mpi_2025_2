//
// Created by 고기호 on 25. 10. 10.
//

#include "koh.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

bool send_dir[4] = {true,true,true,true};//상좌하우
bool recv_dir[4] = {true,true,true,true};//하우상좌

void exchange_bound(int rank, struct MPI_info info, double**T_old) {
    //pex, pey로 끝부분 탐지 모서리 탐지해서 그 부분은 그냥 패스 차피 0 들어있으니까
    // printf("exchange_bound\n");

    //단열
    insulate(info,T_old);

    MPI_Status st;
    double *message;
    //통신 시작
    if ((info.coord_px+info.coord_py)%2==0) {//짝수 먼저 recv 순서 반대로 하우상좌 뒤에 부터 좌상우하-------------------------------
        for (int i=3;i>=0;i--) {//recv
            if (!recv_dir[i]) continue;
            int origin_rank=get_origin(i,rank,info);
            int count;

            if (i%2!=0) {//상하 일떄 메시지 사이즈
                count = info.local_xsize;
                message=(double*)malloc(info.local_xsize*sizeof(double));
            }else {//좌우 일때 메시지 사이즈
                count = info.local_ysize;
                message=(double*)malloc(info.local_ysize*sizeof(double));
            }

            // printf("rank %d : waiting : origin rank %d count %d dir %d\n",rank, origin_rank,count,3-i);


            MPI_Recv(message, count, MPI_DOUBLE, origin_rank, 3-i, MPI_COMM_WORLD, &st);
            // printf("rank %d : received origin rank %d count %d dir %d\n",rank, origin_rank,count,3-i);

            read_message(i,info,T_old,message);
            free(message);
        }
        for (int i=0;i<4;i++) {//send
            if (!send_dir[i]) continue;
            int dest_rank=get_dest(i,rank,info);
            int count;

            if (i%2==0) {//상하 일떄 메시지 사이즈
                count = info.local_xsize;
                message=(double*)malloc(info.local_xsize*sizeof(double));
            }else {//좌우 일때 메시지 사이즈
                count = info.local_ysize;
                message=(double*)malloc(info.local_ysize*sizeof(double));
            }

            make_message(i,info,T_old,message);
            // printf("rank %d sending : dest rank %d count %d dir %d\n",rank,dest_rank,count,i);
            MPI_Send(message,count,MPI_DOUBLE,dest_rank,i,MPI_COMM_WORLD);

            free (message);
        }
    }else {//홀수 먼저 send 상좌우하------------------------------------------------
        for (int i=0;i<4;i++) {//send
            if (!send_dir[i]) continue;
            int dest_rank=get_dest(i,rank,info);
            int count;

            if (i%2==0) {//상하 일떄 메시지 사이즈
                count = info.local_xsize;
                message=(double*)malloc(info.local_xsize*sizeof(double));
            }else {//좌우 일때 메시지 사이즈
                count = info.local_ysize;
                message=(double*)malloc(info.local_ysize*sizeof(double));
            }

            make_message(i,info,T_old,message);
            // printf("rank %d sending : dest rank %d count %d dir %d\n",rank,dest_rank,count,i);
            MPI_Send(message,count,MPI_DOUBLE,dest_rank,i,MPI_COMM_WORLD);

            free (message);
        }
        for (int i=3;i>=0;i--) {//recv
            if (!recv_dir[i]) continue;
            int origin_rank=get_origin(i,rank,info);
            int count;

            if (i%2!=0) {//상하 일떄 메시지 사이즈
                count = info.local_xsize;
                message=(double*)malloc(info.local_xsize*sizeof(double));
            }else {//좌우 일때 메시지 사이즈
                count = info.local_ysize;
                message=(double*)malloc(info.local_ysize*sizeof(double));
            }

            // printf("rank %d : waiting : origin rank %d count %d dir %d\n",rank, origin_rank,count,3-i);


            MPI_Recv(message, count, MPI_DOUBLE, origin_rank, 3-i, MPI_COMM_WORLD, &st);
            // printf("rank %d : received origin rank %d count %d dir %d\n",rank, origin_rank,count,3-i);

            read_message(i,info,T_old,message);
            free(message);
        }
    }

}

void make_message(int i, struct MPI_info info, double**T_old, double*message) {
    switch (i) {
        case 0://상
            for (int i=0;i<info.local_xsize;i++) {
                message[i]=T_old[i+1][1];
            }
            break;
        case 1://좌
            for (int i=0;i<info.local_ysize;i++) {
                message[i]=T_old[1][i+1];
            }
            break;
        case 2://하
            for (int i=0;i<info.local_xsize;i++) {
                message[i]=T_old[i+1][info.local_ysize];
            }
            break;
        case 3://우
            for (int i=0;i<info.local_ysize;i++) {
                message[i]=T_old[info.local_xsize][i+1];
            }
            break;
        default:
            printf("make_message error\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

void read_message(int i, struct MPI_info info, double **T_old, double *message) {//하우상좌
    switch (i) {
        case 0: //좌
            for (int p=1;p<=info.local_ysize;p++) {
                T_old[0][p]=message[p-1];
            }
            break;
        case 1: //상
            for (int p=1;p<=info.local_xsize;p++) {
                T_old[p][0]=message[p-1];
            }
            break;
        case 2: //우
            for (int p=1;p<=info.local_ysize;p++) {
                T_old[info.local_xsize+1][p]=message[p-1];
            }
            break;
        case 3: //하
            for (int p=1;p<=info.local_xsize;p++) {
                T_old[p][info.local_ysize+1]=message[p-1];
            }
            break;
        default:
            printf("get_dest error i %d 프로그램 종료\n",i);
            MPI_Abort(MPI_COMM_WORLD, 1);
            break;

    }
}

// void insulate(struct MPI_info info, double **T_old) {
//     int nx = info.local_xsize;
//     int ny = info.local_ysize;
//
//     if (!send_dir[0]) {//상
//         for (int j = 1; j <= ny; j++) {
//             T_old[0][j] = T_old[1][j];
//         }
//     }
//
//     if (!send_dir[1]) {//좌
//         for (int i = 1; i <= nx; i++) {
//             T_old[i][0] = T_old[i][1];
//         }
//     }
//
//     if (!send_dir[2]) {//하
//         for (int j = 1; j <= ny; j++) {
//             T_old[nx + 1][j] = T_old[nx][j];
//         }
//     }
//
//     if (!send_dir[3]) {//우
//         for (int i = 1; i <= nx; i++) {
//             T_old[i][ny + 1] = T_old[i][ny];
//         }
//     }
// }

void insulate(struct MPI_info info, double **T_old) {
    int nx = info.local_xsize;
    int ny = info.local_ysize;

    if (!send_dir[0]) {//상
        for (int j = 1; j <= nx; j++) {
            T_old[j][0] = T_old[j][1];
        }
    }

    if (!send_dir[1]) {//좌
        for (int i = 1; i <= ny; i++) {
            T_old[0][i] = T_old[1][i];
        }
    }

    if (!send_dir[2]) {//하
        for (int j = 1; j <= nx; j++) {
            T_old[j][nx + 1] = T_old[j][nx];
        }
    }

    if (!send_dir[3]) {//우
        for (int i = 1; i <= ny; i++) {
            T_old[ny + 1][i] = T_old[ny][i];
        }
    }
}

void check_bound(int rank, struct MPI_info info) {
    if (rank<info.px) {//상
        send_dir[0]=false;
        recv_dir[1]=false;
    }
    if (rank%info.px==0) {//좌
        send_dir[1]=false;
        recv_dir[0]=false;
    }
    if (rank+1>info.px*(info.py-1)) {//하
        // printf("이거동작?");
        send_dir[2]=false;
        recv_dir[3]=false;
    }
    if ((rank+1)%info.px==0) {//우
        send_dir[3]=false;
        recv_dir[2]=false;
    }
    // printf("rank %d, %d %d %d %d\n",rank, send_dir[0],send_dir[1],send_dir[2],send_dir[3]);
}

int get_dest(int i, int rank, struct MPI_info info) {
    switch (i) {
        case 0: //상
            return rank-info.px;
        case 1: //좌
            return rank-1;
        case 2: //하
            return rank+info.px;
        case 3: //우
            return rank+1;
        default:
            printf("get_dest error i %d 프로그램 종료\n",i);
            MPI_Abort(MPI_COMM_WORLD, 1);

    }
    return -1;
}

int get_origin(int i, int rank, struct MPI_info info) {
    switch (i) {
        case 0: //좌
            return rank-1;
        case 1: //상
            return rank-info.px;
        case 2: //우
            return rank+1;
        case 3: //하
            return rank+info.px;
        default:
            printf("get_origin error 프로그램 종료\n");
            MPI_Abort(MPI_COMM_WORLD, 1);

    }
    return -1;
}

void print_info(struct MPI_info info) {
    int x_size = info.px;
    int y_size = info.py;
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

void gather_total_heat(struct MPI_info info, double **T_new, int rank) {
    int nx = info.local_xsize;
    int ny = info.local_ysize;
    int size = info.px*info.py;
    MPI_Status st;
    double local_sum = 0.0;
    double max=-1;

    for (int i = 1; i <= nx; i++) {
        for (int j = 1; j <= ny; j++) {
            local_sum += T_new[i][j];
            max=(max>T_new[i][j])?max:T_new[i][j];
        }
    }

    if (rank == 0) {
        double global_sum = local_sum;
        // double recv_sum = 0.0;
        // double recv_max=0.0;

        // for (int src = 1; src < size; src++) {
        //     MPI_Recv(&recv_sum, 1, MPI_DOUBLE, src, 100, MPI_COMM_WORLD, &st);
        //     global_sum += recv_sum;
        // }
        // for (int src = 1; src < size; src++) {
        //     MPI_Recv(&recv_max, 1, MPI_DOUBLE, src, 1000, MPI_COMM_WORLD, &st);
        //     max=(max>recv_max)?max:recv_max;
        // }
        struct RES res;

        for (int src = 1; src < size; src++) {
            MPI_Recv(&res, 2, MPI_DOUBLE, src, 1000, MPI_COMM_WORLD, &st);
            global_sum += res.local_sum;
            max=(max>res.local_max)?max:res.local_max;
        }

        printf("\n rank 0 전체 열 합계 %f\n", global_sum);
        printf(" 전체 max %f\n",max);
    } else {
        // MPI_Send(&local_sum, 1, MPI_DOUBLE, 0, 100, MPI_COMM_WORLD);
        // MPI_Send(&max, 1, MPI_DOUBLE, 0, 1000, MPI_COMM_WORLD);
        struct RES send_res;
        send_res.local_sum = local_sum;
        send_res.local_max = max;
        MPI_Send(&max, 2, MPI_DOUBLE, 0, 1000, MPI_COMM_WORLD);
    }
}