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


    MPI_Status st;
    double *message;
    //통신 시작
    if (info.coord_px+info.coord_py%2==0) {//짝수 먼저 recv 순서 반대로 하우상좌 뒤에 부터 좌상우하-------------------------------
        for (int i=3;i>=0;i--) {//recv
            if (!recv_dir[i]) continue;
            int origin_rank=get_origin(i,rank,info);
            int count;

            if (i%2==0) {//상하 일떄 메시지 사이즈
                count = info.local_xsize;
                message=(double*)malloc(info.local_xsize*sizeof(double));
            }else {//좌우 일때 메시지 사이즈
                count = info.local_ysize;
                message=(double*)malloc(info.local_ysize*sizeof(double));
            }

            // printf("rank %d : waiting : origin rank %d count %d dir %d\n",rank, origin_rank,count,3-i);


            MPI_Recv(message, count, MPI_DOUBLE, origin_rank, 3-i, MPI_COMM_WORLD, &st);
            printf("rank %d : received origin rank %d count %d dir %d\n",rank, origin_rank,count,3-i);

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
            printf("rank %d sending : dest rank %d count %d dir %d\n",rank,dest_rank,count,i);
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
            printf("rank %d sending : dest rank %d count %d dir %d\n",rank,dest_rank,count,i);
            MPI_Send(message,count,MPI_DOUBLE,dest_rank,i,MPI_COMM_WORLD);

            free (message);
        }
        for (int i=3;i>=0;i--) {//recv
            if (!recv_dir[i]) continue;
            int origin_rank=get_origin(i,rank,info);
            int count;

            if (i%2==0) {//상하 일떄 메시지 사이즈
                count = info.local_xsize;
                message=(double*)malloc(info.local_xsize*sizeof(double));
            }else {//좌우 일때 메시지 사이즈
                count = info.local_ysize;
                message=(double*)malloc(info.local_ysize*sizeof(double));
            }

            // printf("rank %d : waiting : origin rank %d count %d dir %d\n",rank, origin_rank,count,3-i);


            MPI_Recv(message, count, MPI_DOUBLE, origin_rank, 3-i, MPI_COMM_WORLD, &st);
            printf("rank %d : received origin rank %d count %d dir %d\n",rank, origin_rank,count,3-i);

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