/*
 * pagerank_skeleton_HW2.c
 * MPI Collective-based Parallel PageRank Skeleton Code for HW2
 *
 * 이 스켈레톤 코드는 여러분이 PageRank 알고리즘의 병렬화에 집중할 수 있도록
 * 파라미터 파싱, 데이터 파티셔닝 등 부가적인 기능들을 미리 구현해두었습니다.
 * 여러분은 main() 함수의 메인 반복문 내부에 있는 TODO 항목들을 채워
 * MPI 집단 통신을 이용한 병렬 PageRank를 완성해야 합니다.
 *
 * Build: mpicc -O2 -Wall -o pagerank pagerank_skeleton_HW2.c
 * Run (example):
 * mpirun -np 4 ./pagerank --graph ring 100000 --damp 0.85 --tol 1e-10 --maxit 100
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Parameters */
typedef struct {
    const char* gtype;    /* ring or star */
    long long N;          /* number of nodes (global) */
    double    damp;       /* damping factor, default 0.85 */
    double    tol;        /* L1 convergence criterion */
    int       maxit;      /* maximum iterations */
} Params;

/* Parameter parsing */
static void parse_args(int rank, int argc, char** argv, Params* p) {
    /* Default values */
    p->gtype = "ring";
    p->N = 100000;
    p->damp = 0.85;
    p->tol = 1e-10;
    p->maxit = 100;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--graph") && i + 2 < argc) {
            if (!strcmp(argv[i+1], "ring")) p->gtype = "ring";
            else if (!strcmp(argv[i+1], "star")) p->gtype = "star";
            p->N = atoll(argv[i+2]);
            i += 2;
        } else if (!strcmp(argv[i], "--damp") && i + 1 < argc) {
            p->damp = atof(argv[++i]);
        } else if (!strcmp(argv[i], "--tol") && i + 1 < argc) {
            p->tol = atof(argv[++i]);
        } else if (!strcmp(argv[i], "--maxit") && i + 1 < argc) {
            p->maxit = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--help")) {
            if (rank == 0) {
                printf("Usage: mpirun -np P ./pagerank --graph (ring|star) N [--damp d] [--tol t] [--maxit m]\n");
            }
            MPI_Abort(MPI_COMM_WORLD, 0);
        }
    }
}

/* Graph outdegree vector generation (ring/star) */
static void build_outdegree(const Params* P, double* outdeg) {
    long long N = P->N;
    if (strcmp(P->gtype, "ring") == 0) {
        for (long long u = 0; u < N; ++u) outdeg[u] = 1.0;
    } else { /* STAR */
        outdeg[0] = 0.0; /* Hub node is dangling */
        for (long long u = 1; u < N; ++u) outdeg[u] = 1.0;
    }
}

/* Partitioning (contiguous block decomposition) */
static void partition(long long N, int size, int rank, long long* locN, long long* off, int** counts, int** displacements) {
    long long base = N / size;
    long long rem = N % size;
    *locN = base + (rank < rem? 1 : 0);
    *off = (rank < rem)? (rank * (base + 1)) : (rem * (base + 1) + (rank - rem) * base);

    if (counts && displacements) {
        *counts = (int*)malloc(sizeof(int) * size);
        *displacements = (int*)malloc(sizeof(int) * size);
        for (int r = 0; r < size; r++) {
            (*counts)[r] = base + (r < rem? 1 : 0);
        }
        int acc = 0;
        for (int r = 0; r < size; r++) {
            (*displacements)[r] = acc;
            acc += (*counts)[r];
        }
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    Params P;
    parse_args(rank, argc, argv, &P);

    /* Partitioning information */
    long long locN = 0, off = 0;
    int *counts = NULL, *displacements = NULL;
    partition(P.N, size, rank, &locN, &off, &counts, &displacements);

    /* PageRank vectors (globally replicated) */
    double *pagerank_global = (double*)malloc(sizeof(double) * P.N);
    double *pagerank_new_local_chunk = (double*)malloc(sizeof(double) * locN);
    /* Local chunk to send with Allgatherv */

    /* Initialization: pr = 1/N */
    for (long long i = 0; i < P.N; i++) {
        pagerank_global[i] = 1.0 / (double)P.N;
    }

    /* Outdegree vector (globally replicated) */
    double *outdeg = (double*)malloc(sizeof(double) * P.N);
    build_outdegree(&P, outdeg);

    /* Iteration */
    const double d = P.damp;
    const double one_over_N = 1.0 / (double)P.N;
    int it = 0;
    int converged = 0;
    double global_l1_diff = 0.0;

    for (it = 0; it < P.maxit &&!converged; ++it) {
        /* (1) Calculate dangling mass = sum_{u where outdeg(u)==0} pagerank(u) */
        double local_dang_mass = 0.0;
        /*
         * TODO: local node들 중 out-degree가 0인 노드를 찾아
         *       해당 node들의 pagerank_global 값을 local_dang_mass와 합산
         */
        for (long long u = off; u < off + locN; ++u){
            if (outdeg[u] == 0.0){
                local_dang_mass += pagerank_global[u];
            }
        }
        double global_dang_mass = 0.0;
        /*
         * TODO: local_dang_mass를 합산 & 그 결과를 공유
         */
        MPI_Allreduce(&local_dang_mass, &global_dang_mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        /* (2) Calculate total sum (for verification and STAR graph optimization) */
        double local_sum = 0.0;
        /*
         * TODO: local_sum 계산
         */
        for (long long u = off; u < off + locN; ++u){
            local_sum += pagerank_global[u];
        }

        double global_sum = 0.0;
        /*
         * TODO: global_sum 계산
         */
        MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        /* (3) Calculate new PageRank values for the local chunk */
        double local_l1_diff = 0.0;
        /*
         * TODO: PageRank 계산 반복
         */
        for (long long i = 0; i < locN; ++i){
            long long u = off + i;
            double new_pr = 0.0;

            if (strcmp(P.gtype, "ring") == 0){
                long long v = (u - 1 + P.N) % P.N;
                new_pr = (1.0 - d) * one_over_N
                       + d * (pagerank_global[v] / outdeg[v])
                       + d * (global_dang_mass / P.N);
            }
            else if (strcmp(P.gtype, "star") == 0) {
                if (u == 0) {
                    double hub_gain = global_sum - pagerank_global[0];
                    double weighted_contrib = d * hub_gain;
                    new_pr = (1.0 - d) * one_over_N + d * (global_dang_mass / P.N) + weighted_contrib;
                } else {
                    // 리프 노드들: 허브로만 보냄 → in-degree 없음
                    new_pr = (1.0 - d) * one_over_N + d * (global_dang_mass / P.N);
                }
            }
            local_l1_diff += fabs(new_pr - pagerank_global[u]);
            pagerank_new_local_chunk[i] = new_pr;
        }

        /* (4) Check for convergence: L1 norm of two vectors (pagerank_new - pagerank_old) */
        MPI_Allreduce(&local_l1_diff, &global_l1_diff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        if (global_l1_diff < P.tol) converged = 1;

        /* (5) Replicate the new PageRank vector to all processes using Allgatherv */
        /*
         * TODO: pagerank_global 벡터 갱신
         */
        MPI_Allgatherv(
            pagerank_new_local_chunk,   // 보낼 데이터 (각 rank의 로컬 PR 부분)
            locN, MPI_DOUBLE,           // 각 rank가 보낼 데이터 크기
            pagerank_global,            // 모든 rank가 받을 버퍼 (전역 벡터)
            counts, displacements,      // 각 rank별 데이터 개수, 시작 위치
            MPI_DOUBLE, MPI_COMM_WORLD  // 자료형, 통신 그룹
        );
    }

    /* Print results (from rank 0) */
    if (rank == 0) {
        printf("[Result] Graph=%s, N=%lld, Procs=%d, Iters=%d, L1_Residual=%.3e, Damp=%.2f\n",
            P.gtype, P.N, size, it, global_l1_diff, P.damp);

        /* Print top-k (for large N, use a top-k selection algorithm) */
        int k = 100;
        // Adjust k if N is smaller than 100
        if (k > (int)P.N) k = (int)P.N;

        int *idx = (int*)malloc(sizeof(int) * P.N);
        for (long long i = 0; i < P.N; i++) idx[i] = (int)i;

        /* Partial selection sort (k passes) */
        for (int t = 0; t < k; t++) {
            int best = t;
            for (long long j = t + 1; j < P.N; j++) {
                if (pagerank_global[idx[j]] > pagerank_global[idx[best]]) best = j;
            }
            int tmp = idx[t];
            idx[t] = idx[best];
            idx[best] = tmp;
        }

        printf("Top %d nodes:\n", k);
        for (int t = 0; t < k; t++) {
            int u = idx[t];
            printf("  rank %3d: node %d, PR=%.6e\n", t + 1, u, pagerank_global[u]);
        }
        free(idx);
    }

    free(outdeg);
    free(pagerank_new_local_chunk);
    free(pagerank_global);
    free(counts);
    free(displacements);
    MPI_Finalize();
    return 0;
}