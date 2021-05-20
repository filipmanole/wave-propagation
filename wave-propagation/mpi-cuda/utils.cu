/*
 * Student:	Trascau Mihai
 * Grupa:	344C4
 * 
 * Lucrare:	Ecuatia undelor pentru acustica 2D
 * Fisier:	utils.c
 * Descriere:	Fisier sursa in care sunt implementate cateva functii utile aplicatiei
 */

#include "acoustics.h"

void load_scenario()
{
    nx = scenario[scn_index].nx;
    ny = scenario[scn_index].ny;
    H = scenario[scn_index].H;
    MAX_TIME = scenario[scn_index].MAX_TIME;
    TIME_STEP = scenario[scn_index].TIME_STEP;
    SAVE_TIME = scenario[scn_index].SAVE_TIME;
}

void unload_scenario_gpu(double *cudaUa, double *cudaUb, double *cudaUc)
{
    nx = ny = H = MAX_TIME = TIME_STEP = SAVE_TIME = 0;

    cudaFree(cudaUa);
    cudaFree(cudaUb);
    cudaFree(cudaUc);
}

void unload_scenario()
{
    nx = ny = H = MAX_TIME = TIME_STEP = SAVE_TIME = 0;
    free(uc);
    free(ub);
    free(ua);
}

__global__ void cudaInitVectors(int rank, int l_ny, int nx, double *cudaUa, double *cudaUb, double *cudaUc)
{
    int cudaI = blockDim.y * blockIdx.y + threadIdx.y;
    int cudaJ = blockDim.x * blockIdx.x + threadIdx.x;
    int cudaNx = blockDim.x * gridDim.x;
    int matrixIdx = cudaNx * cudaI + cudaJ;
    // int matrixIdx = (blockDim.x * gridDim.x) * (blockDim.y * blockIdx.y + threadIdx.y) +
    //                 blockDim.x * blockIdx.x + threadIdx.x;

    if (matrixIdx < l_ny * nx)
    {
        cudaUa[matrixIdx] = 0;
        cudaUb[matrixIdx] = 0;
        cudaUc[matrixIdx] = 0;
        // cudaUa[matrixIdx] = rank;
        // cudaUb[matrixIdx] = rank;
        // cudaUc[matrixIdx] = rank;
    }
}

__global__ void cudaPrinter(int rank, int l_ny, int nx, double *cudaUa, double *cudaUb, double *cudaUc)
{
    printf("RANK: %d\n\n", rank);
    for (int i = 0; i < l_ny; ++i)
    {
        printf("rank: %i, line: %i, length: %i\n", rank, i, l_ny);
        for (int j = 0; j < nx; ++j)
        {
            printf("%.1f ", cudaUa[i * nx + j]);
        }
        printf("\n");
    }
}

void init_scenario_gpu(int rank, int l_ny, double **cudaUa, double **cudaUb, double **cudaUc)
{
    dim3 dimBlock(32, 32);
    dim3 dimGrid0(6, 6);
    dim3 dimGridX(6, 2);

    size_t allocSize = nx * l_ny * sizeof(double *);

    cudaError_t cudaStatus;
    cudaStatus = cudaMalloc((void **)cudaUa, allocSize);
    if (cudaStatus != cudaSuccess)
    {
        EXIT_ERROR(cudaGetErrorString(cudaStatus));
    }

    cudaStatus = cudaMalloc((void **)cudaUb, allocSize);
    if (cudaStatus != cudaSuccess)
    {
        EXIT_ERROR(cudaGetErrorString(cudaStatus));
    }

    cudaStatus = cudaMalloc((void **)cudaUc, allocSize);
    if (cudaStatus != cudaSuccess)
    {
        EXIT_ERROR(cudaGetErrorString(cudaStatus));
    }

    if (rank == 0)
    {
        cudaInitVectors<<<dimGrid0, dimBlock>>>(rank, l_ny, nx, *cudaUa, *cudaUb, *cudaUc);
    }
    else
    {
        cudaInitVectors<<<dimGridX, dimBlock>>>(rank, l_ny, nx, *cudaUa, *cudaUb, *cudaUc);
    }
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess)
    {
        EXIT_ERROR(cudaGetErrorString(cudaStatus));
    }

    // if (rank == 5)
    // {
    //     cudaPrinter<<<1, 1>>>(rank, l_ny, nx, *cudaUa, *cudaUb, *cudaUc);
    //     cudaStatus = cudaDeviceSynchronize();
    //     if (cudaStatus != cudaSuccess)
    //     {
    //         EXIT_ERROR(cudaGetErrorString(cudaStatus));
    //     }
    // }
}

void init_scenario(int l_ny)
{
    int i;

    ua = (double **)malloc(l_ny * sizeof(double *));
    ub = (double **)malloc(l_ny * sizeof(double *));
    uc = (double **)malloc(l_ny * sizeof(double *));
    for (i = 0; i < l_ny; i++)
    {
        ua[i] = (double *)calloc(nx, sizeof(double));
        ub[i] = (double *)calloc(nx, sizeof(double));
        uc[i] = (double *)calloc(nx, sizeof(double));
    }
}

void recalculate_positions(int rank, int numtask)
{
    if (rank == 0)
        return;

    int i, j;

    scenario[scn_index].source.x -= rank * ((ny / numtask) - 1);

    for (i = 0; i < scenario[scn_index].nr_struct; i++)
        for (j = 0; j < 4; j++)
            scenario[scn_index].structure[i].c_points[j][0] -= rank * ((ny / numtask) - 1);
}
