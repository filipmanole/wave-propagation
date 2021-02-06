/*
 * Student:	Trascau Mihai, Dorin Ionita
 * Grupa:	344C4
 * 
 * Lucrare:	Ecuatia undelor pentru acustica 2D
 * Fisier:	acoustics_alg.h
 * Descriere:	Fisier sursa care contine implementarile pentru algoritmul utilizat (in cazul nostru MDF pentru ecuatia propagarii undei)
 */

// TODO: document CUDA version 9

#include "acoustics.h"
#include <math.h>

extern "C" double **ua, **ub, **uc, gain;

double *xchg_gpu;

// Variabile globale pentru rularea unui scenariu
extern "C" int nx;
extern "C" int ny;
extern "C" double H;
extern "C" double MAX_TIME;
extern "C" double TIME_STEP;
extern "C" int SAVE_TIME;

extern "C" scenario_t scenario[MAX_SCENARIOS];
extern "C" int num_scenarios;
extern "C" int scn_index;

double P_0;
double tau;
double P_r[1000][1000];
double chi[200][200];
//__device__ double gpu_chi[1000][1000];
// TODO: take notes of this

//__device__ scenario_t scenario_gpu[MAX_SCENARIOS];

/***************************** EXTERNAL REQUIREMENTS ***************************/

// required by the export_to_vtk function
int in_structure(int x, int y)
{
    int i;
    for (i = 0; i < scenario[scn_index].nr_struct; i++)
    {
        if (x > scenario[scn_index].structure[i].c_points[0][0] && x < scenario[scn_index].structure[i].c_points[3][0])
            if (y > scenario[scn_index].structure[i].c_points[0][1] && y < scenario[scn_index].structure[i].c_points[1][1])
                return 1;
    }
    return 0;
}

/******************************** DUMPING MODEL ***********************************/

__global__ void free_path_loss_kernel(double amp, double *gpu_chi)
{
    double G = 10;
    double sigma = 1;
    double L = 0.051; // Initial 1
    double f = 3 * 1000000;
    double lightspeed = 3 * 1000000; // missing one zero maybe
    double env_pent = 1.001;
    double lambda = lightspeed / f;
    double R;
    double P_r;
    //double gpu_chi[1000][1000];
    double tau = 1;
    double P_0 = tau * (amp * amp);

    R = sqrt((double)(blockIdx.x * blockIdx.x + blockIdx.y * blockIdx.y));
    P_r = P_0 * G * G * sigma * lambda * lambda;
    P_r /= (pow(4 * M_PI, 3) * pow(R, env_pent) * L);
    gpu_chi[blockIdx.x * blockDim.x + blockIdx.y] = sqrt(P_r / P_0);
    gpu_chi[blockIdx.x * blockDim.x + blockIdx.y] /= 10;
    gpu_chi[blockIdx.x * blockDim.x + blockIdx.y] += 0.9;
}

/************************ SOURCE VARIATION AND REFLEXION MODEL ********************/

__device__ int is_source_gpu(int i, int j,
                             int radius, int source_active,
                             int src_x, int src_y)
{
    if (!source_active)
        return 0;
    if (sqrt(pow((float)(src_x - i), 2) + pow((float)(src_y - j), 2)) <= radius)
        return 1;
    return 0;
}

__global__ void set_all_zero_kernel(double *ua_gpu, double *ub_gpu, double *uc_gpu)
{
    ua_gpu[blockIdx.x * blockDim.x + blockIdx.y] = 0;
    ub_gpu[blockIdx.x * blockDim.x + blockIdx.y] = 0;
    uc_gpu[blockIdx.x * blockDim.x + blockIdx.y] = 0;
    // TODO: sync CPU after this -> move to utils.cu file
}

__global__ void wireless_src_pulse_kernel(int step, double amp,
                                          double MAX_TIME, double TIME_STEP,
                                          int radius, int source_active, int src_x, int src_y,
                                          double *ua_gpu, double *ub_gpu, double *uc_gpu)
{
    int i, j;
    int i_start, j_start;
    int i_final, j_final;
    int line_length;
    int global_thread_x, global_thread_y;
    int thread_work = 32;

    line_length = gridDim.y * blockDim.y;

    global_thread_x = blockDim.x * blockIdx.x + threadIdx.x;
    global_thread_y = blockDim.y * blockIdx.y + threadIdx.y;

    i_start = global_thread_x * thread_work;
    j_start = global_thread_y * thread_work;
    i_final = global_thread_x * (thread_work + 1);
    j_final = global_thread_y * (thread_work + 1);

    if (step < (int)(MAX_TIME / TIME_STEP) / 2)
    {
        // Pulse source
        for (i = i_start; i < i_final; i++)
        {
            for (j = j_start; j < j_final; j++)
            {
                if (is_source_gpu(i, j, radius, 1, src_x, src_y))
                    uc_gpu[i * line_length + j] = amp * fabs(sin(step * M_PI / 4));
            }
        }
    }
    else if (source_active)
    {
        for (i = i_start; i < i_final; i++)
        {
            for (j = j_start; j < j_final; j++)
            {
                if (is_source_gpu(i, j, radius, source_active, src_x, src_y))
                {
                    ua_gpu[i * line_length + j] = 0;
                    ub_gpu[i * line_length + j] = 0;
                    uc_gpu[i * line_length + j] = 0;
                }
            }
        }
    }
    // All threads should reach this point before setting source_active.
    // Option 1:  need a thread barrier here -> not done, I chose option 2
    // Option 2:  simply write 2 kernels and syncCPU -> done, I chose this option
    //	 	CPU is setting source_active = 0 after this kernel is done executing.
}
__device__ int on_edge_gpu(int i, int j)
{
    int line_length = blockDim.y * gridDim.y;    // ny
    int collumn_length = blockDim.x * gridDim.x; // nx

    if (0 == i && 0 != j && j != collumn_length - 1)
        return N_EDGE;
    if (i == j - 1 && j != 0 && j != collumn_length - 1)
        return S_EDGE;
    if (0 == j && 0 != i && i != line_length - 1)
        return W_EDGE;
    if (j == collumn_length - 1 && i != 0 && i != line_length - 1)
        return E_EDGE;
    return 0;
}

__device__ int on_corner_gpu(int i, int j)
{
    int line_length = blockDim.y * gridDim.y;    // ny
    int collumn_length = blockDim.x * gridDim.x; // nx

    if (0 == i && j == 0)
        return NW_CORNER;
    if (0 == i && j == collumn_length - 1)
        return NE_CORNER;
    if (line_length - 1 == i && j == 0)
        return SW_CORNER;
    if (line_length - 1 == i && j == collumn_length - 1)
        return SE_CORNER;
    return 0;
}

__device__ int on_structure_edge_gpu(int i, int j, scenario_t scenario)
{
    int k;

    // TODO: generalize to any scenario - see the original C code
    for (k = 0; k < scenario.nr_struct; k++)
    {
        if (j > scenario.structure[k].c_points[0][1] &&
            j < scenario.structure[k].c_points[1][1])
            if (i == scenario.structure[i].c_points[0][0])
                return N_EDGE;
        if (i > scenario.structure[k].c_points[1][0] &&
            i < scenario.structure[k].c_points[2][0])
            if (j == scenario.structure[i].c_points[1][1])
                return E_EDGE;
        if (j > scenario.structure[k].c_points[3][1] &&
            j < scenario.structure[k].c_points[2][1])
            if (i == scenario.structure[i].c_points[3][0])
                return S_EDGE;
        if (i > scenario.structure[k].c_points[0][0] &&
            i < scenario.structure[k].c_points[3][0])
            if (j == scenario.structure[i].c_points[0][1])
                return W_EDGE;
    }

    return 0;
}

__device__ int in_structure_gpu(int i, int j, scenario_t scenario)
{

    int k;

    // TODO: generalize to any scenario - see the original C code
    for (k = 0; k < scenario.nr_struct; k++)
    {
        if (i > scenario.structure[k].c_points[0][0] &&
            i < scenario.structure[k].c_points[3][0])
            if (j > scenario.structure[k].c_points[0][1] &&
                j < scenario.structure[k].c_points[1][1])
                return 1;
    }

    return 0;
}

__device__ int on_structure_corner_gpu(int i, int j, scenario_t scenario)
{
    int k;

    // TODO: generalize to any scenario - see the original C code
    for (k = 0; k < scenario.nr_struct; k++)
    {
        if (i == scenario.structure[k].c_points[0][0] &&
            j == scenario.structure[k].c_points[0][1])
            return NW_CORNER;
        if (i == scenario.structure[k].c_points[1][0] &&
            j == scenario.structure[k].c_points[1][1])
            return NE_CORNER;
        if (i == scenario.structure[k].c_points[2][0] &&
            j == scenario.structure[k].c_points[2][1])
            return SE_CORNER;
        if (i == scenario.structure[k].c_points[3][0] &&
            j == scenario.structure[k].c_points[3][1])
            return SW_CORNER;
    }
    return 0;
}

__device__ double compute_edge_node_gpu(int i, int j, int side, double *ub_gpu)
{
    int line_length = blockDim.y * gridDim.y;    // ny
    int collumn_length = blockDim.x * gridDim.x; // nx

    switch (side)
    {
    case N_EDGE:
        return ub_gpu[(i + 1) * line_length + j];
    case E_EDGE:
        return ub_gpu[i * line_length + j - 1];
    case S_EDGE:
        return ub_gpu[(i - 1) * line_length + j];
    case W_EDGE:
        return ub_gpu[i * line_length + j + 1];
    default:
        return 0;
    }
}

__device__ double compute_corner_node_gpu(int i, int j, int corner, double *ub_gpu)
{
    int line_length = blockDim.y * gridDim.y;    // ny
    int collumn_length = blockDim.x * gridDim.x; // nx

    switch (corner)
    {
    case NW_CORNER:
        return (ub_gpu[i * line_length + j + 1] +
                ub_gpu[(i + 1) * line_length + j]) /
               2;
    case NE_CORNER:
        return (ub_gpu[(i + 1) * line_length + j] +
                ub_gpu[i * line_length + j - 1]) /
               2;
    case SE_CORNER:
        return (ub_gpu[i * line_length + j - 1] +
                ub_gpu[(i - 1) * line_length + j]) /
               2;
    case SW_CORNER:
        return (ub_gpu[(i - 1) * line_length + j] +
                ub_gpu[i * line_length + j + 1]) /
               2;
    default:
        return 0;
    }
}

__device__ double compute_structure_corner_node_gpu(int i, int j, int corner, double *ub_gpu)
{
    int line_length = blockDim.y * gridDim.y;    // ny
    int collumn_length = blockDim.x * gridDim.x; // nx

    switch (corner)
    {
    case NW_CORNER:
        return (ub_gpu[i * line_length + j - 1] +
                ub_gpu[(i - 1) * line_length + j]) /
               2;
    case NE_CORNER:
        return (ub_gpu[(i - 1) * line_length + j] +
                ub_gpu[i * line_length + j + 1]) /
               2;
    case SE_CORNER:
        return (ub_gpu[i * line_length + j + 1] +
                ub_gpu[(i + 1) * line_length + j]) /
               2;
    case SW_CORNER:
        return (ub_gpu[(i + 1) * line_length + j] +
                ub_gpu[i * line_length + j - 1]) /
               2;
    default:
        return 0;
    }
}

__device__ double compute_structure_edge_node_gpu(int i, int j, int side, double *ub_gpu)
{
    int line_length = blockDim.y * gridDim.y;    // ny
    int collumn_length = blockDim.x * gridDim.x; // nx

    switch (side)
    {
    case N_EDGE:
        return ub_gpu[(i - 1) * line_length + j];
    case E_EDGE:
        return ub_gpu[i * line_length + j + 1];
    case S_EDGE:
        return ub_gpu[(i + 1) * line_length + j];
    case W_EDGE:
        return ub_gpu[i * line_length + j - 1];
    default:
        return 0;
    }
}

__device__ double compute_node_gpu(int i, int j, double gain,
                                   double *ua_gpu, double *ub_gpu, double *uc_gpu)
{
    int line_length = blockDim.y * gridDim.y;    // ny
    int collumn_length = blockDim.x * gridDim.x; // nx

    return (2 * ub_gpu[i * line_length + j] -
            ua_gpu[i * line_length + j] +
            gain * (ub_gpu[(i + 1) * line_length + j] -
                    4 * ub_gpu[i * line_length + j] +
                    ub_gpu[(i - 1) * line_length + j] +
                    ub_gpu[i * line_length + j + 1] +
                    ub_gpu[i * line_length + j - 1]));
}

__global__ void wireless_propagate_kernel(double gain, int radius, int source_active,
                                          int src_x, int src_y, scenario_t scenario,
                                          double *ua_gpu, double *ub_gpu, double *uc_gpu, double *chi_gpu)
{
    int i, j;
    int i_start, i_final;
    int j_start, j_final;
    int place;
    int line_length;
    int global_thread_x, global_thread_y;
    int thread_work = 32;

    line_length = gridDim.y * blockDim.y;

    global_thread_x = blockDim.x * blockIdx.x + threadIdx.x;
    global_thread_y = blockDim.y * blockIdx.y + threadIdx.y;

    i_start = global_thread_x * thread_work;
    j_start = global_thread_y * thread_work;
    i_final = global_thread_x * (thread_work + 1);
    j_final = global_thread_y * (thread_work + 1);

    for (i = i_start; i < i_final; i++)
    {
        for (j = j_start; j < j_final; j++)
        {
            if (!on_corner_gpu(i, j) &&                                      //ok
                !on_edge_gpu(i, j) &&                                        //ok
                !is_source_gpu(i, j, radius, source_active, src_x, src_y) && //ok
                !on_structure_edge_gpu(i, j, scenario) &&                    //ok
                !on_structure_corner_gpu(i, j, scenario) &&                  //ok
                !in_structure_gpu(i, j, scenario))                           //OK
                uc_gpu[i * line_length + j] = compute_node_gpu(i, j, gain, ua_gpu, ub_gpu, uc_gpu);
            else if (place = on_edge_gpu(i, j))
                uc_gpu[i * line_length + j] = compute_edge_node_gpu(i, j, place, ub_gpu);
            else if (place = on_corner_gpu(i, j))
                uc_gpu[i * line_length + j] = compute_corner_node_gpu(i, j, place, ub_gpu);
            else if (place = on_structure_edge_gpu(i, j, scenario))
                uc_gpu[i * line_length + j] = compute_structure_edge_node_gpu(i, j, place, ub_gpu);
            else if (place = on_structure_corner_gpu(i, j, scenario))
                uc_gpu[i * line_length + j] = compute_structure_corner_node_gpu(i, j, place, ub_gpu);

            ua_gpu[i * line_length + j] = 0;

            uc_gpu[i * line_length + j] *= chi_gpu[i * line_length + j];
        }
    }
}

void s_compute_acoustics()
{
    int i, j;

    double *gpu_chi;
    double *ua_gpu, *ub_gpu, *uc_gpu;
    cudaError_t cuda_status;
    //dim3 dimBlock(nx, ny);
    dim3 dimBlock(16, 16);
    dim3 dimGrid(4, 4);

    cuda_status = cudaMalloc((void **)&gpu_chi, nx * ny * sizeof(double));
    if (cudaSuccess != cuda_status)
    {
        printf("Failed cudaMalloc gpu_chi with message %s\n", cudaGetErrorString(cuda_status));
    }
    cuda_status = cudaMalloc((void **)&ua_gpu, nx * ny * sizeof(double *));
    if (cudaSuccess != cuda_status)
    {
        printf("Failed cudaMalloc ua_gpu with message %s\n", cudaGetErrorString(cuda_status));
    }
    cuda_status = cudaMalloc((void **)&ub_gpu, nx * ny * sizeof(double *));
    if (cudaSuccess != cuda_status)
    {
        printf("Failed cudaMalloc ub_gpu with message %s\n", cudaGetErrorString(cuda_status));
    }
    cuda_status = cudaMalloc((void **)&uc_gpu, nx * ny * sizeof(double *));
    if (cudaSuccess != cuda_status)
    {
        printf("Failed cudaMalloc uc_gpu with message %s\n", cudaGetErrorString(cuda_status));
    }

    set_all_zero_kernel<<<dimGrid, dimBlock>>>(ua_gpu, ub_gpu, uc_gpu);
    cuda_status = cudaPeekAtLastError();
    if (cudaSuccess != cuda_status)
    {
        printf("Failed launching set_all_zero_kernel  %s\n", cudaGetErrorString(cuda_status));
    }
    // I need to sync GPU with CPU here so that u* vectors are safely zeroized
    cuda_status = cudaDeviceSynchronize();
    if (cudaSuccess != cuda_status)
    {
        printf("Failed synchronizing CPU and GPU after zero kerner launch %s\n", cudaGetErrorString(cuda_status));
    }
    free_path_loss_kernel<<<dimGrid, dimBlock>>>(scenario[scn_index].source.p_amp, gpu_chi);
    cuda_status = cudaPeekAtLastError();
    if (cudaSuccess != cuda_status)
    {
        printf("Failed launching free_path_loss kernel %s\n", cudaGetErrorString(cuda_status));
    }
    cuda_status = cudaDeviceSynchronize();
    if (cudaSuccess != cuda_status)
    {
        printf("Failed synchronizing CPU and GPU after free path loss kernel launch %s\n", cudaGetErrorString(cuda_status));
    }

    // TODO: Document about paraview - How to use, what version

    int step = 0;
    int source_active = 1;
    int radius = scenario[scn_index].source.radius;

    while (step < (int)(MAX_TIME / TIME_STEP))
    {
        // Pulse source
        wireless_src_pulse_kernel<<<dimGrid, dimBlock>>>(
            step,
            scenario[scn_index].source.p_amp,
            MAX_TIME,
            TIME_STEP,
            radius,
            source_active,
            scenario[scn_index].source.x /* src_x */,
            scenario[scn_index].source.y /* src_y */,
            ua_gpu,
            ub_gpu,
            uc_gpu);
        cudaDeviceSynchronize();
        if (step >= (int)(MAX_TIME / TIME_STEP) / 2)
            source_active = 0;

        // Propagate wave
        wireless_propagate_kernel<<<dimGrid, dimBlock>>>(
            gain,
            radius,
            source_active,
            scenario[scn_index].source.x,
            scenario[scn_index].source.y,
            scenario[0], //TODO: generalize this
            ua_gpu,
            ub_gpu,
            uc_gpu,
            gpu_chi);
        cudaDeviceSynchronize();

        // !!!!!! TODO 2: save time should be extremely rare here, maybe just once

        xchg_gpu = ua_gpu;
        ua_gpu = ub_gpu;
        ub_gpu = uc_gpu;
        uc_gpu = xchg_gpu;

        step++;
    }

    cudaFree(gpu_chi);
    cudaFree(ua_gpu);
    cudaFree(ub_gpu);
    cudaFree(uc_gpu);
}
