#include "acoustics.h"

int in_structure(int x, int y)
{
    int i;
    for (i = 0; i < scenario[scn_index].nr_struct; i++)
    {
        if (x >= scenario[scn_index].structure[i].c_points[0][0] && x <= scenario[scn_index].structure[i].c_points[3][0])
            if (y >= scenario[scn_index].structure[i].c_points[0][1] && y <= scenario[scn_index].structure[i].c_points[1][1])
                return 1;
    }
    return 0;
}

__device__ int cudaOnEdge(int rank, int numtask, int x, int y, int local_ny, int nx)
{
    if (x == 0 && y != 0 && y != nx - 1 && rank == 0)
        return N_EDGE;
    if (x == local_ny && y != 0 && y != nx - 1 && rank == numtask - 1)
        return S_EDGE;
    if (y == 0 && x != 0)
        if (rank == numtask - 1 && x == local_ny)
            return 0;
        else
            return W_EDGE;
    if (y == nx - 1 && x != 0)
        if (rank == numtask - 1 && x == local_ny)
            return 0;
        else
            return E_EDGE;
    return 0;
}

__device__ int cudaOnCorner(int rank, int numtask, int x, int y, int local_ny, int nx)
{
    if (x == 0 && y == 0 && rank == 0)
        return NW_CORNER;
    if (x == 0 && y == nx - 1 && rank == 0)
        return NE_CORNER;
    if (x == local_ny && y == 0 && rank == numtask - 1)
        return SW_CORNER;
    if (x == local_ny && y == nx - 1 && rank == numtask - 1)
        return SE_CORNER;
    return 0;
}

__device__ int cudaOnStructureEdge(int x, int y, scenario_t scenario)
{
    int i;
    for (i = 0; i < scenario.nr_struct; i++)
    {
        if (y > scenario.structure[i].c_points[0][1] && y < scenario.structure[i].c_points[1][1])
            if (x == scenario.structure[i].c_points[0][0])
                return N_EDGE;
        if (x > scenario.structure[i].c_points[1][0] && x < scenario.structure[i].c_points[2][0])
            if (y == scenario.structure[i].c_points[1][1])
                return E_EDGE;
        if (y > scenario.structure[i].c_points[0][1] && y < scenario.structure[i].c_points[1][1])
            if (x == scenario.structure[i].c_points[3][0])
                return S_EDGE;
        if (x > scenario.structure[i].c_points[0][0] && x < scenario.structure[i].c_points[3][0])
            if (y == scenario.structure[i].c_points[0][1])
                return W_EDGE;
    }
    return 0;
}

__device__ int cudaOnStructureCorner(int x, int y, scenario_t scenario)
{
    int i;
    for (i = 0; i < scenario.nr_struct; i++)
    {
        if (x == scenario.structure[i].c_points[0][0] && y == scenario.structure[i].c_points[0][1])
            return NW_CORNER;
        if (x == scenario.structure[i].c_points[1][0] && y == scenario.structure[i].c_points[1][1])
            return NE_CORNER;
        if (x == scenario.structure[i].c_points[2][0] && y == scenario.structure[i].c_points[2][1])
            return SE_CORNER;
        if (x == scenario.structure[i].c_points[3][0] && y == scenario.structure[i].c_points[3][1])
            return SW_CORNER;
    }
    return 0;
}

__device__ int cudaInStructure(int x, int y, scenario_t scenario)
{
    int i;
    for (i = 0; i < scenario.nr_struct; i++)
    {
        if (x >= scenario.structure[i].c_points[0][0] && x <= scenario.structure[i].c_points[3][0])
            if (y >= scenario.structure[i].c_points[0][1] && y <= scenario.structure[i].c_points[1][1])
                return 1;
    }
    return 0;
}

__device__ double cudaComputeNode(int x, int y, double TIME_STEP, double H, double *cudaUa, double *cudaUb)
{
    int cudaNx = blockDim.x * gridDim.x;
    return (2 * cudaUb[x * cudaNx + y] - cudaUa[x * cudaNx + y] + pow(TIME_STEP, 2) / pow(H, 2) * (cudaUb[(x + 1) * cudaNx + y] - 4 * cudaUb[x * cudaNx + y] + cudaUb[(x - 1) * cudaNx + y] + cudaUb[x * cudaNx + (y + 1)] + cudaUb[x * cudaNx + (y - 1)]));
}

__device__ double cudaComputeEdgeNode(int i, int j, int side, double *cudaUb)
{
    int cudaNx = blockDim.x * gridDim.x;

    switch (side)
    {
    case N_EDGE:
        return cudaUb[(i + 1) * cudaNx + j];
    case E_EDGE:
        return cudaUb[i * cudaNx + (j - 1)];
    case S_EDGE:
        return cudaUb[(i - 1) * cudaNx + j];
    case W_EDGE:
        return cudaUb[i * cudaNx + (j + 1)];
    default:
        return 0;
    }
}

__device__ double cudaComputeCornerNode(int i, int j, int corner, double *cudaUb)
{
    int cudaNx = blockDim.x * gridDim.x;

    switch (corner)
    {
    case NW_CORNER:
        return (cudaUb[i * cudaNx + (j + 1)] + cudaUb[(i + 1) * cudaNx + j]) / 2;
    case NE_CORNER:
        return (cudaUb[(i + 1) * cudaNx + j] + cudaUb[i * cudaNx + (j - 1)]) / 2;
    case SE_CORNER:
        return (cudaUb[i * cudaNx + (j - 1)] + cudaUb[(i - 1) * cudaNx + j]) / 2;
    case SW_CORNER:
        return (cudaUb[(i - 1) * cudaNx + j] + cudaUb[i * cudaNx + (j + 1)]) / 2;
    default:
        return 0;
    }
}

__device__ double cudaComputeStructureCornerNode(int i, int j, int corner, double *cudaUb)
{
    int cudaNx = blockDim.x * gridDim.x;

    switch (corner)
    {
    case NW_CORNER:
        return (cudaUb[i * cudaNx + (j - 1)] + cudaUb[(i - 1) * cudaNx + j]) / 2;
    case NE_CORNER:
        return (cudaUb[(i - 1) * cudaNx + j] + cudaUb[i * cudaNx + (j + 1)]) / 2;
    case SE_CORNER:
        return (cudaUb[i * cudaNx + (j + 1)] + cudaUb[(i + 1) * cudaNx + j]) / 2;
    case SW_CORNER:
        return (cudaUb[(i + 1) * cudaNx + j] + cudaUb[i * cudaNx + (j - 1)]) / 2;
    default:
        return 0;
    }
}

__device__ double cudaComputeStructureEdgeNode(int i, int j, int side, double *cudaUb)
{
    int cudaNx = blockDim.x * gridDim.x;

    switch (side)
    {
    case N_EDGE:
        return cudaUb[(i - 1) * cudaNx + j];
    case E_EDGE:
        return cudaUb[i * cudaNx + (j + 1)];
    case S_EDGE:
        return cudaUb[(i + 1) * cudaNx + j];
    case W_EDGE:
        return cudaUb[i * cudaNx + (j - 1)];
    default:
        return 0;
    }
}

__device__ int cudaIsSource(int x, int y, int source_active, int radius, int src_x, int src_y)
{
    int length = ceil(sqrt(pow((double)(src_x - x), 2) + pow((double)(src_y - y), 2)));
    if (source_active && length <= radius)
        return 1;

    return 0;
}

__global__ void cudaPulseSouce(
    int rank,
    int step,
    double *cudaUa,
    double *cudaUb,
    double *cudaUc,
    int source_active,
    double amp,
    double MAX_TIME,
    double TIME_STEP,
    int radius,
    int src_x,
    int src_y,
    int local_ny,
    int numtask,
    scenario_t scenario)
{
    int cudaI = blockDim.y * blockIdx.y + threadIdx.y;
    int cudaJ = blockDim.x * blockIdx.x + threadIdx.x;
    int cudaNx = blockDim.x * gridDim.x;

    int matrixIdx = cudaNx * cudaI + cudaJ;

    int start, stop;
    if (rank == 0)
        start = 0;
    else
        start = 1;
    if (rank == numtask - 1 || rank == 0)
        stop = local_ny + 1;
    else
        stop = local_ny + 2;

    if (step < (int)(MAX_TIME / TIME_STEP) / 2)
    {
        if (cudaI >= 1 && cudaI < local_ny + 1)
        {
            if (cudaIsSource(cudaI, cudaJ, 1, radius, src_x, src_y))
            {
                cudaUc[matrixIdx] = amp * fabs(sin(step * M_PI / 4));
            }
        }
    }
    else if (source_active)
    {
        if (cudaI >= start && cudaI < stop)
        {
            if (cudaIsSource(cudaI, cudaJ, source_active, radius, src_x, src_y))
            {
                cudaUa[matrixIdx] = 0;
                cudaUb[matrixIdx] = 0;
                cudaUc[matrixIdx] = 0;
            }
        }
    }
}

void pulse_source_gpu(int rank, int step, double *cudaUa, double *cudaUb, double *cudaUc, int *source_active, int local_ny, int numtask)
{
    dim3 dimBlock(32, 32);
    dim3 dimGrid(6, 2);

    double amp = scenario[scn_index].amp;
    int radius = scenario[scn_index].source.radius;
    int src_x = scenario[scn_index].source.x;
    int src_y = scenario[scn_index].source.y;

    cudaPulseSouce<<<dimGrid, dimBlock>>>(rank, step, cudaUa, cudaUb, cudaUc, *source_active, amp, MAX_TIME, TIME_STEP, radius, src_x, src_y, local_ny, numtask, scenario[scn_index]);
    if (step >= (int)(MAX_TIME / TIME_STEP) / 2)
    {
        *source_active = 0;
    }
}

__global__ void compute(
    int rank,
    int numtask,
    double *cudaUa,
    double *cudaUb,
    double *cudaUc,
    int source_active,
    double TIME_STEP,
    double H,
    int radius,
    int src_x,
    int src_y,
    scenario_t scenario,
    int local_ny)
{
    int cudaI = blockDim.y * blockIdx.y + threadIdx.y;
    int cudaJ = blockDim.x * blockIdx.x + threadIdx.x;
    int cudaNx = blockDim.x * gridDim.x;
    int matrixIdx = cudaNx * cudaI + cudaJ;
    int place, start, stop;

    if (rank == 0)
        start = 0;
    else
        start = 1;
    stop = local_ny + 1;

    if (cudaI >= start && cudaI < stop)
    {
        if (!cudaOnCorner(rank, numtask, cudaI, cudaJ, local_ny, cudaNx) &&
            !cudaOnEdge(rank, numtask, cudaI, cudaJ, local_ny, cudaNx) &&
            !cudaIsSource(cudaI, cudaJ, source_active, radius, src_x, src_y) &&
            !cudaOnStructureEdge(cudaI, cudaJ, scenario) &&
            !cudaOnStructureCorner(cudaI, cudaJ, scenario) &&
            !cudaInStructure(cudaI, cudaJ, scenario))
            cudaUc[matrixIdx] = cudaComputeNode(cudaI, cudaJ, TIME_STEP, H, cudaUa, cudaUb);
        else if ((place = cudaOnEdge(rank, numtask, cudaI, cudaJ, local_ny, cudaNx)))
            cudaUc[matrixIdx] = cudaComputeEdgeNode(cudaI, cudaJ, place, cudaUb);
        else if ((place = cudaOnCorner(rank, numtask, cudaI, cudaJ, local_ny, cudaNx)))
            cudaUc[matrixIdx] = cudaComputeCornerNode(cudaI, cudaJ, place, cudaUb);
        else if ((place = cudaOnStructureEdge(cudaI, cudaJ, scenario)))
            cudaUc[matrixIdx] = cudaComputeStructureEdgeNode(cudaI, cudaJ, place, cudaUb);
        else if ((place = cudaOnStructureCorner(cudaI, cudaJ, scenario)))
            cudaUc[matrixIdx] = cudaComputeStructureCornerNode(cudaI, cudaJ, place, cudaUb);

        cudaUa[matrixIdx] = 0;
    }
}

void compute_acoustics_gpu(int rank, int numtask, double *cudaUa, double *cudaUb, double *cudaUc, int source_active, int local_ny)
{
    dim3 dimBlock(32, 32);
    dim3 dimGrid(6, 2);

    scenario_t scn = scenario[scn_index];
    int radius = scenario[scn_index].source.radius;
    int src_x = scenario[scn_index].source.x;
    int src_y = scenario[scn_index].source.y;

    compute<<<dimGrid, dimBlock>>>(rank,
                                   numtask,
                                   cudaUa,
                                   cudaUb,
                                   cudaUc,
                                   source_active,
                                   TIME_STEP,
                                   H,
                                   radius,
                                   src_x,
                                   src_y,
                                   scn,
                                   local_ny);
}
