/*
 * Student:	Trascau Mihai
 * Grupa:	344C4
 * 
 * Lucrare:	Ecuatia undelor pentru acustica 2D
 * Fisier:	acoustics.c
 * Descriere:	Fisier sursa in care este implementat tot programul folosind functiile descrise in celelalte fisiere sursa.
 */

#include "acoustics.h"

scenario_t scenario[MAX_SCENARIOS];
int num_scenarios;
int scn_index;
int nx;
int ny;
int local_ny;
double H;
double MAX_TIME;
double TIME_STEP;
int SAVE_TIME;
double **ua;
double **ub;
double **uc;
double **xchg;

// double *cudaUa, *cudaUb, *cudaUc, *cudaXchg;

__global__ void cudaPrinterX(int rank, int l_ny, int nx, double *cudaUa, double *cudaUb, double *cudaUc)
{
    printf("RANK: %d\n\n", rank);
    for (int i = 0; i < 2; ++i)
    {
        printf("rank: %i, line: %i, length: %i\n", rank, i, l_ny);
        for (int j = 0; j < nx; ++j)
        {
            printf("%.1f ", cudaUc[i * nx + j]);
        }
        printf("\n");
    }

    for (int i = l_ny; i < l_ny + 2; ++i)
    {
        printf("rank: %i, line: %i, length: %i\n", rank, i, l_ny);
        for (int j = 0; j < nx; ++j)
        {
            printf("%.1f ", cudaUc[i * nx + j]);
        }
        printf("\n");
    }
}

__global__ void cudaPrinterMatrix(int rank, int ny, int nx, double *cudaUc)
{
    printf("RANK: %d\n\n", rank);
    for (int i = 0; i < ny; ++i)
    {
        printf("rank: %i, line: %i, length: %i\n", rank, i, ny);
        for (int j = 0; j < nx; ++j)
        {
            printf("%.1f ", cudaUc[i * nx + j]);
        }
        printf("\n");
    }
}

int main(int argc, char *argv[])
{
    double *cudaUa, *cudaUb, *cudaUc;
    double *cudaMpiBuf = NULL;
    double *cudaXchg;
    int i, j;
    int err;
    int numtask, rank;
    time_t start_time;

    double **exportMatrix;

    MPI_Status status;
    cudaError_t cudaStatus;

    MPI_Datatype MPI_STRUCTURE, str_types[1];
    MPI_Datatype MPI_SOURCE, src_types[1];
    MPI_Datatype MPI_SCENARIO, scn_types[4];
    MPI_Aint str_offsets[1], src_offsets[1], scn_offsets[4], extent;
    int str_blockcounts[1], src_blockcounts[2], scn_blockcounts[4];

    err = MPI_Init(&argc, &argv);
    if (err != MPI_SUCCESS)
        EXIT_ERROR("[EROARE] Initializare MPI esuata\n\n");

    MPI_Comm_size(MPI_COMM_WORLD, &numtask);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    str_offsets[0] = 0;
    str_types[0] = MPI_INT;
    str_blockcounts[0] = 8;

    MPI_Type_struct(1, str_blockcounts, str_offsets, str_types, &MPI_STRUCTURE);
    MPI_Type_commit(&MPI_STRUCTURE);

    src_offsets[0] = 0;
    src_types[0] = MPI_INT;
    src_blockcounts[0] = 3;

    MPI_Type_struct(1, src_blockcounts, src_offsets, src_types, &MPI_SOURCE);
    MPI_Type_commit(&MPI_SOURCE);

    scn_offsets[0] = 0;
    scn_types[0] = MPI_INT;
    scn_blockcounts[0] = 5;

    MPI_Type_extent(MPI_INT, &extent);
    scn_offsets[1] = scn_offsets[0] + scn_blockcounts[0] * extent;
    scn_types[1] = MPI_DOUBLE;
    scn_blockcounts[1] = 4;

    MPI_Type_extent(MPI_DOUBLE, &extent);
    scn_offsets[2] = scn_offsets[1] + scn_blockcounts[1] * extent;
    scn_types[2] = MPI_SOURCE;
    scn_blockcounts[2] = 1;

    MPI_Type_extent(MPI_SOURCE, &extent);
    scn_offsets[3] = scn_offsets[2] + scn_blockcounts[2] * extent;
    scn_types[3] = MPI_STRUCTURE;
    scn_blockcounts[3] = MAX_STRUCTURES;

    MPI_Type_struct(4, scn_blockcounts, scn_offsets, scn_types, &MPI_SCENARIO);
    MPI_Type_commit(&MPI_SCENARIO);

    if (rank == 0)
    {
        if (argc != 2)
            EXIT_ERROR("[EROARE] Numar incorect de argumente. Folosire: ./acoustics <input_file_name>\n\n");

        if (import_data(argv[1]))
            EXIT_ERROR("[EROARE] Datele de intrare nu au putut fi incarcate\n\n");

        for (i = 1; i < numtask; i++)
            MPI_Send(&num_scenarios, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
    }
    else
        MPI_Recv(&num_scenarios, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);

    if (rank == 0)
        for (i = 1; i < numtask; i++)
            MPI_Send(&scenario, num_scenarios, MPI_SCENARIO, i, 1, MPI_COMM_WORLD);
    else
        MPI_Recv(&scenario, num_scenarios, MPI_SCENARIO, 0, 1, MPI_COMM_WORLD, &status);

    //print_import_data(rank);

    int scn_index = 0;
    while (scn_index < num_scenarios)
    {
        int step = 0;
        int source_active = 1;

        start_time = time(NULL);
        load_scenario();

        if (rank == 0)
        {
            local_ny = ny / numtask;
            init_scenario(ny);
            init_scenario_gpu(rank, ny, &cudaUa, &cudaUb, &cudaUc);

            exportMatrix = (double **)malloc(ny * sizeof(double *));
            for (int ii = 0; ii < ny; ii++)
            {
                exportMatrix[ii] = (double *)calloc(nx, sizeof(double));
            }
        }
        else if (rank != numtask - 1)
        {
            local_ny = ny / numtask;
            init_scenario(local_ny + 2);
            init_scenario_gpu(rank, local_ny + 2, &cudaUa, &cudaUb, &cudaUc);
        }
        else
        {
            local_ny = ny - (numtask - 1) * (ny / numtask);
            init_scenario(local_ny + 1);
            init_scenario_gpu(rank, local_ny + 1, &cudaUa, &cudaUb, &cudaUc);
        }

        recalculate_positions(rank, numtask);

        int radius = scenario[scn_index].source.radius;
        cudaMpiBuf = (double *)calloc(nx, sizeof(double));
        if (cudaMpiBuf == NULL)
        {
            EXIT_ERROR("Error calloc");
        }

        int source_active_gpu = source_active;
        while (step < (int)(MAX_TIME / TIME_STEP))
        {

            pulse_source(rank, numtask, radius, step, scenario[scn_index].amp, &source_active);
            pulse_source_gpu(rank, step, cudaUa, cudaUb, cudaUc, &source_active_gpu, local_ny, numtask);
            cudaDeviceSynchronize();

            // if (source_active != source_active_gpu)
            // {
            //     printf("There is a problem\n");
            // }

            // m_compute_acoustics(rank, numtask, source_active, radius);
            // compute_acoustics_gpu(rank, numtask, cudaUa, cudaUb, cudaUc, source_active_gpu, local_ny);
            // cudaDeviceSynchronize();

            // MPI_Barrier(MPI_COMM_WORLD);

            // if (rank == 0)
            // {
            //     MPI_Recv(cudaMpiBuf, nx, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &status);
            //     cudaStatus = cudaMemcpy(cudaUc + (local_ny * nx), (const double *)cudaMpiBuf, nx * sizeof(double), cudaMemcpyHostToDevice);
            //     cudaStatus = cudaMemcpy(cudaMpiBuf, (const double *)cudaUc + ((local_ny - 1) * nx), nx * sizeof(double), cudaMemcpyDeviceToHost);
            //     MPI_Send(cudaMpiBuf, nx, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);
            //     ////
            //     MPI_Recv(uc[local_ny], nx, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &status);
            //     MPI_Send(uc[local_ny - 1], nx, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);
            // }
            // else if (rank % 2 == 1 && rank != numtask - 1)
            // {
            //     cudaStatus = cudaMemcpy(cudaMpiBuf, cudaUc + nx, nx * sizeof(double), cudaMemcpyDeviceToHost);
            //     MPI_Send(cudaMpiBuf, nx, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD);
            //     MPI_Recv(cudaMpiBuf, nx, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &status);
            //     cudaStatus = cudaMemcpy(cudaUc, (const double *)cudaMpiBuf, nx * sizeof(double), cudaMemcpyHostToDevice);

            //     cudaStatus = cudaMemcpy(cudaMpiBuf, (const double *)cudaUc + (local_ny * nx), nx * sizeof(double), cudaMemcpyDeviceToHost);
            //     MPI_Send(cudaMpiBuf, nx, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
            //     MPI_Recv(cudaMpiBuf, nx, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &status);
            //     cudaStatus = cudaMemcpy(cudaUc + ((local_ny + 1) * nx), (const double *)cudaMpiBuf, nx * sizeof(double), cudaMemcpyHostToDevice);
            //     ////
            //     MPI_Send(uc[1], nx, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD);
            //     MPI_Recv(uc[0], nx, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &status);
            //     MPI_Send(uc[local_ny], nx, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
            //     MPI_Recv(uc[local_ny + 1], nx, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &status);
            // }
            // else if (rank % 2 == 0 && rank != numtask - 1)
            // {
            //     MPI_Recv(cudaMpiBuf, nx, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &status);
            //     cudaStatus = cudaMemcpy(cudaUc + ((local_ny + 1) * nx), (const double *)cudaMpiBuf, nx * sizeof(double), cudaMemcpyHostToDevice);
            //     cudaStatus = cudaMemcpy(cudaMpiBuf, (const double *)cudaUc + (local_ny * nx), nx * sizeof(double), cudaMemcpyDeviceToHost);
            //     MPI_Send(cudaMpiBuf, nx, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);

            //     MPI_Recv(cudaMpiBuf, nx, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &status);
            //     cudaStatus = cudaMemcpy(cudaUc, (const double *)cudaMpiBuf, nx * sizeof(double), cudaMemcpyHostToDevice);
            //     cudaStatus = cudaMemcpy(cudaMpiBuf, cudaUc + nx, nx * sizeof(double), cudaMemcpyDeviceToHost);
            //     MPI_Send(cudaMpiBuf, nx, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD);
            //     ////
            //     MPI_Recv(uc[local_ny + 1], nx, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &status);
            //     MPI_Send(uc[local_ny], nx, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
            //     MPI_Recv(uc[0], nx, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &status);
            //     MPI_Send(uc[1], nx, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD);
            // }
            // else if (rank == numtask - 1)
            // {
            //     cudaStatus = cudaMemcpy(cudaMpiBuf, (const double *)cudaUc + nx, nx * sizeof(double), cudaMemcpyDeviceToHost);
            //     MPI_Send(cudaMpiBuf, nx, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD);
            //     MPI_Recv(cudaMpiBuf, nx, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &status);
            //     cudaStatus = cudaMemcpy(cudaUc, (const double *)cudaMpiBuf, nx * sizeof(double), cudaMemcpyHostToDevice);
            //     ////
            //     MPI_Send(uc[1], nx, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD);
            //     MPI_Recv(uc[0], nx, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &status);
            // }

            if (step % SAVE_TIME == 1)
            {
                if (rank == 0)
                {
                    for (j = 0; j < ny; j++)
                    {
                        cudaStatus = cudaMemcpy(exportMatrix[j], (const double *)cudaUc + (j * nx), nx * sizeof(double), cudaMemcpyDeviceToHost);
                    }

                    for (i = 1; i < numtask; i++)
                    {
                        if (i != numtask - 1)
                        {
                            for (j = 0; j < ny / numtask; j++)
                            {
                                MPI_Recv(exportMatrix[i * (ny / numtask) + j], nx, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &status);
                                //
                                MPI_Recv(uc[i * (ny / numtask) + j], nx, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &status);
                            }
                        }
                        if (i == numtask - 1)
                        {
                            for (j = (numtask - 1) * (ny / numtask); j < ny; j++)
                            {
                                MPI_Recv(exportMatrix[j], nx, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &status);
                                //
                                MPI_Recv(uc[j], nx, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &status);
                            }
                        }
                    }
                    export_to_vtk(step, exportMatrix);
                    // export_to_vtk(step, uc);
                }
                else
                {
                    for (i = 1; i < local_ny + 1; i++)
                    {
                        cudaStatus = cudaMemcpy(cudaMpiBuf, (const double *)cudaUc + (i * nx), nx * sizeof(double), cudaMemcpyDeviceToHost);
                        MPI_Send(cudaMpiBuf, nx, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
                        //
                        MPI_Send(uc[i], nx, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
                    }
                }
            }

            cudaXchg = cudaUa;
            cudaUa = cudaUb;
            cudaUb = cudaUc;
            cudaUc = cudaXchg;

            xchg = ua;
            ua = ub;
            ub = uc;
            uc = xchg;

            MPI_Barrier(MPI_COMM_WORLD);

            step++;
        }

        unload_scenario();
        unload_scenario_gpu(cudaUa, cudaUb, cudaUc);

        if (rank == 0)
        {
            time_t stop = time(NULL);
            double compute_time = difftime(stop, start_time);
            export_to_gnuplot(scn_index, compute_time);
        }

        scn_index++;
    }

    MPI_Type_free(&MPI_SCENARIO);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}
