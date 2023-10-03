#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    //=======================================================
    printf("id = %d \n",world_rank);
    int data = 0;
    for (int i = 0; 10>i;i++){
        //printf("%d", i);
        if(world_rank == 0){
                data=data+1;
                MPI_Send(&data,1, MPI_INT,1,0,MPI_COMM_WORLD);
                MPI_Recv(&data,1,MPI_INT,1,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                printf("data in node %d is %d\n", world_rank, data);

        }else if(world_rank==1){
                MPI_Recv(&data,1,MPI_INT,0,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                data=data+1;
                printf("data in node %d is %d\n", world_rank, data);
                MPI_Send(&data,1,MPI_INT,0,0,MPI_COMM_WORLD);
        }

    }
    /*
    if (world_rank==0) {
         //--- calculation -> data
         data = -1;
         MPI_Send(&data,1,MPI_INT,1,0,MPI_COMM_WORLD);
    } else if (world_rank==1) {
         MPI_Recv(&data,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
         printf("data in node %d is %d \n",world_rank,data);
    }*/

    //=======================================================

    // Finalize the MPI environment.
    MPI_Finalize();
}