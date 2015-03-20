//
//  Main.c
//  Multi-Core Assignment 3
//
//  Created by Myles Ringle on 18/12/2014.
//  Copyright (c) 2014 Myles Ringle. All rights reserved.
//

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include "Network.h"

int main(void) {
    
    // Initialising MPI, getting comm size and each process rank
    int comm_sz;
    int rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // Other initialization
    int N, L, K, numAttackers, M;
    clock_t begin, end;
    double time_taken;
    
    // Seed rand() function to make numbers random.
    // Rank must be included since processes run in parallel they have the same time.
    srand (time(NULL)+rank);
    
    if (rank == 0) {
        // Taking user input
        printf("Please specify the number of inputs: ");
        fflush(stdout);
        scanf("%d", &N);

        printf("Please specify the numbers of perceptrons per Network: ");
        fflush(stdout);
        scanf("%d", &K);

        printf("Please specify the boundaries of the weights: ");
        fflush(stdout);
        scanf("%d", &L);

        printf("Please specify the number of attackers: ");
        fflush(stdout);
        scanf("%d", &numAttackers);


        // Fix the number of inputs (increase until divisible by no. of perceptrons)
        while (N % K != 0) {
            N = N + 1;
        }
        printf("\nCorrected number of inputs: %d\n", N);

        // Calculate num of inputs per perceptron
        M = N/K;
    }
    
    // Broadcasting the user inputs to all other processes
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&K, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&L, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numAttackers, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    //------------------------------------- SETUP --------------------------------------
    
    int **weightsA;
    int **weightsB;
    
    if (rank == 0) {
        printf("\n------------------------------------------------------\n");
        printf("\nINITIAL NETWORK STATE\n\n");

        // Generate weights for A and B
        weightsA = generateWeights(N, K, L, M);
        weightsB = generateWeights(N, K, L, M);
    }
    
    /*
     * If numAttackers is not divisible by the number of processes, then some
     * processes will need to do more work than other.
     * Taking this into account, the following calculates the local size each
     * processes can handle.
     */
    int localSize = 0;
    int remainder = numAttackers % comm_sz;
    if (remainder == 0) {
        localSize = numAttackers/comm_sz;
    }
    else {
        localSize = numAttackers/comm_sz;
        if (rank < remainder) {
            localSize++;
        }
    }
    
    // ATTACKER: Generate weights for attackers
    int ***attackerWeights = malloc(localSize*sizeof(int**));
    for (int i = 0; i < localSize; i++) {
         attackerWeights[i] = generateWeights(N, K, L, M);
    }

    //----------------------------------- END SETUP ------------------------------------
    
    
    
    //-------------------------------- INITIALIZATION ----------------------------------
    
    // Synchronisation steps, total step and input array
    int steps = 0;
    int S = 0;
    int* S_attacker = calloc(localSize, sizeof(int));
    int* inputsArray = malloc(N*sizeof(int));
    
    // Single output figure for each network
    int netAoutput = 0;
    int netBoutput = 0;

    // Arrays with individual perceptron outputs for A and B
    int* pOutputsA = malloc(K*sizeof(int));
    int* pOutputsB = malloc(K*sizeof(int));
    
    // ATTACKER: Output struct for attackers
    Outputs outputsAttackers[localSize];
    
    // ATTACKER: Arrays with individual perceptron outputs for attackers
    int **pOutputsAttackers = malloc(localSize*sizeof(int*));
    for (int i = 0; i < localSize; i++) {
        pOutputsAttackers[i] = malloc(K*sizeof(int));
    }
    

    if (rank == 0) {
        // Store the starting clock time
        begin = clock();
    }
    //------------------------------ END INITIALIZATION --------------------------------

    
    
    //--------------------------------- CALCULATION ------------------------------------
    
    // Loop until either network has been synchronised for 50 consecutive steps or
    // a total of 1,000,000 steps have been carried out.
    while ((S < 50) && (steps < 1000000)) {
        
        if (rank == 0) {
            // Generate inputs
            inputsArray = generateInputs(N, inputsArray);

            // Overall outputs of each network
            Outputs outputsA = generateOutputs(K, M, pOutputsA, weightsA, inputsArray);
            netAoutput = outputsA.netOutput;
            pOutputsA = outputsA.pOutputs;

            Outputs outputsB = generateOutputs(K, M, pOutputsB, weightsB, inputsArray);
            netBoutput = outputsB.netOutput;
            pOutputsB = outputsB.pOutputs;
        }
        // Broadcasting the inputsArray to all other processes
        MPI_Bcast(&(inputsArray[0]), N, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&netAoutput, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&netBoutput, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        // ATTACKER: Output for each attacker network
        for (int i = 0; i < localSize; i++) {
            outputsAttackers[i] = generateOutputs(K, M, pOutputsAttackers[i], attackerWeights[i], inputsArray);
        }
        
        if (rank == 0) {
            // Printing the initial perceptron weights and network outputs
            if (steps == 0) {
                printOutputs(K, M, weightsA, weightsB, netAoutput, netBoutput);
            }
        }
                
        if (rank == 0) {
            // Check network outputs and apply learning rule if neccessary
            if (netAoutput == netBoutput) {

                // Change A & B weights
                weightsA = changeWeights(K, L, M, netAoutput, pOutputsA, weightsA, inputsArray);
                weightsB = changeWeights(K, L, M, netBoutput, pOutputsB, weightsB, inputsArray);

                S++;
            }
            else {
                S = 0;
            }
        }
                
        // ATTACKER: Check outputs for each attacker and apply learning rule if equal
        for (int i = 0; i < localSize; i++) {
            
            if (outputsAttackers[i].netOutput == netAoutput) {
                
                if (netAoutput == netBoutput) {
                    attackerWeights[i] = changeWeights(K, L, M, outputsAttackers[i].netOutput, outputsAttackers[i].pOutputs, attackerWeights[i], inputsArray);
                }
                
                S_attacker[i] = S_attacker[i] + 1;
            }
            else {
                if (netAoutput == netBoutput) {
                
                    int w = 0;
                    int wMin = 0;
                    int chosenPerceptron = 0;
                    
                    for (int j = 0; j < K; j++) {
                        
                        for (int k = 0; k < M; k++) {
                            
                            // Getting absolute minimum perceptron weight * input
                            w = abs(attackerWeights[i][j][k] * inputsArray[(j*M)+k]);
                            
                            if ((w < wMin) || ((j == 0) && (k == 0))) {
                                wMin = w;
                                chosenPerceptron = j;
                            }
                            
                        }
                    }
                    
                    // Flipping the sign and assuming output of A
                    outputsAttackers[i].pOutputs[chosenPerceptron] = outputsAttackers[i].pOutputs[chosenPerceptron] * -1;
                    outputsAttackers[i].netOutput = netAoutput;
                    
                    // Applying learning rule
                    attackerWeights[i] = changeWeights(K, L, M, outputsAttackers[i].netOutput, outputsAttackers[i].pOutputs, attackerWeights[i], inputsArray);
                }
                
                S_attacker[i] = 0;
            }
        }
        
        steps++;
        
        // Broadcast 'S' so other processes know when to break out of loop
        MPI_Bcast(&S, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    //------------------------------- END CALCULATION ----------------------------------

    
    
    //----------------------------------- OUTPUTS --------------------------------------
    
    if (rank == 0) {
        printf("\n------------------------------------------------------\n");
        printf("\nFINAL NETWORK STATE\n\n");

        // Printing the final network state and outputs
        printOutputs(K, M, weightsA, weightsB, netAoutput, netBoutput);

        printf("\n------------------------------------------------------");

        if (S == 50) {
            printf("\n*** NETWORKS SYNCHRONISED ***\n");
            printf("Synchronised steps: %d\n", S);
            printf("Total steps: %d", steps);

            // Take ending clock time and print time taken
            end = clock();
            time_taken = ((double)(end - begin) / CLOCKS_PER_SEC)*1000;
            printf("\nTime Taken: %.3f ms\n", time_taken);
        }
        else {
            printf("\n*** NETWORKS FAILED TO SYNCHRONISE ***\n");
            printf("Total steps: %d\n", steps);
        }

        printf("------------------------------------------------------\n");
        
        printf("ATTACKERS THAT SYNCHRONISED WITH A:\n");
    }
        
    // Barrier to prevent other processes going ahead and printing erratically
    MPI_Barrier(MPI_COMM_WORLD);
    
    /*
     * Since each process will start iterating at i = 0, the following calculates an offset
     * from 0 for each process, otherwise they would all print attacker numbers starting from 0.
     * This also takes into account again that numAttackers may not be divisible by the
     * number of processes.
     */
    int offset = 0;
    if (rank < remainder) {
        offset = rank * localSize;
    }
    else {
        offset = (rank * localSize) + remainder;
    }
    
    int localSynced = 0;
    int attackSynced = 0;
    
    // ATTACKER: Print out any attackers that synchronised
    for (int i = 0; i < localSize; i++) {
        
        if (S_attacker[i] >= 50) {
            printf("- Attacker %d (Synchronised for '%d' steps)\n", i+1+offset, S_attacker[i]);
            localSynced++;
        }
    }
        
    // Reducing all the local temp values of localSynced into 1 overall value for attackSynced.
    // If this value is > 1, we know that atleast one attacker synchronised.
    MPI_Reduce(&localSynced, &attackSynced, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    // Barrier to prevent other processes going ahead and printing erratically
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (rank == 0) {
        // If no attackers synchronised then print 'None'.
        if (attackSynced == 0) {
            printf("None\n");
        }
    }
    
    //--------------------------------- END OUTPUTS ------------------------------------
    
    
    // Releasing memory
    if (rank == 0) {
        free(*weightsA);
        free(weightsA);
        free(*weightsB);
        free(weightsB);
    }
    
    free(inputsArray);
    free(pOutputsA);
    free(pOutputsB);
    free(**attackerWeights);
    free(*attackerWeights);
    free(attackerWeights);
    free(S_attacker);
    
    MPI_Finalize();
    return 0;
}