//
//  Network.c
//  Multi-Core Assignment 4
//
//  Created by Myles Ringle on 18/12/2014.
//  Copyright (c) 2014 Myles Ringle. All rights reserved.
//

#include "Network.h"

// Method to generate weights for a network
int **generateWeights(int N, int K, int L, int M) {
    
    int **weights = malloc(K*sizeof(int*));
    
    // Outer loop handles num of perceptrons
    for(int i = 0; i < K; i++){
        
        int *values = malloc(M*sizeof(int));
        
        // Inner loop handles num of inputs per perceptron
        for(int j = 0; j < M; j++){
            
            // Generate random weight for the perceptron i
            values[j] = (rand() % (2 * L + 1) - L);
            
        }
        
        weights[i] = values;
    }
    
    return weights;
}


// Method generates an array of random inputs (either -1 or +1)
int *generateInputs(int N, int *inputsArray) {
        
    // Generating the Inputs
    for (int i = 0; i < N; i++) {
        
        int r = rand() % 2;
        
        if (r == 0) {
            r = -1;
        }
        inputsArray[i] = r;
    }
    
    return inputsArray;
}


// Method generates outputs for each perceptron in the network and then overall network output
struct Outputs generateOutputs(int K, int M, int *pOutputs, int **weights, int *inputsArray) {
    
    int netOutput = 1;
    
    // Outer loop handles num of perceptrons
    for (int i = 0; i < K; i++) {
        
        pOutputs[i] = 0;
        
        // Inner loop handles num of inputs per perceptron
        for (int j = 0; j < M; j++) {
            
            // Getting perceptron output
            pOutputs[i] = pOutputs[i] + (weights[i][j] * inputsArray[(i*M)+j]);
            
        }
        
        // Sign the output of perceptron in A
        if (pOutputs[i] <= 0) {
            pOutputs[i] = -1;
        }
        else {
            pOutputs[i] = 1;
        }
        
        netOutput = netOutput * pOutputs[i];
    }
    
    struct Outputs outputs;
    outputs.netOutput = netOutput;
    outputs.pOutputs = pOutputs;
    
    return outputs;
}


// Method prints out the weights for each perceptron in networks A & B
void printOutputs(int K, int M, int **weightsA, int **weightsB, int netAoutput, int netBoutput) {
    
    // Printing the perceptron weights
    for (int i = 0; i < K; i++) {
        
        // Prints network A weights
        printf("A.%d.weights: ", i);
        printf("[");
        for (int j = 0; j < M; j++) {
            
            printf("%d", weightsA[i][j]);
            
            if (j < M-1) {
                printf(", ");
            }
        }
        printf("]\n");
        
        // Print network B weights
        printf("B.%d.weights: ", i);
        printf("[");
        for (int j = 0; j < M; j++) {
            
            printf("%d", weightsB[i][j]);
            
            if (j < M-1) {
                printf(", ");
            }
        }
        printf("]\n");
    }
    
    // Printng the overall network outputs
    printf("\nOutput A: %d", netAoutput);
    printf("\nOutput B: %d\n", netBoutput);
}


// Method prints out weights for an attacker network
// *** Currently un-used but can be implemented to assist with debugging
void printAttackerOutputs(int attacker, int K, int M, int **weights, int netOutput) {
    
    printf("\n------------------------------------------------------\n");
    
    // Printing the perceptron weights
    for (int i = 0; i < K; i++) {
        
        printf("Attacker %d.%d.weights: ", attacker, i);
        printf("[");
        for (int j = 0; j < M; j++) {
            
            printf("%d", weights[i][j]);
            
            if (j < M-1) {
                printf(", ");
            }
        }
        printf("]\n");
    }
    
    // Printng the overall network output
    printf("\nOutput Attacker %d: %d", attacker, netOutput);
    
}


// Method to change weights in network.
// This happens only where individual perceptron output = overall network output
int **changeWeights(int K, int L, int M, int netOutput, int *pOutputs, int **weights, int *inputsArray) {
    
    for (int i = 0; i < K; i++) {
        
        if (pOutputs[i] == netOutput) {
            
            for (int j = 0; j < M; j++) {
                
                int w = (weights[i][j] - (pOutputs[i] * inputsArray[(i*M)+j]));
                
                if (w < (0-L)) {
                    w = (0-L);
                }
                else if (w > L) {
                    w = L;
                }
                
                weights[i][j] = w;
            }
        }
        
    }
    
    return weights;
}