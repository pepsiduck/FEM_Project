/*
 *  main.c
 *  Library for EPL1110 : Finite Elements for dummies
 *  Elasticite lineaire plane
 *  Calcul des densités de force aux noeuds contraints
 *
 *  Copyright (C) 2024 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */
 
#include "fem.h"
#include <time.h>

int main(int argc, char *argv[])
{  
    //parsing
    if(argc > 3)
    {
        printf("Incorrect numbert of arguments.\n");
        exit(EXIT_FAILURE);
    }

    char fileDataName[MAXNAME] = "";
    char fileProblemName[MAXNAME] = "";
    

    femGeo *geometry = geoGetGeometry();

    if(geoMeshRead(fileDataName) == -1)
    {
        printf("Unable to data problem file.\n");
        exit(EXIT_FAILURE);
    }

    femMeshRenumber(geometry->theElements, FEM_CUTHILL_MCKEE);

    printf("\n\n %d\n",femMeshComputeBand(geometry->theElements));

    geoMeshPrint();
    
    femProblem *problem = femElasticityRead(geometry,fileProblemName);
    if(problem == NULL)
    {
        printf("Unable to open problem file.\n");
        exit(EXIT_FAILURE);
    }

    femElasticityPrint(problem);
    
    clock_t start = clock();
    
    double *soluce = femElasticitySolve(problem);
    if(soluce == NULL)
    {
        printf("Aaah eto, bleh.\n");
        exit(EXIT_FAILURE);
    }

    clock_t finish = clock();
    
    if(femSolutionWrite(geometry->theNodes->nNodes, 2, soluce, "../data/result.txt") == -1)
    {
        printf("Unable to open result file.\n");
        exit(EXIT_FAILURE);
    }

    double time = (double)(finish - start) / CLOCKS_PER_SEC;
    printf("Execution time : %f sec\n",time);
    
    femElasticityFree(problem);
    geoFinalize();
    
    exit(EXIT_SUCCESS);
    return 0;  
}

 
