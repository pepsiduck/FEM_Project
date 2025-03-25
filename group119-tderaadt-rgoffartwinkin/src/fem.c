/*
 *  fem.c
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2021 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#include "fem.h"

femGeo theGeometry;

femGeo *geoGetGeometry()                        { return &theGeometry; }

double geoSizeDefault(double x, double y)       { return theGeometry.h; }

double geoGmshSize(int dim, int tag, double x, double y, double z, double lc, void *data)
                                                { return theGeometry.geoSize(x,y);    }
void geoInitialize() 
{
    int ierr;
    theGeometry.geoSize = geoSizeDefault;
    theGeometry.theNodes = NULL;
    theGeometry.theElements = NULL;
    theGeometry.theEdges = NULL;
    theGeometry.nDomains = 0;
    theGeometry.theDomains = NULL;
}

void geoFinalize() 
{
    int ierr;
    
    if (theGeometry.theNodes) {
        free(theGeometry.theNodes->X);
        free(theGeometry.theNodes->Y);
        free(theGeometry.theNodes->number);
        free(theGeometry.theNodes); }
    if (theGeometry.theElements) {
        free(theGeometry.theElements->elem);
        free(theGeometry.theElements); }
    if (theGeometry.theEdges) {
        free(theGeometry.theEdges->elem);
        free(theGeometry.theEdges); }
    for (int i=0; i < theGeometry.nDomains; i++) {
        free(theGeometry.theDomains[i]->elem);
        free(theGeometry.theDomains[i]);  }
    free(theGeometry.theDomains);
}


void geoSetSizeCallback(double (*geoSize)(double x, double y)) 
{
    theGeometry.geoSize = geoSize; 
}

void geoMeshPrint() 
{
   femNodes *theNodes = theGeometry.theNodes;
   if (theNodes != NULL) {
      printf("Number of nodes %d \n", theNodes->nNodes);
      for (int i = 0; i < theNodes->nNodes; i++) {
        printf("%6d : %14.7e %14.7e \n",i,theNodes->X[i],theNodes->Y[i]); }}
   femMesh *theEdges = theGeometry.theEdges;
   if (theEdges != NULL) {
     printf("Number of edges %d \n", theEdges->nElem);
     int *elem = theEdges->elem;
     for (int i = 0; i < theEdges->nElem; i++) {
        printf("%6d : %6d %6d \n",i,elem[2*i],elem[2*i+1]); }}
   femMesh *theElements = theGeometry.theElements;
   if (theElements != NULL) {
     if (theElements->nLocalNode == 3) {
        printf("Number of triangles %d \n", theElements->nElem);
        int *elem = theElements->elem;
        for (int i = 0; i < theElements->nElem; i++) {
            printf("%6d : %6d %6d %6d\n",i,elem[3*i],elem[3*i+1],elem[3*i+2]); }}
     if (theElements->nLocalNode == 4) {
        printf("Number of quads %d \n", theElements->nElem);
        int *elem = theElements->elem;
        for (int i = 0; i < theElements->nElem; i++) {
            printf("%6d : %6d %6d %6d %6d\n",i,elem[4*i],elem[4*i+1],elem[4*i+2],elem[4*i+3]); }}}
   int nDomains = theGeometry.nDomains;
   printf("Number of domains %d\n", nDomains);
   for (int iDomain = 0; iDomain < nDomains; iDomain++) {
      femDomain *theDomain = theGeometry.theDomains[iDomain];
      printf("  Domain : %6d \n", iDomain);
      printf("  Name : %s\n", theDomain->name);
      printf("  Number of elements : %6d\n", theDomain->nElem);
      for (int i=0; i < theDomain->nElem; i++){
 //         if (i != theDomain->nElem  && (i % 10) != 0)  printf(" - ");
          printf("%6d",theDomain->elem[i]);
          if ((i+1) != theDomain->nElem  && (i+1) % 10 == 0) printf("\n"); }
      printf("\n"); }
  
  
}


void geoMeshWrite(const char *filename) 
{
   FILE* file = fopen(filename,"w");
 
   femNodes *theNodes = theGeometry.theNodes;
   fprintf(file, "Number of nodes %d \n", theNodes->nNodes);
   for (int i = 0; i < theNodes->nNodes; i++) {
      fprintf(file,"%6d : %14.7e %14.7e \n",i,theNodes->X[i],theNodes->Y[i]); }
      
   femMesh *theEdges = theGeometry.theEdges;
   fprintf(file,"Number of edges %d \n", theEdges->nElem);
   int *elem = theEdges->elem;
   for (int i = 0; i < theEdges->nElem; i++) {
      fprintf(file,"%6d : %6d %6d \n",i,elem[2*i],elem[2*i+1]); }
      
   femMesh *theElements = theGeometry.theElements;
   if (theElements->nLocalNode == 3) {
      fprintf(file,"Number of triangles %d \n", theElements->nElem);
      elem = theElements->elem;
      for (int i = 0; i < theElements->nElem; i++) {
          fprintf(file,"%6d : %6d %6d %6d\n",i,elem[3*i],elem[3*i+1],elem[3*i+2]); }}
   if (theElements->nLocalNode == 4) {
      fprintf(file,"Number of quads %d \n", theElements->nElem);
      elem = theElements->elem;
      for (int i = 0; i < theElements->nElem; i++) {
          fprintf(file,"%6d : %6d %6d %6d %6d\n",i,elem[4*i],elem[4*i+1],elem[4*i+2],elem[4*i+3]); }}
     
   int nDomains = theGeometry.nDomains;
   fprintf(file,"Number of domains %d\n", nDomains);
   for (int iDomain = 0; iDomain < nDomains; iDomain++) {
      femDomain *theDomain = theGeometry.theDomains[iDomain];
      fprintf(file,"  Domain : %6d \n", iDomain);
      fprintf(file,"  Name : %s\n", theDomain->name);
      fprintf(file,"  Number of elements : %6d\n", theDomain->nElem);
      for (int i=0; i < theDomain->nElem; i++){
          fprintf(file,"%6d",theDomain->elem[i]);
          if ((i+1) != theDomain->nElem  && (i+1) % 10 == 0) fprintf(file,"\n"); }
      fprintf(file,"\n"); }
    
   fclose(file);
}

int geoMeshRead(const char *filename) 
{
   FILE* file = NULL;
   if(strcmp(filename,"") == 0)
     file = fopen("../data/mesh.txt","r");
   else
     file = fopen(filename,"r");
   if(!file)
     return -1;
   
   int trash, *elem;
   
   femNodes *theNodes = malloc(sizeof(femNodes));
   theGeometry.theNodes = theNodes;
   ErrorScan(fscanf(file, "Number of nodes %d \n", &theNodes->nNodes));
   theNodes->X = malloc(sizeof(double)*(theNodes->nNodes));
   theNodes->Y = malloc(sizeof(double)*(theNodes->nNodes));
   theNodes->number = malloc(sizeof(int)*(theNodes->nNodes));
   for (int i = 0; i < theNodes->nNodes; i++) {
       ErrorScan(fscanf(file,"%d : %le %le \n",&theNodes->number[i],&theNodes->X[i],&theNodes->Y[i]));} 

   femMesh *theEdges = malloc(sizeof(femMesh));
   theGeometry.theEdges = theEdges;
   theEdges->nLocalNode = 2;
   theEdges->nodes = theNodes;
   ErrorScan(fscanf(file, "Number of edges %d \n", &theEdges->nElem));
   theEdges->elem = malloc(sizeof(int)*theEdges->nLocalNode*theEdges->nElem);
   for(int i=0; i < theEdges->nElem; ++i) {
        elem = theEdges->elem;
        ErrorScan(fscanf(file, "%6d : %6d %6d \n", &trash,&elem[2*i],&elem[2*i+1])); }
  
   femMesh *theElements = malloc(sizeof(femMesh));
   theGeometry.theElements = theElements;
   theElements->nLocalNode = 0;
   theElements->nodes = theNodes;
   char elementType[MAXNAME];  
   ErrorScan(fscanf(file, "Number of %s %d \n",elementType,&theElements->nElem));  
   if (strncasecmp(elementType,"triangles",MAXNAME) == 0) {
      theElements->nLocalNode = 3;
      theElements->elem = malloc(sizeof(int)*theElements->nLocalNode*theElements->nElem);
      for(int i=0; i < theElements->nElem; ++i) {
          elem = theElements->elem;
          ErrorScan(fscanf(file, "%6d : %6d %6d %6d \n", 
                    &trash,&elem[3*i],&elem[3*i+1],&elem[3*i+2])); }}
   if (strncasecmp(elementType,"quads",MAXNAME) == 0) {
      theElements->nLocalNode = 4;
      theElements->elem = malloc(sizeof(int)*theElements->nLocalNode*theElements->nElem);
      for(int i=0; i < theElements->nElem; ++i) {
          elem = theElements->elem;
          ErrorScan(fscanf(file, "%6d : %6d %6d %6d %6d \n", 
                    &trash,&elem[4*i],&elem[4*i+1],&elem[4*i+2],&elem[4*i+3])); }}
           
   ErrorScan(fscanf(file, "Number of domains %d\n", &theGeometry.nDomains));
   int nDomains = theGeometry.nDomains;
   theGeometry.theDomains = malloc(sizeof(femDomain*)*nDomains);
   for (int iDomain = 0; iDomain < nDomains; iDomain++) {
      femDomain *theDomain = malloc(sizeof(femDomain)); 
      theGeometry.theDomains[iDomain] = theDomain;
      theDomain->mesh = theEdges; 
      ErrorScan(fscanf(file,"  Domain : %6d \n", &trash));
      ErrorScan(fscanf(file,"  Name : %[^\n]s \n", (char*)&theDomain->name));
      ErrorScan(fscanf(file,"  Number of elements : %6d\n", &theDomain->nElem));
      theDomain->elem = malloc(sizeof(int)*2*theDomain->nElem); 
      for (int i=0; i < theDomain->nElem; i++){
          ErrorScan(fscanf(file,"%6d",&theDomain->elem[i]));
          if ((i+1) != theDomain->nElem  && (i+1) % 10 == 0) ErrorScan(fscanf(file,"\n")); }}
    
   fclose(file);
}

void geoSetDomainName(int iDomain, char *name) 
{
    if (iDomain >= theGeometry.nDomains)  Error("Illegal domain number");
    if (geoGetDomain(name) != -1)         Error("Cannot use the same name for two domains");
    sprintf(theGeometry.theDomains[iDomain]->name,"%s",name);
} 

int geoGetDomain(char *name)
{
    int theIndex = -1;
    int nDomains = theGeometry.nDomains;
    for (int iDomain = 0; iDomain < nDomains; iDomain++) {
        femDomain *theDomain = theGeometry.theDomains[iDomain];
        if (strncasecmp(name,theDomain->name,MAXNAME) == 0)
            theIndex = iDomain;  }
    return theIndex;
            
}


int femMinIndex(int *tab, int n)
{
    if(n <= 0)
        return FALSE;
    int minimal = tab[0];
    int index = 0;
    for(int i = 1; i < n; ++i)
    {
        if(tab[i] < minimal)
        {
            index = i;
            minimal = tab[i];
        }
    }
    return index;
}

int femMaxIndex(int *tab, int n)
{
    if(n <= 0)
        return FALSE;
    int maximal = tab[0];
    int index = 0;
    for(int i = 1; i < n; ++i)
    {
        if(tab[i] > maximal)
        {
            index = i;
            maximal = tab[i];
        }
    }
    return index;
}

int abs_value(int arg)
{
    if(arg > 0)
        return arg;
    return -1 * arg;
}

double *global_xy;

int cmp_xy(const void *a, const void *b)
{
    int *A = (int *)a;
    int *B = (int *)b;

    if(global_xy[*A] < global_xy[*B])
        return -1;
    else if(global_xy[*A] > global_xy[*B])
        return 1;
    return 0; 
}

void reverse_tab(int *tab, int n)
{
    int buffer;
    for(int i = 0; i < n / 2; ++i)
    {
        buffer = tab[i];
        tab[i] = tab[n - 1 - i];
        tab[n - 1 - i] = buffer;
    }
}

void femMeshDegrees(femMesh *theMesh, int *degrees, int *redundancy, int n)
{
    for(int d = 0; d < n; ++d)
    {
        degrees[d] = 0;
        for(int d2 = 0; d2 < n; ++d2)
            redundancy[d2] = -1;
        for(int elem = 0; elem < theMesh->nElem; ++elem)
        {
            for(int i = 0; i < theMesh->nLocalNode; ++i)
            {
                if(theMesh->elem[elem * theMesh->nLocalNode + i] == d)
                {
                    if(redundancy[theMesh->elem[elem * theMesh->nLocalNode + ((i + 1) % theMesh->nLocalNode)]] == -1)
                    {
                        degrees[d] += 1;
                        redundancy[theMesh->elem[elem * theMesh->nLocalNode + ((i + 1) % theMesh->nLocalNode)]] = 1;
                    }
                    if(redundancy[theMesh->elem[elem * theMesh->nLocalNode + ((i + theMesh->nLocalNode - 1) % theMesh->nLocalNode)]] == -1)
                    {
                        degrees[d] += 1;
                        redundancy[theMesh->elem[elem * theMesh->nLocalNode + ((i + theMesh->nLocalNode - 1) % theMesh->nLocalNode)]] = 1;
                    }   
                }
            }
        }
    }
}

int *degrees_global;

int cmp_degrees(const void *a, const void *b)
{
    int *A = (int *)a;
    int *B = (int *)b;
    
    if(degrees_global[*A] < degrees_global[*B])
        return -1;
    if(degrees_global[*A] > degrees_global[*B])
        return 1;
    return 0;
}

void femMeshNeighboors(femMesh *theMesh, int **neighboors, int *degrees, int *redundancy, int n)
{
    int counter;
    degrees_global = degrees;
    for(int d = 0; d < n; ++d)
    {
        counter = 0;
        for(int d2 = 0; d2 < n; ++d2)
            redundancy[d2] = -1;
        for(int elem = 0; elem < theMesh->nElem; ++elem)
        {
            for(int i = 0; i < theMesh->nLocalNode; ++i)
            {
                if(theMesh->elem[elem * theMesh->nLocalNode + i] == d)
                {
                    if(redundancy[theMesh->elem[elem * theMesh->nLocalNode + ((i + 1) % theMesh->nLocalNode)]] == -1)
                    {
                        neighboors[d][counter] = theMesh->elem[elem * theMesh->nLocalNode + ((i + 1) % theMesh->nLocalNode)];
                        redundancy[theMesh->elem[elem * theMesh->nLocalNode + ((i + 1) % theMesh->nLocalNode)]] = 1;
                        counter++;
                    }
                    if(redundancy[theMesh->elem[elem * theMesh->nLocalNode + ((i + theMesh->nLocalNode - 1) % theMesh->nLocalNode)]] == -1)
                    {
                        neighboors[d][counter] = theMesh->elem[elem * theMesh->nLocalNode + ((i  + theMesh->nLocalNode - 1) % theMesh->nLocalNode)];
                        redundancy[theMesh->elem[elem * theMesh->nLocalNode + ((i + theMesh->nLocalNode - 1) % theMesh->nLocalNode)]] = 1;
                        counter++;
                    }   
                }
            }
        }
        
        qsort(neighboors[d],degrees[d],sizeof(int),cmp_degrees);

    }
}

void advance_queue(int *queue, int n) //optimized for -1 values
{
    for(int i = 0; i < n - 1; ++i)
    {
        if(queue[i] != -1)
            queue[i] = queue[i + 1];
        else
            break;
    }
    queue[n - 1] = -1;
}

int isInArray(int *tab, int n, int arg) //optimized for -1 values
{
    for(int i = 0; i < n; ++i)
    {
        if(tab[i] == -1)
            return FALSE;
        if(tab[i] == arg)
            return TRUE;
    }
    return FALSE;
}

int femMeshRenumber(femMesh *theMesh, femRenumType renumType)
{
    
    double *array = NULL;
    double buffer;
    int swapped;
    int *tab;
    switch (renumType) {
        case FEM_NO :
            for (int i = 0; i < theMesh->nodes->nNodes; i++) 
                theMesh->nodes->number[i] = i;
            break;
    case FEM_XNUM :
            tab = (int *) malloc(sizeof(int)*theMesh->nodes->nNodes);
            for (int i = 0; i < theMesh->nodes->nNodes; i++) 
                tab[i] = i;

            global_xy = theMesh->nodes->X;

            qsort(tab, theMesh->nodes->nNodes, sizeof(int) , cmp_xy);

            for (int i = 0; i < theMesh->nodes->nNodes; i++) 
                theMesh->nodes->number[tab[i]] = i;
            
            free(tab);
            
            break;
        case FEM_YNUM : 
            tab = (int *) malloc(sizeof(int)*theMesh->nodes->nNodes);
            for (int i = 0; i < theMesh->nodes->nNodes; i++) 
                tab[i] = i;

            global_xy = theMesh->nodes->Y;
            
            qsort(tab, theMesh->nodes->nNodes, sizeof(int) , cmp_xy);

            for (int i = 0; i < theMesh->nodes->nNodes; i++) 
                theMesh->nodes->number[tab[i]] = i;
            
            free(tab);
            
            break;
        case FEM_CUTHILL_MCKEE :
            for (int i = 0; i < theMesh->nodes->nNodes; i++) 
                theMesh->nodes->number[i] = i;

            int *degrees = (int*) malloc(sizeof(int)*theMesh->nodes->nNodes);
            int **neighboors = (int**) malloc(sizeof(int*)*theMesh->nodes->nNodes);
            int *R = (int*) malloc(sizeof(int)*theMesh->nodes->nNodes);
        
            femMeshDegrees(theMesh, degrees, R ,theMesh->nodes->nNodes);

            for (int i = 0; i < theMesh->nodes->nNodes; i++)
                neighboors[i] = (int*) malloc(sizeof(int)*degrees[i]);

            femMeshNeighboors(theMesh, neighboors, degrees, R, theMesh->nodes->nNodes);

            for(int i = 0; i < theMesh->nodes->nNodes; ++i)
                R[i] = -1;

            int *Q = (int*) malloc(sizeof(int)*theMesh->nodes->nNodes);
            for(int i = 0; i < theMesh->nodes->nNodes; ++i)
                Q[i] = -1;

            R[0] = femMinIndex(degrees,theMesh->nodes->nNodes);
            for(int d = 0; d < degrees[R[0]]; ++d)
                Q[d] = neighboors[R[0]][d];

            int counter_queue = degrees[R[0]];  // index du premier élément égal de -1 dans Q. En principe c'est jamais égal à 0
            int counter_R = 1; //index du premier élément égal à -1 dans R
            int buffer;            

            while(counter_R < theMesh->nodes->nNodes)
            {
                //En principe, je peux pas avoir ici Q[0] = -1
                buffer = Q[0]; 
                advance_queue(Q, theMesh->nodes->nNodes);
                counter_queue--; 
                if(isInArray(R, theMesh->nodes->nNodes, buffer) == FALSE)
                {
                    R[counter_R] = buffer;
                    counter_R++;
                    for(int i = 0; i < degrees[buffer]; ++i)
                    {
                        if(isInArray(Q, theMesh->nodes->nNodes, neighboors[buffer][i]) == FALSE && isInArray(R, theMesh->nodes->nNodes, neighboors[buffer][i]) == FALSE)
                        {
                            Q[counter_queue] = neighboors[buffer][i];
                            counter_queue++;
                        }
                    }
                }
            }

            reverse_tab(R, theMesh->nodes->nNodes);
            for (int i = 0; i < theMesh->nodes->nNodes; i++)
                theMesh->nodes->number[R[i]] = i;
                        
            for (int i = 0; i < theMesh->nodes->nNodes; i++)
                free(neighboors[i]);
            free(Q);
            free(neighboors);
            free(degrees);
            free(R);

            break;
        default : 
            Error("Unexpected renumbering option"); 
            return -1;
            break;
    }

    return 0;
}

int femMeshComputeBand(femMesh *theMesh)
{
    int myBand = 0;
    for(int n = 0; n < theMesh->nElem; ++n)
    {
        for(int i = 0; i < theMesh->nLocalNode; ++i)
        {
            for(int j = i + 1; j < theMesh->nLocalNode; ++j)
            {
                if(abs_value(theMesh->nodes->number[theMesh->elem[theMesh->nLocalNode * n + i]] - theMesh->nodes->number[theMesh->elem[theMesh->nLocalNode * n + j]]) > myBand)
                    myBand = abs_value(theMesh->nodes->number[theMesh->elem[theMesh->nLocalNode * n + i]] - theMesh->nodes->number[theMesh->elem[theMesh->nLocalNode * n + j]]);    
            }
        }
    }
    return myBand + 1;
}

femBandSystem *femBandSystemCreate(int size, int band)
{
    femBandSystem *myBandSystem = malloc(sizeof(femBandSystem));
    myBandSystem->B = malloc(sizeof(double)*size*(band+1));
    myBandSystem->A = malloc(sizeof(double*)*size);        
    myBandSystem->size = size;
    myBandSystem->band = band;
    myBandSystem->A[0] = myBandSystem->B + size;
    int i;
    for (i=1 ; i < size ; i++) 
        myBandSystem->A[i] = myBandSystem->A[i-1] + band - 1;
    femBandSystemInit(myBandSystem);
    return(myBandSystem);
}
 
void femBandSystemFree(femBandSystem *myBandSystem)
{
    free(myBandSystem->B);
    free(myBandSystem->A); 
    free(myBandSystem);
}
 
void femBandSystemInit(femBandSystem *myBandSystem)
{
    int i;
    int size = myBandSystem->size;
    int band = myBandSystem->band;
    for (i=0 ; i < size*(band+1) ; i++) 
        myBandSystem->B[i] = 0;        
}
 
void femBandSystemPrint(femBandSystem *myBand)
{
    double  **A, *B;
    int     i, j, band, size;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;

    for (i=0; i < size; i++) {
        for (j=i; j < i+band; j++)
            if (A[i][j] == 0) printf("         ");   
            else              printf(" %+.1e",A[i][j]);
        printf(" :  %+.1e \n",B[i]); }
}
  
void femBandSystemPrintInfos(femBandSystem *myBand)
{
    int size = myBand->size;
    int band = myBand->band;
    printf(" \n");
    printf("    Banded Gaussian elimination \n");
    printf("    Storage informations \n");
    printf("    Matrix size      : %8d\n",size);
    printf("    Matrix band      : %8d\n",band);
    printf("    Bytes required   : %8d\n",(int)sizeof(double)*size*(band+1));     
}


double femBandSystemGet(femBandSystem* myBandSystem, int myRow, int myCol)
{
    double value = 0;
    if (myCol >= myRow && myCol < myRow+myBandSystem->band)  value = myBandSystem->A[myRow][myCol]; 
    return(value);
}

static const double _gaussQuad4Xsi[4]    = {-0.577350269189626,-0.577350269189626, 0.577350269189626, 0.577350269189626};
static const double _gaussQuad4Eta[4]    = { 0.577350269189626,-0.577350269189626,-0.577350269189626, 0.577350269189626};
static const double _gaussQuad4Weight[4] = { 1.000000000000000, 1.000000000000000, 1.000000000000000, 1.000000000000000};
static const double _gaussTri3Xsi[3]     = { 0.166666666666667, 0.666666666666667, 0.166666666666667};
static const double _gaussTri3Eta[3]     = { 0.166666666666667, 0.166666666666667, 0.666666666666667};
static const double _gaussTri3Weight[3]  = { 0.166666666666667, 0.166666666666667, 0.166666666666667};
static const double _gaussEdge2Xsi[2]    = { 0.577350269189626,-0.577350269189626};
static const double _gaussEdge2Weight[2] = { 1.000000000000000, 1.000000000000000};



femIntegration *femIntegrationCreate(int n, femElementType type)
{
    femIntegration *theRule = malloc(sizeof(femIntegration));
    if (type == FEM_QUAD && n == 4) {
        theRule->n      = 4;
        theRule->xsi    = _gaussQuad4Xsi;
        theRule->eta    = _gaussQuad4Eta;
        theRule->weight = _gaussQuad4Weight; }
    else if (type == FEM_TRIANGLE && n == 3) {
        theRule->n      = 3;
        theRule->xsi    = _gaussTri3Xsi;
        theRule->eta    = _gaussTri3Eta;
        theRule->weight = _gaussTri3Weight; }
    else if (type == FEM_EDGE && n == 2) {
        theRule->n      = 2;
        theRule->xsi    = _gaussEdge2Xsi;
        theRule->eta    = NULL;
        theRule->weight = _gaussEdge2Weight; }
    else Error("Cannot create such an integration rule !");
    return theRule; 
}

void femIntegrationFree(femIntegration *theRule)
{
    free(theRule);
}

void _q1c0_x(double *xsi, double *eta) 
{
    xsi[0] =  1.0;  eta[0] =  1.0;
    xsi[1] = -1.0;  eta[1] =  1.0;
    xsi[2] = -1.0;  eta[2] = -1.0;
    xsi[3] =  1.0;  eta[3] = -1.0;
}

void _q1c0_phi(double xsi, double eta, double *phi)
{
    phi[0] = (1.0 + xsi) * (1.0 + eta) / 4.0;  
    phi[1] = (1.0 - xsi) * (1.0 + eta) / 4.0;
    phi[2] = (1.0 - xsi) * (1.0 - eta) / 4.0;
    phi[3] = (1.0 + xsi) * (1.0 - eta) / 4.0;
}

void _q1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    dphidxsi[0] =   (1.0 + eta) / 4.0;  
    dphidxsi[1] = - (1.0 + eta) / 4.0;
    dphidxsi[2] = - (1.0 - eta) / 4.0;
    dphidxsi[3] =   (1.0 - eta) / 4.0;
    dphideta[0] =   (1.0 + xsi) / 4.0;  
    dphideta[1] =   (1.0 - xsi) / 4.0;
    dphideta[2] = - (1.0 - xsi) / 4.0;
    dphideta[3] = - (1.0 + xsi) / 4.0;

}

void _p1c0_x(double *xsi, double *eta) 
{
    xsi[0] =  0.0;  eta[0] =  0.0;
    xsi[1] =  1.0;  eta[1] =  0.0;
    xsi[2] =  0.0;  eta[2] =  1.0;
}

void _p1c0_phi(double xsi, double eta, double *phi)
{
    phi[0] = 1 - xsi - eta;  
    phi[1] = xsi;
    phi[2] = eta;
}

void _p1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    dphidxsi[0] = -1.0;  
    dphidxsi[1] =  1.0;
    dphidxsi[2] =  0.0;
    dphideta[0] = -1.0;  
    dphideta[1] =  0.0;
    dphideta[2] =  1.0;
}

void _e1c0_x(double *xsi) 
{
    xsi[0] = -1.0;  
    xsi[1] =  1.0;  
}

void _e1c0_phi(double xsi,  double *phi)
{
    phi[0] = (1 - xsi) / 2.0;  
    phi[1] = (1 + xsi) / 2.0;
}

void _e1c0_dphidx(double xsi, double *dphidxsi)
{
    dphidxsi[0] = -0.5;  
    dphidxsi[1] =  0.5;
}



femDiscrete *femDiscreteCreate(int n, femElementType type)
{
    femDiscrete *theSpace = malloc(sizeof(femDiscrete));
    theSpace->type = type;  
    theSpace->n = 0;
    theSpace->x = NULL;
    theSpace->phi = NULL;   
    theSpace->dphidx = NULL;    
    theSpace->x2 = NULL;    
    theSpace->phi2 = NULL;
    theSpace->dphi2dx = NULL;
 
    if (type == FEM_QUAD && n == 4) {
        theSpace->n       = 4;
        theSpace->x2      = _q1c0_x;
        theSpace->phi2    = _q1c0_phi;
        theSpace->dphi2dx = _q1c0_dphidx; }
    else if (type == FEM_TRIANGLE && n == 3) {
        theSpace->n       = 3;
        theSpace->x2      = _p1c0_x;
        theSpace->phi2    = _p1c0_phi;
        theSpace->dphi2dx = _p1c0_dphidx; }
    else if (type == FEM_EDGE && n == 2) {
        theSpace->n       = 2;
        theSpace->x       = _e1c0_x;
        theSpace->phi     = _e1c0_phi;
        theSpace->dphidx  = _e1c0_dphidx; }
    else Error("Cannot create such a discrete space !");
    return theSpace; 
}

void femDiscreteFree(femDiscrete *theSpace)
{
    free(theSpace);
}

void femDiscreteXsi2(femDiscrete* mySpace, double *xsi, double *eta)
{
    mySpace->x2(xsi,eta);
}

void femDiscretePhi2(femDiscrete* mySpace, double xsi, double eta, double *phi)
{
    mySpace->phi2(xsi,eta,phi);
}

void femDiscreteDphi2(femDiscrete* mySpace, double xsi, double eta, double *dphidxsi, double *dphideta)
{
    mySpace->dphi2dx(xsi,eta,dphidxsi,dphideta);
}

void femDiscreteXsi(femDiscrete* mySpace, double *xsi)
{
    mySpace->x(xsi);
}

void femDiscretePhi(femDiscrete* mySpace, double xsi, double *phi)
{
    mySpace->phi(xsi,phi);
}

void femDiscreteDphi(femDiscrete* mySpace, double xsi, double *dphidxsi)
{
    mySpace->dphidx(xsi,dphidxsi);
}

void femDiscretePrint(femDiscrete *mySpace)
{
    int i,j;
    int n = mySpace->n;
    double xsi[4], eta[4], phi[4], dphidxsi[4], dphideta[4];

    if (mySpace->type == FEM_EDGE) {
        femDiscreteXsi(mySpace,xsi);
        for (i=0; i < n; i++) {           
            femDiscretePhi(mySpace,xsi[i],phi);
            femDiscreteDphi(mySpace,xsi[i],dphidxsi);
            for (j=0; j < n; j++)  {
                printf("(xsi=%+.1f) : ",xsi[i]);
                printf(" phi(%d)=%+.1f",j,phi[j]);  
                printf("   dphidxsi(%d)=%+.1f \n",j,dphidxsi[j]); }
            printf(" \n"); }}
    
    if (mySpace->type == FEM_QUAD || mySpace->type == FEM_TRIANGLE) {
        femDiscreteXsi2(mySpace, xsi, eta);
        for (i = 0; i < n; i++)  {    
            femDiscretePhi2(mySpace, xsi[i], eta[i], phi);
            femDiscreteDphi2(mySpace, xsi[i], eta[i], dphidxsi, dphideta);
            for (j = 0; j < n; j++) {  
                printf("(xsi=%+.1f,eta=%+.1f) : ", xsi[i], eta[i]);  
                printf(" phi(%d)=%+.1f", j, phi[j]);
                printf("   dphidxsi(%d)=%+.1f", j, dphidxsi[j]);
                printf("   dphideta(%d)=%+.1f \n", j, dphideta[j]); }
            printf(" \n"); }}   
}

femFullSystem *femFullSystemCreate(int size)
{
    femFullSystem *theSystem = malloc(sizeof(femFullSystem));
    femFullSystemAlloc(theSystem, size);
    femFullSystemInit(theSystem);

    return theSystem; 
}

void femFullSystemFree(femFullSystem *theSystem)
{
    free(theSystem->A);
    free(theSystem->B);
    free(theSystem);
}

void femFullSystemAlloc(femFullSystem *mySystem, int size)
{
    int i;  
    double *elem = malloc(sizeof(double) * size * (size+1)); 
    mySystem->A = malloc(sizeof(double*) * size); 
    mySystem->B = elem;
    mySystem->A[0] = elem + size;  
    mySystem->size = size;
    for (i=1 ; i < size ; i++) 
        mySystem->A[i] = mySystem->A[i-1] + size;
}

void femFullSystemInit(femFullSystem *mySystem)
{
    int i,size = mySystem->size;
    for (i=0 ; i < size*(size+1) ; i++) 
        mySystem->B[i] = 0;}


void femFullSystemPrint(femFullSystem *mySystem)
{
    double  **A, *B;
    int     i, j, size;
    
    A    = mySystem->A;
    B    = mySystem->B;
    size = mySystem->size;
    
    for (i=0; i < size; i++) {
        for (j=0; j < size; j++)
            if (A[i][j] == 0)  printf("         ");   
            else               printf(" %+.1e",A[i][j]);
        printf(" :  %+.1e \n",B[i]); }
}

double* femFullSystemEliminate(femFullSystem *mySystem)
{
    double  **A, *B, factor;
    int     i, j, k, size;
    
    A    = mySystem->A;
    B    = mySystem->B;
    size = mySystem->size;
    
    /* Gauss elimination */
    
    for (k=0; k < size; k++) {
        if ( fabs(A[k][k]) <= 1e-16 ) {
            printf("Pivot index %d  ",k);
            printf("Pivot value %e  ",A[k][k]);
            Error("Cannot eliminate with such a pivot"); }
        for (i = k+1 ; i <  size; i++) {
            factor = A[i][k] / A[k][k];
            for (j = k+1 ; j < size; j++) 
                A[i][j] = A[i][j] - A[k][j] * factor;
            B[i] = B[i] - B[k] * factor; }}
    
    /* Back-substitution */
    
    for (i = size-1; i >= 0 ; i--) {
        factor = 0;
        for (j = i+1 ; j < size; j++)
            factor += A[i][j] * B[j];
        B[i] = ( B[i] - factor)/A[i][i]; }
    
    return(mySystem->B);    
}

void  femFullSystemConstrain(femFullSystem *mySystem, 
                             int myNode, double myValue) 
{
    double  **A, *B;
    int     i, size;
    
    A    = mySystem->A;
    B    = mySystem->B;
    size = mySystem->size;
    
    for (i=0; i < size; i++) {
        B[i] -= myValue * A[i][myNode];
        A[i][myNode] = 0; }
    
    for (i=0; i < size; i++) 
        A[myNode][i] = 0; 
    
    A[myNode][myNode] = 1;
    B[myNode] = myValue;
}


femProblem *femElasticityCreate(femGeo* theGeometry, 
                  double E, double nu, double rho, double gx, double gy, femElasticCase iCase)
{
    femProblem *theProblem = malloc(sizeof(femProblem));
    theProblem->E   = E;
    theProblem->nu  = nu;
    theProblem->gx   = gx;
    theProblem->gy   = gy;
    theProblem->rho = rho;
    
    if (iCase == PLANAR_STRESS) {
        theProblem->A = E/(1-nu*nu);
        theProblem->B = E*nu/(1-nu*nu);
        theProblem->C = E/(2*(1+nu)); }
    else if (iCase == PLANAR_STRAIN || iCase == AXISYM) {
        theProblem->A = E*(1-nu)/((1+nu)*(1-2*nu));
        theProblem->B = E*nu/((1+nu)*(1-2*nu));
        theProblem->C = E/(2*(1+nu)); }

    theProblem->planarStrainStress = iCase;
    theProblem->nBoundaryConditions = 0;
    theProblem->conditions = NULL;
    
    int size = 2*theGeometry->theNodes->nNodes;
    theProblem->constrainedNodes = malloc(size*sizeof(int));
    theProblem->soluce = malloc(size*sizeof(double));
    theProblem->residuals = malloc(size*sizeof(double));
    for (int i=0; i < size; i++) {
        theProblem->constrainedNodes[i] = -1;
        theProblem->soluce[i] = 0.0;
        theProblem->residuals[i] = 0.0;}


    
    theProblem->geometry = theGeometry;  
    if (theGeometry->theElements->nLocalNode == 3) {
        theProblem->space    = femDiscreteCreate(3,FEM_TRIANGLE);
        theProblem->rule     = femIntegrationCreate(3,FEM_TRIANGLE); }
    if (theGeometry->theElements->nLocalNode == 4) {
        theProblem->space    = femDiscreteCreate(4,FEM_QUAD);
        theProblem->rule     = femIntegrationCreate(4,FEM_QUAD); }
    theProblem->spaceEdge    = femDiscreteCreate(2,FEM_EDGE);
    theProblem->ruleEdge     = femIntegrationCreate(2,FEM_EDGE); 
    theProblem->system       = femFullSystemCreate(size); 

    femDiscretePrint(theProblem->space);   
    femDiscretePrint(theProblem->spaceEdge);  
  
    return theProblem;
}

femProblem *femElasticityRead(femGeo *theGeometry, const char *filename) 
{
  FILE *file = NULL;
  if(strcmp(filename,"") == 0)
    file = fopen("../data/problem.txt","r");
  else
    file = fopen(filename,"r");
  if (!file) {
    printf("Error at %s:%d\nUnable to open file %s\n", __FILE__, __LINE__, filename);
    return NULL;
  }
  femProblem *theProblem = malloc(sizeof(femProblem));
  theProblem->nBoundaryConditions = 0;
  theProblem->conditions = NULL;

  int nNodes = theGeometry->theNodes->nNodes;
  int size = 2 * nNodes;
  theProblem->soluce = malloc(size * sizeof(double));
  theProblem->residuals = malloc(size * sizeof(double));
  for (int i = 0; i < size; i++) {
    theProblem->soluce[i] = 0.0;
    theProblem->residuals[i] = 0.0;
  }

  theProblem->constrainedNodes = malloc(size * sizeof(int));
  for (int i = 0; i < size; i++) 
    theProblem->constrainedNodes[i] = -1;

  theProblem->geometry = theGeometry;
  if (theGeometry->theElements->nLocalNode == 3) {
    theProblem->space = femDiscreteCreate(3, FEM_TRIANGLE);
    theProblem->rule = femIntegrationCreate(3, FEM_TRIANGLE);
  }
  if (theGeometry->theElements->nLocalNode == 4) {
    theProblem->space = femDiscreteCreate(4, FEM_QUAD);
    theProblem->rule = femIntegrationCreate(4, FEM_QUAD);
  }
  theProblem->spaceEdge = femDiscreteCreate(2, FEM_EDGE);
  theProblem->ruleEdge = femIntegrationCreate(2, FEM_EDGE);
  theProblem->system = femFullSystemCreate(size);

  char theLine[MAXNAME];
  char theDomain[MAXNAME];
  char theArgument[MAXNAME];
  double value1, value2;
  femBoundaryType typeCondition;

  theProblem->gx = 0;
  theProblem->gy = 0;

  while (!feof(file)) {
    ErrorScan(fscanf(file, "%19[^\n]s \n", (char *)&theLine));
    if (strncasecmp(theLine, "Type of problem     ", 19) == 0) {
      ErrorScan(fscanf(file, ":  %[^\n]s \n", (char *)&theArgument));
      if (strncasecmp(theArgument, "Planar stresses", 13) == 0)
        theProblem->planarStrainStress = PLANAR_STRESS;
      if (strncasecmp(theArgument, "Planar strains", 13) == 0)
        theProblem->planarStrainStress = PLANAR_STRAIN;
      if (strncasecmp(theArgument, "Axi-symetric problem", 13) == 0)
        theProblem->planarStrainStress = AXISYM;
    }
    if (strncasecmp(theLine, "Young modulus       ", 19) == 0) {
      ErrorScan(fscanf(file, ":  %le\n", &theProblem->E));
    }
    if (strncasecmp(theLine, "Poisson ratio       ", 19) == 0) {
      ErrorScan(fscanf(file, ":  %le\n", &theProblem->nu));
    }
    if (strncasecmp(theLine, "Mass density        ", 19) == 0) {
      ErrorScan(fscanf(file, ":  %le\n", &theProblem->rho));
    }
    if (strncasecmp(theLine, "Gravity-X           ", 19) == 0) {
      ErrorScan(fscanf(file, ":  %le\n", &theProblem->gx));
    }
    if (strncasecmp(theLine, "Gravity-Y           ", 19) == 0 || strncasecmp(theLine, "Gravity             ", 19) == 0) {
      ErrorScan(fscanf(file, ":  %le\n", &theProblem->gy));
    }
    if (strncasecmp(theLine, "Boundary condition  ", 19) == 0) {
      ErrorScan(fscanf(file, ":  %19s = %le, %le : %[^\n]s\n", (char *)&theArgument, &value1, &value2, (char *)&theDomain));
      if (strncasecmp(theArgument, "Dirichlet-X", 19) == 0)
      {
        typeCondition = DIRICHLET_X;
        femElasticityAddBoundaryCondition(theProblem, theDomain, typeCondition, value1);
      }
      if (strncasecmp(theArgument, "Dirichlet-Y", 19) == 0)
      {
        typeCondition = DIRICHLET_Y;
        femElasticityAddBoundaryCondition(theProblem, theDomain, typeCondition, value1);
      }
      if (strncasecmp(theArgument, "Dirichlet-XY", 19) == 0)
      {
        typeCondition = DIRICHLET_XY;
        femElasticityAddBoundaryCondition(theProblem, theDomain, DIRICHLET_X, value1);
        femElasticityAddBoundaryCondition(theProblem, theDomain, DIRICHLET_Y, value2);
      }
      if (strncasecmp(theArgument, "Neumann-X", 19) == 0)
      {
        typeCondition = NEUMANN_X;
        femElasticityAddBoundaryCondition(theProblem, theDomain, typeCondition, value1);
      }
      if (strncasecmp(theArgument, "Neumann-Y", 19) == 0)
      {
        typeCondition = NEUMANN_Y;
        femElasticityAddBoundaryCondition(theProblem, theDomain, typeCondition, value1);
      }
      
    }
    ErrorScan(fscanf(file, "\n"));
  }

  int iCase = theProblem->planarStrainStress;
  double E = theProblem->E;
  double nu = theProblem->nu;

  if (iCase == PLANAR_STRESS) {
    theProblem->A = E / (1 - nu * nu);
    theProblem->B = E * nu / (1 - nu * nu);
    theProblem->C = E / (2 * (1 + nu));
  } else if (iCase == PLANAR_STRAIN || iCase == AXISYM) {
    theProblem->A = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));
    theProblem->B = E * nu / ((1 + nu) * (1 - 2 * nu));
    theProblem->C = E / (2 * (1 + nu));
  }

  fclose(file);
  return theProblem;
}

void femElasticityFree(femProblem *theProblem)
{
    femFullSystemFree(theProblem->system);
    femIntegrationFree(theProblem->rule);
    femDiscreteFree(theProblem->space);
    femIntegrationFree(theProblem->ruleEdge);
    femDiscreteFree(theProblem->spaceEdge);
    free(theProblem->conditions);
    free(theProblem->constrainedNodes);
    free(theProblem->soluce);
    free(theProblem->residuals);
    free(theProblem);
}
    
void femElasticityAddBoundaryCondition(femProblem *theProblem, char *nameDomain, femBoundaryType type, double value) // celle là
{
    int iDomain = geoGetDomain(nameDomain);
    if (iDomain == -1)  Error("Undefined domain :-(");

    femBoundaryCondition* theBoundary = malloc(sizeof(femBoundaryCondition));
    theBoundary->domain = theProblem->geometry->theDomains[iDomain];
    theBoundary->value = value;
    theBoundary->type = type;
    theProblem->nBoundaryConditions++;
    int size = theProblem->nBoundaryConditions;
    
    if (theProblem->conditions == NULL)
        theProblem->conditions = malloc(size*sizeof(femBoundaryCondition*));
    else 
        theProblem->conditions = realloc(theProblem->conditions, size*sizeof(femBoundaryCondition*));
    theProblem->conditions[size-1] = theBoundary;
    
        
    if(type == DIRICHLET_X || type == DIRICHLET_Y)
    {
            int shift=-1;
            if (type == DIRICHLET_X)  shift = 0;      
            if (type == DIRICHLET_Y)  shift = 1;  
            if (shift == -1) return; 
            int *elem = theBoundary->domain->elem;
            int nElem = theBoundary->domain->nElem;
            for (int e=0; e<nElem; e++) {
                for (int i=0; i<2; i++) {
                    int node = theBoundary->domain->mesh->elem[2*elem[e]+i];
                    theProblem->constrainedNodes[2*node+shift] = size-1; }}   
    } 
    else if(type == DIRICHLET_XY)
        printf("Warning : DIRICHLET_XY case handled wrong.\n");
}

void femElasticityPrint(femProblem *theProblem)  //celle là
{    
    printf("\n\n ======================================================================================= \n\n");
    printf(" Linear elasticity problem \n");
    printf("   Young modulus   E   = %14.7e [N/m2]\n",theProblem->E);
    printf("   Poisson's ratio nu  = %14.7e [-]\n",theProblem->nu);
    printf("   Density         rho = %14.7e [kg/m3]\n",theProblem->rho);
    printf("   Gravity-x         gx   = %14.7e [m/s2]\n",theProblem->gx);
    printf("   Gravity-y         gy   = %14.7e [m/s2]\n",theProblem->gy);
    
    if (theProblem->planarStrainStress == PLANAR_STRAIN)  printf("   Planar strains formulation \n");
    if (theProblem->planarStrainStress == PLANAR_STRESS)  printf("   Planar stresses formulation \n");
    if (theProblem->planarStrainStress == AXISYM)         printf("   Axisymmetric formulation \n");

    printf("   Boundary conditions : \n");
    for(int i=0; i < theProblem->nBoundaryConditions; i++) {
          femBoundaryCondition *theCondition = theProblem->conditions[i];
          double value = theCondition->value;
          printf("  %20s :",theCondition->domain->name);
          if (theCondition->type==DIRICHLET_X)  printf(" imposing %9.2e as the horizontal displacement  \n",value);
          if (theCondition->type==DIRICHLET_Y)  printf(" imposing %9.2e as the vertical displacement  \n",value); 
          //if (theCondition->type==DIRICHLET_XY) printf( "imposing %9.2e as the horizontal displacement and %9.2e as the vertical displacement  \n",value1,theCondition->value2);
          //Si tout est bien fait il n'y a pas  de truct condition DIRICHLET_XY : on les décompose en DIRICHLET_X et DIRICHLET_Y
          if (theCondition->type==NEUMANN_X)    printf(" imposing %9.2e as the horizontal force desnity \n",value); 
          if (theCondition->type==NEUMANN_Y)    printf(" imposing %9.2e as the vertical force density \n",value);}
    printf(" ======================================================================================= \n\n");
}

double femElasticityIntegrate(femProblem *theProblem, double (*f)(double x, double y)){
    femIntegration *theRule = theProblem->rule;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theMesh = theGeometry->theElements;
    femDiscrete    *theSpace = theProblem->space;

    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,i,map[4];
    int nLocal = theMesh->nLocalNode;
    double value = 0.0;
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (i=0; i < nLocal; i++) {
            map[i]  = theMesh->elem[iElem*nLocal+i];
            x[i]    = theNodes->X[map[i]];
            y[i]    = theNodes->Y[map[i]];} 
        for (iInteg=0; iInteg < theRule->n; iInteg++) {    
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];  
            femDiscretePhi2(theProblem->space,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0; 
            double dydeta = 0.0;
            for (i = 0; i < theSpace->n; i++) {  
                dxdxsi += x[i]*dphidxsi[i];       
                dxdeta += x[i]*dphideta[i];   
                dydxsi += y[i]*dphidxsi[i];   
                dydeta += y[i]*dphideta[i]; }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            for (i = 0; i < theProblem->space->n; i++) {    
                value += phi[i] * f(x[i],y[i]) * jac * weight; }}}
    return value;

}


int femSolutionWrite(int nNodes, int nfields, double *data, const char *filename) {
  FILE *file = fopen(filename, "w");
  if (!file) {
    printf("Error at %s:%d\nUnable to open file %s\n", __FILE__, __LINE__, filename);
    return -1;
  }
  fprintf(file, "Size %d,%d\n", nNodes, nfields);
  for (int i = 0; i < nNodes; i++) {
    for (int j = 0; j < nfields - 1; j++) {
      fprintf(file, "%.18le,", data[i * nfields + j]);
    }
    fprintf(file, "%.18le", data[i * nfields + nfields - 1]);
    fprintf(file, "\n");
  }
  fclose(file);
}

double femMin(double *x, int n) 
{
    double myMin = x[0];
    int i;
    for (i=1 ;i < n; i++) 
        myMin = MIN(myMin,x[i]);
    return myMin;
}

double femMax(double *x, int n) 
{
    double myMax = x[0];
    int i;
    for (i=1 ;i < n; i++) 
        myMax = MAX(myMax,x[i]);
    return myMax;
}


void femError(char *text, int line, char *file)                                  
{ 
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
    exit(69);                                                 
}

void femErrorScan(int test, int line, char *file)                                  
{ 
    if (test >= 0)  return;
    
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in fscanf or fgets in %s at line %d : \n", file, line);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");   
    exit(69);                                       
}

void femWarning(char *text, int line, char *file)                                  
{ 
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Warning in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");                                              
}
