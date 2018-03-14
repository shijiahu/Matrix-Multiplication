//
//  main.cpp
//  matrix
//
//  Created by shijia hu on 3/13/17.
//  Copyright © 2017 shijia hu. All rights reserved.
//

#include <iostream>
#include <time.h>
#include <math.h>
//#include <limits.h>

#define MAX_INTEGER 32767

//To make you not handle the potential overflow problem in the integer range, I want to restrain the array size n up to 1024 (n<=1024). Please make the maximum random integer for each array size n as floor(sqrt(MAX_INTEGER/n)).  	32767
using namespace std;



//standard matrix multiplication
void StanMatrixMul(int **matrixA, int **matrixB, int **matrixC, int n)
{
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            matrixC[i][j] = 0;
            for(int k = 0; k < n; ++k)
            {
                matrixC[i][j] += matrixA[i][k] * matrixB[k][j];
            }
        }
    }
}

void AddMatrix(int **matrixA, int **matrixB, int **matrixC, int n){
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            matrixC[i][j] = matrixA[i][j] + matrixB[i][j];
        }
    }
}

void SubMatrix(int **matrixA, int **matrixB, int **matrixC, int n){
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            matrixC[i][j] = matrixA[i][j] - matrixB[i][j];
        }
    }
}

void StrassenMatrixMul(int **matrixA, int **matrixB, int **matrixC, int n){
    
    if(n == 1)
    {
        matrixC[0][0] = matrixA[0][0] * matrixB[0][0];
    }else{
    
    
        int **A11,**A12,**A21,**A22,**B11,**B12,**B21,**B22,**C11,**C12,**C21,**C22,**M1,**M2,**M3,**M4,**M5,**M6,**M7;
        int **matrixTemp1, **matrixTemp2;
        
        A11 = new int*[n/2];
        A12 = new int*[n/2];
        A21 = new int*[n/2];
        A22 = new int*[n/2];
        
        B11 = new int*[n/2];
        B12 = new int*[n/2];
        B21 = new int*[n/2];
        B22 = new int*[n/2];
        
        C11 = new int*[n/2];
        C12 = new int*[n/2];
        C21 = new int*[n/2];
        C22 = new int*[n/2];
        
        M1 = new int*[n/2];
        M2 = new int*[n/2];
        M3 = new int*[n/2];
        M4 = new int*[n/2];
        M5 = new int*[n/2];
        M6 = new int*[n/2];
        M7 = new int*[n/2];
        
        matrixTemp1 = new int*[n/2];
        matrixTemp2 = new int*[n/2];
        
        for(int i=0;i<n/2;i++)
        {
            A11[i] = new int[n/2];
            A12[i] = new int[n/2];
            A21[i] = new int[n/2];
            A22[i] = new int[n/2];
            B11[i] = new int[n/2];
            B12[i] = new int[n/2];
            B21[i] = new int[n/2];
            B22[i] = new int[n/2];
            C11[i] = new int[n/2];
            C12[i] = new int[n/2];
            C21[i] = new int[n/2];
            C22[i] = new int[n/2];
            
            M1[i] = new int[n/2];
            M2[i] = new int[n/2];
            M3[i] = new int[n/2];
            M4[i] = new int[n/2];
            M5[i] = new int[n/2];
            M6[i] = new int[n/2];
            M7[i] = new int[n/2];
            
            matrixTemp1[i] = new int[n/2];
            matrixTemp2[i] = new int[n/2];
            
            memset(A11[i],0,n/2*sizeof(int));
            memset(A12[i],0,n/2*sizeof(int));
            memset(A21[i],0,n/2*sizeof(int));
            memset(A22[i],0,n/2*sizeof(int));
            memset(B11[i],0,n/2*sizeof(int));
            memset(B12[i],0,n/2*sizeof(int));
            memset(B21[i],0,n/2*sizeof(int));
            memset(B22[i],0,n/2*sizeof(int));
            memset(C11[i],0,n/2*sizeof(int));
            memset(C12[i],0,n/2*sizeof(int));
            memset(C22[i],0,n/2*sizeof(int));
            memset(M1[i],0,n/2*sizeof(int));
            memset(M2[i],0,n/2*sizeof(int));
            memset(M3[i],0,n/2*sizeof(int));
            memset(M4[i],0,n/2*sizeof(int));
            memset(M5[i],0,n/2*sizeof(int));
            memset(M6[i],0,n/2*sizeof(int));
            memset(M7[i],0,n/2*sizeof(int));
            memset(matrixTemp1[i],0,n/2*sizeof(int));
            memset(matrixTemp2[i],0,n/2*sizeof(int));
        }
        
        for(int i=0;i<n/2;i++)
        {
            for(int j=0;j<n/2;j++)
            {
                A11[i][j] = matrixA[i][j];
                A12[i][j] = matrixA[i][j+n/2];
                A21[i][j] = matrixA[i+n/2][j];
                A22[i][j] = matrixA[i+n/2][j+n/2];
                
                B11[i][j] = matrixB[i][j];
                B12[i][j] = matrixB[i][j+n/2];
                B21[i][j] = matrixB[i+n/2][j];
                B22[i][j] = matrixB[i+n/2][j+n/2];
            }
        }

        //calculate M1
        AddMatrix(A11,A22,matrixTemp1,n/2);
        AddMatrix(B11,B22,matrixTemp2,n/2);
        StrassenMatrixMul(matrixTemp1,matrixTemp2,M1,n/2);
        
        //calculate M2
        AddMatrix(A21,A22,matrixTemp1,n/2);
        StrassenMatrixMul(matrixTemp1,B11,M2,n/2);
        
        //calculate M3
        SubMatrix(B12,B22,matrixTemp1,n/2);
        StrassenMatrixMul(A11,matrixTemp1,M3,n/2);
        
        //calculate M4
        SubMatrix(B21,B11,matrixTemp1,n/2);
        StrassenMatrixMul(A22,matrixTemp1,M4,n/2);
        
        //calculate M5
        AddMatrix(A11,A12,matrixTemp1,n/2);
        StrassenMatrixMul(matrixTemp1,B22,M5,n/2);
        
        //calculate M6
        SubMatrix(A21,A11,matrixTemp1,n/2);
        AddMatrix(B11,B12,matrixTemp2,n/2);
        StrassenMatrixMul(matrixTemp1,matrixTemp2,M6,n/2);
        
        //calculate M7
        SubMatrix(A12,A22,matrixTemp1,n/2);
        AddMatrix(B21,B22,matrixTemp2,n/2);
        StrassenMatrixMul(matrixTemp1,matrixTemp2,M7,n/2);
        
        //C11
        AddMatrix(M1,M4,C11,n/2);
        SubMatrix(C11,M5,C11,n/2);
        AddMatrix(C11,M7,C11,n/2);
        
        //C12
        AddMatrix(M3,M5,C12,n/2);
        
        //C21
        AddMatrix(M2,M4,C21,n/2);
        
        //C22
        SubMatrix(M1,M2,C22,n/2);
        AddMatrix(C22,M3,C22,n/2);
        AddMatrix(C22,M6,C22,n/2);
        
        for(int i=0;i<n/2;i++)
        {
            for(int j=0;j<n/2;j++)
            {
                matrixC[i][j] = C11[i][j];
                matrixC[i][j+n/2] = C12[i][j];
                matrixC[i+n/2][j] = C21[i][j];
                matrixC[i+n/2][j+n/2] = C22[i][j];
            }
        }
        
        for(int i=0;i<n/2;i++)
        {
            delete[] A11[i];    delete[] A12[i];    delete[] A21[i];    delete[] A22[i];
            delete[] B11[i];    delete[] B12[i];    delete[] B21[i];    delete[] B22[i];
            delete[] C11[i];    delete[] C12[i];    delete[] C21[i];    delete[] C22[i];
            
            delete[] M1[i]; delete[] M2[i]; delete[] M3[i]; delete[] M4[i]; delete[] M5[i]; delete[] M6[i]; delete[] M7[i];
            delete[] matrixTemp1[i];    delete[] matrixTemp2[i];
        }
        delete[] A11;   delete[] A12;   delete[] A21;   delete[] A22;
        delete[] B11;   delete[] B12;   delete[] B21;   delete[] B22;
        delete[] C11;   delete[] C12;   delete[] C21;   delete[] C22;
        
        delete[] M1;    delete[] M2;    delete[] M3;    delete[] M4;    delete[] M5;    delete[] M6;
        delete[] M7;
        delete[] matrixTemp1;   delete[] matrixTemp2;
    }
}



int main(int argc, const char * argv[]) {
    cout<<"Please enter the size of matrix n:"<<endl;
    int n,temp;
    int **matrixA, **matrixB, **matrixC, **matrixD;
    
    
    cin>>n;
    temp = (int)sqrt(MAX_INTEGER/n);
    
    matrixA=new int *[n];
    matrixB=new int *[n];
    matrixC=new int *[n];
    matrixD=new int *[n];
    
    
    //matrixTemp1 = new int*[n/2];
    //matrixTemp2 = new int*[n/2];
    
    for(int i=0;i<n;i++)
    {
        matrixA[i]=new int [n];
        memset(matrixA[i],0,n*sizeof(int));
        matrixB[i]=new int [n];
        memset(matrixB[i],0,n*sizeof(int));
        matrixC[i]=new int [n];
        memset(matrixC[i],0,n*sizeof(int));
        matrixD[i]=new int [n];
        memset(matrixD[i],0,n*sizeof(int));
    }
    
    srand(time(0));
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            matrixA[i][j] = rand() % temp;
            matrixB[i][j] = rand() % temp;
        }
    }
    
    cout<<endl<<"matrixA:"<<endl;
    for (int i = 0; i<n; i++) {
        for (int j = 0; j<n; j++) {
            cout<<matrixA[i][j]<<"  ";
        }
        cout<<endl;
    }
    cout<<endl<<"matrixB:"<<endl;
    for (int i = 0; i<n; i++) {
        for (int j = 0; j<n; j++) {
            cout<<matrixB[i][j]<<"  ";
        }
        cout<<endl;
    }
    
    StanMatrixMul(matrixA, matrixB, matrixC, n);
    StrassenMatrixMul(matrixA, matrixB, matrixD, n);
    //compare two method
    int flag = 0;
    
    for (int i = 0; i<n; i++) {
        for (int j = 0; j<n; j++) {
            if(matrixC[i][j]!=matrixD[i][j]){
                flag = 1;
            }
        }
    }
    if (flag == 1) {
        cout<<endl<<"Wrong Answer!!!"<<endl;
    }else{
        cout<<endl<<"Compared. Correct Answer!!!"<<endl;
        cout<<endl<<"matrixC:(Standard Matrix Multiplication)"<<endl;
        for (int i = 0; i<n; i++) {
            for (int j = 0; j<n; j++) {
                cout<<matrixC[i][j]<<"  ";
            }
            cout<<endl;
        }
        cout<<endl<<"matrixD:(Strassen's Matrix Multiplication)"<<endl;
        for (int i = 0; i<n; i++) {
            for (int j = 0; j<n; j++) {
                cout<<matrixD[i][j]<<"  ";
            }
            cout<<endl;
        }
        
    }
    
    
    for (int i = 0; i < n; i ++)
    {
        delete[] matrixA[i];
        delete[] matrixB[i];
        delete[] matrixC[i];
        delete[] matrixD[i];
    }
    delete [] matrixA;
    delete [] matrixB;
    delete [] matrixC;
    delete [] matrixD;
    
    return 0;
}












