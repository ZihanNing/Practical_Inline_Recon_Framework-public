//--------------------------------------------------------------
// file: mf2.cpp - Matlab Mex function application
//
// Debugging notes:
//   (1) In the Project Settings, under the Debug tab, set
//       the "Executable for debug session" to:
//         C:\MATLAB6p1\bin\win32\matlab.exe
//       and set the working directory to:
//         debug
//       and (optionally) set the "Program arguments" to:
//          /nosplash
//   (2) From the menu bar, select Debug->Exceptions to bring up 
//       the exception handling dialog. Change the trapping of 
//       "access violation" from "Stop in not handled" to "Stop
//       always." This will allow the lines of code causing
//       exceptions to be readily identified.
//   (3) When you begin a debug session, Matlab is invoked.
//       Your mex code will be executed when the mex function
//       is called from Matlab.
//   (4) You may use a simple startup.m file to call your mex 
//       function as soon as Matlab starts
//   (5) A virus scanner may cause long delays when Matlab is
//       starting. It may be worthwhile to temporarily disable
//       the virus scanner while you are developing.
// Development Notes:
//   (1) All global variables in your mex function remain in
//       scope after the mex function is executed the first time
//       in a Matlab session.
//   (2) Likewise, all static variables are persistent between
//       function calls.
//   (3) The "clear functions" command may be used to unload
//       your mex function DLL, and force all globals to be
//       reinitialized the next time the mex function is called.
//   (4) Do not use malloc() or free(). Use mxCalloc() and 
//       mxFree() instead. 
//   (5) If you are writing code to be portable to other
//       environments, you may use the MATLAB_MEX_FILE macro
//       to determine if the code is targeting a mex file.
//--------------------------------------------------------------
#include <iostream>
//#include <iostream.h>//Changed
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "graph.h"
#include <string.h>

extern "C" {
#include "mex.h"
}

// TODO: Add you supporting functions here


//--------------------------------------------------------------
// function: mf2 - Entry point from Matlab environment (via 
//   mexFucntion(), below)
// INPUTS:
//   nlhs - number of left hand side arguments (outputs)
//   plhs[] - pointer to table where created matrix pointers are
//            to be placed
//   nrhs - number of right hand side arguments (inputs)
//   prhs[] - pointer to table of input matrices
//--------------------------------------------------------------
void mf2( int nlhs, mxArray *plhs[], int nrhs, const mxArray  *prhs[] )
{
  
  //Declarations
  //int *sourcesinkVal;
  //int *remainVal;
  //int *outArrayFlow,*outArrayAssign;      
  float *sourcesinkVal;
  float *remainVal;
  float *outArrayFlow,*outArrayAssign;
  //double *sourcesinkVal;
  //double *remainVal;
  //double *outArrayFlow,*outArrayAssign;
    
  int dims[1] = {1};
  
  int sourcerowLen, sourcecolLen;
  int remainrowLen, remaincolLen;
  
  //Get matrix sourcesink
  //sourcesinkVal = (int *) mxGetPr(prhs[0]);
  sourcesinkVal = (float *) mxGetPr(prhs[0]);
  //sourcesinkVal = (double *) mxGetPr(prhs[0]);
  sourcerowLen = (int) mxGetN(prhs[0]);
  sourcecolLen = (int) mxGetM(prhs[0]);
	
  //Get matrix remain
  //remainVal    = (int *) mxGetPr(prhs[1]);
  remainVal    = (float *) mxGetPr(prhs[1]);
  //remainVal    = (double *) mxGetPr(prhs[1]);
  remainrowLen = (int) mxGetN(prhs[1]);
  remaincolLen = (int) mxGetM(prhs[1]);	
  
  
  //typedef Graph<int,int,int> GraphType;
  typedef Graph<float,float,float> GraphType;
  //typedef Graph<double,double,double> GraphType;
	GraphType *g = new GraphType(/*estimated # of nodes*/ sourcecolLen, /*estimated # of edges*/ remaincolLen); 

    g->add_node(sourcecolLen);
    
    for(int k=0; k<(sourcecolLen); k++) // sourcesink
      g -> add_tweights(sourcesinkVal[k], sourcesinkVal[sourcecolLen + k], sourcesinkVal[2*sourcecolLen + k]); 
  	
    for(int h=0; h<remaincolLen; h++)	
		g -> add_edge(remainVal[h],remainVal[remaincolLen + h], remainVal[2*remaincolLen + h],remainVal[3*remaincolLen + h]);

    //int flow = g-> maxflow();
	float flow = g -> maxflow();
    //double flow = g->maxflow();


  

  
  //Allocate memory and assign output pointer

  //plhs[0] = mxCreateNumericArray(1,dims, mxINT32_CLASS, mxREAL);
  plhs[0] = mxCreateNumericArray(1,dims, mxSINGLE_CLASS, mxREAL);
  //plhs[0] = mxCreateNumericArray(1,dims, mxDOUBLE_CLASS, mxREAL);
  int dims2[] = {sourcecolLen,2};
  //plhs[1] = mxCreateNumericArray(2,dims2, mxINT32_CLASS, mxREAL);
  plhs[1] = mxCreateNumericArray(2,dims2, mxSINGLE_CLASS, mxREAL);
  //plhs[1] = mxCreateNumericArray(2,dims2, mxDOUBLE_CLASS, mxREAL);
  
  
  //Get a pointer to the data space in our newly allocated memory
  //outArrayFlow = (int*) mxGetPr(plhs[0]);
  outArrayFlow = (float*) mxGetPr(plhs[0]);
  //outArrayFlow = (double*) mxGetPr(plhs[0]);
  
  outArrayFlow[0] = flow;
  
  //Get a pointer to the data space in our newly allocated memory
  //outArrayAssign = (int*) mxGetPr(plhs[1]);
  outArrayAssign = (float*) mxGetPr(plhs[1]);
  //outArrayAssign = (double*) mxGetPr(plhs[1]);

  for(int i=0;i<sourcecolLen;i++)
{
    
	//assigns the number of the nodes
	outArrayAssign[i] = sourcesinkVal[i];
	
    //If o n� becomes assigned to source we assign the value 0
	//Otherwise we assign the value 1.
	
    //if (g->what_segment(nodes[(int)sourcesinkVal[i]]) == Graph::SOURCE)
    
    if (g->what_segment((int)sourcesinkVal[i],GraphType::SINK) == GraphType::SOURCE)//The opposite than in previous version...
		
		outArrayAssign[sourcecolLen + i] = 0;
			
	else
		outArrayAssign[sourcecolLen + i] = 1;
	
}

  
  delete g;

} // end mf2()

extern "C" {
  //--------------------------------------------------------------
  // mexFunction - Entry point from Matlab. From this C function,
  //   simply call the C++ application function, above.
  //--------------------------------------------------------------
  void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray  *prhs[] )
  {
    mf2(nlhs, plhs, nrhs, prhs);
  }



}


