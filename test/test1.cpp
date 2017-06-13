/*
This is the test file for the Brocokett problem defined in StieBrockett.h and StieBrockett.cpp.

---- WH
*/

#ifndef TESTSIMPLEEXAMPLE_CPP
#define TESTSIMPLEEXAMPLE_CPP

/*Output to console*/
#include <iostream>
/*Generate random number*/
#include "Others/randgen.h"
/*Computational time*/
#include <ctime>

/*Problem related classes*/
#include "Problems/Problem.h"
#include "Problems/estimateRotation/estimaterotation.h"

/*Manifold related classes*/
#include "Manifolds/Manifold.h"
#include "Manifolds/Stiefel/StieVector.h"
#include "Manifolds/Stiefel/StieVariable.h"
#include "Manifolds/Stiefel/Stiefel.h"


/*Trust-region based solvers*/
#include "Solvers/SolversTR.h"
#include "Solvers/RTRNewton.h"

/*The global head file*/
#include "Others/def.h"

//opencv
#include <opencv2/core.hpp>

using namespace ROPTLIB;

int main(void)
{
        // choose a random seed based on current clock
//        unsigned tt = (unsigned)time(NULL);
//        genrandseed(0);

        // size of the Stiefel manifold
        integer n = 3, p = 2;

        double data0[6] = { 0.1858,0.9329,0.1930,0.3907,0.3416,0.2732 };
        cv::Mat A = cv::Mat(3, 2, CV_64F, data0);

        A=A.t();
        std::cout << "matrix A"<<std::endl;
        std::cout << A <<std::endl;

        double data1[4] = { 0.1519,0.3747,0.3971,0.1311 };
        cv::Mat B = cv::Mat(2, 2, CV_64F, data1);

        B=B.t();

        std::cout << "matrix B"<<std::endl;
        std::cout << B <<std::endl;


        // Obtain an initial iterate
        StieVariable StieX(n, p);
        StieX.RandInManifold();

//        double *xxM1 = StieX.ObtainWriteEntireData();

//        std::cout<< "StieX:"<<std::endl;

////        for (integer i = 0; i < StieX.Getlength(); i++)
////        {
////            std::cout<< xxM1[i]<<std::endl;
//////            xxM1[i]=0.1;
////        }

//        xxM1[1]=2;

//        StieX.Print();


        // Define the Stiefel manifold
        Stiefel Domain(n, p);
//        Domain.ChooseStieParamsSet4();

        // Define the Brockett problem
        estimateRotation Prob(A, B, n, p);

        // Set the domain of the problem to be the Stiefel manifold
        Prob.SetDomain(&Domain);

        // output the parameters of the manifold of domain
//        Domain.CheckParams();


        //// test RSD
        //printf("********************************Check all line search algorithm in RSD*****************************************\n");
        //for (integer i = 0; i < INPUTFUN; i++)
        //{
        //	RSD *RSDsolver = new RSD(&Prob, &StieX);
        //	RSDsolver->Debug = FINALRESULT;
        //	RSDsolver->CheckParams();
        //	RSDsolver->Run();
        //	delete RSDsolver;
        //}

        //// test RNewton
        //printf("********************************Check all line search algorithm in RNewton*************************************\n");
        //for (integer i = 0; i < INPUTFUN; i++)
        //{
        //	RNewton *RNewtonsolver = new RNewton(&Prob, &StieX);
        //	RNewtonsolver->LineSearch_LS = static_cast<LSAlgo> (i);
        //	RNewtonsolver->Debug = FINALRESULT;
        //	RNewtonsolver->CheckParams();
        //	RNewtonsolver->Run();
        //	delete RNewtonsolver;
        //}

        //// test RCG
        //printf("********************************Check all Formulas in RCG*************************************\n");
        //for (integer i = 0; i < RCGMETHODSLENGTH; i++)
        //{
        //	RCG *RCGsolver = new RCG(&Prob, &StieX);
        //	RCGsolver->RCGmethod = static_cast<RCGmethods> (i);
        //	RCGsolver->LineSearch_LS = STRONGWOLFE;
        //	RCGsolver->LS_beta = 0.1;
        //	RCGsolver->Debug = FINALRESULT;
        //	RCGsolver->CheckParams();
        //	RCGsolver->Run();
        //	delete RCGsolver;
        //}

        //// test RBroydenFamily
        //printf("********************************Check all Formulas in RBroydenFamily*************************************\n");
        //for (integer i = 0; i < INPUTFUN; i++)
        //{
        //	RBroydenFamily *RBroydenFamilysolver = new RBroydenFamily(&Prob, &StieX);
        //	RBroydenFamilysolver->LineSearch_LS = static_cast<LSAlgo> (i);
        //	RBroydenFamilysolver->Debug = FINALRESULT;
        //	RBroydenFamilysolver->CheckParams();
        //	RBroydenFamilysolver->Run();
        //	delete RBroydenFamilysolver;
        //}

        //// test RWRBFGS
        //printf("********************************Check all line search algorithm in RWRBFGS*************************************\n");
        //for (integer i = 0; i < INPUTFUN; i++)
        //{
        //	RWRBFGS *RWRBFGSsolver = new RWRBFGS(&Prob, &StieX);
        //	RWRBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
        //	RWRBFGSsolver->Debug = FINALRESULT; //ITERRESULT;//
        //	RWRBFGSsolver->CheckParams();
        //	RWRBFGSsolver->Run();
        //	delete RWRBFGSsolver;
        //}

        //// test RBFGS
        //printf("********************************Check all line search algorithm in RBFGS*************************************\n");
        //for (integer i = 0; i < INPUTFUN; i++)
        //{
        //	RBFGS *RBFGSsolver = new RBFGS(&Prob, &StieX);
        //	RBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
        //	RBFGSsolver->Debug = FINALRESULT;
        //	RBFGSsolver->CheckParams();
        //	RBFGSsolver->Run();
        //	delete RBFGSsolver;
        //}

        //// test LRBFGS
        //printf("********************************Check all line search algorithm in LRBFGS*************************************\n");
        //for (integer i = 0; i < INPUTFUN; i++)//LSALGOLENGTH
        //{
        //	LRBFGS *LRBFGSsolver = new LRBFGS(&Prob, &StieX);
        //	LRBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
        //	LRBFGSsolver->Debug = FINALRESULT; //ITERRESULT;//
        //	LRBFGSsolver->CheckParams();
        //	LRBFGSsolver->Run();
        //	delete LRBFGSsolver;
        //}

        //// test RTRSD
        //printf("********************************Check RTRSD*************************************\n");
        //RTRSD RTRSDsolver(&Prob, &StieX);
        //RTRSDsolver.Debug = FINALRESULT;
        //RTRSDsolver.CheckParams();
        //RTRSDsolver.Run();

        // test RTRNewton
        printf("********************************Check RTRNewton*************************************\n");
        RTRNewton RTRNewtonsolver(&Prob, &StieX);
        RTRNewtonsolver.Debug = FINALRESULT;
        RTRNewtonsolver.Max_Iteration=100;
        RTRNewtonsolver.Tolerance=0.001;
//        RTRNewtonsolver.CheckParams();
        RTRNewtonsolver.Run();



        //// test RTRSR1
        //printf("********************************Check RTRSR1*************************************\n");
        //RTRSR1 RTRSR1solver(&Prob, &StieX);
        //RTRSR1solver.Debug = FINALRESULT;
        //RTRSR1solver.CheckParams();
        //RTRSR1solver.Run();

        //// test LRTRSR1
        //printf("********************************Check LRTRSR1*************************************\n");
        //LRTRSR1 LRTRSR1solver(&Prob, &StieX);
        //LRTRSR1solver.Debug = FINALRESULT;
        //LRTRSR1solver.CheckParams();
        //LRTRSR1solver.Run();

        // Check gradient and Hessian
        Prob.CheckGradHessian(&StieX);
        const Variable *xopt = RTRNewtonsolver.GetXopt();
        const double *xxM = xopt->ObtainReadData();

        for (integer i = 0; i < xopt->Getlength(); i++)
        {
            std::cout<< xxM[i]<<std::endl;
        }


        Prob.CheckGradHessian(xopt);

        return 0;
}

#endif // end of TESTSIMPLEEXAMPLE_CPP

