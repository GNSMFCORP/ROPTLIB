#ifndef ESTIMATEROTATION_H
#define ESTIMATEROTATION_H

#include "Manifolds/Stiefel/Stiefel.h"
#include "Manifolds/Stiefel/StieVariable.h"
#include "Manifolds/Stiefel/StieVector.h"
#include "Problems/Problem.h"
#include "Manifolds/SharedSpace.h"
#include "Others/def.h"

#include <opencv2/core.hpp>

using namespace cv;

/*Define the namespace*/
namespace ROPTLIB{
class estimateRotation:public Problem
{
public:
    estimateRotation(Mat, Mat, int, int);

    //functions to be overloaded
    virtual double f(Variable *x) const;
    virtual void EucGrad(Variable *x, Vector *egf) const;
    virtual void EucHessianEta(Variable *x, Vector *etax, Vector *exix) const;

    Mat A_,B_;
    mutable Mat E_; //error matrix

    int n_;
    int p_;
};
};

#endif // ESTIMATEROTATION_H
