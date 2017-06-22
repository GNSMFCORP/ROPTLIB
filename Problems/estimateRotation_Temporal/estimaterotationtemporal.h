#ifndef ESTIMATEROTATIONTEMPORAL_H
#define ESTIMATEROTATIONTEMPORAL_H

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
class estimateRotationTemporal:public Problem
{
public:
    estimateRotationTemporal(const vector<Mat>&,const vector<Mat>&,int F, double lambda);

    //functions to be overloaded
    virtual double f(Variable *x) const;
    virtual void EucGrad(Variable *x, Vector *egf) const;
    virtual void EucHessianEta(Variable *x, Vector *etax, Vector *exix) const;

    Mat P2_;
    int F_;
    int lambda_;
    std::vector<Mat> S_list_;
    std::vector<Mat> W_list_;
};
};

#endif // ESTIMATEROTATIONTEMPORAL_H
