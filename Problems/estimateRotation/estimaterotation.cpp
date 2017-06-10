#include "estimaterotation.h"

/*Define the namespace*/
namespace ROPTLIB{

estimateRotation::estimateRotation(Mat iA, Mat iB, int in, int ip){

    A_ = iA;
    B_ = iB;
    n_ = in;
    p_ = ip;

    E_ = Mat::zeros(iB.rows,iB.cols,CV_64F);

}

//cost function
double estimateRotation::f(Variable *x) const{

    const double *xxM = x->ObtainReadData();

    //convertinf xxM to Mat
    Mat X = Mat::zeros(n_,p_,CV_64F);

    for(int c=0;c<p_;c++){
        for(int r=0;r<n_;r++){
            X.at<double>(r,c) = xxM[c*n_+r];
        }
    }

    E_ = A_*X-B_;

    Mat E_stack = Mat::zeros(E_.rows*E_.cols,1,CV_64F);
    for(int i=0;i<E_.cols;i++){
        E_.col(i).copyTo(E_stack(Range(i*E_.rows,(i+1)*E_.rows),Range::all()));
    }

    Mat EtE = (E_stack.t()*E_stack);
    double cost = EtE.at<double>(0,0)/2;

    return cost;
}

void estimateRotation::EucGrad(Variable *x, Vector *egf) const{


    const double *xxM = x->ObtainReadData();

    //convertinf xxM to Mat
    Mat X = Mat::zeros(n_,p_,CV_64F);

    for(int c=0;c<p_;c++){
        for(int r=0;r<n_;r++){
            X.at<double>(r,c) = xxM[c*n_+r];
        }
    }

    E_ = A_*X-B_;

    Mat Egrad = A_.t()*E_ ;

    double *egfPtr = egf->ObtainWriteEntireData();

    for (integer i = 0; i < egf->Getlength(); i++)
        egfPtr[i] = Egrad.at<double>(i%n_,i/n_);


}

void estimateRotation::EucHessianEta(Variable *x, Vector *etax, Vector *exix) const{


    const double *etaxPtr= etax->ObtainReadData();
    double *exixPtr = exix->ObtainWriteEntireData();

    Mat Eta = Mat::zeros(n_,p_,CV_64F);

    for(int c=0;c<p_;c++){
        for(int r=0;r<n_;r++){
            Eta.at<double>(r,c) = etaxPtr[c*n_+r];
        }
    }

    Mat hessian_action = A_.t()*A_*Eta;

    for (integer i = 0; i < exix->Getlength(); i++)
    {
        exixPtr[i] = hessian_action.at<double>(i%n_,i/n_);
    }

}

}
