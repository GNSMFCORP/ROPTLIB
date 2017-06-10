#include "estimaterotationtemporal.h"

/*Define the namespace*/
namespace ROPTLIB{

estimateRotationTemporal::estimateRotationTemporal(const Mat& S,const Mat& W,int F, double lambda){

    S_ =S;
    W_=W;
    lambda_=lambda;
    F_=F;

    split(S_,S_list_);
    split(W_,W_list_);

    P2_ = Mat::zeros(2,3,CV_64F);
    P2_.at<double>(0,0)=1.0;
    P2_.at<double>(1,1)=1.0;


    std::cout << "P2:"<<P2_<<std::endl;
    std::cout << "P2 size:"<<P2_.size()<<std::endl;

}

//cost function
double estimateRotationTemporal::f(Variable *x) const{
//    std::cout << "c2"<<std::endl;
    const double *xxM = x->ObtainReadData();

    double totalcost=0;

    for(int i =0;i<F_;i++){

        Mat Xi = Mat::zeros(3,3,CV_64F);
        Mat Xi1 = Mat::zeros(3,3,CV_64F);

        int starting_index_Xi = 3*i*(F_+1);
        int starting_index_Xi1 = 3*(i+1)*(F_+1);

        for(int c=0;c<3;c++){
            for(int r=0;r<3;r++){
                Xi.at<double>(r,c) = xxM[starting_index_Xi+r+3*c];

                if(i<F_-1)
                    Xi1.at<double>(r,c) = xxM[starting_index_Xi1+r+3*c];

            }
        }

        totalcost += 0.5*pow(norm(W_list_[i] - P2_*Xi*S_list_[i]),2);

        if(i<F_-1)
            totalcost += 0.5*lambda_*pow(norm(Xi - Xi1),2);

    }

    return totalcost;

}

void estimateRotationTemporal::EucGrad(Variable *x, Vector *egf) const{

//    std::cout << "c3"<<std::endl;
    const double *xxM = x->ObtainReadData();
    double *egfPtr = egf->ObtainWriteEntireData();

    Mat g = Mat::zeros(3*F_,3*F_,CV_64F);

    for(int i =0;i<F_;i++){

        Mat Xi = Mat::zeros(3,3,CV_64F);
        Mat Xip1 = Mat::zeros(3,3,CV_64F); // i+1
        Mat Xim1 = Mat::zeros(3,3,CV_64F); // i-1

        int starting_index_Xi = 3*i*(F_+1);
        int starting_index_Xip1 = 3*(i+1)*(F_+1);
        int starting_index_Xim1 = 3*(i-1)*(F_+1);

        for(int c=0;c<3;c++){
            for(int r=0;r<3;r++){
                Xi.at<double>(r,c) = xxM[starting_index_Xi+r+3*c];

                if(i<F_-1)
                    Xip1.at<double>(r,c) = xxM[starting_index_Xip1+r+3*c];

                if(i>0)
                    Xim1.at<double>(r,c) = xxM[starting_index_Xim1+r+3*c];

            }
        }
//        std::cout << "c31"<<std::endl;

        Mat Gi = Mat::zeros(3,3,CV_64F);

        Gi = g(Range(3*i,3*(i+1)),Range(3*i,3*(i+1))).clone();
//        int starting_index_Gi = 3*i*(F_+1);
//        for(int c=0;c<3;c++){
//            for(int r=0;r<3;r++){
//                Gi.at<double>(r,c) = xxM[starting_index_Gi+r+3*c];
//            }
//        }
//        std::cout << "c32"<<std::endl;
        Gi = Gi- P2_.t()*W_list_[i]*S_list_[i].t() +
                (P2_.t()*P2_)*Xi*S_list_[i]*S_list_[i].t();

        if(i<F_-1)
            Gi = Gi - lambda_*Xip1;
        if(i>0)
            Gi = Gi - lambda_*Xim1;

        Gi.copyTo(g(Range(3*i,3*(i+1)),Range(3*i,3*(i+1))));
//        std::cout << "c33"<<std::endl;
    }

    for (integer i = 0; i < egf->Getlength(); i++)
        egfPtr[i] = g.at<double>(i%(3*F_),i/(3*F_));

}

void estimateRotationTemporal::EucHessianEta(Variable *x, Vector *etax, Vector *exix) const{


    const double *etaxPtr= etax->ObtainReadData();
    double *exixPtr = exix->ObtainWriteEntireData();

    Mat h = Mat::zeros(3*F_,3*F_,CV_64F); //entire hessian action matrix

    for(int i =0;i<F_;i++){
//        std::cout << "c36"<<std::endl;

        Mat Eta = Mat::zeros(3,3,CV_64F);
        int starting_index = 3*i*(F_+1);

        for(int c=0;c<3;c++){
            for(int r=0;r<3;r++){
                Eta.at<double>(r,c) = etaxPtr[starting_index+r+3*c];
            }
        }

        Mat hessian_action = Mat::zeros(3,3,CV_64F);
        hessian_action = (P2_.t()*P2_)*Eta*(S_list_[i]*S_list_[i].t());
        hessian_action.copyTo(h(Range(3*i,3*(i+1)),Range(3*i,3*(i+1))));

    }
//    std::cout << "c38"<<std::endl;
    for (integer i = 0; i < exix->Getlength(); i++)
        exixPtr[i] = h.at<double>(i%(3*F_),i/(3*F_));

}

}
