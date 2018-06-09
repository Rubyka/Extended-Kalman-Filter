#include <iostream>
#include "tools.h"
#define PI 3.14159265

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector <VectorXd> &estimations,
                              const vector <VectorXd> &ground_truth) {
    /**
      * Calculate the RMSE here.
    */
    VectorXd rmse(4);
    rmse << 0,0,0,0;

    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    //  * the estimation vector size should equal ground truth vector size
    if(estimations.size() != ground_truth.size()
       || estimations.size() == 0){
        cout << "Invalid estimation or ground_truth data" << endl;
        return rmse;
    }

    //accumulate squared residuals
    for(unsigned int i=0; i < estimations.size(); ++i){

        VectorXd residual = estimations[i] - ground_truth[i];

        //coefficient-wise multiplication
        residual = residual.array()*residual.array();
        rmse += residual;
    }

    //calculate the mean
    rmse = rmse/estimations.size();

    //calculate the squared root
    rmse = rmse.array().sqrt();

    //return the result
    return rmse;
}


MatrixXd Tools::CalculateJacobian(const VectorXd &x_state) {
    /**
      * Calculate a Jacobian here.
    */
    MatrixXd Hj(3,4);

    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    //check division by zero
    if( px == 0 && py == 0 )
    {
        cout << "Error:  division by zero in CalculateJacobian" << endl;
        return Hj;
    }

    //compute the Jacobian
    float rho = sqrt( px*px + py* py );
    float rho2 = rho*rho;
    float rho3 = rho2*rho;
    Hj <<                 px/rho,                    py/rho,      0,      0,
            -py/rho2,                   px/rho2,      0,      0,
            py*( vx*py - vy*px )/rho3, px*( vy*px - vx*py )/rho3, px/rho, py/rho;

    return Hj;
}

VectorXd Tools::CalculateCoordinates(const VectorXd &x) {
    /**
     * Calculate the change of coordinates (cartesian -> polar)
     */

    // Extracting the state info
    float px = x[0];
    float py = x[1];
    float vx = x[2];
    float vy = x[3];

    if (fabs(x[0]+x[1]) < 1e-4) {
        px = 1e-4;
        py = 1e-4;
    }

    // Compute polar
    float rho = sqrt(px * px + py * py);
    float phi = atan2(py, px);
    float rho_dot = (px * vx + py * vy) / rho;

    VectorXd x_polar(3);
    x_polar << rho, phi, rho_dot;

    return x_polar;

}


void Tools::NormalizePhi(VectorXd &y) {
    /**
     * Normalize the angle between -PI and PI
     */
    if( y[1] > PI/2 )
        //y[1] -= 2.f*PI;
        y[1] = y[1] - PI;
    if( y[1] < -PI/2 )
        //y[1] += 2.f*PI;
        y[1] += PI;
}