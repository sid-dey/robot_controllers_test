#include <iostream>
#include "fetch_rnea.cpp"
#include<cmath>
#include<vector>

class FetchArmTorqueController {
    public:        
        Eigen::VectorXd q_, qd_, q_des_, qd_des_, qdd_des_;
        unsigned int num_joints_;
        Eigen::MatrixXd Kr_;        

    FetchArmTorqueController() {
        num_joints_ = 7;
        Kr_ = Eigen::MatrixXd::Identity(num_joints_, num_joints_) * 10;

        //Resize input vectors to the controller
        q_.resize(num_joints_);
        qd_.resize(num_joints_);
        q_des_.resize(num_joints_);
        qd_des_.resize(num_joints_);
        qdd_des_.resize(num_joints_);  
    }

    FetchArmTorqueController(int num_joints, Eigen::MatrixXd Kr) {
        num_joints_ = num_joints;
        Kr_ = Kr;        

        //Resize input vectors to the controller
        q_.resize(num_joints_);
        qd_.resize(num_joints_);
        q_des_.resize(num_joints_);
        qd_des_.resize(num_joints_);
        qdd_des_.resize(num_joints_);  
    }

    /* Torque Control
    * 
    * Inputs
    *   q - Nx1 joint angle state vector
    *   q_d - Nx1 joint velocity state vector
    *   qd - Nx1 reference joint angle
    *   qd_d - Nx1 reference joint velocity
    *   qd_dd - Nx1 reference joint acceleration
    *
    * Outputs
    *   u - Nx1 joint torque vector
    *   
    * Steps
    * 1. Calculate reference terms qd_aux, qdd_aux, r
    * 2. u1 = rnea(q, qd, qd_aux, qdd_aux)
    *    This gives the nominal torques required to follow the trajectory
    * 3. u = u1
    */
    Eigen::VectorXd update(Eigen::VectorXd q, Eigen::VectorXd qd, Eigen::VectorXd q_des, 
                    Eigen::VectorXd qd_des, Eigen::VectorXd qdd_des) {
        //Step 1, calculate reference terms
        Eigen::VectorXd qd_aux = qd_des + Kr_*(q_des - q);
        Eigen::VectorXd qdd_aux = qdd_des + Kr_*(qd_des - qd);        

        //Step 2, calculate nominal torque from passivity RNEA
        Eigen::VectorXd u1 = m_rnea(q, qd, qd_aux, qdd_aux, true);                

        // Step 3. u = u1
        return u1;
    }

    std::vector<double> update(std::vector<double> q_in, std::vector<double> qd_in, std::vector<double> q_des_in, 
                    std::vector<double> qd_des_in, std::vector<double> qdd_des_in) {

        //Copy over into eigen arrays    
        for (size_t i = 0; i < q_in.size(); i++) {
            q_[i] = q_in.at(i);
            qd_[i] = qd_in.at(i);
            q_des_[i] = q_des_in.at(i);
            qd_des_[i] = qd_des_in.at(i);
            qdd_des_[i] = qdd_des_in.at(i);
        }
        
        Eigen::VectorXd tau = update(q_, qd_, q_des_, qd_des_, qdd_des_);        
        
        std::vector<double> control_input ;
        for (size_t i = 0; i < q_in.size(); i++) {
            control_input.push_back(tau[i]) ;
        }

        return control_input ;
    }
};