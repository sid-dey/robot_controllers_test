#include<iostream>
#include<eigen3/Eigen/Dense>
#include<eigen3/Eigen/Geometry>


Eigen::Matrix<double, 7, 1> m_rnea(Eigen::VectorXd q, Eigen::VectorXd qd, Eigen::VectorXd q_aux_d, Eigen::VectorXd qdd, bool use_gravity) {

    /****************************** manually adding Fetch robot params **************************/
    // need to setup struct for it
    // link masses
    Eigen::VectorXd mass(7) ;     
    mass << 2.558700000000000,   2.661500000000000,   2.331100000000000,   2.129900000000000,   1.656300000000000,   1.725000000000000,   1.812500000000000;
    
    // center of mass for each link
    Eigen::Matrix<double, 3, 7> com ; 
    com << 
    0.092700000000000,    0.143200000000000,    0.116500000000000,    0.127900000000000,    0.109700000000000,    0.088200000000000,    0.078493017931034,
   -0.005600000000000,    0.007200000000000,    0.001400000000000,    0.007300000000000,   -0.026600000000000,    0.000900000000000,   -0.000053842758621,
    0.056400000000000,   -0.000100000000000,    0.000000000000000,    0.000000000000000,    0.000000000000000,   -0.000100000000000,   -0.001438251034483 ;
    
    // inertia wrt to center of mass frame
    Eigen::Matrix<double, 21, 3> I_t;
    I_t <<
    0.004300000000000,  -0.000100000000000,   0.001000000000000,
  -0.000100000000000,   0.008700000000000,  -0.000100000000000,
   0.001000000000000,  -0.000100000000000,   0.008700000000000,
   0.002800000000000,  -0.002100000000000,                   0,
  -0.002100000000000,   0.011100000000000,                   0,
                   0,                   0,   0.011200000000000,
   0.001900000000000,  -0.000100000000000,                   0,
  -0.000100000000000,   0.004500000000000,                   0,
                   0,                   0,   0.004700000000000,
   0.002400000000000,  -0.001600000000000,                   0,
  -0.001600000000000,   0.008200000000000,                   0,
                   0,                   0,   0.008400000000000,
   0.001600000000000,  -0.000300000000000,                   0,
  -0.000300000000000,   0.003000000000000,                   0,
                   0,                   0,   0.003500000000000,
   0.001800000000000,  -0.000100000000000,                   0,
  -0.000100000000000,   0.004200000000000,                   0,
                   0,                   0,   0.004200000000000,
   0.005438647027291,   0.000003426633880,  -0.000007138806433,
   0.000003426633880,   0.003621420239448,  -0.000000106784082,
  -0.000007138806433,  -0.000000106784082,   0.004158783836627;

    Eigen::Matrix<double, 3, 21> I = I_t.transpose();

//     I << 
//     0.0043,   -0.0001,    0.0010,    0.0028,   -0.0021,         0,    0.0019,   -0.0001,         0,    0.0024,   -0.0016,         0,    0.0016,   -0.0003,         0,   0.0018,   -0.0001,         0,    0.0033,    0.0000,    0.0000,
//    -0.0001,    0.0088,   -0.0001,   -0.0021,    0.0112,         0,   -0.0001,    0.0045,         0,   -0.0016,    0.0083,         0,   -0.0003,    0.0030,         0,  -0.0001,    0.0042,         0,    0.0000,    0.0126,   -0.0000,
//     0.0010,   -0.0001,    0.0088,         0,         0,    0.0113,         0,         0,    0.0047,         0,         0,    0.0085,         0,         0,    0.0035,   0     ,    0     ,    0.0042,    0.0000,   -0.0000,    0.0132 ;  
    
    // number of active joints
    const int num_joints = 7; // needs to be const to appease the compiler
    const int num_bodies = 7; // needs to be const to appease the compiler
    
    // use interval arithmetic
    // bool use_interval = false;
    
    // fixed transforms
    Eigen::Matrix<double, 4, 28> H ;
    H << 
    1.0000,         0,         0,   0.03265,    1.0000,         0,         0,    0.1170,    1.0000,         0,         0,    0.2190,    1.0000,         0,         0,    0.1330,    1.0000,         0,         0,    0.1970,    1.0000,         0,         0,    0.1245,    1.0000,         0,         0,    0.1385,
         0,    1.0000,         0,         0,         0,    1.0000,         0,         0,         0,    1.0000,         0,         0,         0,    1.0000,         0,         0,         0,    1.0000,         0,         0,         0,    1.0000,         0,         0,         0,    1.0000,         0,         0,
         0,         0,    1.0000,   0.72601,         0,         0,    1.0000,    0.0600,         0,         0,    1.0000,         0,         0,         0,    1.0000,         0,         0,         0,    1.0000,         0,         0,         0,    1.0000,         0,         0,         0,    1.0000,         0,
         0,         0,         0,    1.0000,         0,         0,         0,    1.0000,         0,         0,         0,    1.0000,         0,         0,         0,    1.0000,         0,         0,         0,    1.0000,         0,         0,         0,    1.0000,         0,         0,         0,    1.0000 ;
    // H << 
    // 1.0000,         0,         0,         0,    1.0000,         0,         0,    0.1170,    1.0000,         0,         0,    0.2190,    1.0000,         0,         0,    0.1330,    1.0000,         0,         0,    0.1970,    1.0000,         0,         0,    0.1245,    1.0000,         0,         0,    0.1385,
    //      0,    1.0000,         0,         0,         0,    1.0000,         0,         0,         0,    1.0000,         0,         0,         0,    1.0000,         0,         0,         0,    1.0000,         0,         0,         0,    1.0000,         0,         0,         0,    1.0000,         0,         0,
    //      0,         0,    1.0000,         0,         0,         0,    1.0000,    0.0600,         0,         0,    1.0000,         0,         0,         0,    1.0000,         0,         0,         0,    1.0000,         0,         0,         0,    1.0000,         0,         0,         0,    1.0000,         0,
    //      0,         0,         0,    1.0000,         0,         0,         0,    1.0000,         0,         0,         0,    1.0000,         0,         0,         0,    1.0000,         0,         0,         0,    1.0000,         0,         0,         0,    1.0000,         0,         0,         0,    1.0000 ;
    
    Eigen::Matrix<double, 3, 7> joint_axes ;
    joint_axes << 
    0,     0,     1,     0,     1,     0,     1,
    0,     1,     0,     1,     0,     1,     0,
    1,     0,     0,     0,     0,     0,     0 ;


    Eigen::Vector3d gravity_vector = {0, 0, -9.81} ;

    // Eigen::Matrix<double, 7, 1> joint_pos = q; 
    Eigen::Matrix<double, 7, 1> joint_vel = qd;
    Eigen::Matrix<double, 7, 1> joint_acc = qdd;
    Eigen::Matrix<double, 7, 1> joint_vel_aux = q_aux_d;

    /* setup reference frames */
    // rotation axis of base frame
    // Eigen::Vector3d z0 = {0, 0, 1};

    // orientation of frame i with respect to frame i-1
    // R = repmat(eye(3), [1, 1, num_joints+1]); // Need to make 3 dimensional matrix for this (or write it more cleverly as a block 2D matrix)
    Eigen::Matrix<double, 3 * (num_joints+1), 3> R ;
    for (size_t i = 0; i < num_joints + 1 ; i++) {
        R.block(3*i, 0, 3, 3) << 1, 0, 0,
                                0, 1, 0,
                                0, 0, 1 ;
    }
    // std::cout << "R matrix\n" << R << std::endl ;


    // position of the origin of frame i with respect to frame i-1
    Eigen::Matrix<double, 3, num_joints+1> P ;
    P.setZero() ;

    // orientation of frame i-1 with respect to frame i
    // R_t = repmat(eye(3), [1, 1, num_joints]); // Need to make 3 dimensional matrix for this (or write it more cleverly as a block 2D matrix)
    Eigen::Matrix<double, 3 * (num_joints), 3> R_t ;
    for (size_t i = 0; i < num_joints ; i++) {
        R_t.block(3*i, 0, 3, 3) << 1, 0, 0,
                                0, 1, 0,
                                0, 0, 1 ;
    }    

    // Frame {i} axis of rotation expressed in Frame {i}
    Eigen::Matrix<double, 3, num_joints> z ;
    z.setZero() ;
    Eigen::Matrix4d T ;

    // calculate frame-to-frame transformations based on DH table
    for (int i = 0; i < num_joints; i++) {
        // get rotation axis and matrix
        // int idx = find(joint_axes(:,i)); // there m,ust be a better way, but oh well 
        int idx = 3;
        Eigen::Matrix4d Rot ;
        // equivalent to find function in MATLAB
        for(size_t j=0; j<=2; j++) {
            if (joint_axes(j,i) == 1) {
                idx = j + 1 ;
                break ;
            }
        }
        if (idx == 1) {
            // Rot = rx(q(i));
            Rot << 1, 0, 0, 0, 
                    0, cos(q(i)), -sin(q(i)), 0, 
                    0, sin(q(i)), cos(q(i)), 0, 
                    0, 0, 0, 1 ;
        }            
        else if (idx == 2) {
            // Rot = ry(q(i));
            Rot << cos(q(i)), 0, sin(q(i)), 0, 
                    0, 1, 0, 0, 
                    -sin(q(i)), 0, cos(q(i)), 0, 
                    0, 0, 0, 1 ;
        }            
        else {
            // Rot = rz(q(i));
            Rot << cos(q(i)), -sin(q(i)), 0, 0,
                    sin(q(i)), cos(q(i)), 0, 0,
                    0, 0, 1, 0, 
                    0, 0, 0, 1 ;
        }

        
        T = H.block(0, 4*i, 4, 4) * Rot; // need to find out how H is defined/represented
        
        // orientation of Frame {i} with respect to Frame {i-1}
        R.block(3*i, 0, 3, 3) = T.block(0, 0, 3, 3);    // consider using Eigen::seq, seqN, last, lastN
        R_t.block(3*i, 0, 3, 3) = R.block(3*i, 0, 3, 3).transpose(); // line 7
        
        // position of Frame {i} with respect to Frame {i-1}
        P.col(i) = T.block(0, 3, 3, 1);
        
        // orientation of joint i axis of rotation with respect to Frame {i}
        // z.col(i) << 0, 0, 1;
        z.col(i) = joint_axes.col(i);
        // z(:,i) = robot.Bodies{i}.Joint.JointAxis' ;
    }
 
    // get transform to end-effector
    if (num_bodies > num_joints) {
        R.block(R.rows()-3, 0, 3, 3) = T.block(0, 0, 3, 3) ; //what is T? a cell array of 4x4 matrices?
        P.col(P.cols()-1) = T.block(0, 3, 3, 1) ;
    }
    
    /* INITIALIZE */
    // base link/frame
    Eigen::Vector3d w0 = Eigen::Vector3d::Zero();
    Eigen::Vector3d w0dot = Eigen::Vector3d::Zero();
    Eigen::Vector3d linear_acc0 = Eigen::Vector3d::Zero();
    Eigen::Vector3d w0_aux = Eigen::Vector3d::Zero(); // auxilliary

    // set gravity
    if (use_gravity == true) {
        linear_acc0 = -gravity_vector;
    }

    // angular velocity/acceleration
    Eigen::Matrix<double, 3, num_joints> w ; 
    Eigen::Matrix<double, 3, num_joints> wdot ;
    Eigen::Matrix<double, 3, num_joints> w_aux ; 
    w.setZero() ;
    wdot.setZero();
    w_aux.setZero() ;

    // linear acceleration of frame
    Eigen::Matrix<double, 3, num_joints> linear_acc ;
    linear_acc.setZero() ;
    
    
    // linear acceleration of com
    Eigen::Matrix<double, 3, num_joints> linear_acc_com ;

    // link forces/torques
    Eigen::Matrix<double, 3, num_joints> F ;
    Eigen::Matrix<double, 3, num_joints> N ;
    
    // intialize f, n, u
    Eigen::Matrix<double, 3, num_joints + 1> f ;
    Eigen::Matrix<double, 3, num_joints + 1> n ;
    Eigen::Matrix<double, num_joints, 1> u ;

    linear_acc_com.setZero() ;
    F.setZero() ;
    N.setZero() ;
    f.setZero() ;
    n.setZero() ;
    u.setZero() ;
    
    /* RNEA forward recursion */
    for (int i = 0; i < num_joints; i++) {
        // MATLAB has no zero indexing, but C++ does...
        if (i == 0) {
            // (6.45) angular velocity
            // w(:,i) = R_t(:,:,i) * w0 + joint_vel(i)*z(:,i); // line 13
            // Eigen::Matrix3d vbn = R_t.block(3*i, 0, 3, 3);
            w.col(i) = R_t.block(3*i, 0, 3, 3) * w0 + joint_vel(i)*z.col(i) ;

            // auxillary angular velocity
            // w_aux(:,i) = R_t(:,:,i) * w0_aux + joint_vel_aux(i)*z(:,i); // line 13
            w_aux.col(i) = R_t.block(3*i, 0, 3, 3) * w0_aux + joint_vel_aux(i)*z.col(i); 

            // (6.46) angular acceleration
            // wdot.col(i) = R_t.block(3*i, 0, 3, 3) * w0dot  // line 15
            //             + (R_t.block(3*i, 0, 3, 3) * w0_aux).cross(joint_vel(i)*z.col(i)) 
            //             + joint_acc(i)*z.col(i);

            Eigen::Vector3d temp_vect1, temp_vect2 ;
            temp_vect1 = R_t.block(3*i, 0, 3, 3) * w0_aux ;
            temp_vect2 = joint_vel(i)*z.col(i) ;
            wdot.col(i) = R_t.block(3*i, 0, 3, 3) * w0dot  // line 15
                        + temp_vect1.cross(temp_vect2) 
                        + joint_acc(i)*z.col(i);

                    
            // (6.47) linear acceleration        
            linear_acc.col(i) = R_t.block(3*i, 0, 3, 3) * (linear_acc0 + w0dot.cross(P.col(i)) 
                                + w0.cross(w0.cross(P.col(i)))) ;    // line 16 (TYPO IN PAPER)
        }
        else {
            // (6.45) angular velocity
            w.col(i) = R_t.block(3*i, 0, 3, 3) * w.col(i-1) + joint_vel(i)*z.col(i); // line 13


            // auxillar angular velocity
            w_aux.col(i) = R_t.block(3*i, 0, 3, 3) * w_aux.col(i-1) + joint_vel_aux(i)*z.col(i); // line 14


            // // (6.46) angular acceleration
            // wdot.col(i) = R_t.block(3*i, 0, 3, 3) * wdot.col(i-1)  // line 15
            //                     + (R_t.block(3*i, 0, 3, 3) * w_aux.col(i-1)).cross(joint_vel(i)*z.col(i))
            //                     + joint_acc(i)*z.col(i);
        
            // (6.46) angular acceleration
            Eigen::Vector3d temp_vect3, temp_vect4 ;
            temp_vect3 = R_t.block(3*i, 0, 3, 3) * w_aux.col(i-1) ;
            temp_vect4 = joint_vel(i)*z.col(i) ;
            wdot.col(i) = R_t.block(3*i, 0, 3, 3) * wdot.col(i-1)  // line 15
                                + (temp_vect3).cross(temp_vect4)
                                + joint_acc(i)*z.col(i);

                            
            // (6.47) linear acceleration
            linear_acc.col(i) = R_t.block(3*i, 0, 3, 3)*(linear_acc.col(i-1) 
                                    + wdot.col(i-1).cross(P.col(i)) // line 16 (TYPO IN PAPER)
                                    + w.col(i-1).cross(w_aux.col(i-1).cross(P.col(i))));
        }

        // (6.48) linear acceleration of CoM auxilliary
        linear_acc_com.col(i) = linear_acc.col(i)  // line 23 (modified for standard RNEA)
                        + wdot.col(i).cross(com.col(i))  
                        + w.col(i).cross(w_aux.col(i).cross(com.col(i)));                

        // (6.49) calculate forces
        F.col(i) = mass[i] * linear_acc_com.col(i); // line 27

        // (6.50) calculate torques
        // N.col(i) = I{i}(:,:)*wdot.col(i)  // calculated in line 29
        //          + w_aux.col(i).cross( (I{i}(:,:)*w.col(i)));
        
        // N.col(i) = I.block(0, 3*i, 3, 3)*wdot.col(i)  // calculated in line 29
        //          + w_aux.col(i).cross( (I.block(0, 3*i, 3, 3)*w.col(i)));

        Eigen::Vector3d temp_vect5 ;
        temp_vect5 = I.block(0, 3*i, 3, 3) * w.col(i) ;
        N.col(i) = I.block(0, 3*i, 3, 3)*wdot.col(i)  // calculated in line 29
            + w_aux.col(i).cross(temp_vect5);

    }    

    /*  RNEA reverse recursion  */
    for (int i = num_joints-1; i >= 0; i--) {
                // (6.51)
        // f.col(i) = R(:,:,i+1) * f.col(i+1) + F.col(i); // line 28
        f.col(i) = R.block(3*i+3, 0, 3, 3) * f.col(i+1) + F.col(i); // line 28

        // // (6.52)
        // n.col(i) = N.col(i) 
        //        + R.block(3*i+3, 0, 3, 3) * n.col(i+1)  // line 29
        //        + com.col(i).cross(F.col(i))  // P(:,i) might not be right
        //        + P.col(i+1).cross(R.block(3*i+3, 0, 3, 3) * f.col(i+1)); // line 29 (TYPO IN PAPER)

        // (6.52)
        Eigen::Vector3d temp_vect6 ;
        temp_vect6 = R.block(3*i+3, 0, 3, 3) * f.col(i+1) ;
        n.col(i) = N.col(i) 
               + R.block(3*i+3, 0, 3, 3) * n.col(i+1)  // line 29
               + com.col(i).cross(F.col(i))  // P(:,i) might not be right
               + P.col(i+1).cross(temp_vect6); // line 29 (TYPO IN PAPER)

    }

    // std::cout << "u before : " << u << std::endl ;

    // calculate joint torques
    for (int i = 0; i < num_joints; i++) {
        // (6.53)
        u[i] = n.col(i).transpose() * z.col(i); // line 31
    }
    return u;
}