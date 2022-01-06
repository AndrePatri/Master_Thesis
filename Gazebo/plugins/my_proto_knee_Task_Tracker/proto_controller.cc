#include <gazebo/gazebo.hh>
#include <gazebo/physics/physics.hh>
#include <gazebo/common/common.hh>
#include <ignition/math/Vector3.hh>
#include <ignition/math/Pose3.hh>
#include <ignition/math/Quaternion.hh>

#include <functional>
#include <cmath>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>

#include <Eigen/Core>
#include <chrono>
#include <limits>
#include "include/QPSolverInterface.cpp" //courtesy of F.Smaldone

#include <boost/numeric/odeint.hpp> // numerical integration (needed)-->use runge_kutta4 class or runge_kutta54_cash_karp with integrate()

using namespace boost::numeric::odeint;

using namespace std; 
using namespace Eigen;

namespace gazebo
{
  class ModelController : public ModelPlugin
  {
    public: ModelController(){//class constructor
    
      this->mu=0.5;//rubber concrete 0.6-0.8-->conservatively 0.5 chosen

      // this->Kp<<92,0,0,100;
      // this->Kd<<42,0,0,4;

      this->Kp<<150,0,0,300;
      this->Kd<<60,0,0,30;
    
      // this->Kp<<100,0.0000,0.0000,200;
      // this->Kd<<15,0.0000,0.0000,20;
      
      this->planar_compensation_Kp<<0,0,0,0;
      this->planar_compensation_Kd<<0,0,0,0;

      this->u_lim<<-2.0,2.0,-20.0,20.0;//input limits
      double beta_lb=M_PI/180.0*62.0;
      double beta_ub=M_PI/180.0*100.0;
      double phi_r_lb=-M_PI/180.0*45.0;
      double phi_r_ub=M_PI/180.0*45.0;      
      this->beta_lim<<beta_lb,beta_ub;// state limits (not working for now)
      this->phi_r_lim<<phi_r_lb,phi_r_ub;
      this->cost_threshold=1;

      this->delta_lambda=1;//used to maintain the vertical reaction force strictly > 0 (>= delta_lambda), [N]
      this->max_step_size=0.001;
      this->sim_steps_over_control_step=20;//number of simulation steps per sample step->fixed the control dt as a multiple of the max_step_size
      this->delta_control=max_step_size*this->sim_steps_over_control_step;
      this->sim_time=20; //amount of time to be simulated, in seconds
      this->N=round(this->sim_time/this->delta_control); //total number of control iterations --> necessary to assign solutions for post-processing
      // Note that,for convenience, only the solution @ the control instants is assigned 
      this->forward_prop_state_meas=1;//if 1 simulate the system as if the controller forward propagated the meas. to try
      // // to compensate for the measurement delay
      this->model_properties_are_manually_set=0;
      this->use_equiv_central_knee=0;

      this->Q_sol=MatrixXf::Zero(this->N+1,10);
      this->U_left=MatrixXf::Zero(this->N,2);
      this->U_right=MatrixXf::Zero(this->N,2);
      this->U_middle=MatrixXf::Zero(this->N,2);
      this->dt_control_calc=VectorXf::Zero(this->N);
      this->Time=VectorXf::Zero(this->N+1);
      this->Chi_ref=MatrixXf::Zero(this->N-1,Chi_ref_k.rows());
      this->Chi_dot_ref=MatrixXf::Zero(this->N-1,Chi_ref_k.rows());
      this->Chi_ddot_ref=MatrixXf::Zero(this->N-1,Chi_ref_k.rows());
      this->Chi=MatrixXf::Zero(this->N-1,Chi_ref_k.rows());
      this->Chi_dot=MatrixXf::Zero(this->N-1,Chi_ref_k.rows());
      this->Chi_ddot=MatrixXf::Zero(this->N-1,Chi_ref_k.rows());
      this->QP_cost_value=VectorXf::Zero(this->N);
      this->X_k=MatrixXf::Zero(this->N,this->A_eq_k.cols());
      
      //this->phi_r0=-0.00032025;
      this->phi_r0=0.0*M_PI/180.0;
    }
    
    public: void Load(physics::ModelPtr _parent, sdf::ElementPtr /*_sdf*/)
    {
      /////////////////PUT HERE EVERY INITIALIZATION WHICH CAN/NEEDS TO BE DONE BEFORE ENTERING THE SIMULATION LOOP/////////////////
      // Store the pointer to the model
      this->model = _parent;
      // Store the pointer to the world where the model is in and gets the grav. acceleration modulus
      this->world=this->model->GetWorld();
      this->gravity=this->world->Gravity();
      this->g=abs(this->gravity[2]);
      // Store the pointer to the used physics engine-->useful for getting simulation parameters which were set in the .world file
      physics_engine_ptr=this->world->Physics();
      physics_engine_ptr->SetMaxStepSize(this->max_step_size); // sets the max_step_size     
      physics_engine_ptr->SetTargetRealTimeFactor(1.0);
      physics_engine_ptr->SetRealTimeUpdateRate(1.0/this->max_step_size);
      this->PointerGetter();// gets all necessary link and joint pointers

      //
       //Important parameters (check the sdf.erb file) for state calculation--> they NEED to be modified if changed in the sdf.erb file
      // These are necessary because joint states are given so that the 0 position is wrt initial relative pose between parent and child link frame
      this->beta_spawn=M_PI/180.0*75.0;
      this->beta_ref_el=M_PI/180.0*90.0;// not easy to retrieve from sdf file, since there is the linkage and those relationships are not easily invertible
      this->phi_w_spawn=0;
      this->phi_r_spawn=0;

      this->ModelPropGetter();// gets necessary model properties from pointed links and joints

      // Listen to the update event. This event is broadcast every
      // simulation iteration-->crucial for the controller
      this->updateConnection = event::Events::ConnectWorldUpdateBegin(
          std::bind(&ModelController::OnUpdate, this));// method to be called before sim update
      this->end_updateConnection = event::Events::ConnectWorldUpdateEnd(
          std::bind(&ModelController::OnUpdateEnd, this));// method to be called after sim update
      this->before_physics_updateConnection = event::Events::ConnectBeforePhysicsUpdate(
          std::bind(&ModelController::BeforePhysicsUpdate, this));// method to be called before physics update
      
          
      this->initial_q_setter();
      this->Q_sol_assigner();// reads the initial state and assignes the initial condition to the first row of Q_sol matrix
      // // this works because sim_iteration_count is initialized to 0      
      }
    
    public:void BeforePhysicsUpdate(){
                /////////////////////////////////////////////////////////////////////////////////////////
      //Reading current time (referred to the end of simulation iteration)
      this->current_time=(this->world->SimTime()).Double();//gets current sim_time as double
      this->sim_iteration_count=this->world->Iterations();//gets current simulation iteration count
      //(it is always the end time of the control step)
      if (current_time>this->sim_time) {
        this->world->SetPaused(true);
        this->WriteToTxT();//write simulation results to .txt file	
        // this->world->Reset();
        //this->world->Stop();

      }//if current time is greater than the required simulation horizon, stop/pause the simulation 
      // (the simulation will never exceed sim_time, due to how sim_time is output)
     
      else{ //simulate
        
        if (((this->sim_iteration_count-1)%this->sim_steps_over_control_step)==0){//Corresponds to a control interval-->compute and apply a new control
         
          this->control_iteration_count=((this->sim_iteration_count-1)/this->sim_steps_over_control_step)+1;
          this->Time(this->control_iteration_count,0)=this->control_iteration_count*this->delta_control;

          if (this->control_iteration_count==1){// first control sample-->the state is not available yet, so no control input

            this->left_wheel_hinge->SetForce( 0,  0.0);
            this->left_ah_d_hinge->SetForce( 0,  0.0);
            this->right_wheel_hinge->SetForce( 0,  0.0); 
            this->right_ah_d_hinge->SetForce( 0,  0.0);

            this->U_left.row(this->control_iteration_count-1)<<0,0;//assigning input
            this->U_right.row(this->control_iteration_count-1)<<0,0;//assigning input
          }

          else{ // state measurement available (one control sample time delayed) --> compute and apply controls 
        
            this->q_k_getter();//gets the state from Q_sol and assignes it to q_k; to be more realistic, this state is delayed by a control sample time
            // by getting the Q_sol row relative to previous simulation iteration
            // In a real system, it might be possible to forward propagate the measurement by integrating the dynamics
            // using the last input as constant forcing term. The time of propagation would be dt_control_sample+dt_control_computation
            // This way it is possible to approximately compensate for the delay, since 
            ///////////////////////////////////////////////////////////////////////////////////////
            chrono::time_point<chrono::high_resolution_clock> start, end;
            start = chrono::high_resolution_clock::now();

            this->Chi_traj_getter();//sets the current reference trajectory
            this->Chi_Chi_dot_k_getter();// gets current value of the Task and its derivative
            this->Chi_ddot_command_calc();// computes the Task command (acceleration) using a PID on the current Task error
            this->QPMatrixCalc();// computes all necessary QP matrices
            this->x_k_calc();
            this->cost_fnctn_val_calc();// computes the cost function-->necessary to check if the QP goes crazy; if it does, do not apply any input or apply the previous one
            this->QP_cost_value(this->control_iteration_count-1,0)=this->cost_fun_value_k;
            this->TaskController_input_calc();

            end= chrono::high_resolution_clock::now();
            std::chrono::duration<double> dt_control = end - start;

            this->dt_control_calc(this->control_iteration_count-1,0)=dt_control.count();

            this->PrintStuff();// prints quantities for debug
          }

        }
        else{//not in a control interval-->use previous input

          this->left_wheel_hinge->SetForce( 0,  this->U_left(this->control_iteration_count-1,0));
          this->left_ah_d_hinge->SetForce( 0,   this->U_left(this->control_iteration_count-1,1));
          this->right_wheel_hinge->SetForce( 0, this->U_right(this->control_iteration_count-1,0)); 
          this->right_ah_d_hinge->SetForce( 0,  this->U_right(this->control_iteration_count-1,1));

        }

      }

      // Trying to keep the motion planar kinematically
      ignition::math::Vector3d current_body_angular_vel=this->body->WorldAngularVel();
      ignition::math::Vector3d angular_vel_correction=ignition::math::Vector3d(0, current_body_angular_vel.Y(),0);
      this->body->SetAngularVel(angular_vel_correction);	//sets to zero vertical angular velocity of body link
      // // ignition::math::Vector3d current_body_linear_vel=this->body->WorldCoGLinearVel();
      // ignition::math::Vector3d linear_vel_correction=ignition::math::Vector3d(current_body_linear_vel.X(), 0,current_body_linear_vel.Z());
      // this->body->SetLinearVel(linear_vel_correction);	//sets planar motion


      
    }


    public: void OnUpdate(){// Called by the world update START event (@ each sim update --> control sample)

     
      
    }
    
    public: void OnUpdateEnd(){// Called by the world update end event (@ each sim update --> control sample)
      // Put here code to be run after the integration update
      
      if (((this->sim_iteration_count)%this->sim_steps_over_control_step)==0){//Corresponds to the end of a control instant-->assign state after integration
        this->Q_sol_assigner();// assignes current solution to Q_sol (after integration)
        // This needs to be put here because otherwise the very last system state would be missed
      }
      
    }
    ///////////////////////////////////////////////////////////////////////////////////////////
    //private method 
    private:void PrintStuff(){
       
    
      // printf("q_k:\n");
      // for(int i=0;i<(this->q_k).rows();i++)  // loop 3 times for three lines
      //   {
      //     for(int j=0;j<(this->q_k).cols();j++)  // loop for the three elements on the line
      //     {
      //         cout<<this->q_k(i,j);  // display the current element out of the array
      //         printf("\t");
      //     }
      //     cout<<endl;  // when the inner loop is done, go to a new line
      //   }
      //   printf("\n");



      // printf("alpha_spawn:\t%f",this->alpha_spawn);
      // printf("\n");

  
      // printf("gamma_spawn:\t%f",this->gamma_spawn);
      // printf("\n");


      printf("p_cm_tot_rel_wheel:\n");
      for(int i=0;i<(this->p_cm_tot_rel_wheel_k).rows();i++)  // loop 3 times for three lines
        {
          for(int j=0;j<(this->p_cm_tot_rel_wheel_k).cols();j++)  // loop for the three elements on the line
          {
              cout<<this->p_cm_tot_rel_wheel_k(i,j);  // display the current element out of the array
              printf("\t");
          }
          cout<<endl;  // when the inner loop is done, go to a new line
        }
        printf("\n");

      printf("v_cm_tot_rel_wheel:\n");
      for(int i=0;i<(this->v_cm_tot_rel_wheel_k).rows();i++)  // loop 3 times for three lines
        {
          for(int j=0;j<(this->v_cm_tot_rel_wheel_k).cols();j++)  // loop for the three elements on the line
          {
              cout<<this->v_cm_tot_rel_wheel_k(i,j);  // display the current element out of the array
              printf("\t");
          }
          cout<<endl;  // when the inner loop is done, go to a new line
        }
        printf("\n");

      printf(" Chi_k:\n");
      for(int i=0;i<(this-> Chi_k).rows();i++)  // loop 3 times for three lines
        {
          for(int j=0;j<(this-> Chi_k).cols();j++)  // loop for the three elements on the line
          {
              cout<<this-> Chi_k(i,j);  // display the current element out of the array
              printf("\t");
          }
          cout<<endl;  // when the inner loop is done, go to a new line
        }
        printf("\n");

      printf(" Chi_dot_k:\n");
      for(int i=0;i<(this-> Chi_dot_k).rows();i++)  // loop 3 times for three lines
        {
          for(int j=0;j<(this-> Chi_dot_k).cols();j++)  // loop for the three elements on the line
          {
              cout<<this-> Chi_dot_k(i,j);  // display the current element out of the array
              printf("\t");
          }
          cout<<endl;  // when the inner loop is done, go to a new line
        }
        printf("\n");

      printf(" Chi_ddot_command:\n");
      for(int i=0;i<(this-> Chi_ddot_command_k).rows();i++)  // loop 3 times for three lines
        {
          for(int j=0;j<(this-> Chi_ddot_command_k).cols();j++)  // loop for the three elements on the line
          {
              cout<<this-> Chi_ddot_command_k(i,j);  // display the current element out of the array
              printf("\t");
          }
          cout<<endl;  // when the inner loop is done, go to a new line
        }
        printf("\n");

    

      // printf("robot_inertial_par:\n");
      // for(int i=0;i<(this->robot_inertial_par).rows();i++)  // loop 3 times for three lines
      //   {
      //     for(int j=0;j<(this->robot_inertial_par).cols();j++)  // loop for the three elements on the line
      //     {
      //         cout<<this->robot_inertial_par(i,j);  // display the current element out of the array
      //         printf("\t");
      //     }
      //     cout<<endl;  // when the inner loop is done, go to a new line
      //   }
      //   printf("\n");

      // printf("robot_el_par:\n");
      // for(int i=0;i<(this->robot_knee_el_par).rows();i++)  // loop 3 times for three lines
      //   {
      //     for(int j=0;j<(this->robot_knee_el_par).cols();j++)  // loop for the three elements on the line
      //     {
      //         cout<<this->robot_knee_el_par(i,j);  // display the current element out of the array
      //         printf("\t");
      //     }
      //     cout<<endl;  // when the inner loop is done, go to a new line
      //   }
      //   printf("\n");

      //   printf("robot dim:\n");
      // for(int i=0;i<(this->robot_dim).rows();i++)  // loop 3 times for three lines
      //   {
      //     for(int j=0;j<(this->robot_dim).cols();j++)  // loop for the three elements on the line
      //     {
      //         cout<<this->robot_dim(i,j);  // display the current element out of the array
      //         printf("\t");
      //     }
      //     cout<<endl;  // when the inner loop is done, go to a new line
      //   }
      //   printf("\n");


      // printf("A_nonholonomic_k:\n");
      // for(int i=0;i<(this->A_nonholonomic_k).rows();i++)  // loop 3 times for three lines
      //   {
      //     for(int j=0;j<(this->A_nonholonomic_k).cols();j++)  // loop for the three elements on the line
      //     {
      //         cout<<this->A_nonholonomic_k(i,j);  // display the current element out of the array
      //         printf("\t");
      //     }
      //     cout<<endl;  // when the inner loop is done, go to a new line
      //   }
      //   printf("\n");

      // printf("A_lambda_dyn_k:\n");
      // for(int i=0;i<(this->A_lambda_dyn_k).rows();i++)  // loop 3 times for three lines
      //   {
      //     for(int j=0;j<(this->A_lambda_dyn_k).cols();j++)  // loop for the three elements on the line
      //     {
      //         cout<<this->A_lambda_dyn_k(i,j);  // display the current element out of the array
      //         printf("\t");
      //     }
      //     cout<<endl;  // when the inner loop is done, go to a new line
      //   }
      //   printf("\n");

      //   printf("A_feas_k:\n");
      //   for(int i=0;i<(this->A_feas_k).rows();i++)  // loop 3 times for three lines
      //   {
      //     for(int j=0;j<(this->A_feas_k).cols();j++)  // loop for the three elements on the line
      //     {
      //         cout<<this->A_feas_k(i,j);  // display the current element out of the array
      //         printf("\t");
      //     }
      //     cout<<endl;  // when the inner loop is done, go to a new line
      //   }
      //   printf("\n");

      //   printf("A_lambda_con_k:\n");
      // for(int i=0;i<(this->A_lambda_con_k).rows();i++)  // loop 3 times for three lines
      //   {
      //     for(int j=0;j<(this->A_lambda_con_k).cols();j++)  // loop for the three elements on the line
      //     {
      //         cout<<this->A_lambda_con_k(i,j);  // display the current element out of the array
      //         printf("\t");
      //     }
      //     cout<<endl;  // when the inner loop is done, go to a new line
      //   }
      //   printf("\n");
        
      //   printf("A_u_check_k:\n");

      // for(int i=0;i<(this->A_u_check_k).rows();i++)  // loop 3 times for three lines
      //   {
      //     for(int j=0;j<(this->A_u_check_k).cols();j++)  // loop for the three elements on the line
      //     {
      //         cout<<this->A_u_check_k(i,j);  // display the current element out of the array
      //         printf("\t");
      //     }
      //     cout<<endl;  // when the inner loop is done, go to a new line
      //   }
      //   printf("\n");
      //   printf("B_check_inv_k:\n");

      //   for(int i=0;i<(this->B_check_inv_k).rows();i++)  // loop 3 times for three lines
      //   {
      //     for(int j=0;j<(this->B_check_inv_k).cols();j++)  // loop for the three elements on the line
      //     {
      //         cout<<this->B_check_inv_k(i,j);  // display the current element out of the array
      //         printf("\t");
      //     }
      //     cout<<endl;  // when the inner loop is done, go to a new line
      //   }
      //   printf("\n");
      //   printf("A_u_lim_k:\n");
      // for(int i=0;i<(this->A_u_lim_k).rows();i++)  // loop 3 times for three lines
      //   {
      //     for(int j=0;j<(this->A_u_lim_k).cols();j++)  // loop for the three elements on the line
      //     {
      //         cout<<this->A_u_lim_k(i,j);  // display the current element out of the array
      //         printf("\t");
      //     }
      //     cout<<endl;  // when the inner loop is done, go to a new line
      //   }
      //   printf("\n");
      //   printf("A_state_lim_k:\n");
      //   for(int i=0;i<(this->A_state_lim_k).rows();i++)  // loop 3 times for three lines
      //   {
      //     for(int j=0;j<(this->A_state_lim_k).cols();j++)  // loop for the three elements on the line
      //     {
      //         cout<<this->A_state_lim_k(i,j);  // display the current element out of the array
      //         printf("\t");
      //     }
      //     cout<<endl;  // when the inner loop is done, go to a new line
      //   }
      //   printf("\n");
      //   printf("b_nonholonomic_k:\n");
      //   for(int i=0;i<(this->b_nonholonomic_k).rows();i++)  // loop 3 times for three lines
      //   {
      //     for(int j=0;j<(this->b_nonholonomic_k).cols();j++)  // loop for the three elements on the line
      //     {
      //         cout<<this->b_nonholonomic_k(i,j);  // display the current element out of the array
      //         printf("\t");
      //     }
      //     cout<<endl;  // when the inner loop is done, go to a new line
      //   }
      //   printf("\n");
      //   printf("b_lambda_dyn_k:\n");
      //   for(int i=0;i<(this->b_lambda_dyn_k).rows();i++)  // loop 3 times for three lines
      //   {
      //     for(int j=0;j<(this->b_lambda_dyn_k).cols();j++)  // loop for the three elements on the line
      //     {
      //         cout<<this->b_lambda_dyn_k(i,j);  // display the current element out of the array
      //         printf("\t");
      //     }
      //     cout<<endl;  // when the inner loop is done, go to a new line
      //   }
      //   printf("\n");
      //   printf("b_feas_k:\n");
      //   for(int i=0;i<(this->b_feas_k).rows();i++)  // loop 3 times for three lines
      //   {
      //     for(int j=0;j<(this->b_feas_k).cols();j++)  // loop for the three elements on the line
      //     {
      //         cout<<this->b_feas_k(i,j);  // display the current element out of the array
      //         printf("\t");
      //     }
      //     cout<<endl;  // when the inner loop is done, go to a new line
      //   }
      //   printf("\n");
      //   printf("b_lambda_con_k:\n");
      //   for(int i=0;i<(this->b_lambda_con_k).rows();i++)  // loop 3 times for three lines
      //   {
      //     for(int j=0;j<(this->b_lambda_con_k).cols();j++)  // loop for the three elements on the line
      //     {
      //         cout<<this->b_lambda_con_k(i,j);  // display the current element out of the array
      //         printf("\t");
      //     }
      //     cout<<endl;  // when the inner loop is done, go to a new line
      //   }
      //   printf("\n");
      //   printf("r_u_check_k:\n");
      //   for(int i=0;i<(this->r_u_check_k).rows();i++)  // loop 3 times for three lines
      //   {
      //     for(int j=0;j<(this->r_u_check_k).cols();j++)  // loop for the three elements on the line
      //     {
      //         cout<<this->r_u_check_k(i,j);  // display the current element out of the array
      //         printf("\t");
      //     }
      //     cout<<endl;  // when the inner loop is done, go to a new line
      //   }
      //   printf("\n");
      //   printf("b_u_lim_k:\n");
      //   for(int i=0;i<(this->b_u_lim_k).rows();i++)  // loop 3 times for three lines
      //   {
      //     for(int j=0;j<(this->b_u_lim_k).cols();j++)  // loop for the three elements on the line
      //     {
      //         cout<<this->b_u_lim_k(i,j);  // display the current element out of the array
      //         printf("\t");
      //     }
      //     cout<<endl;  // when the inner loop is done, go to a new line
      //   }
      //   printf("\n");
      //   printf("b_state_lim_k:\n");
      //   for(int i=0;i<(this->b_state_lim_k).rows();i++)  // loop 3 times for three lines
      //   {
      //     for(int j=0;j<(this->b_state_lim_k).cols();j++)  // loop for the three elements on the line
      //     {
      //         cout<<this->b_state_lim_k(i,j);  // display the current element out of the array
      //         printf("\t");
      //     }
      //     cout<<endl;  // when the inner loop is done, go to a new line
      //   }
      //   printf("\n");
      //   printf("A_eq_k:\n");
      //   for(int i=0;i<(this->A_eq_k).rows();i++)  // loop 3 times for three lines
      //   {
      //     for(int j=0;j<(this->A_eq_k).cols();j++)  // loop for the three elements on the line
      //     {
      //         cout<<this->A_eq_k(i,j);  // display the current element out of the array
      //         printf("\t");
      //     }
      //     cout<<endl;  // when the inner loop is done, go to a new line
      //   }
      //   printf("\n");
      //   printf("A_ineq_k:\n");
      //   for(int i=0;i<(this->A_ineq_k).rows();i++)  // loop 3 times for three lines
      //   {
      //     for(int j=0;j<(this->A_ineq_k).cols();j++)  // loop for the three elements on the line
      //     {
      //         cout<<this->A_ineq_k(i,j);  // display the current element out of the array
      //         printf("\t");
      //     }
      //     cout<<endl;  // when the inner loop is done, go to a new line
      //   }
      //   printf("\n");
      //   printf("b_eq_k:\n");
      //   for(int i=0;i<(this->b_eq_k).rows();i++)  // loop 3 times for three lines
      //   {
      //     for(int j=0;j<(this->b_eq_k).cols();j++)  // loop for the three elements on the line
      //     {
      //         cout<<this->b_eq_k(i,j);  // display the current element out of the array
      //         printf("\t");
      //     }
      //     cout<<endl;  // when the inner loop is done, go to a new line
      //   }
      //   printf("\n");
      //   printf("b_ineq_k:\n");
      //   for(int i=0;i<(this->b_ineq_k).rows();i++)  // loop 3 times for three lines
      //   {
      //     for(int j=0;j<(this->b_ineq_k).cols();j++)  // loop for the three elements on the line
      //     {
      //         cout<<this->b_ineq_k(i,j);  // display the current element out of the array
      //         printf("\t");
      //     }
      //     cout<<endl;  // when the inner loop is done, go to a new line
      //   }
      //   printf("\n");

        printf("x_k:\n");
        for(int i=0;i<(this->x_k).rows();i++)  // loop 3 times for three lines
        {
          for(int j=0;j<(this->x_k).cols();j++)  // loop for the three elements on the line
          {
              cout<<this->x_k(i,j);  // display the current element out of the array
              printf("\t");
          }
          cout<<endl;  // when the inner loop is done, go to a new line
        }
        printf("\n");

        printf("u_left_k:\n");
        for(int i=0;i<(this->u_left_k).rows();i++)  // loop 3 times for three lines
        {
          for(int j=0;j<(this->u_left_k).cols();j++)  // loop for the three elements on the line
          {
              cout<<this->u_left_k(i,j);  // display the current element out of the array
              printf("\t");
          }
          cout<<endl;  // when the inner loop is done, go to a new line
        }
        printf("\n");

        
        printf("%d\n\n\n",this->control_iteration_count);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //private method 

    private:void WriteToTxT(){

      ofstream myfile;

      myfile.open ("~/IIT_CRIS_HHCM/Gazebo/sim_results/task_tracker/Time.txt");
      myfile << this->Time; 
      myfile.close();
      myfile.open ("~/IIT_CRIS_HHCM/Gazebo/sim_results/task_tracker/Q_sol.txt");
      myfile << this->Q_sol; 
      myfile.close();
      myfile.open ("~/IIT_CRIS_HHCM/Gazebo/sim_results/task_tracker/U_left.txt");
      myfile << this->U_left; 
      myfile.close();
      myfile.open ("~/IIT_CRIS_HHCM/Gazebo/sim_results/task_tracker/U_right.txt");
      myfile << this->U_right; 
      myfile.close();
      myfile.open ("~/IIT_CRIS_HHCM/Gazebo/sim_results/task_tracker/Chi_ref.txt");
      myfile << this->Chi_ref; 
      myfile.close();
      myfile.open ("~/IIT_CRIS_HHCM/Gazebo/sim_results/task_tracker/Chi_dot_ref.txt");
      myfile <<this->Chi_dot_ref; 
      myfile.close();
      myfile.open ("~/IIT_CRIS_HHCM/Gazebo/sim_results/task_tracker/Chi_ddot_ref.txt");
      myfile << this->Chi_ddot_ref; 
      myfile.close();
      myfile.open ("~/IIT_CRIS_HHCM/Gazebo/sim_results/task_tracker/Chi.txt");
      myfile << this->Chi; 
      myfile.close();
      myfile.open ("~/IIT_CRIS_HHCM/Gazebo/sim_results/task_tracker/Chi_dot.txt");
      myfile <<this->Chi_dot; 
      myfile.close();
      myfile.open ("~/IIT_CRIS_HHCM/Gazebo/sim_results/task_tracker/dt_control.txt");
      myfile << this->dt_control_calc; 
      myfile.close();
      myfile.open ("~/IIT_CRIS_HHCM/Gazebo/sim_results/task_tracker/cost.txt");
      myfile << this->QP_cost_value; 
      myfile.close();
      myfile.open ("~/IIT_CRIS_HHCM/Gazebo/sim_results/task_tracker/U_middle.txt");
      myfile << this->U_middle; 
      myfile.close();
      myfile.open ("~/IIT_CRIS_HHCM/Gazebo/sim_results/task_tracker/X_k.txt");
      myfile << this->X_k; 
      myfile.close();

    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //private method 
    private: void PointerGetter(){
      // Pointers to joints and links
      //Joints:
      this->left_wheel_hinge=this->model->GetJoint("left_wheel_hinge");
      this->left_ah_d_hinge=this->model->GetJoint("left_ah_d_hinge");
      this->left_ah_c_hinge=this->model->GetJoint("left_ah_c_hinge");
      this->left_c_body_hinge=this->model->GetJoint("left_c_body_hinge");
      this->left_d_body_hinge=this->model->GetJoint("left_d_body_hinge");

      this->right_wheel_hinge=this->model->GetJoint("right_wheel_hinge");
      this->right_ah_d_hinge=this->model->GetJoint("right_ah_d_hinge");
      this->right_ah_c_hinge=this->model->GetJoint("right_ah_c_hinge");
      this->right_c_body_hinge=this->model->GetJoint("right_c_body_hinge");
      this->right_d_body_hinge=this->model->GetJoint("right_d_body_hinge");

      //Links:
      this->left_wheel=this->model->GetLink("left_wheel");
      this->left_ah_link=this->model->GetLink("left_ah_link");
      this->left_d_link=this->model->GetLink("left_d_link");
      this->left_c_link=this->model->GetLink("left_c_link");

      this->right_wheel=this->model->GetLink("right_wheel");
      this->right_ah_link=this->model->GetLink("right_ah_link");
      this->right_d_link=this->model->GetLink("right_d_link");
      this->right_c_link=this->model->GetLink("right_c_link");
      this->body=this->model->GetLink("body");
    }
  
    ///////////////////////////////////////////////////////////////////////////////////////////
    //private method for the actual controller-->to be called inside Load and Onupdate (similar to Matlab)
    private: void TaskController_input_calc(){

      VectorXf Delta_left=VectorXf::Zero(2);
      VectorXf Delta_dot_left=VectorXf::Zero(2);
      VectorXf Delta_right=VectorXf::Zero(2);
      VectorXf Delta_dot_right=VectorXf::Zero(2);

      this->u_middle_k=this->B_check_inv_k*(this->A_u_check_k*this->x_k-this->r_u_check_k);//input computed with the "equivalent" central knee
      
      Delta_left<<this->x_w_left_k-this->q_k(0,0),this->beta_left_k-this->q_k(4,0);
      Delta_dot_left<<this->x_w_left_dot_k-this->q_k(5,0),this->beta_left_dot_k-this->q_k(9,0);
      Delta_right<<this->x_w_right_k-this->q_k(0,0),this->beta_right_k-this->q_k(4,0);
      Delta_dot_right<<this->x_w_right_dot_k-this->q_k(5,0),this->beta_right_dot_k-this->q_k(9,0);

      this->u_left_k=this->u_middle_k-this->planar_compensation_Kp*Delta_left-this->planar_compensation_Kd*Delta_dot_left;
      this->u_right_k=this->u_middle_k-this->planar_compensation_Kp*Delta_right-this->planar_compensation_Kd*Delta_dot_right;// compensation for maintaining the motion planar
      
      

      // this->u_left_k<<this->u_left_k(0,0),0; // uncomment to apply only wheel torques
      // this->u_right_k<<this->u_right_k(0,0),0;

      // this->u_left_k<<0,0; // uncomment to apply no input
      // this->u_right_k<<0,0;

      this->U_middle.row(this->control_iteration_count-1)=this->u_middle_k.transpose();//assigning inputs
    
      if ( this->cost_fun_value_k>this->cost_threshold || isnan(this->u_middle_k(0,0)) || isnan(this->u_middle_k(1,0)) || abs(this->q_k(3,0))>=M_PI/180.0*70.0){ 
        //if the cost threshold is exceeded or if the QP cannot compute a solution or the robot fell, do not apply any input
        this->u_left_k<<0,0;
        this->u_right_k<<0,0;
      }

      // MatrixXf K_LQR= MatrixXf::Zero(2,6); //crude LQR test
      // VectorXf q_hat_tilde= VectorXf::Zero(6,1);
      // VectorXf u_ref= VectorXf::Zero(2,1);
      // u_ref<<0,2.325270695412623;
      // K_LQR <<-0.1363,   16.4941,    4.2691,   -0.2126,    3.2229,     0.3421,
      //         -0.5611,   25.4323,   20.4953,   -0.6929,   5.6745,     4.3798;
      // q_hat_tilde<<this->q_k(2,0),this->q_k(3,0),this->q_k(4,0)-75.0*M_PI/180.0,this->q_k(7,0),this->q_k(8,0),this->q_k(9,0) ;
      // this->u_left_k=-K_LQR*q_hat_tilde+u_ref;
      // this->u_right_k=-K_LQR*q_hat_tilde+u_ref;
      
      // this->u_left_k<<0,0;//zeros th inputs
      // this->u_right_k<<0,0;

      this->U_left.row(this->control_iteration_count-1)=(this->u_left_k).transpose();
      this->U_right.row(this->control_iteration_count-1)=(this->u_right_k).transpose();

      this->left_wheel_hinge->SetForce(0,  this->u_left_k(0,0));
      this->left_ah_d_hinge->SetForce(0,   this->u_left_k(1,0));
      this->right_wheel_hinge->SetForce(0, this->u_right_k(0,0)); 
      this->right_ah_d_hinge->SetForce(0,  this->u_right_k(1,0));
  
  
    }
    ///////////////////////////////////////////////////////////////////////////////////////////
    //PRIVATE method for reading and assigning to the controller all necessary model properties (dimensions, physical properties, inertial properties) 

    private: void ModelPropGetter(){
      
      double alpha_aux;
      double b_gamma_aux;
      double h1_aux;
      double h2_aux;

      if (this->model_properties_are_manually_set==1){ 
        
        double l_a=0.1000;
        double l_b=0.1450;
        double l_c=0.2360;
        double l_d=0.2560;
        double l_h=0.3550;
        double theta_st=M_PI/4.0;
        double theta_st_CoM_body=4.0/5.0*M_PI;

        double r_w=0.06;
        double l_cm_w=r_w;
        double l_cm_d=l_d/2.0;
        double l_cm_ah=(l_a+l_h)/2.0;
        double l_cm_c=l_c/2.0;
        double l_cm_body=0.0178;

        this->robot_dim<<l_a,l_b,l_c,l_d,l_h,theta_st,theta_st_CoM_body;

        this->alpha_gamma_calc(this->beta_spawn, &alpha_aux, &b_gamma_aux);
        this->alpha_spawn=alpha_aux;
        this->gamma_spawn=b_gamma_aux;

        double rho_alluminium=2700;
        double profilate_fullness=0.6;
        double profilate_area=profilate_fullness*0.02*0.02;

        double m_d=rho_alluminium*profilate_area*l_d;
        double m_ah=rho_alluminium*profilate_area*(l_a+l_h);
        double m_c=rho_alluminium*profilate_area*l_c;
        double m_body=10.0;
        double m_w=0.3;

        double J_c=m_c*l_c*l_c/12.0;
        double J_d=m_d*l_d*l_d/12.0;
        double J_ah=m_ah*(l_a+l_h)*(l_a+l_h)/12.0;
        double J_body=0.1;
        double J_w=0.001;

        double k1t=30;
        double k2t=0;
        double k3t=0;
        double k4t=0;
        this->alpha_gamma_calc(this->beta_ref_el, &alpha_aux, &b_gamma_aux);
        double alpha0_el=alpha_aux;
        double gamma0_el=b_gamma_aux;
        
        this->robot_inertial_par<<l_cm_w,l_cm_d,l_cm_ah,l_cm_c,l_cm_body,m_w,m_ah,m_d,m_c,m_body,J_w,J_ah,J_d,J_c,J_body;
        this->robot_knee_el_par<<k1t,k2t,k3t,k4t,this->beta_ref_el,alpha0_el,gamma0_el;

      } else{
        // Using only left links, being the system symmetric

        physics::InertialPtr inertial_wheel_link=this->left_wheel->GetInertial();
        physics::InertialPtr inertial_ah_link=this->left_ah_link->GetInertial();
        physics::InertialPtr inertial_d_link=this->left_d_link->GetInertial();
        physics::InertialPtr inertial_c_link=this->left_c_link->GetInertial();
        physics::InertialPtr inertial_body_link=this->body->GetInertial();

        double l_a=0.1000;
        double l_b=0.1450;
        double l_c=0.2360;
        double l_d=0.2560;
        double l_h=0.3550;
        double theta_st=M_PI/4.0;
        double theta_st_CoM_body=4.0/5.0*M_PI;

        this->robot_dim<<l_a,l_b,l_c,l_d,l_h,theta_st,theta_st_CoM_body;

        this->alpha_gamma_calc(this->beta_spawn, &alpha_aux, &b_gamma_aux);
        this->alpha_spawn=alpha_aux;
        this->gamma_spawn=b_gamma_aux;
        
        double r_w=0.06;
        double l_cm_w=r_w;
        double l_cm_d=l_d/2.0;
        double l_cm_ah=(l_a+l_h)/2.0;
        double l_cm_c=l_c/2.0;
        double l_cm_body=0.0178;

        double m_d=inertial_d_link->Mass();
        double m_ah=inertial_ah_link->Mass();
        double m_c=inertial_c_link->Mass();
        double m_body=inertial_body_link->Mass();
        double m_w=inertial_wheel_link->Mass();

        double J_d=inertial_d_link->IYY();
        double J_ah=inertial_ah_link->IYY();
        double J_c=inertial_c_link->IYY();
        double J_body=inertial_body_link->IYY();
        double J_w=inertial_wheel_link->IYY();
       
        double k1t=this->left_ah_c_hinge->GetStiffness(0);
        double k2t=this->left_ah_d_hinge->GetStiffness(0);
        double k3t=this->left_c_body_hinge->GetStiffness(0);
        double k4t=this->left_d_body_hinge->GetStiffness(0);

        this->alpha_gamma_calc(this->beta_ref_el, &alpha_aux, &b_gamma_aux);
        double alpha0_el=alpha_aux;
        double gamma0_el=b_gamma_aux;

        this->robot_inertial_par<<l_cm_w,l_cm_d,l_cm_ah,l_cm_c,l_cm_body,m_w,m_ah,m_d,m_c,m_body,J_w,J_ah,J_d,J_c,J_body;
        this->robot_knee_el_par<<k1t,k2t,k3t,k4t,this->beta_ref_el,alpha0_el,gamma0_el;
      }
    }
    ///////////////////////////////////////////////////////////////////////////////////////////
    //PRIVATE method for computing all necessary QP matrices

    private: void QPMatrixCalc(){
      
      double alpha;
      double b_gamma;
      double h1;
      double h2;
      this->knee_jacobians_h1_h2_calc(this->q_k(4,0), &h1, &h2);
      this->alpha_gamma_calc(this->q_k(4,0), &alpha, &b_gamma);

      this->H_rel_CoM_calc(alpha, b_gamma, h1, h2);
      this->f_T_rel_CoM_calc(alpha, b_gamma, h1, h2);

      this->A_nonholonomic_calc();
      this->A_lambda_dyn_calc(alpha, b_gamma, h1, h2);
      this->A_feas_calc(alpha, b_gamma, h1, h2);
      this->A_lambda_con_calc();
      this->A_u_check_k_calc(alpha, b_gamma, h1, h2);
      this->B_check_inv_calc(h1);
      this->A_u_lim_calc(alpha, b_gamma, h1, h2);
      //this->A_state_lim_calc();//A_state_lim_k not assigned since does not work (left initialized to zeros)

      this->b_nonholonomic_calc();
      this->b_lambda_dyn_calc(alpha, b_gamma, h1, h2);
      this->b_feas_calc(alpha, b_gamma, h1, h2);
      this->b_lambda_con_calc();
      this->r_u_check_calc(alpha, b_gamma, h1, h2);
      this->b_u_lim_calc(alpha, b_gamma, h1, h2);
      //this->b_state_lim_calc();//b_state_lim_k not assigned since does not work (left initialized to zeros)

      this->A_eq_k<<this->A_nonholonomic_k,this->A_lambda_dyn_k,this->A_feas_k;
      this->b_eq_k<<this->b_nonholonomic_k,this->b_lambda_dyn_k,this->b_feas_k;
      this->A_ineq_k<<this->A_lambda_con_k,this->A_u_lim_k,this->A_state_lim_k; 
      this->b_ineq_k<<this->b_lambda_con_k,this->b_u_lim_k,this->b_state_lim_k; 
    }
    ///////////////////////////////////////////////////////////////////////////////////////////
    //PRIVATE method 

    private: void initial_q_setter(){
        
        ignition::math::Pose3d left_wheel_pose;
        left_wheel_pose.Set(0, 0, this->robot_inertial_par(0,0), 0, -phi_r0, 0);
        this->model->SetLinkWorldPose(left_wheel_pose,this->left_wheel); 

        // ignition::math::Pose3d right_wheel_pose;
        // right_wheel_pose.Set(0, 0, this->robot_inertial_par(0,0), 0, -phi_r0, 0);
        // this->model->SetLinkWorldPose(right_wheel_pose,this->right_wheel); 

        //  double beta0=M_PI/180.0*90;
        //  this->model->SetJointPosition(this->left_ah_d_hinge->GetName(),this->beta_spawn-beta0);
        //  this->model->SetJointPosition(this->right_ah_d_hinge->GetName(),this->beta_spawn-beta0);

         //this->left_d_body_hinge->SetPosition(0,this->beta_spawn-beta0,false);
         //this->left_d_body_hinge->SetVelocity(0,1000);
         //this->right_d_body_hinge->SetPosition(0,this->beta_spawn-beta0,false);
         //this->right_d_body_hinge->SetVelocity(0,1000);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////
    //PRIVATE method for assigning the current state solution

    private: void Q_sol_assigner(){
      // Calling state constructors for necessary joints and links; necessary to access state at the current sim iteration
      // left side joints' states
      this->left_wheel_hinge_state=physics::JointState(this->left_wheel_hinge);
      this->left_ah_d_hinge_state=physics::JointState(this->left_ah_d_hinge);
      this->left_ah_c_hinge_state=physics::JointState(this->left_ah_c_hinge);
      this->left_c_body_hinge_state=physics::JointState(this->left_c_body_hinge);
      this->left_d_body_hinge_state=physics::JointState(this->left_d_body_hinge);
      // right side joints' states
      this->right_wheel_hinge_state=physics::JointState(this->right_wheel_hinge); 
      this->right_ah_d_hinge_state=physics::JointState(this->right_ah_d_hinge);
      this->right_ah_c_hinge_state=physics::JointState(this->right_ah_c_hinge);
      this->right_c_body_hinge_state=physics::JointState(this->right_c_body_hinge);
      this->right_d_body_hinge_state=physics::JointState(this->right_d_body_hinge);     
      // Link states
      this->left_wheel_state=physics::LinkState(this->left_wheel);
      this->left_ah_link_state=physics::LinkState(this->left_ah_link);
      this->left_d_link_state=physics::LinkState(this->left_d_link);
      this->left_c_link_state=physics::LinkState(this->left_c_link);
      this->right_wheel_state=physics::LinkState(this->right_wheel);
      this->right_ah_link_state=physics::LinkState(this->right_ah_link);
      this->right_d_link_state=physics::LinkState(this->right_d_link);
      this->right_c_link_state=physics::LinkState(this->right_c_link);
      this->body_state=physics::LinkState(this->body);
      // Getting useful state/velocity info
      ignition::math::Pose3d left_wheel_pose=this->left_wheel_state.Pose();
      ignition::math::Pose3d left_wheel_vel=this->left_wheel_state.Velocity();
      ignition::math::Pose3d right_wheel_pose=this->right_wheel_state.Pose();
      ignition::math::Pose3d right_wheel_vel=this->right_wheel_state.Velocity();
      ignition::math::Pose3d body_pose=this->body_state.Pose();
      ignition::math::Pose3d body_vel=this->body_state.Velocity();

    
      this->x_w_left_k=left_wheel_pose.Pos()[0];
      this->x_w_left_dot_k=left_wheel_vel.Pos()[0];
      this->x_w_right_k=right_wheel_pose.Pos()[0];
      this->x_w_right_dot_k=right_wheel_vel.Pos()[0];

      this->z_w_left_k=left_wheel_pose.Pos()[2]-this->robot_inertial_par(0,0);
      this->z_w_left_dot_k=left_wheel_vel.Pos()[2];
      this->z_w_right_k=right_wheel_pose.Pos()[2]-this->robot_inertial_par(0,0);
      this->z_w_right_dot_k=right_wheel_vel.Pos()[2];

      this->phi_r_k=-(body_pose.Rot()).Pitch();
      this->phi_r_dot_k=-(body_vel.Rot()).Pitch();

      double alpha_aux;
      double b_gamma_aux;
      double h1_aux;
      double h2_aux;

      this->beta_left_k=this->beta_spawn-this->left_d_body_hinge->Position(0);
      this->knee_jacobians_h1_h2_calc(this->beta_left_k, &h1_aux, &h2_aux);
      this->alpha_gamma_calc(this->beta_left_k, &alpha_aux, &b_gamma_aux);
      this->alpha_left_k=alpha_aux;
      this->h1_left_k= h1_aux;
      this->beta_left_dot_k=-this->left_d_body_hinge->GetVelocity(0);

      this->beta_right_k=this->beta_spawn-this->right_d_body_hinge->Position(0);
      this->knee_jacobians_h1_h2_calc(this->beta_right_k, &h1_aux, &h2_aux);
      this->alpha_gamma_calc(this->beta_right_k, &alpha_aux, &b_gamma_aux);
      this->alpha_right_k=alpha_aux;
      this->h1_right_k= h1_aux;
      this->beta_right_dot_k=-this->right_d_body_hinge->GetVelocity(0);

      this->phi_w_left_k=this->left_wheel_hinge->Position(0)+this->alpha_left_k-this->alpha_spawn-this->phi_r_k+this->phi_r_spawn;
      this->phi_w_left_dot_k=this->left_wheel_hinge->GetVelocity(0)+this->h1_left_k*this->beta_left_dot_k-this->phi_r_dot_k;
      this->phi_w_right_k=this->right_wheel_hinge->Position(0)+this->alpha_right_k-this->alpha_spawn-this->phi_r_k+this->phi_r_spawn;
      this->phi_w_right_dot_k=this->right_wheel_hinge->GetVelocity(0)+this->h1_right_k*this->beta_right_dot_k-this->phi_r_dot_k;
      
      double beta_k=(this->beta_left_k+this->beta_right_k)/2.0;

      double beta_dot_k=(this->beta_left_dot_k+this->beta_right_dot_k)/2.0;
      this->knee_jacobians_h1_h2_calc(beta_k, &h1_aux, &h2_aux);
      this->alpha_gamma_calc(beta_k, &alpha_aux, &b_gamma_aux);
      double phi_w_k=(this->left_wheel_hinge->Position(0)+this->right_wheel_hinge->Position(0))/2.0+
                     alpha_aux-this->alpha_spawn-this->phi_r_k+this->phi_r_spawn;
      double phi_w_dot_k=(this->left_wheel_hinge->GetVelocity(0)+this->right_wheel_hinge->GetVelocity(0))/2.0+
                          h1_aux*beta_dot_k-this->phi_r_dot_k;

      if (this->use_equiv_central_knee==1){
        this->Q_sol.row(this->control_iteration_count)<<(this->x_w_left_k+this->x_w_right_k)/2.0,
                                                      (this->z_w_left_k+this->z_w_right_k)/2.0,
                                                      phi_w_k,
                                                      this->phi_r_k,
                                                      beta_k,
                                                      (this->x_w_left_dot_k+this->x_w_right_dot_k)/2.0,
                                                      (this->z_w_left_dot_k+this->z_w_right_dot_k)/2.0,
                                                      phi_w_dot_k,
                                                      this->phi_r_dot_k,
                                                      beta_dot_k;
      }
      else{//use one side as a reference for the state
        this->Q_sol.row(this->control_iteration_count)<<this->x_w_left_k,
                                                        this->z_w_left_k,
                                                        this->phi_w_left_k,
                                                        this->phi_r_k,
                                                        this->beta_left_k,
                                                        this->x_w_left_dot_k,
                                                        this->z_w_left_dot_k,
                                                        this->phi_w_left_dot_k,
                                                        this->phi_r_dot_k,
                                                        this->beta_left_dot_k;


      }

      // OBTAINING beta and phi_r derivatives numerically to correct for strange velocity errors
      if (this->control_iteration_count==0){//before first sample time
      }
      
      else{
        this->x_w_left_dot_k=(this->Q_sol((this->control_iteration_count),0)-this->Q_sol((this->control_iteration_count-1),0))/(this->delta_control);
        this->z_w_left_dot_k=(this->Q_sol((this->control_iteration_count),1)-this->Q_sol((this->control_iteration_count-1),1))/(this->delta_control);
        this->phi_w_left_dot_k=(this->Q_sol((this->control_iteration_count),2)-this->Q_sol((this->control_iteration_count-1),2))/(this->delta_control);
        this->phi_r_dot_k=(this->Q_sol((this->control_iteration_count),3)-this->Q_sol((this->control_iteration_count-1),3))/(this->delta_control);
        this->beta_left_dot_k=(this->Q_sol((this->control_iteration_count),4)-this->Q_sol((this->control_iteration_count-1),4))/(this->delta_control);
        // this->beta_right_dot_k=(this->Q_sol((this->control_iteration_count),4)-this->Q_sol((this->control_iteration_count-1),4))/(this->delta_control);
        this->Q_sol.row(this->control_iteration_count)<<this->x_w_left_k,
                                                        this->z_w_left_k,
                                                        this->phi_w_left_k,
                                                        this->phi_r_k,
                                                        this->beta_left_k,
                                                        this->x_w_left_dot_k,
                                                        this->z_w_left_dot_k,
                                                        this->phi_w_left_dot_k,
                                                        this->phi_r_dot_k,
                                                        this->beta_left_dot_k;
      

      }

    } 
      
    ///////////////////////////////////////////////////////////////////////////////////////////
    //PRIVATE method for getting the current state 

    private: void q_k_getter(){
      if (this->forward_prop_state_meas==0){//state delayed by one control sample

      double x_w_k=this->Q_sol((this->control_iteration_count-2),0);//-2 so that the assigned state is always one control sample time behind
      double z_w_k=this->Q_sol((this->control_iteration_count-2),1);
      double phi_w_k=this->Q_sol((this->control_iteration_count-2),2);
      double phi_r_k=this->Q_sol((this->control_iteration_count-2),3);
      double beta_k=this->Q_sol((this->control_iteration_count-2),4);
      double x_w_dot_k=this->Q_sol((this->control_iteration_count-2),5);
      double z_w_dot_k=this->Q_sol((this->control_iteration_count-2),6);
      double phi_w_dot_k=this->Q_sol((this->control_iteration_count-2),7);
      double phi_r_dot_k=this->Q_sol((this->control_iteration_count-2),8);
      double beta_dot_k=this->Q_sol((this->control_iteration_count-2),9);

      this->q_k<<x_w_k,    z_w_k,    phi_w_k,    phi_r_k,    beta_k,
                 x_w_dot_k,z_w_dot_k,phi_w_dot_k,phi_r_dot_k,beta_dot_k;}

      else{//state not delayed
      double x_w_k=this->Q_sol((this->control_iteration_count-1),0);
      double z_w_k=this->Q_sol((this->control_iteration_count-1),1);
      double phi_w_k=this->Q_sol((this->control_iteration_count-1),2);
      double phi_r_k=this->Q_sol((this->control_iteration_count-1),3);
      double beta_k=this->Q_sol((this->control_iteration_count-1),4);
      double x_w_dot_k=this->Q_sol((this->control_iteration_count-1),5);
      double z_w_dot_k=this->Q_sol((this->control_iteration_count-1),6);
      double phi_w_dot_k=this->Q_sol((this->control_iteration_count-1),7);
      double phi_r_dot_k=this->Q_sol((this->control_iteration_count-1),8);
      double beta_dot_k=this->Q_sol((this->control_iteration_count-1),9);

      this->q_k<<x_w_k,    z_w_k,    phi_w_k,    phi_r_k,    beta_k,
                 x_w_dot_k,z_w_dot_k,phi_w_dot_k,phi_r_dot_k,beta_dot_k;}
    
    }
    //////////////////////////////////////////////////////////////////////////////////////////
    //PRIVATE method for getting the current value of the task (and its derivative) 

    private: void Chi_traj_getter(){
      //////////////////////////////////////
      double x_freq=1;
      double y_freq=1;
      double x_amplitute=0.05;
      double y_amplitute=0.04;
      double y_offset=0.4;
      double phase_lag_x=0.0*M_PI/180.0;
      double phase_lag_y=90.0*M_PI/180.0;

      
      // // ATAN-like y trajectory
      // this->Chi_ref_k<<x_amplitute*sin(2.0*M_PI*x_freq*this->Time(this->control_iteration_count-1,0)-phase_lag_x),
      //                  y_offset+2.0*((atan(this->Time(this->control_iteration_count-1,0)-this->sim_time/2))-(atan(-this->sim_time/2)))*0.01;
      
      // this->Chi_dot_ref_k<<2.0*M_PI*x_freq*x_amplitute*cos(2.0*M_PI*x_freq*this->Time(this->control_iteration_count-1,0)-phase_lag_x),
      //                            1/(50.0*(pow(this->sim_time/2.0 - this->Time(this->control_iteration_count-1,0),2) + 1.0));

      // this->Chi_ddot_ref_k<<-4.0*pow(M_PI,2.0)*pow(x_freq,2.0)*x_amplitute*sin(2.0*M_PI*x_freq*this->Time(this->control_iteration_count-1,0)-phase_lag_x),
      //                             (50.0*this->sim_time - 100.0*this->Time(this->control_iteration_count-1,0))/pow(50.0*pow(this->sim_time/2.0 - this->Time(this->control_iteration_count-1,0),2) + 50.0,2);

      // Sinusoidal trajectories
      this->Chi_ref_k<<x_amplitute*sin(2.0*M_PI*x_freq*this->Time(this->control_iteration_count-1,0)-phase_lag_x),
                       y_offset+y_amplitute*sin(2.0*M_PI*y_freq*this->Time(this->control_iteration_count-1,0)-phase_lag_y);
      
      this->Chi_dot_ref_k<<2.0*M_PI*x_freq*x_amplitute*cos(2.0*M_PI*x_freq*this->Time(this->control_iteration_count-1,0)-phase_lag_x),
                           2.0*M_PI*y_freq*y_amplitute*cos(2.0*M_PI*y_freq*this->Time(this->control_iteration_count-1,0)-phase_lag_y);

      this->Chi_ddot_ref_k<<-4.0*pow(M_PI,2.0)*pow(x_freq,2.0)*x_amplitute*sin(2.0*M_PI*x_freq*this->Time(this->control_iteration_count-1,0)-phase_lag_x),
                            -4.0*pow(M_PI,2.0)*pow(y_freq,2.0)*y_amplitute*sin(2.0*M_PI*y_freq*this->Time(this->control_iteration_count-1,0)-phase_lag_y);

      this->Chi_ref.row(this->control_iteration_count-2)<<this->Chi_ref_k.transpose();
      this->Chi_dot_ref.row(this->control_iteration_count-2)<<this->Chi_dot_ref_k.transpose();
      this->Chi_ddot_ref.row(this->control_iteration_count-2)<<this->Chi_ddot_ref_k.transpose();

      
    }
    ///////////////////////////////////////////////////////////////////////////////////////////
    //PRIVATE method for getting the current value of the task (and its derivative) 

    private: void Chi_Chi_dot_k_getter(){
      double alpha;
      double b_gamma;
      double h1;
      double h2;
      this->knee_jacobians_h1_h2_calc(this->q_k(4,0), &h1, &h2);
      this->alpha_gamma_calc(this->q_k(4,0), &alpha, &b_gamma);
      
      this->p_cm_tot_rel_wheel_calc(alpha,b_gamma);
      this->v_cm_tot_rel_wheel_calc(alpha, b_gamma, h1, h2);

      this->Chi_k=this->p_cm_tot_rel_wheel_k;
      this->Chi_dot_k=this->v_cm_tot_rel_wheel_k;

      this->Chi.row(this->control_iteration_count-2)<<this->Chi_k.transpose();
      this->Chi_dot.row(this->control_iteration_count-2)<<this->Chi_dot_k.transpose();

    }
    //////////////////////////////////////////////////////////////////////////////////////////
    //PRIVATE method for getting the current state 

    private: void Chi_ddot_command_calc(){
      this->Chi_ddot_command_k=this->Chi_ddot_ref_k+this->Kd*(this->Chi_dot_ref_k-this->Chi_dot_k)+this->Kp*(this->Chi_ref_k-this->Chi_k);
    }
    //////////////////////////////////////////////////////////////////////////////////////////
    //PRIVATE method for computing the current value of the cost function

    private: void cost_fnctn_val_calc(){
      double alpha;
      double b_gamma;
      double h1;
      double h2;
      this->knee_jacobians_h1_h2_calc(this->q_k(4,0), &h1, &h2);
      this->alpha_gamma_calc(this->q_k(4,0), &alpha, &b_gamma);

      this->J_CoMv_tot_rel_wheel_calc(alpha, b_gamma, h1, h2);
      this->J_CoMv_dot_tot_rel_wheel_calc(alpha,b_gamma, h1,h2);

      VectorXf q_p_dot=VectorXf::Zero(5);
      VectorXf q_p_ddot=VectorXf::Zero(5);
      q_p_dot<<this->q_k(5,0),this->q_k(6,0),this->q_k(7,0),this->q_k(8,0),this->q_k(9,0);
      q_p_ddot<<this->x_k(0,0),this->x_k(1,0),this->x_k(2,0),this->x_k(3,0),this->x_k(4,0);
      MatrixXf vector_cost=this->Chi_ddot_command_k-this->J_CoMv_tot_rel_wheel_k*q_p_ddot-this->J_CoMv_dot_tot_rel_wheel_k*q_p_dot;
      this->cost_fun_value_k=sqrt(vector_cost(0,0)*vector_cost(0,0)+vector_cost(1,0)*vector_cost(1,0));
    }
    ///////////////////////////////////////////////////////////////////////////////////////////
    //PRIVATE method for computing the QP solution x_k ( [\ddot{q}_p;lambda] )

    private: void x_k_calc(){

      this->x_k=solveQP_hpipm(this->H_k, this->f_k, this->A_eq_k, this->b_eq_k, this->A_ineq_k,this->b_neg_inf_QP, this->b_ineq_k); 
      (this->X_k).row(this->control_iteration_count-1)<<(this->x_k).transpose();//assigning inputs

    }
    //////////////////////////////////////////////////////////////////////////////////////////

    //Private "low-level" methods
    private: 

      void H_rel_CoM_calc(double alpha,double b_gamma,double h1,double h2){
        
        double phi_r=this->q_k(3,0);
        double beta=this->q_k(4,0);

        double t11;
        double t15;
        double t16;
        double t19;
        double t2;
        double t20;
        double t21;
        double t22;
        double t23;
        double t24;
        double t3;
        double t35;
        double t39;
        double t4;
        double t5;
        double t50;
        double t7;
        double t8;
        double t80;
        //     This function was generated by the Symbolic Math Toolbox version 8.7.
        //     29-Jul-2021 16:00:58
        t2 = this->robot_dim(0,0) + this->robot_dim(4,0);
        t3 = this->robot_inertial_par(6,0) * 2.0;
        t4 = this->robot_inertial_par(8,0) * 2.0;
        t5 = this->robot_inertial_par(7,0) * 2.0;
        t7 = (beta + phi_r) + this->robot_dim(5,0);
        t8 = (b_gamma + phi_r) + this->robot_dim(5,0);
        t11 = std::cos(t7);
        t7 = std::sin(t7);
        t15 = phi_r + -this->robot_dim(6,0);
        t16 = (phi_r + -alpha) + this->robot_dim(5,0);
        t80 = (((this->robot_inertial_par(9,0) + t3) + t4) + t5) +
              this->robot_inertial_par(5,0) * 2.0;
        t19 = this->robot_dim(3,0) * t11;
        t20 = this->robot_inertial_par(1,0) * t11;
        t21 = this->robot_inertial_par(3,0) * std::cos(t8);
        t22 = this->robot_dim(3,0) * t7;
        t23 = this->robot_inertial_par(1,0) * t7;
        t24 = this->robot_inertial_par(3,0) * std::sin(t8);
        t11 = std::cos(t16);
        t8 = std::sin(t16);
        t50 = 1.0 / (t80 * t80);
        t35 = this->robot_dim(4,0) * t11;
        t16 = this->robot_dim(4,0) * t8;
        t39 = t2 * t11;
        t2 *= t8;
        t7 = this->robot_inertial_par(2,0) * t3;
        t3 = t7 * t11;
        t7 *= t8;
        t8 = h1 * t35;
        t11 = h1 * t16;
        t80 = ((t7 + -t5 * (t23 - t16)) + -t4 * (t21 - t2)) +
              this->robot_inertial_par(9,0) *
                  ((this->robot_inertial_par(4,0) * std::cos(t15) + -t22) + t16);
        t16 = ((h1 * t7 + this->robot_inertial_par(9,0) * (t22 + t11)) +
              t5 * (t23 + t11)) +
              t4 * (h2 * t21 + h1 * t2);
        t8 =
            ((h1 * t3 + this->robot_inertial_par(9,0) * (t19 + t8)) + t5 * (t20 + t8)) +
            this->robot_inertial_par(8,0) * (h2 * t24 + -(h1 * t39)) * -2.0;
        t7 = ((t3 + t4 * (t24 + t39)) + -t5 * (t20 + -t35)) +
            -(this->robot_inertial_par(9,0) *
              ((t19 + this->robot_inertial_par(4,0) * std::sin(t15)) + -t35));
        t11 = -(t50 * t16 * t80 * 2.0) + -(t50 * t8 * t7 * 2.0);
        this->H_k(3,3) = t50 * (t7 * t7) * 2.0 + t50 * (t80 * t80) * 2.0;
        this->H_k(4,3) = t11;
        this->H_k(3,4) = t11;
        this->H_k(4,4) = t50 * (t16 * t16) * 2.0 + t50 * (t8 * t8) * 2.0;
      }
      
      void f_T_rel_CoM_calc(double alpha,double b_gamma,double h1,double h2){
        
        double phi_r=this->q_k(3,0);
        double beta=this->q_k(4,0);
        double phi_r_dot=this->q_k(8,0);
        double beta_dot=this->q_k(9,0);

        double t11;
        double t15;
        double t16;
        double t19;
        double t2;
        double t20;
        double t21;
        double t22;
        double t23;
        double t24;
        double t3;
        double t30;
        double t39;
        double t4;
        double t40;
        double t43;
        double t46;
        double t5;
        double t54;
        double t7;
        double t8;
        double t83;
        double t84;
        //     This function was generated by the Symbolic Math Toolbox version 8.7.
        //     29-Jul-2021 16:00:59
        t2 = this->robot_dim(0,0) + this->robot_dim(4,0);
        t3 = this->robot_inertial_par(6,0) * 2.0;
        t4 = this->robot_inertial_par(8,0) * 2.0;
        t5 = this->robot_inertial_par(7,0) * 2.0;
        t7 = (beta + phi_r) + this->robot_dim(5,0);
        t8 = (b_gamma + phi_r) + this->robot_dim(5,0);
        t11 = std::cos(t7);
        t7 = std::sin(t7);
        t15 = phi_r + -this->robot_dim(6,0);
        t16 = (phi_r + -alpha) + this->robot_dim(5,0);
        t19 = this->robot_dim(3,0) * t11;
        t20 = this->robot_inertial_par(1,0) * t11;
        t21 = this->robot_inertial_par(3,0) * std::cos(t8);
        t22 = this->robot_dim(3,0) * t7;
        t23 = this->robot_inertial_par(1,0) * t7;
        t24 = this->robot_inertial_par(3,0) * std::sin(t8);
        t11 = std::cos(t16);
        t30 = std::sin(t16);
        t54 = 1.0 / ((((this->robot_inertial_par(9,0) + t3) + t4) + t5) +
                    this->robot_inertial_par(5,0) * 2.0);
        t39 = this->robot_dim(4,0) * t11;
        t40 = this->robot_dim(4,0) * t30;
        t43 = t2 * t11;
        t2 *= t30;
        t7 = this->robot_inertial_par(2,0) * t3;
        t46 = t7 * t11;
        t11 = t7 * t30;
        t16 = h1 * t39;
        t8 = h1 * t40;
        t7 = beta_dot * t54;
        t83 = t7 * (this->robot_inertial_par(9,0) * t19 + t5 * t20);
        t84 = t7 * (this->robot_inertial_par(9,0) * t22 + t5 * t23);
        t40 = ((t11 + -t5 * (t23 - t40)) + -t4 * (t21 - t2)) +
              this->robot_inertial_par(9,0) *
                  ((this->robot_inertial_par(4,0) * std::cos(t15) + -t22) + t40);
        t3 = ((h1 * t11 + this->robot_inertial_par(9,0) * (t22 + t8)) +
              t5 * (t23 + t8)) +
            t4 * (h2 * t21 + h1 * t2);
        t2 = ((h1 * t46 + this->robot_inertial_par(9,0) * (t19 + t16)) +
              t5 * (t20 + t16)) +
            this->robot_inertial_par(8,0) * (h2 * t24 + -(h1 * t43)) * -2.0;
        t30 = ((t46 + t4 * (t24 + t43)) + -t5 * (t20 + -t39)) +
              -(this->robot_inertial_par(9,0) *
                ((t19 + this->robot_inertial_par(4,0) * std::sin(t15)) + -t39));
        t7 = phi_r_dot * t54;
        t16 = beta_dot * (t84 + t7 * t3) + phi_r_dot * (t84 + -(t7 * t40));
        t7 = beta_dot * (t83 + t7 * t2) + phi_r_dot * (t83 + -(t7 * t30));
        t11 = this->Chi_ddot_command_k(0,0) * t54;
        t8 = this->Chi_ddot_command_k(1,0) * t54;
        this->f_k(3,0) =
            ((t11 * t30 * 2.0 + t8 * t40 * 2.0) + t54 * t30 * t16 * 2.0) -
            t54 * t40 * t7 * 2.0;
        this->f_k(4,0) = ((t11 * t2 * -2.0 - t8 * t3 * 2.0) + t54 * t3 * t7 * 2.0) -
                        t54 * t2 * t16 * 2.0;
}
          
      void A_nonholonomic_calc(){
         this->A_nonholonomic_k<<1.0,0,-this->robot_inertial_par(0,0),0,0,0,0,0,1.0,0,0,0,0,0;
      }
      
      void A_lambda_dyn_calc(double alpha,double b_gamma,double h1,double h2){

        double phi_r=this->q_k(3,0);
        double beta=this->q_k(4,0);

        double A_lambda_dyn_tmp;
        double b_A_lambda_dyn_tmp;
        double c_A_lambda_dyn_tmp;
        double d_A_lambda_dyn_tmp;
        double e_A_lambda_dyn_tmp;
        double f_A_lambda_dyn_tmp;
        double t10;
        double t11;
        double t12;
        double t13;
        double t14;
        double t16_tmp;
        double t18;
        double t19;
        double t2;
        double t3;
        double t4;
        double t6;
        double t7;
        //     This function was generated by the Symbolic Math Toolbox version 8.7.
        //     29-Jul-2021 16:01:00
        t2 = this->robot_inertial_par(6,0) * 2.0;
        t3 = this->robot_inertial_par(8,0) * 2.0;
        t4 = this->robot_inertial_par(7,0) * 2.0;
        t6 = (beta + phi_r) + this->robot_dim(5,0);
        t7 = (b_gamma + phi_r) + this->robot_dim(5,0);
        t10 = std::cos(t6);
        t11 = std::cos(t7);
        t12 = std::sin(t6);
        t13 = std::sin(t7);
        t14 = phi_r + -this->robot_dim(6,0);
        t6 = (phi_r + -alpha) + this->robot_dim(5,0);
        t7 = (((this->robot_inertial_par(9,0) + t2) + t3) + t4) +
            this->robot_inertial_par(5,0) * 2.0;
        t16_tmp = this->robot_dim(3,0) * this->robot_inertial_par(9,0);
        t18 = std::cos(t6);
        t19 = std::sin(t6);
        t6 = this->robot_inertial_par(1,0) * t4;
        this->A_lambda_dyn_k(0,0) = t7;
        this->A_lambda_dyn_k(1,1) = t7;
        A_lambda_dyn_tmp = this->robot_inertial_par(4,0) * this->robot_inertial_par(9,0);
        b_A_lambda_dyn_tmp = this->robot_dim(0,0) * this->robot_inertial_par(8,0);
        c_A_lambda_dyn_tmp =
            this->robot_inertial_par(2,0) * this->robot_inertial_par(6,0);
        d_A_lambda_dyn_tmp = this->robot_dim(4,0) * this->robot_inertial_par(8,0);
        e_A_lambda_dyn_tmp = this->robot_dim(4,0) * this->robot_inertial_par(9,0);
        f_A_lambda_dyn_tmp = this->robot_dim(4,0) * this->robot_inertial_par(7,0);
        t7 = t16_tmp * t10 + t6 * t10;
        this->A_lambda_dyn_k(0,3) = ((((((t7 + A_lambda_dyn_tmp * std::sin(t14)) -
                              b_A_lambda_dyn_tmp * t18 * 2.0) -
                              c_A_lambda_dyn_tmp * t18 * 2.0) -
                            this->robot_inertial_par(3,0) *
                                this->robot_inertial_par(8,0) * t13 * 2.0) -
                            d_A_lambda_dyn_tmp * t18 * 2.0) -
                          e_A_lambda_dyn_tmp * t18) -
                          f_A_lambda_dyn_tmp * t18 * 2.0;
        t6 = t16_tmp * t12 + t6 * t12;
        this->A_lambda_dyn_k(1,3) = ((((((t6 - A_lambda_dyn_tmp * std::cos(t14)) -
                              b_A_lambda_dyn_tmp * t19 * 2.0) -
                              c_A_lambda_dyn_tmp * t19 * 2.0) -
                            d_A_lambda_dyn_tmp * t19 * 2.0) -
                            e_A_lambda_dyn_tmp * t19) -
                          f_A_lambda_dyn_tmp * t19 * 2.0) +
                          this->robot_inertial_par(3,0) * t3 * t11;
        A_lambda_dyn_tmp = h1 * this->robot_dim(4,0);
        b_A_lambda_dyn_tmp = A_lambda_dyn_tmp * this->robot_inertial_par(9,0);
        c_A_lambda_dyn_tmp = h1 * this->robot_dim(0,0) * t3;
        d_A_lambda_dyn_tmp = h2 * this->robot_inertial_par(3,0);
        e_A_lambda_dyn_tmp = h1 * this->robot_inertial_par(2,0) * t2;
        f_A_lambda_dyn_tmp = A_lambda_dyn_tmp * t3;
        A_lambda_dyn_tmp *= t4;
        this->A_lambda_dyn_k(0,4) =
            (((((t7 - d_A_lambda_dyn_tmp * this->robot_inertial_par(8,0) * t13 * 2.0) +
                b_A_lambda_dyn_tmp * t18) +
              c_A_lambda_dyn_tmp * t18) +
              e_A_lambda_dyn_tmp * t18) +
            f_A_lambda_dyn_tmp * t18) +
            A_lambda_dyn_tmp * t18;
        this->A_lambda_dyn_k(1,4) =
            (((((t6 + b_A_lambda_dyn_tmp * t19) + c_A_lambda_dyn_tmp * t19) +
              d_A_lambda_dyn_tmp * t3 * t11) +
              e_A_lambda_dyn_tmp * t19) +
            f_A_lambda_dyn_tmp * t19) +
            A_lambda_dyn_tmp * t19;
        this->A_lambda_dyn_k(0,5) = -1.0;
        this->A_lambda_dyn_k(1,6) = -1.0;
      }

      void A_feas_calc(double alpha,double b_gamma,double h1,double h2){
        
        double phi_r=this->q_k(3,0);
        double beta=this->q_k(4,0);
      
        double A_feas_tmp;
        double A_feas_tmp_tmp;
        double A_feas_tmp_tmp_tmp;
        double b_A_feas_tmp;
        double b_A_feas_tmp_tmp;
        double c_A_feas_tmp;
        double c_A_feas_tmp_tmp;
        double d_A_feas_tmp;
        double d_A_feas_tmp_tmp;
        double e_A_feas_tmp;
        double f_A_feas_tmp;
        double g_A_feas_tmp;
        double h_A_feas_tmp;
        double i_A_feas_tmp;
        double j_A_feas_tmp;
        double k_A_feas_tmp;
        double l_A_feas_tmp;
        double m_A_feas_tmp;
        double n_A_feas_tmp;
        double o_A_feas_tmp;
        double p_A_feas_tmp;
        double t10;
        double t11;
        double t12;
        double t13;
        double t15;
        double t17;
        double t19;
        double t2;
        double t20;
        double t22;
        double t25;
        double t26;
        double t3;
        double t30;
        double t31;
        double t32;
        double t32_tmp;
        double t33;
        double t34;
        double t34_tmp;
        double t4;
        double t5;
        double t6;
        double t7;
        double t8;
        double t9;
        //     This function was generated by the Symbolic Math Toolbox version 8.7.
        //     29-Jul-2021 17:21:19
        t2 = std::cos(alpha);
        t3 = std::cos(beta);
        t4 = std::cos(b_gamma);
        t5 = std::cos(phi_r);
        t6 = std::sin(alpha);
        t7 = std::sin(beta);
        t8 = std::cos(this->robot_dim(5,0));
        t9 = std::cos(this->robot_dim(6,0));
        t10 = std::sin(b_gamma);
        t11 = std::sin(phi_r);
        t12 = std::sin(this->robot_dim(5,0));
        t13 = std::sin(this->robot_dim(6,0));
        t15 = this->robot_dim(0,0) * this->robot_dim(0,0);
        t17 = this->robot_inertial_par(3,0) * this->robot_inertial_par(3,0);
        t19 = this->robot_inertial_par(2,0) * this->robot_inertial_par(2,0);
        t20 = this->robot_dim(4,0) * this->robot_dim(4,0);
        t22 = -(this->robot_inertial_par(12,0) * 2.0);
        t25 = std::cos((beta + phi_r) + this->robot_dim(5,0));
        t26 = std::sin((b_gamma + phi_r) + this->robot_dim(5,0));
        t30 = -(this->robot_inertial_par(9,0) *
                (this->robot_dim(3,0) * this->robot_dim(3,0)));
        t31 = -(this->robot_inertial_par(7,0) *
                (this->robot_inertial_par(1,0) * this->robot_inertial_par(1,0)) * 2.0);
        t32_tmp = this->robot_dim(3,0) * this->robot_inertial_par(9,0);
        t32 = t32_tmp * t25;
        t33 = std::cos((phi_r + -alpha) + this->robot_dim(5,0));
        t34_tmp = this->robot_inertial_par(1,0) * this->robot_inertial_par(7,0);
        t34 = t34_tmp * t25 * 2.0;
        A_feas_tmp_tmp = this->robot_inertial_par(4,0) * this->robot_inertial_par(9,0);
        A_feas_tmp = A_feas_tmp_tmp * this->robot_inertial_par(0,0);
        b_A_feas_tmp_tmp = this->robot_dim(0,0) * this->robot_inertial_par(8,0);
        b_A_feas_tmp = b_A_feas_tmp_tmp * this->robot_inertial_par(0,0);
        c_A_feas_tmp = b_A_feas_tmp * t2;
        c_A_feas_tmp_tmp = this->robot_inertial_par(2,0) * this->robot_inertial_par(6,0);
        d_A_feas_tmp = c_A_feas_tmp_tmp * this->robot_inertial_par(0,0);
        e_A_feas_tmp = d_A_feas_tmp * t2;
        A_feas_tmp_tmp_tmp =
            this->robot_inertial_par(3,0) * this->robot_inertial_par(8,0);
        d_A_feas_tmp_tmp = A_feas_tmp_tmp_tmp * this->robot_inertial_par(0,0);
        f_A_feas_tmp = d_A_feas_tmp_tmp * t4;
        g_A_feas_tmp = t32_tmp * this->robot_inertial_par(0,0);
        h_A_feas_tmp = g_A_feas_tmp * t3;
        i_A_feas_tmp = t34_tmp * this->robot_inertial_par(0,0);
        j_A_feas_tmp = i_A_feas_tmp * t3;
        t25 = this->robot_dim(4,0) * this->robot_inertial_par(8,0);
        k_A_feas_tmp = t25 * this->robot_inertial_par(0,0);
        l_A_feas_tmp = k_A_feas_tmp * t2;
        t32_tmp = this->robot_dim(4,0) * this->robot_inertial_par(9,0);
        m_A_feas_tmp = t32_tmp * this->robot_inertial_par(0,0);
        n_A_feas_tmp = m_A_feas_tmp * t2;
        t34_tmp = this->robot_dim(4,0) * this->robot_inertial_par(7,0);
        o_A_feas_tmp = t34_tmp * this->robot_inertial_par(0,0);
        p_A_feas_tmp = o_A_feas_tmp * t2;
        this->A_feas_k(0,2) =
            (((((((((((((((((((((((this->robot_inertial_par(10,0) * 2.0 +
                                  this->robot_inertial_par(0,0) *
                                      this->robot_inertial_par(0,0) *
                                      ((((this->robot_inertial_par(6,0) * 2.0 +
                                          this->robot_inertial_par(8,0) * 2.0) +
                                          this->robot_inertial_par(9,0)) +
                                        this->robot_inertial_par(7,0) * 2.0) +
                                        this->robot_inertial_par(5,0) * 2.0)) +
                                  A_feas_tmp * t5 * t13) -
                                A_feas_tmp * t9 * t11) +
                                c_A_feas_tmp * t5 * t8 * 2.0) +
                              b_A_feas_tmp * t5 * t6 * t12 * 2.0) -
                              c_A_feas_tmp * t11 * t12 * 2.0) +
                            b_A_feas_tmp * t6 * t8 * t11 * 2.0) +
                            e_A_feas_tmp * t5 * t8 * 2.0) +
                          d_A_feas_tmp * t5 * t6 * t12 * 2.0) -
                          e_A_feas_tmp * t11 * t12 * 2.0) +
                        d_A_feas_tmp * t6 * t8 * t11 * 2.0) +
                        f_A_feas_tmp * t5 * t12 * 2.0) +
                      f_A_feas_tmp * t8 * t11 * 2.0) +
                      d_A_feas_tmp_tmp * t5 * t8 * t10 * 2.0) -
                    d_A_feas_tmp_tmp * t10 * t11 * t12 * 2.0) -
                    h_A_feas_tmp * t5 * t8) +
                  g_A_feas_tmp * t5 * t7 * t12) +
                  h_A_feas_tmp * t11 * t12) +
                g_A_feas_tmp * t7 * t8 * t11) -
                j_A_feas_tmp * t5 * t8 * 2.0) +
              i_A_feas_tmp * t5 * t7 * t12 * 2.0) +
              j_A_feas_tmp * t11 * t12 * 2.0) +
            i_A_feas_tmp * t7 * t8 * t11 * 2.0) +
            (((((((((((l_A_feas_tmp * t5 * t8 * 2.0 +
                      k_A_feas_tmp * t5 * t6 * t12 * 2.0) -
                      l_A_feas_tmp * t11 * t12 * 2.0) +
                    k_A_feas_tmp * t6 * t8 * t11 * 2.0) +
                    n_A_feas_tmp * t5 * t8) +
                  m_A_feas_tmp * t5 * t6 * t12) -
                  n_A_feas_tmp * t11 * t12) +
                m_A_feas_tmp * t6 * t8 * t11) +
                p_A_feas_tmp * t5 * t8 * 2.0) +
              o_A_feas_tmp * t5 * t6 * t12 * 2.0) -
              p_A_feas_tmp * t11 * t12 * 2.0) +
            o_A_feas_tmp * t6 * t8 * t11 * 2.0);
        A_feas_tmp = this->robot_dim(0,0) * this->robot_inertial_par(3,0) *
                    this->robot_inertial_par(8,0);
        b_A_feas_tmp = this->robot_inertial_par(3,0) * this->robot_dim(4,0) *
                      this->robot_inertial_par(8,0);
        c_A_feas_tmp =
            this->robot_dim(3,0) * this->robot_dim(4,0) * this->robot_inertial_par(9,0);
        d_A_feas_tmp = this->robot_inertial_par(1,0) * this->robot_dim(4,0) *
                      this->robot_inertial_par(7,0);
        d_A_feas_tmp_tmp = this->robot_dim(3,0) * this->robot_inertial_par(4,0) *
                          this->robot_inertial_par(9,0);
        e_A_feas_tmp = d_A_feas_tmp_tmp * t3;
        f_A_feas_tmp = d_A_feas_tmp_tmp * t7;
        d_A_feas_tmp_tmp = this->robot_inertial_par(4,0) * this->robot_dim(4,0) *
                          this->robot_inertial_par(9,0);
        g_A_feas_tmp = d_A_feas_tmp_tmp * t2;
        h_A_feas_tmp = d_A_feas_tmp_tmp * t6;
        i_A_feas_tmp = c_A_feas_tmp * t2 * t3;
        c_A_feas_tmp = c_A_feas_tmp * t6 * t7;
        j_A_feas_tmp = d_A_feas_tmp * t2 * t3;
        d_A_feas_tmp = d_A_feas_tmp * t6 * t7;
        k_A_feas_tmp = e_A_feas_tmp * t8 * t13;
        e_A_feas_tmp = e_A_feas_tmp * t9 * t12;
        l_A_feas_tmp = f_A_feas_tmp * t8 * t9;
        f_A_feas_tmp = f_A_feas_tmp * t12 * t13;
        this->A_feas_k(0,3) =
            (((((((((((((((((((((((((this->robot_inertial_par(11,0) * -2.0 -
                                    this->robot_inertial_par(13,0) * 2.0) -
                                    this->robot_inertial_par(14,0)) +
                                  t22) +
                                  t30) +
                                t31) -
                                this->robot_inertial_par(0,0) *
                                    ((((((((-t32 - t34) -
                                          A_feas_tmp_tmp *
                                              std::sin(phi_r - this->robot_dim(6,0))) +
                                          b_A_feas_tmp_tmp * t33 * 2.0) +
                                        c_A_feas_tmp_tmp * t33 * 2.0) +
                                        A_feas_tmp_tmp_tmp * t26 * 2.0) +
                                      t25 * t33 * 2.0) +
                                      t32_tmp * t33) +
                                    t34_tmp * t33 * 2.0)) -
                              this->robot_inertial_par(6,0) * t19 * 2.0) -
                              this->robot_inertial_par(8,0) * t15 * 2.0) -
                            this->robot_inertial_par(8,0) * t17 * 2.0) -
                            this->robot_inertial_par(8,0) * t20 * 2.0) -
                          this->robot_inertial_par(9,0) * t20) -
                          this->robot_inertial_par(7,0) * t20 * 2.0) -
                        this->robot_inertial_par(4,0) * this->robot_inertial_par(4,0) *
                            this->robot_inertial_par(9,0)) -
                        this->robot_dim(0,0) * this->robot_dim(4,0) *
                            this->robot_inertial_par(8,0) * 4.0) -
                      A_feas_tmp * t4 * t6 * 4.0) -
                      A_feas_tmp * t2 * t10 * 4.0) -
                    b_A_feas_tmp * t4 * t6 * 4.0) -
                    b_A_feas_tmp * t2 * t10 * 4.0) +
                  i_A_feas_tmp * 2.0) -
                  c_A_feas_tmp * 2.0) +
                j_A_feas_tmp * 4.0) -
                d_A_feas_tmp * 4.0) +
              k_A_feas_tmp * 2.0) +
              e_A_feas_tmp * 2.0) +
            l_A_feas_tmp * 2.0) +
            ((((f_A_feas_tmp * -2.0 - g_A_feas_tmp * t8 * t13 * 2.0) -
              g_A_feas_tmp * t9 * t12 * 2.0) +
              h_A_feas_tmp * t8 * t9 * 2.0) -
            h_A_feas_tmp * t12 * t13 * 2.0);
        A_feas_tmp = h1 * this->robot_dim(4,0);
        b_A_feas_tmp = h1 * this->robot_inertial_par(8,0);
        g_A_feas_tmp = h1 * this->robot_dim(0,0);
        h_A_feas_tmp =
            g_A_feas_tmp * this->robot_inertial_par(3,0) * this->robot_inertial_par(8,0);
        m_A_feas_tmp = h2 * this->robot_dim(0,0) * this->robot_inertial_par(3,0) *
                      this->robot_inertial_par(8,0);
        n_A_feas_tmp = h2 * this->robot_inertial_par(3,0);
        o_A_feas_tmp = h1 * this->robot_inertial_par(3,0) * this->robot_dim(4,0) *
                      this->robot_inertial_par(8,0);
        p_A_feas_tmp =
            n_A_feas_tmp * this->robot_dim(4,0) * this->robot_inertial_par(8,0);
        t32_tmp = h1 * this->robot_dim(3,0) * this->robot_dim(4,0) *
                  this->robot_inertial_par(9,0);
        t34_tmp = h1 * this->robot_inertial_par(1,0) * this->robot_dim(4,0) *
                  this->robot_inertial_par(7,0);
        A_feas_tmp_tmp = h1 * this->robot_inertial_par(4,0) * this->robot_dim(4,0) *
                        this->robot_inertial_par(9,0);
        b_A_feas_tmp_tmp = A_feas_tmp_tmp * t2;
        t25 = A_feas_tmp_tmp * t6;
        this->A_feas_k(0,4) =
            ((((((((((((((((((((((((((t22 + t30) + t31) +
                                    this->robot_inertial_par(11,0) * h1 * 2.0) -
                                  this->robot_inertial_par(13,0) * h2 * 2.0) +
                                  this->robot_inertial_par(0,0) *
                                      (((((((t32 + t34) +
                                            g_A_feas_tmp *
                                                this->robot_inertial_par(8,0) * t33 *
                                                2.0) +
                                          h1 * this->robot_inertial_par(2,0) *
                                              this->robot_inertial_par(6,0) * t33 *
                                              2.0) -
                                          n_A_feas_tmp * this->robot_inertial_par(8,0) *
                                              t26 * 2.0) +
                                        A_feas_tmp * this->robot_inertial_par(8,0) *
                                            t33 * 2.0) +
                                        A_feas_tmp * this->robot_inertial_par(9,0) *
                                            t33) +
                                      A_feas_tmp * this->robot_inertial_par(7,0) *
                                          t33 * 2.0)) +
                                h1 * this->robot_inertial_par(6,0) * t19 * 2.0) +
                                b_A_feas_tmp * t15 * 2.0) -
                              h2 * this->robot_inertial_par(8,0) * t17 * 2.0) +
                              b_A_feas_tmp * t20 * 2.0) +
                            h1 * this->robot_inertial_par(9,0) * t20) +
                            h1 * this->robot_inertial_par(7,0) * t20 * 2.0) +
                          g_A_feas_tmp * this->robot_dim(4,0) *
                              this->robot_inertial_par(8,0) * 4.0) +
                          i_A_feas_tmp) -
                        c_A_feas_tmp) +
                        j_A_feas_tmp * 2.0) -
                      d_A_feas_tmp * 2.0) +
                      k_A_feas_tmp) +
                    e_A_feas_tmp) +
                    l_A_feas_tmp) -
                  f_A_feas_tmp) +
                  h_A_feas_tmp * t4 * t6 * 2.0) -
                m_A_feas_tmp * t4 * t6 * 2.0) +
                h_A_feas_tmp * t2 * t10 * 2.0) -
              m_A_feas_tmp * t2 * t10 * 2.0) +
              o_A_feas_tmp * t4 * t6 * 2.0) -
            p_A_feas_tmp * t4 * t6 * 2.0) +
            (((((((((o_A_feas_tmp * t2 * t10 * 2.0 - p_A_feas_tmp * t2 * t10 * 2.0) -
                    t32_tmp * t2 * t3) +
                  t32_tmp * t6 * t7) -
                  t34_tmp * t2 * t3 * 2.0) +
                t34_tmp * t6 * t7 * 2.0) +
                b_A_feas_tmp_tmp * t8 * t13) +
              b_A_feas_tmp_tmp * t9 * t12) -
              t25 * t8 * t9) +
            t25 * t12 * t13);
      }

      void A_lambda_con_calc(){
        this->A_lambda_con_k(0,6)=-1;
        this->A_lambda_con_k(1,5)=1;
        this->A_lambda_con_k(1,6)=-(this->mu);
        this->A_lambda_con_k(2,5)=-1;
        this->A_lambda_con_k(2,6)=-(this->mu);
      }

      void A_u_check_k_calc(double alpha,double b_gamma,double h1,double h2){

        double phi_r=this->q_k(3,0);
        double beta=this->q_k(4,0);

        double b_t47_tmp;
        double b_t69_tmp;
        double c_t47_tmp;
        double d_t47_tmp;
        double e_t47_tmp;
        double f_t47_tmp;
        double g_t47_tmp;
        double h_t47_tmp;
        double i_t47_tmp;
        double j_t47_tmp;
        double k_t47_tmp;
        double t11;
        double t12;
        double t15;
        double t16;
        double t17;
        double t19;
        double t22;
        double t23;
        double t25;
        double t27;
        double t28;
        double t31;
        double t32;
        double t4;
        double t46;
        double t47_tmp;
        double t49;
        double t5;
        double t50;
        double t51;
        double t52_tmp;
        double t53_tmp;
        double t55;
        double t56;
        double t6;
        double t69;
        double t69_tmp;
        double t7;
        double t9;
        //     This function was generated by the Symbolic Math Toolbox version 8.7.
        //     29-Jul-2021 16:01:02
        t4 = this->robot_inertial_par(12,0) * 2.0;
        t5 = h1 * h1;
        t6 = h2 * h2;
        t7 = this->robot_dim(0,0) * this->robot_dim(0,0);
        t9 = this->robot_inertial_par(3,0) * this->robot_inertial_par(3,0);
        t11 = this->robot_inertial_par(2,0) * this->robot_inertial_par(2,0);
        t12 = this->robot_dim(4,0) * this->robot_dim(4,0);
        t17 = (beta + phi_r) + this->robot_dim(5,0);
        t19 = (b_gamma + phi_r) + this->robot_dim(5,0);
        t15 = std::cos(alpha + beta);
        t16 = std::sin(alpha + b_gamma);
        t22 = std::cos(t17);
        t23 = std::cos(t19);
        t25 = std::sin(t17);
        t27 = std::sin(t19);
        t28 = this->robot_inertial_par(9,0) * (this->robot_dim(3,0) * this->robot_dim(3,0));
        t31 = this->robot_inertial_par(7,0) *
              (this->robot_inertial_par(1,0) * this->robot_inertial_par(1,0)) * 2.0;
        t32 = phi_r + -this->robot_dim(6,0);
        t17 = (phi_r + -alpha) + this->robot_dim(5,0);
        t46 = this->robot_dim(3,0) * this->robot_dim(4,0) * this->robot_inertial_par(9,0) *
              t15;
        t47_tmp = this->robot_dim(3,0) * this->robot_inertial_par(9,0);
        t49 = std::cos(t17);
        t50 = std::sin(t17);
        t51 = std::sin((-alpha + this->robot_dim(5,0)) + this->robot_dim(6,0));
        t52_tmp = this->robot_inertial_par(1,0) * this->robot_dim(4,0) *
                  this->robot_inertial_par(7,0) * t15;
        t17 = t52_tmp * 2.0;
        t53_tmp = this->robot_inertial_par(1,0) * this->robot_inertial_par(7,0);
        t56 = this->robot_dim(3,0) * this->robot_inertial_par(4,0) *
              this->robot_inertial_par(9,0) *
              std::sin((beta + this->robot_dim(5,0)) + this->robot_dim(6,0));
        t55 = h1 * t46;
        t19 = h1 * this->robot_inertial_par(8,0);
        t69_tmp = h1 * this->robot_dim(0,0);
        b_t69_tmp = h2 * this->robot_inertial_par(3,0);
        t69 = ((((((((((((((((((((t4 + this->robot_inertial_par(13,0) * h2 * 2.0) +
                                -(this->robot_inertial_par(11,0) * h1 * 2.0)) +
                                t28) +
                              t31) +
                              h2 * this->robot_inertial_par(8,0) * t9 * 2.0) +
                            -(t69_tmp * this->robot_dim(4,0) *
                              this->robot_inertial_par(8,0) * 4.0)) +
                            -(t19 * t7 * 2.0)) +
                          -(h1 * this->robot_inertial_par(6,0) * t11 * 2.0)) +
                          -(t19 * t12 * 2.0)) +
                        -(h1 * this->robot_inertial_par(9,0) * t12)) +
                        -(h1 * this->robot_inertial_par(7,0) * t12 * 2.0)) +
                      t55) +
                      -t46) +
                    -t17) +
                    h1 * t17) +
                  h2 * this->robot_dim(0,0) * this->robot_inertial_par(3,0) *
                      this->robot_inertial_par(8,0) * t16 * 2.0) +
                  b_t69_tmp * this->robot_dim(4,0) * this->robot_inertial_par(8,0) * t16 *
                      2.0) +
                -(t69_tmp * this->robot_inertial_par(3,0) *
                  this->robot_inertial_par(8,0) * t16 * 2.0)) +
                -(h1 * this->robot_inertial_par(3,0) * this->robot_dim(4,0) *
                  this->robot_inertial_par(8,0) * t16 * 2.0)) +
              -t56) +
              -(h1 * this->robot_inertial_par(4,0) * this->robot_dim(4,0) *
                this->robot_inertial_par(9,0) * t51);
        b_t47_tmp = t47_tmp * t22 + t53_tmp * t22 * 2.0;
        c_t47_tmp = this->robot_inertial_par(4,0) * this->robot_inertial_par(9,0);
        d_t47_tmp = this->robot_dim(0,0) * this->robot_inertial_par(8,0);
        e_t47_tmp = this->robot_inertial_par(2,0) * this->robot_inertial_par(6,0);
        f_t47_tmp = this->robot_inertial_par(3,0) * this->robot_inertial_par(8,0);
        g_t47_tmp = this->robot_dim(4,0) * this->robot_inertial_par(8,0);
        h_t47_tmp = this->robot_dim(4,0) * this->robot_inertial_par(9,0);
        i_t47_tmp = this->robot_dim(4,0) * this->robot_inertial_par(7,0);
        this->A_u_check_k(0,0) =
            ((((((b_t47_tmp + c_t47_tmp * std::sin(t32)) - d_t47_tmp * t49 * 2.0) -
                e_t47_tmp * t49 * 2.0) -
              f_t47_tmp * t27 * 2.0) -
              g_t47_tmp * t49 * 2.0) -
            h_t47_tmp * t49) -
            i_t47_tmp * t49 * 2.0;
        j_t47_tmp = h1 * this->robot_dim(4,0);
        t69_tmp *= this->robot_inertial_par(8,0);
        k_t47_tmp = h1 * this->robot_inertial_par(2,0) * this->robot_inertial_par(6,0);
        t17 = b_t69_tmp * this->robot_inertial_par(8,0);
        t19 = j_t47_tmp * this->robot_inertial_par(8,0);
        t22 = j_t47_tmp * this->robot_inertial_par(9,0);
        j_t47_tmp *= this->robot_inertial_par(7,0);
        this->A_u_check_k(1,0) =
            (((((b_t47_tmp + t69_tmp * t49 * 2.0) + k_t47_tmp * t49 * 2.0) -
              t17 * t27 * 2.0) +
              t19 * t49 * 2.0) +
            t22 * t49) +
            j_t47_tmp * t49 * 2.0;
        t47_tmp = t47_tmp * t25 + t53_tmp * t25 * 2.0;
        this->A_u_check_k(0,1) =
            ((((((t47_tmp - c_t47_tmp * std::cos(t32)) - d_t47_tmp * t50 * 2.0) -
                e_t47_tmp * t50 * 2.0) +
              f_t47_tmp * t23 * 2.0) -
              g_t47_tmp * t50 * 2.0) -
            h_t47_tmp * t50) -
            i_t47_tmp * t50 * 2.0;
        this->A_u_check_k(1,1) = (((((t47_tmp + t69_tmp * t50 * 2.0) + k_t47_tmp * t50 * 2.0) +
                          t17 * t23 * 2.0) +
                        t19 * t50 * 2.0) +
                        t22 * t50) +
                      j_t47_tmp * t50 * 2.0;
        t47_tmp =
            this->robot_dim(0,0) * this->robot_dim(4,0) * this->robot_inertial_par(8,0);
        this->A_u_check_k(0,3) =
            ((((((((((((((((((this->robot_inertial_par(11,0) * 2.0 +
                              this->robot_inertial_par(13,0) * 2.0) +
                            this->robot_inertial_par(14,0)) +
                            t4) +
                          t28) +
                          t31) -
                        t46 * 2.0) -
                        t56 * 2.0) +
                      this->robot_inertial_par(6,0) * t11 * 2.0) +
                      this->robot_inertial_par(8,0) * t7 * 2.0) +
                    this->robot_inertial_par(8,0) * t9 * 2.0) +
                    this->robot_inertial_par(8,0) * t12 * 2.0) +
                  this->robot_inertial_par(9,0) * t12) +
                  this->robot_inertial_par(7,0) * t12 * 2.0) +
                this->robot_inertial_par(4,0) * this->robot_inertial_par(4,0) *
                    this->robot_inertial_par(9,0)) +
                t47_tmp * 4.0) +
              this->robot_dim(0,0) * this->robot_inertial_par(3,0) *
                  this->robot_inertial_par(8,0) * t16 * 4.0) +
              this->robot_inertial_par(3,0) * this->robot_dim(4,0) *
                  this->robot_inertial_par(8,0) * t16 * 4.0) -
            t52_tmp * 4.0) +
            this->robot_inertial_par(4,0) * this->robot_dim(4,0) *
                this->robot_inertial_par(9,0) * t51 * 2.0;
        this->A_u_check_k(1,3) = t69;
        this->A_u_check_k(0,4) = t69;
        b_t47_tmp = this->robot_inertial_par(8,0) * t5;
        c_t47_tmp = h1 * h2;
        this->A_u_check_k(1,4) = ((((((((((((((t4 + t28) + t31) + t55 * 2.0) +
                                  this->robot_inertial_par(11,0) * t5 * 2.0) +
                                this->robot_inertial_par(13,0) * t6 * 2.0) +
                                this->robot_inertial_par(6,0) * t5 * t11 * 2.0) +
                              b_t47_tmp * t7 * 2.0) +
                              this->robot_inertial_par(8,0) * t6 * t9 * 2.0) +
                            b_t47_tmp * t12 * 2.0) +
                            this->robot_inertial_par(9,0) * t5 * t12) +
                          this->robot_inertial_par(7,0) * t5 * t12 * 2.0) +
                          t47_tmp * t5 * 4.0) +
                        h1 * this->robot_inertial_par(1,0) * this->robot_dim(4,0) *
                            this->robot_inertial_par(7,0) * t15 * 4.0) -
                        c_t47_tmp * this->robot_dim(0,0) * this->robot_inertial_par(3,0) *
                            this->robot_inertial_par(8,0) * t16 * 4.0) -
                      c_t47_tmp * this->robot_inertial_par(3,0) * this->robot_dim(4,0) *
                          this->robot_inertial_par(8,0) * t16 * 4.0;
      }

      void B_check_inv_calc(double h1){
        
        double beta=this->q_k(4,0);

        this->B_check_inv_k(0,0) = 0.5;
        this->B_check_inv_k(1,0) = (h1) / 2.0;
        this->B_check_inv_k(0,1) = 0.0;
        this->B_check_inv_k(1,1) = 0.5;
      }

      void A_u_lim_calc(double alpha,double b_gamma,double h1,double h2){

        double phi_r=this->q_k(3,0);
        double beta=this->q_k(4,0);
 
        double t10;
        double t100;
        double t101;
        double t102;
        double t103;
        double t113;
        double t114;
        double t115;
        double t116;
        double t118;
        double t120;
        double t120_tmp;
        double t122;
        double t122_tmp;
        double t129;
        double t130;
        double t134;
        double t134_tmp;
        double t135;
        double t136;
        double t137;
        double t138;
        double t139;
        double t140;
        double t147;
        double t155;
        double t166;
        double t168;
        double t17;
        double t173;
        double t174;
        double t175;
        double t179;
        double t180;
        double t181;
        double t182;
        double t183;
        double t184;
        double t192;
        double t195;
        double t2;
        double t20;
        double t200;
        double t201;
        double t202;
        double t203;
        double t206;
        double t207;
        double t208;
        double t209;
        double t21;
        double t22;
        double t251;
        double t252;
        double t259;
        double t263;
        double t264;
        double t272;
        double t28;
        double t3;
        double t31;
        double t36;
        double t41;
        double t42;
        double t43;
        double t44;
        double t45;
        double t46;
        double t47;
        double t48;
        double t49;
        double t50;
        double t50_tmp;
        double t51;
        double t52;
        double t53;
        double t54;
        double t55;
        double t56;
        double t57;
        double t58;
        double t64;
        double t67;
        double t75;
        double t8;
        double t9;
        double t96;
        double t97;
        double t98;
        double t99;
        //     This function was generated by the Symbolic Math Toolbox version 8.7.
        //     29-Jul-2021 16:01:03
        t2 = this->robot_inertial_par(11,0) * h1;
        t3 = this->robot_inertial_par(13,0) * h2;
        t8 = this->robot_inertial_par(12,0) * 2.0;
        t9 = h1 * h1;
        t10 = h2 * h2;
        t17 = this->robot_dim(4,0) * this->robot_dim(4,0);
        t22 = (beta + phi_r) + this->robot_dim(5,0);
        t28 = (b_gamma + phi_r) + this->robot_dim(5,0);
        t50_tmp =
            this->robot_dim(0,0) * this->robot_dim(4,0) * this->robot_inertial_par(8,0);
        t50 = t50_tmp * 2.0;
        t51 = t50_tmp * 4.0;
        t52 = this->robot_inertial_par(14,0) / 2.0;
        t20 = std::cos(alpha + beta);
        t21 = std::sin(alpha + b_gamma);
        t31 = std::cos(t22);
        t36 = std::sin(t22);
        t41 = this->robot_inertial_par(8,0) * (this->robot_dim(0,0) * this->robot_dim(0,0));
        t42 = this->robot_inertial_par(6,0) *
              (this->robot_inertial_par(2,0) * this->robot_inertial_par(2,0));
        t43 = this->robot_inertial_par(8,0) *
              (this->robot_inertial_par(3,0) * this->robot_inertial_par(3,0));
        t44 = this->robot_inertial_par(9,0) * (this->robot_dim(3,0) * this->robot_dim(3,0));
        t45 = this->robot_inertial_par(7,0) *
              (this->robot_inertial_par(1,0) * this->robot_inertial_par(1,0));
        t46 = this->robot_inertial_par(9,0) *
              (this->robot_inertial_par(4,0) * this->robot_inertial_par(4,0));
        t47 = this->robot_inertial_par(8,0) * t17;
        t48 = this->robot_inertial_par(9,0) * t17;
        t49 = this->robot_inertial_par(7,0) * t17;
        t67 = phi_r + -this->robot_dim(6,0);
        t75 = (phi_r + -alpha) + this->robot_dim(5,0);
        t53 = h1 * t41;
        t54 = h1 * t42;
        t55 = h2 * t43;
        t56 = h1 * t47;
        t57 = h1 * t48;
        t58 = h1 * t49;
        t64 = t45 * 2.0;
        t96 = this->robot_dim(3,0) * this->robot_dim(4,0) * this->robot_inertial_par(9,0) *
              t20;
        t97 = this->robot_inertial_par(1,0) * this->robot_dim(4,0) *
              this->robot_inertial_par(7,0) * t20;
        t22 = this->robot_dim(3,0) * this->robot_inertial_par(9,0);
        t98 = t22 * t31;
        t17 = this->robot_inertial_par(1,0) * this->robot_inertial_par(7,0);
        t99 = t17 * t31;
        t20 = this->robot_inertial_par(3,0) * this->robot_inertial_par(8,0);
        t100 = t20 * std::cos(t28);
        t101 = t22 * t36;
        t102 = t17 * t36;
        t103 = t20 * std::sin(t28);
        t31 = std::cos(t75);
        t75 = std::sin(t75);
        t113 = t44 / 2.0;
        t114 = t46 / 2.0;
        t115 = t48 / 2.0;
        t120_tmp = this->robot_dim(0,0) * this->robot_inertial_par(3,0) *
                  this->robot_inertial_par(8,0) * t21;
        t120 = t120_tmp * 2.0;
        t122_tmp = this->robot_inertial_par(3,0) * this->robot_dim(4,0) *
                  this->robot_inertial_par(8,0) * t21;
        t122 = t122_tmp * 2.0;
        t134_tmp = h1 * this->robot_dim(0,0);
        t134 = t134_tmp * this->robot_inertial_par(3,0) * this->robot_inertial_par(8,0) *
              t21;
        t135 = h2 * this->robot_dim(0,0) * this->robot_inertial_par(3,0) *
              this->robot_inertial_par(8,0) * t21;
        t136 = h1 * this->robot_inertial_par(3,0) * this->robot_dim(4,0) *
              this->robot_inertial_par(8,0) * t21;
        t137 = h2 * this->robot_inertial_par(3,0) * this->robot_dim(4,0) *
              this->robot_inertial_par(8,0) * t21;
        t140 = this->robot_dim(3,0) * this->robot_inertial_par(4,0) *
              this->robot_inertial_par(9,0) *
              std::sin((beta + this->robot_dim(5,0)) + this->robot_dim(6,0));
        t116 = t57 / 2.0;
        t118 = t97 * 2.0;
        t129 = h1 * t96;
        t130 = h1 * t97;
        t138 = h2 * t100;
        t139 = h2 * t103;
        t22 = this->robot_inertial_par(4,0) * this->robot_inertial_par(9,0);
        t147 = t22 * std::cos(t67);
        t155 = t22 * std::sin(t67);
        t166 = h1 * t120;
        t168 = h1 * t122;
        t173 = t96 / 2.0;
        t174 = t98 / 2.0;
        t175 = t101 / 2.0;
        t22 = this->robot_dim(0,0) * this->robot_inertial_par(8,0);
        t180 = t22 * t31;
        t17 = this->robot_inertial_par(2,0) * this->robot_inertial_par(6,0);
        t181 = t17 * t31;
        t20 = this->robot_dim(4,0) * this->robot_inertial_par(8,0);
        t182 = t20 * t31;
        t36 = this->robot_dim(4,0) * this->robot_inertial_par(9,0);
        t183 = t36 * t31;
        t28 = this->robot_dim(4,0) * this->robot_inertial_par(7,0);
        t184 = t28 * t31;
        t192 = t22 * t75;
        t67 = t17 * t75;
        t21 = t20 * t75;
        t195 = t36 * t75;
        t22 = t28 * t75;
        t31 = this->robot_inertial_par(4,0) * this->robot_dim(4,0) *
              this->robot_inertial_par(9,0) *
              std::sin((-alpha + this->robot_dim(5,0)) + this->robot_dim(6,0));
        t28 = t140 / 2.0;
        t179 = t130 * -2.0;
        t200 = h1 * t180;
        t201 = h1 * t181;
        t202 = h1 * t182;
        t203 = h1 * t184;
        t206 = h1 * t192;
        t207 = h1 * t67;
        t208 = h1 * t21;
        t209 = h1 * t22;
        t75 = t129 / 2.0;
        t17 = t147 / 2.0;
        t20 = t155 / 2.0;
        t36 = h1 * t31;
        t251 = t183 / 2.0;
        t252 = t195 / 2.0;
        t259 = t36 / 2.0;
        t264 =
            (((((((-t100 + -t102) + t192) + t67) + t21) + t22) + -t175) + t17) + t252;
        t272 = ((((((((((((((((((this->robot_inertial_par(11,0) +
                                this->robot_inertial_par(13,0)) +
                                this->robot_inertial_par(12,0)) +
                              t41) +
                              t42) +
                            t43) +
                            t45) +
                          t47) +
                          t49) +
                        t50) +
                        t52) +
                      t113) +
                      t114) +
                    t115) +
                    t120) +
                  t122) +
                  -t96) +
                -t118) +
                -t140) +
              t31;
        t52 = ((((((((((((((((((-this->robot_inertial_par(11,0) +
                                -this->robot_inertial_par(13,0)) +
                              -this->robot_inertial_par(12,0)) +
                              -t50) +
                            -t52) +
                            -t41) +
                          -t42) +
                          -t43) +
                        -t45) +
                        -t47) +
                      -t49) +
                      t96) +
                    t118) +
                    -t113) +
                  -t114) +
                  -t115) +
                t140) +
                -t120) +
              -t122) +
              -t31;
        t263 = (((((((t103 + -t99) + t180) + t181) + t182) + t184) + -t174) + -t20) +
              t251;
        t147 = h1 *
              ((((((((t147 + -(t100 * 2.0)) + -t101) + -(t102 * 2.0)) + t195) +
                  t192 * 2.0) +
                  t67 * 2.0) +
                t21 * 2.0) +
                t22 * 2.0) /
              2.0;
        t155 = h1 *
              ((((((((t103 * 2.0 + -t98) + -(t99 * 2.0)) + t183) + -t155) +
                  t180 * 2.0) +
                  t181 * 2.0) +
                t182 * 2.0) +
                t184 * 2.0) /
              2.0;
        t122 = (((((((t100 + t102) + t175) + -t192) + -t67) + -t21) + -t22) + -t17) +
              -t252;
        t120 =
            (((((((t99 + -t103) + t174) + t20) + -t180) + -t181) + -t182) + -t184) +
            -t251;
        t114 = h1 *
              (((((((((((((((((((this->robot_inertial_par(14,0) +
                                  this->robot_inertial_par(11,0) * 2.0) +
                                this->robot_inertial_par(13,0) * 2.0) +
                                t8) +
                              t44) +
                              t46) +
                            t48) +
                            t51) +
                          t41 * 2.0) +
                          t42 * 2.0) +
                        t43 * 2.0) +
                        t64) +
                      t47 * 2.0) +
                      t49 * 2.0) +
                    t120_tmp * 4.0) +
                    t122_tmp * 4.0) +
                  -(t96 * 2.0)) +
                  -(t97 * 4.0)) +
                -(t140 * 2.0)) +
                t31 * 2.0) /
              2.0;
        t22 = h1 *
              (((((((((((((((((((((t2 * 2.0 + -t8) + -(t3 * 2.0)) + t57) + h1 * t51) +
                              t53 * 2.0) +
                              t54 * 2.0) +
                            t56 * 2.0) +
                            t58 * 2.0) +
                          -t44) +
                          -t64) +
                        -(t55 * 2.0)) +
                        t96) +
                      t118) +
                      t140) +
                    t166) +
                    t168) +
                  -t129) +
                  t179) +
                t135 * -2.0) +
                t137 * -2.0) +
              t36) /
              2.0;
        t192 = ((((((((((((((((((((t2 + -this->robot_inertial_par(12,0)) + -t3) + t53) +
                                t54) +
                              t56) +
                              t58) +
                            h1 * t50) +
                            -t45) +
                          -t55) +
                          t97) +
                        t116) +
                        -t113) +
                      t134) +
                      t136) +
                    t173) +
                    -t130) +
                  -t135) +
                  -t137) +
                t28) +
                -t75) +
              t259;
        t67 = ((((((((((((((((((((this->robot_inertial_par(12,0) + t3) + -t2) + t45) +
                              t55) +
                              t134_tmp * this->robot_dim(4,0) *
                                  this->robot_inertial_par(8,0) * -2.0) +
                            -t53) +
                            -t54) +
                          -t56) +
                          -t58) +
                        t113) +
                        t130) +
                      t135) +
                      t137) +
                    -t116) +
                    -t97) +
                  -t134) +
                  -t136) +
                -t173) +
                t75) +
              -t28) +
              -t259;
        t21 = (((((((t139 + -t99) + -t174) + -t200) + -t201) + -t202) + -t203) +
              h1 * t183 * -0.5) +
              t155;
        t75 = (((((((-t102 + -t138) + -t175) + -t206) + -t207) + -t208) + -t209) +
              h1 * t195 * -0.5) +
              t147;
        t28 = (((((((((((((((-this->robot_inertial_par(12,0) + h1 * -t2) + h2 * -t3) +
                          -t45) +
                          t50_tmp * t9 * -2.0) +
                        -t113) +
                        t9 * -t41) +
                      t9 * -t42) +
                      t10 * -t43) +
                    t9 * -t47) +
                    t9 * -t49) +
                  h2 * t166) +
                  h2 * t168) +
                -t129) +
                t179) +
              t9 * t48 * -0.5) +
              t22;
        t36 =
            (((((((t102 + t138) + t175) + t206) + t207) + t208) + t209) + h1 * t252) +
            -t147;
        t31 =
            (((((((t99 + t174) + -t139) + t200) + t201) + t202) + t203) + h1 * t251) +
            -t155;
        t17 =
            (((((((((((((((this->robot_inertial_par(12,0) + h1 * t2) + h2 * t3) + t45) +
                        t9 * t41) +
                      t9 * t42) +
                      t10 * t43) +
                    t9 * t47) +
                    t9 * t49) +
                  t9 * t50) +
                  t113) +
                t129) +
                h1 * t118) +
              t9 * t115) +
              h2 * t134 * -2.0) +
            h2 * t136 * -2.0) +
            -t22;
        t20 = -t114 + t192;
        t22 = t114 + t67;
        this->A_u_lim_k(0,0) = t120;
        this->A_u_lim_k(1,0) = t31;
        this->A_u_lim_k(2,0) = t263;
        this->A_u_lim_k(3,0) = t21;
        this->A_u_lim_k(4,0) = t263;
        this->A_u_lim_k(5,0) = t21;
        this->A_u_lim_k(6,0) = t120;
        this->A_u_lim_k(7,0) = t31;
        this->A_u_lim_k(0,1) = t122;
        this->A_u_lim_k(1,1) = t36;
        this->A_u_lim_k(2,1) = t264;
        this->A_u_lim_k(3,1) = t75;
        this->A_u_lim_k(4,1) = t264;
        this->A_u_lim_k(5,1) = t75;
        this->A_u_lim_k(6,1) = t122;
        this->A_u_lim_k(7,1) = t36;
        this->A_u_lim_k(0,3) = t272;
        this->A_u_lim_k(1,3) = t22;
        this->A_u_lim_k(2,3) = t52;
        this->A_u_lim_k(3,3) = t20;
        this->A_u_lim_k(4,3) = t52;
        this->A_u_lim_k(5,3) = t20;
        this->A_u_lim_k(6,3) = t272;
        this->A_u_lim_k(7,3) = t22;
        this->A_u_lim_k(0,4) = t67;
        this->A_u_lim_k(1,4) = t17;
        this->A_u_lim_k(2,4) = t192;
        this->A_u_lim_k(3,4) = t28;
        this->A_u_lim_k(4,4) = t192;
        this->A_u_lim_k(5,4) = t28;
        this->A_u_lim_k(6,4) = t67;
        this->A_u_lim_k(7,4) = t17;
      }
     
      void A_state_lim_calc(){
        double t3;
        //     This function was generated by the Symbolic Math Toolbox version 8.7.
        //     29-Jul-2021 16:01:06
        t3 = this->delta_control * this->delta_control / 2.0;
        this->A_state_lim_k(0,3)= t3;
        this->A_state_lim_k(2,3) = -t3;
        this->A_state_lim_k(1,4) = t3;
        this->A_state_lim_k(3,4) = -t3;
}

      void b_nonholonomic_calc(){
        (this->b_nonholonomic_k)<<0,0;//wouldn't be necessary to reassign it for the planar model
      } 
 
      void b_lambda_dyn_calc(double alpha,double b_gamma,double h1,double h2){

        double phi_r=this->q_k(3,0);
        double beta=this->q_k(4,0);
        double phi_r_dot=this->q_k(8,0);
        double beta_dot=this->q_k(9,0);

        double b_b_lambda_dyn_tmp;
        double b_lambda_dyn_tmp;
        double c_b_lambda_dyn_tmp;
        double d_b_lambda_dyn_tmp;
        double e_b_lambda_dyn_tmp;
        double f_b_lambda_dyn_tmp;
        double g_b_lambda_dyn_tmp;
        double t10;
        double t11;
        double t12;
        double t13;
        double t14;
        double t16;
        double t16_tmp;
        double t17;
        double t18;
        double t19_tmp_tmp;
        double t2;
        double t3;
        double t4;
        double t5;
        double t6;
        double t7;
        //     This function was generated by the Symbolic Math Toolbox version 8.7.
        //     29-Jul-2021 16:01:06
        t2 = beta_dot * beta_dot;
        t3 = h1 * h1;
        t4 = h2 * h2;
        
        t5 = phi_r_dot * phi_r_dot;
        t6 = (beta + phi_r) + this->robot_dim(5,0);
        t7 = (b_gamma + phi_r) + this->robot_dim(5,0);
        t10 = std::cos(t6);
        t11 = std::cos(t7);
        t12 = std::sin(t6);
        t13 = std::sin(t7);
        t14 = phi_r + -this->robot_dim(6,0);
        t6 = (phi_r + -alpha) + this->robot_dim(5,0);
        t16_tmp = this->robot_dim(3,0) * this->robot_inertial_par(9,0);
        t16 = t16_tmp * t12;
        t17 = std::cos(t6);
        t18 = std::sin(t6);
        t19_tmp_tmp = this->robot_inertial_par(1,0) * this->robot_inertial_par(7,0);
        t6 = t19_tmp_tmp * t12;
        t7 = this->robot_dim(0,0) * this->robot_inertial_par(8,0);
        t12 = -t16 + -(t6 * 2.0);
        b_lambda_dyn_tmp = this->robot_inertial_par(2,0) * this->robot_inertial_par(6,0);
        b_b_lambda_dyn_tmp =
            this->robot_inertial_par(3,0) * this->robot_inertial_par(8,0);
        c_b_lambda_dyn_tmp = this->robot_dim(4,0) * this->robot_inertial_par(8,0);
        d_b_lambda_dyn_tmp = this->robot_dim(4,0) * this->robot_inertial_par(9,0);
        e_b_lambda_dyn_tmp = this->robot_dim(4,0) * this->robot_inertial_par(7,0);
        f_b_lambda_dyn_tmp = h1 * this->robot_dim(4,0);
        g_b_lambda_dyn_tmp =
            this->robot_inertial_par(4,0) * this->robot_inertial_par(9,0);
        this->b_lambda_dyn_k(0,0) =
            (-t5 *
                (((((((t12 + g_b_lambda_dyn_tmp * std::cos(t14)) + t7 * t18 * 2.0) +
                      b_lambda_dyn_tmp * t18 * 2.0) -
                    b_b_lambda_dyn_tmp * t11 * 2.0) +
                    c_b_lambda_dyn_tmp * t18 * 2.0) +
                  d_b_lambda_dyn_tmp * t18) +
                  e_b_lambda_dyn_tmp * t18 * 2.0) -
            t2 * ((((((t12 + t7 * t3 * t18 * 2.0) +
                      b_lambda_dyn_tmp * t3 * t18 * 2.0) -
                      b_b_lambda_dyn_tmp * t4 * t11 * 2.0) +
                    c_b_lambda_dyn_tmp * t3 * t18 * 2.0) +
                    d_b_lambda_dyn_tmp * t3 * t18) +
                  e_b_lambda_dyn_tmp * t3 * t18 * 2.0)) +
            beta_dot * phi_r_dot *
                (((((((t16 * 2.0 + t6 * 4.0) + h1 * this->robot_dim(0,0) *
                                                  this->robot_inertial_par(8,0) * t18 *
                                                  4.0) +
                    h1 * this->robot_inertial_par(2,0) * this->robot_inertial_par(6,0) *
                        t18 * 4.0) +
                    h2 * this->robot_inertial_par(3,0) * this->robot_inertial_par(8,0) *
                        t11 * 4.0) +
                  f_b_lambda_dyn_tmp * this->robot_inertial_par(8,0) * t18 * 4.0) +
                  f_b_lambda_dyn_tmp * this->robot_inertial_par(9,0) * t18 * 2.0) +
                f_b_lambda_dyn_tmp * this->robot_inertial_par(7,0) * t18 * 4.0);
        t12 = beta_dot * h1;
        f_b_lambda_dyn_tmp = t12 * this->robot_dim(4,0);
        this->b_lambda_dyn_k(1,0) =
            ((((((((((((((((((((((((((((this->g * this->robot_inertial_par(6,0) * -2.0 -
                                        this->g * this->robot_inertial_par(8,0) * 2.0) -
                                      this->g * this->robot_inertial_par(9,0)) -
                                      this->g * this->robot_inertial_par(7,0) * 2.0) -
                                    this->g * this->robot_inertial_par(5,0) * 2.0) +
                                    t7 * t5 * t17 * 2.0) +
                                  b_lambda_dyn_tmp * t5 * t17 * 2.0) +
                                  b_b_lambda_dyn_tmp * t5 * t13 * 2.0) -
                                t16_tmp * t2 * t10) -
                                t16_tmp * t5 * t10) -
                              t19_tmp_tmp * t2 * t10 * 2.0) -
                              t19_tmp_tmp * t5 * t10 * 2.0) +
                            c_b_lambda_dyn_tmp * t5 * t17 * 2.0) +
                            d_b_lambda_dyn_tmp * t5 * t17) +
                          e_b_lambda_dyn_tmp * t5 * t17 * 2.0) -
                          g_b_lambda_dyn_tmp * t5 * std::sin(t14)) -
                        beta_dot * this->robot_dim(3,0) * this->robot_inertial_par(9,0) *
                            phi_r_dot * t10 * 2.0) -
                        beta_dot * this->robot_inertial_par(1,0) *
                            this->robot_inertial_par(7,0) * phi_r_dot * t10 * 4.0) +
                      t7 * t2 * t3 * t17 * 2.0) +
                      b_lambda_dyn_tmp * t2 * t3 * t17 * 2.0) +
                    b_b_lambda_dyn_tmp * t2 * t4 * t13 * 2.0) +
                    c_b_lambda_dyn_tmp * t2 * t3 * t17 * 2.0) +
                  d_b_lambda_dyn_tmp * t2 * t3 * t17) +
                  e_b_lambda_dyn_tmp * t2 * t3 * t17 * 2.0) -
                t12 * this->robot_dim(0,0) * this->robot_inertial_par(8,0) * phi_r_dot *
                    t17 * 4.0) -
                t12 * this->robot_inertial_par(2,0) * this->robot_inertial_par(6,0) *
                    phi_r_dot * t17 * 4.0) +
              beta_dot * h2 * this->robot_inertial_par(3,0) *
                  this->robot_inertial_par(8,0) * phi_r_dot * t13 * 4.0) -
              f_b_lambda_dyn_tmp * this->robot_inertial_par(8,0) * phi_r_dot * t17 *
                  4.0) -
            f_b_lambda_dyn_tmp * this->robot_inertial_par(9,0) * phi_r_dot * t17 *
                2.0) +
            f_b_lambda_dyn_tmp * this->robot_inertial_par(7,0) * phi_r_dot * t17 * -4.0;
      }

      void b_feas_calc(double alpha,double b_gamma,double h1,double h2){

        double phi_r=this->q_k(3,0);
        double beta=this->q_k(4,0);
        double phi_r_dot=this->q_k(8,0);
        double beta_dot=this->q_k(9,0);

        double b_b_feas_tmp;
        double b_feas_tmp;
        double c_b_feas_tmp;
        double d_b_feas_tmp;
        double e_b_feas_tmp;
        double f_b_feas_tmp;
        double g_b_feas_tmp;
        double h_b_feas_tmp;
        double i_b_feas_tmp;
        double j_b_feas_tmp;
        double k_b_feas_tmp;
        double l_b_feas_tmp;
        double m_b_feas_tmp;
        double n_b_feas_tmp;
        double o_b_feas_tmp;
        double p_b_feas_tmp;
        double q_b_feas_tmp;
        double r_b_feas_tmp;
        double s_b_feas_tmp;
        double t15;
        double t16;
        double t17;
        double t20;
        double t22;
        double t23;
        double t4;
        double t5;
        double t6;
        double t7;
        double t8;
        double t9;
        //     This function was generated by the Symbolic Math Toolbox version 8.7.
        //     29-Jul-2021 17:21:23
        t4 = beta_dot * beta_dot;
        t5 = h1 * h1;
        t6 = h2 * h2;
        t7 = phi_r_dot * phi_r_dot;
        t8 = std::cos(alpha + b_gamma);
        t9 = std::sin(alpha + beta);
        t15 = std::cos((beta + this->robot_dim(5,0)) + this->robot_dim(6,0));
        t16 = std::cos((b_gamma + phi_r) + this->robot_dim(5,0));
        t17 = std::sin((beta + phi_r) + this->robot_dim(5,0));
        t20 = std::cos(phi_r + -this->robot_dim(6,0));
        t22 = std::cos((-alpha + this->robot_dim(5,0)) + this->robot_dim(6,0));
        t23 = std::sin((phi_r + -alpha) + this->robot_dim(5,0));
        b_feas_tmp = this->g * this->robot_dim(4,0);
        b_b_feas_tmp = this->robot_dim(3,0) * this->robot_inertial_par(9,0) *
                      this->robot_inertial_par(0,0);
        c_b_feas_tmp = this->robot_inertial_par(1,0) * this->robot_inertial_par(7,0) *
                      this->robot_inertial_par(0,0);
        d_b_feas_tmp = this->robot_dim(0,0) * this->robot_inertial_par(3,0) *
                      this->robot_inertial_par(8,0) * t4;
        e_b_feas_tmp = this->robot_inertial_par(3,0) * this->robot_dim(4,0) *
                      this->robot_inertial_par(8,0) * t4;
        f_b_feas_tmp = this->robot_inertial_par(1,0) * this->robot_dim(4,0) *
                      this->robot_inertial_par(7,0) * t4;
        g_b_feas_tmp = this->robot_dim(0,0) * this->robot_inertial_par(8,0) *
                      this->robot_inertial_par(0,0);
        h_b_feas_tmp = this->robot_inertial_par(2,0) * this->robot_inertial_par(6,0) *
                      this->robot_inertial_par(0,0);
        i_b_feas_tmp = this->robot_inertial_par(3,0) * this->robot_inertial_par(8,0) *
                      this->robot_inertial_par(0,0);
        j_b_feas_tmp = this->robot_dim(4,0) * this->robot_inertial_par(8,0) *
                      this->robot_inertial_par(0,0);
        k_b_feas_tmp = this->robot_dim(4,0) * this->robot_inertial_par(9,0) *
                      this->robot_inertial_par(0,0);
        l_b_feas_tmp = this->robot_dim(4,0) * this->robot_inertial_par(7,0) *
                      this->robot_inertial_par(0,0);
        m_b_feas_tmp = beta_dot * this->robot_dim(3,0);
        n_b_feas_tmp = beta_dot * this->robot_inertial_par(1,0);
        o_b_feas_tmp = beta_dot * h1;
        p_b_feas_tmp = beta_dot * h2;
        q_b_feas_tmp = o_b_feas_tmp * this->robot_dim(0,0);
        r_b_feas_tmp = p_b_feas_tmp * this->robot_inertial_par(3,0);
        s_b_feas_tmp = o_b_feas_tmp * this->robot_dim(4,0);
        this->b_feas_k(0,0)=(((((((((((((((((((((((((((this->g * this->robot_dim(0,0) *
                                              this->robot_inertial_par(8,0) * t23 *
                                              -2.0 -
                                          this->g * this->robot_inertial_par(2,0) *
                                              this->robot_inertial_par(6,0) * t23 *
                                              2.0) +
                                        this->g * this->robot_inertial_par(3,0) *
                                            this->robot_inertial_par(8,0) * t16 *
                                            2.0) +
                                        this->g * this->robot_dim(3,0) *
                                            this->robot_inertial_par(9,0) * t17) +
                                      this->g * this->robot_inertial_par(1,0) *
                                          this->robot_inertial_par(7,0) * t17 * 2.0) -
                                      this->g * this->robot_inertial_par(4,0) *
                                          this->robot_inertial_par(9,0) * t20) -
                                    b_feas_tmp * this->robot_inertial_par(8,0) * t23 *
                                        2.0) -
                                    b_feas_tmp * this->robot_inertial_par(9,0) * t23) -
                                  b_feas_tmp * this->robot_inertial_par(7,0) * t23 *
                                      2.0) -
                                  this->robot_dim(3,0) * this->robot_inertial_par(4,0) *
                                      this->robot_inertial_par(9,0) * t4 * t15) +
                                this->robot_dim(3,0) * this->robot_dim(4,0) *
                                    this->robot_inertial_par(9,0) * t4 * t9) +
                                f_b_feas_tmp * t9 * 2.0) -
                              g_b_feas_tmp * t7 * t23 * 2.0) -
                              h_b_feas_tmp * t7 * t23 * 2.0) +
                            i_b_feas_tmp * t7 * t16 * 2.0) +
                            b_b_feas_tmp * t4 * t17) +
                          b_b_feas_tmp * t7 * t17) +
                          c_b_feas_tmp * t4 * t17 * 2.0) +
                        c_b_feas_tmp * t7 * t17 * 2.0) -
                        this->robot_inertial_par(4,0) * this->robot_inertial_par(9,0) *
                            this->robot_inertial_par(0,0) * t7 * t20) -
                      j_b_feas_tmp * t7 * t23 * 2.0) -
                      k_b_feas_tmp * t7 * t23) -
                    l_b_feas_tmp * t7 * t23 * 2.0) -
                    d_b_feas_tmp * t5 * t8 * 2.0) +
                  d_b_feas_tmp * t6 * t8 * 2.0) -
                  e_b_feas_tmp * t5 * t8 * 2.0) +
                e_b_feas_tmp * t6 * t8 * 2.0) +
                (((((((((((((((((((((-this->robot_dim(3,0) * this->robot_dim(4,0) *
                                        this->robot_inertial_par(9,0) * t4 * t5 * t9 -
                                    f_b_feas_tmp * t5 * t9 * 2.0) +
                                    this->robot_inertial_par(4,0) * this->robot_dim(4,0) *
                                        this->robot_inertial_par(9,0) * t4 * t5 * t22) -
                                  g_b_feas_tmp * t4 * t5 * t23 * 2.0) -
                                  h_b_feas_tmp * t4 * t5 * t23 * 2.0) +
                                i_b_feas_tmp * t4 * t6 * t16 * 2.0) -
                                j_b_feas_tmp * t4 * t5 * t23 * 2.0) -
                              k_b_feas_tmp * t4 * t5 * t23) -
                              l_b_feas_tmp * t4 * t5 * t23 * 2.0) -
                            m_b_feas_tmp * this->robot_inertial_par(4,0) *
                                this->robot_inertial_par(9,0) * phi_r_dot * t15 *
                                2.0) +
                            m_b_feas_tmp * this->robot_dim(4,0) *
                                this->robot_inertial_par(9,0) * phi_r_dot * t9 * 2.0) +
                          n_b_feas_tmp * this->robot_dim(4,0) *
                              this->robot_inertial_par(7,0) * phi_r_dot * t9 * 4.0) +
                          m_b_feas_tmp * this->robot_inertial_par(9,0) * phi_r_dot *
                              this->robot_inertial_par(0,0) * t17 * 2.0) +
                        n_b_feas_tmp * this->robot_inertial_par(7,0) * phi_r_dot *
                            this->robot_inertial_par(0,0) * t17 * 4.0) +
                        q_b_feas_tmp * this->robot_inertial_par(3,0) *
                            this->robot_inertial_par(8,0) * phi_r_dot * t8 * 4.0) +
                      p_b_feas_tmp * this->robot_dim(0,0) *
                          this->robot_inertial_par(3,0) * this->robot_inertial_par(8,0) *
                          phi_r_dot * t8 * 4.0) +
                      o_b_feas_tmp * this->robot_inertial_par(3,0) *
                          this->robot_dim(4,0) * this->robot_inertial_par(8,0) *
                          phi_r_dot * t8 * 4.0) +
                    r_b_feas_tmp * this->robot_dim(4,0) * this->robot_inertial_par(8,0) *
                        phi_r_dot * t8 * 4.0) +
                    o_b_feas_tmp * this->robot_dim(3,0) * this->robot_dim(4,0) *
                        this->robot_inertial_par(9,0) * phi_r_dot * t9 * 2.0) +
                  o_b_feas_tmp * this->robot_inertial_par(1,0) * this->robot_dim(4,0) *
                      this->robot_inertial_par(7,0) * phi_r_dot * t9 * 4.0) -
                  o_b_feas_tmp * this->robot_inertial_par(4,0) * this->robot_dim(4,0) *
                      this->robot_inertial_par(9,0) * phi_r_dot * t22 * 2.0) +
                q_b_feas_tmp * this->robot_inertial_par(8,0) * phi_r_dot *
                    this->robot_inertial_par(0,0) * t23 * 4.0)) +
              ((((o_b_feas_tmp * this->robot_inertial_par(2,0) *
                      this->robot_inertial_par(6,0) * phi_r_dot *
                      this->robot_inertial_par(0,0) * t23 * 4.0 +
                  r_b_feas_tmp * this->robot_inertial_par(8,0) * phi_r_dot *
                      this->robot_inertial_par(0,0) * t16 * 4.0) +
                  s_b_feas_tmp * this->robot_inertial_par(8,0) * phi_r_dot *
                      this->robot_inertial_par(0,0) * t23 * 4.0) +
                s_b_feas_tmp * this->robot_inertial_par(9,0) * phi_r_dot *
                    this->robot_inertial_par(0,0) * t23 * 2.0) +
                s_b_feas_tmp * this->robot_inertial_par(7,0) * phi_r_dot *
                    this->robot_inertial_par(0,0) * t23 * 4.0);
      }
      
      void b_lambda_con_calc(){

        this->b_lambda_con_k << -(this->delta_lambda),0,0;
    
      }
     
      void r_u_check_calc(double alpha,double b_gamma,double h1,double h2){

        double phi_r=this->q_k(3,0);
        double beta=this->q_k(4,0);
        double phi_r_dot=this->q_k(8,0);
        double beta_dot=this->q_k(9,0);

        double b_r_u_check_tmp;
        double c_r_u_check_tmp;
        double d_r_u_check_tmp;
        double e_r_u_check_tmp;
        double f_r_u_check_tmp;
        double g_r_u_check_tmp;
        double h_r_u_check_tmp;
        double i_r_u_check_tmp;
        double j_r_u_check_tmp;
        double k_r_u_check_tmp;
        double r_u_check_tmp;
        double t14;
        double t15;
        double t16;
        double t19;
        double t20;
        double t25_tmp;
        double t25_tmp_tmp;
        double t26_tmp;
        double t26_tmp_tmp;
        double t4;
        double t5;
        double t6;
        double t7;
        double t8;
        double t9;
        //     This function was generated by the Symbolic Math Toolbox version 8.7.
        //     29-Jul-2021 16:01:07
        t4 = beta_dot * beta_dot;
        t5 = h1 * h1;
        t6 = h2 * h2;
        t7 = phi_r_dot * phi_r_dot;
        t8 = std::cos(alpha + b_gamma);
        t9 = std::sin(alpha + beta);
        t14 = std::cos((beta + this->robot_dim(5,0)) + this->robot_dim(6,0));
        t15 = std::cos((b_gamma + phi_r) + this->robot_dim(5,0));
        t16 = std::sin((beta + phi_r) + this->robot_dim(5,0));
        t19 = std::cos((-alpha + this->robot_dim(5,0)) + this->robot_dim(6,0));
        t20 = std::sin((phi_r + -alpha) + this->robot_dim(5,0));
        t25_tmp_tmp =
            this->robot_dim(3,0) * this->robot_dim(4,0) * this->robot_inertial_par(9,0);
        t25_tmp = t25_tmp_tmp * t4;
        t26_tmp_tmp = this->robot_inertial_par(1,0) * this->robot_dim(4,0) *
                      this->robot_inertial_par(7,0);
        t26_tmp = t26_tmp_tmp * t4;
        r_u_check_tmp = this->g * this->robot_dim(4,0);
        b_r_u_check_tmp = this->robot_dim(0,0) * this->robot_inertial_par(3,0) *
                          this->robot_inertial_par(8,0) * t4;
        c_r_u_check_tmp = this->robot_inertial_par(3,0) * this->robot_dim(4,0) *
                          this->robot_inertial_par(8,0) * t4;
        d_r_u_check_tmp = beta_dot * this->robot_dim(3,0);
        e_r_u_check_tmp = beta_dot * h1;
        f_r_u_check_tmp = beta_dot * h2;
        g_r_u_check_tmp =
            ((-(this->g * this->robot_dim(3,0) * this->robot_inertial_par(9,0) * t16) +
              -(this->g * this->robot_inertial_par(1,0) * this->robot_inertial_par(7,0) *
                t16 * 2.0)) +
            t25_tmp * t5 * t9) +
            t26_tmp * t5 * t9 * 2.0;
        h_r_u_check_tmp = this->robot_dim(3,0) * this->robot_inertial_par(4,0) *
                          this->robot_inertial_par(9,0);
        this->r_u_check_k(0,0) =
            (((((((((((((((((((((((g_r_u_check_tmp + this->g * this->robot_dim(0,0) *
                                                        this->robot_inertial_par(8,0) *
                                                        t20 * 2.0) +
                                  this->g * this->robot_inertial_par(2,0) *
                                      this->robot_inertial_par(6,0) * t20 * 2.0) -
                                this->g * this->robot_inertial_par(3,0) *
                                    this->robot_inertial_par(8,0) * t15 * 2.0) +
                                r_u_check_tmp * this->robot_inertial_par(8,0) * t20 *
                                    2.0) +
                              r_u_check_tmp * this->robot_inertial_par(9,0) * t20) +
                              r_u_check_tmp * this->robot_inertial_par(7,0) * t20 *
                                  2.0) +
                            this->g * this->robot_inertial_par(4,0) *
                                this->robot_inertial_par(9,0) *
                                std::cos(phi_r - this->robot_dim(6,0))) +
                            h_r_u_check_tmp * t4 * t14) -
                          t25_tmp * t9) -
                          t26_tmp * t9 * 2.0) +
                        b_r_u_check_tmp * t5 * t8 * 2.0) -
                        b_r_u_check_tmp * t6 * t8 * 2.0) +
                      c_r_u_check_tmp * t5 * t8 * 2.0) -
                      c_r_u_check_tmp * t6 * t8 * 2.0) -
                    this->robot_inertial_par(4,0) * this->robot_dim(4,0) *
                        this->robot_inertial_par(9,0) * t4 * t5 * t19) +
                    d_r_u_check_tmp * this->robot_inertial_par(4,0) *
                        this->robot_inertial_par(9,0) * phi_r_dot * t14 * 2.0) -
                  d_r_u_check_tmp * this->robot_dim(4,0) *
                      this->robot_inertial_par(9,0) * phi_r_dot * t9 * 2.0) -
                  beta_dot * this->robot_inertial_par(1,0) * this->robot_dim(4,0) *
                      this->robot_inertial_par(7,0) * phi_r_dot * t9 * 4.0) -
                e_r_u_check_tmp * this->robot_dim(0,0) * this->robot_inertial_par(3,0) *
                    this->robot_inertial_par(8,0) * phi_r_dot * t8 * 4.0) -
                f_r_u_check_tmp * this->robot_dim(0,0) * this->robot_inertial_par(3,0) *
                    this->robot_inertial_par(8,0) * phi_r_dot * t8 * 4.0) -
              e_r_u_check_tmp * this->robot_inertial_par(3,0) * this->robot_dim(4,0) *
                  this->robot_inertial_par(8,0) * phi_r_dot * t8 * 4.0) -
              f_r_u_check_tmp * this->robot_inertial_par(3,0) * this->robot_dim(4,0) *
                  this->robot_inertial_par(8,0) * phi_r_dot * t8 * 4.0) -
            e_r_u_check_tmp * this->robot_dim(3,0) * this->robot_dim(4,0) *
                this->robot_inertial_par(9,0) * phi_r_dot * t9 * 2.0) +
            (e_r_u_check_tmp * this->robot_inertial_par(1,0) * this->robot_dim(4,0) *
                this->robot_inertial_par(7,0) * phi_r_dot * t9 * -4.0 +
            e_r_u_check_tmp * this->robot_inertial_par(4,0) * this->robot_dim(4,0) *
                this->robot_inertial_par(9,0) * phi_r_dot * t19 * 2.0);
        r_u_check_tmp = alpha * h1;
        b_r_u_check_tmp = this->robot_knee_el_par(5,0) * h1;
        c_r_u_check_tmp = this->robot_knee_el_par(6,0) * h2;
        d_r_u_check_tmp = b_gamma * h2;
        e_r_u_check_tmp = this->g * h1;
        f_r_u_check_tmp = e_r_u_check_tmp * this->robot_dim(4,0);
        t16 = h1 * this->robot_dim(3,0) * this->robot_dim(4,0) *
              this->robot_inertial_par(9,0);
        t25_tmp = h1 * this->robot_inertial_par(1,0) * this->robot_dim(4,0) *
                  this->robot_inertial_par(7,0);
        t26_tmp = h1 * this->robot_dim(0,0) * this->robot_inertial_par(3,0) *
                  this->robot_inertial_par(8,0);
        i_r_u_check_tmp = h2 * this->robot_dim(0,0) * this->robot_inertial_par(3,0) *
                          this->robot_inertial_par(8,0);
        j_r_u_check_tmp = h1 * this->robot_inertial_par(3,0) * this->robot_dim(4,0) *
                          this->robot_inertial_par(8,0);
        k_r_u_check_tmp = h2 * this->robot_inertial_par(3,0) * this->robot_dim(4,0) *
                          this->robot_inertial_par(8,0);
        this->r_u_check_k(1,0) =
            (((((((((((((((((((((((((((((((g_r_u_check_tmp -
                                          alpha * this->robot_knee_el_par(1,0) * 2.0) +
                                          this->robot_knee_el_par(1,0) *
                                              this->robot_knee_el_par(5,0) * 2.0) +
                                        this->robot_knee_el_par(1,0) *
                                            this->robot_knee_el_par(4,0) * 2.0) +
                                        this->robot_knee_el_par(3,0) *
                                            this->robot_knee_el_par(4,0) * 2.0) -
                                      beta * this->robot_knee_el_par(1,0) * 2.0) -
                                      beta * this->robot_knee_el_par(3,0) * 2.0) -
                                    r_u_check_tmp * this->robot_knee_el_par(0,0) *
                                        2.0) +
                                    b_r_u_check_tmp * this->robot_knee_el_par(0,0) *
                                        2.0) -
                                  r_u_check_tmp * this->robot_knee_el_par(1,0) * 2.0) -
                                  alpha * h2 * this->robot_knee_el_par(0,0) * 2.0) +
                                b_r_u_check_tmp * this->robot_knee_el_par(1,0) * 2.0) +
                                this->robot_knee_el_par(5,0) * h2 *
                                    this->robot_knee_el_par(0,0) * 2.0) +
                              this->robot_knee_el_par(4,0) * h1 *
                                  this->robot_knee_el_par(1,0) * 2.0) -
                              beta * h1 * this->robot_knee_el_par(1,0) * 2.0) +
                            this->robot_knee_el_par(6,0) * h1 *
                                this->robot_knee_el_par(0,0) * 2.0) +
                            c_r_u_check_tmp * this->robot_knee_el_par(0,0) * 2.0) +
                          c_r_u_check_tmp * this->robot_knee_el_par(2,0) * 2.0) -
                          b_gamma * h1 * this->robot_knee_el_par(0,0) * 2.0) -
                        d_r_u_check_tmp * this->robot_knee_el_par(0,0) * 2.0) -
                        d_r_u_check_tmp * this->robot_knee_el_par(2,0) * 2.0) -
                      e_r_u_check_tmp * this->robot_dim(0,0) *
                          this->robot_inertial_par(8,0) * t20 * 2.0) -
                      e_r_u_check_tmp * this->robot_inertial_par(2,0) *
                          this->robot_inertial_par(6,0) * t20 * 2.0) -
                    this->g * h2 * this->robot_inertial_par(3,0) *
                        this->robot_inertial_par(8,0) * t15 * 2.0) -
                    f_r_u_check_tmp * this->robot_inertial_par(8,0) * t20 * 2.0) -
                  f_r_u_check_tmp * this->robot_inertial_par(9,0) * t20) -
                  f_r_u_check_tmp * this->robot_inertial_par(7,0) * t20 * 2.0) -
                h_r_u_check_tmp * t7 * t14) +
                t25_tmp_tmp * t7 * t9) +
              t26_tmp_tmp * t7 * t9 * 2.0) +
              t26_tmp * t7 * t8 * 2.0) +
            i_r_u_check_tmp * t7 * t8 * 2.0) +
            ((((((((((j_r_u_check_tmp * t7 * t8 * 2.0 +
                      k_r_u_check_tmp * t7 * t8 * 2.0) +
                    t16 * t4 * t9) +
                    t16 * t7 * t9) +
                  t25_tmp * t4 * t9 * 2.0) +
                  t25_tmp * t7 * t9 * 2.0) -
                h1 * this->robot_inertial_par(4,0) * this->robot_dim(4,0) *
                    this->robot_inertial_par(9,0) * t7 * t19) +
                t26_tmp * t4 * t6 * t8 * 2.0) +
              i_r_u_check_tmp * t4 * t5 * t8 * 2.0) +
              j_r_u_check_tmp * t4 * t6 * t8 * 2.0) +
            k_r_u_check_tmp * t4 * t5 * t8 * 2.0);
      }
    
      void b_u_lim_calc(double alpha,double b_gamma,double h1,double h2){

        double phi_r=this->q_k(3,0);
        double beta=this->q_k(4,0);
        double phi_r_dot=this->q_k(8,0);
        double beta_dot=this->q_k(9,0);

        double t10;
        double t101;
        double t101_tmp;
        double t102;
        double t102_tmp;
        double t105;
        double t106;
        double t108;
        double t11;
        double t111;
        double t12;
        double t123;
        double t125;
        double t13;
        double t133;
        double t14;
        double t15;
        double t151;
        double t152;
        double t153;
        double t154;
        double t155;
        double t156;
        double t157;
        double t159;
        double t16;
        double t160;
        double t161;
        double t162;
        double t163;
        double t164;
        double t165;
        double t166;
        double t167;
        double t17;
        double t172;
        double t18;
        double t182;
        double t183;
        double t19;
        double t194;
        double t20;
        double t202;
        double t21;
        double t22;
        double t23;
        double t24;
        double t25;
        double t26;
        double t27;
        double t28;
        double t29;
        double t37;
        double t39;
        double t4;
        double t5;
        double t6;
        double t64;
        double t65;
        double t66;
        double t67;
        double t68;
        double t69;
        double t7;
        double t70;
        double t73;
        double t74;
        double t75;
        double t76;
        double t77;
        double t78;
        double t79;
        double t8;
        double t82;
        double t83;
        double t87;
        double t88;
        double t89;
        double t9;
        double t90;
        double t91;
        double t92;
        double t93;
        double t94;
        double t95;
        double t96;
        double t97;
        double t97_tmp;
        double t98;
        double t98_tmp;
        //     This function was generated by the Symbolic Math Toolbox version 8.7.
        //     29-Jul-2021 16:01:07
        t4 = alpha * this->robot_knee_el_par(1,0);
        t5 = this->robot_knee_el_par(1,0) * this->robot_knee_el_par(5,0);
        t6 = this->robot_knee_el_par(1,0) * this->robot_knee_el_par(4,0);
        t7 = this->robot_knee_el_par(3,0) * this->robot_knee_el_par(4,0);
        t8 = beta * this->robot_knee_el_par(1,0);
        t9 = beta * this->robot_knee_el_par(3,0);
        t10 = beta_dot * beta_dot;
        t11 = h1 * h1;
        t12 = h2 * h2;
        t13 = phi_r_dot * phi_r_dot;
        t15 = alpha * h1 * this->robot_knee_el_par(0,0);
        t16 = this->robot_knee_el_par(5,0) * h1 * this->robot_knee_el_par(0,0);
        t18 = alpha * h2 * this->robot_knee_el_par(0,0);
        t20 = this->robot_knee_el_par(5,0) * h2 * this->robot_knee_el_par(0,0);
        t24 = this->robot_knee_el_par(6,0) * h1 * this->robot_knee_el_par(0,0);
        t39 = this->robot_knee_el_par(6,0) * h2;
        t25 = t39 * this->robot_knee_el_par(0,0);
        t26 = t39 * this->robot_knee_el_par(2,0);
        t27 = b_gamma * h1 * this->robot_knee_el_par(0,0);
        t39 = b_gamma * h2;
        t28 = t39 * this->robot_knee_el_par(0,0);
        t29 = t39 * this->robot_knee_el_par(2,0);
        t14 = std::cos(alpha + b_gamma);
        t17 = h1 * t4;
        t19 = h1 * t5;
        t21 = h1 * t6;
        t22 = h1 * t8;
        t23 = std::sin(alpha + beta);
        t37 = std::cos((beta + this->robot_dim(5,0)) + this->robot_dim(6,0));
        t39 = std::sin((beta + phi_r) + this->robot_dim(5,0));
        t64 = std::cos((-alpha + this->robot_dim(5,0)) + this->robot_dim(6,0));
        t65 = std::sin((phi_r + -alpha) + this->robot_dim(5,0));
        t66 = this->g * this->robot_inertial_par(3,0) * this->robot_inertial_par(8,0) *
              std::cos((b_gamma + phi_r) + this->robot_dim(5,0));
        t67 = this->g * this->robot_dim(3,0) * this->robot_inertial_par(9,0) * t39;
        t68 =
            this->g * this->robot_inertial_par(1,0) * this->robot_inertial_par(7,0) * t39;
        t194 = beta_dot * this->robot_dim(3,0);
        t69 =
            t194 * this->robot_dim(4,0) * this->robot_inertial_par(9,0) * phi_r_dot * t23;
        t202 = this->robot_inertial_par(1,0) * this->robot_dim(4,0) *
              this->robot_inertial_par(7,0);
        t73 = t202 * t13 * t23;
        t39 = beta_dot * this->robot_inertial_par(1,0) * this->robot_dim(4,0) *
              this->robot_inertial_par(7,0) * phi_r_dot * t23;
        t75 = t39 * 2.0;
        t76 = t39 * 4.0;
        t78 = t194 * this->robot_inertial_par(4,0) * this->robot_inertial_par(9,0) *
              phi_r_dot * t37;
        t194 = this->robot_dim(3,0) * this->robot_dim(4,0) * this->robot_inertial_par(9,0);
        t82 = t194 * t10 * t23;
        t83 = t202 * t10 * t23;
        t93 = h1 * this->robot_dim(0,0) * this->robot_inertial_par(3,0) *
              this->robot_inertial_par(8,0) * t13 * t14;
        t94 = h2 * this->robot_dim(0,0) * this->robot_inertial_par(3,0) *
              this->robot_inertial_par(8,0) * t13 * t14;
        t95 = h1 * this->robot_inertial_par(3,0) * this->robot_dim(4,0) *
              this->robot_inertial_par(8,0) * t13 * t14;
        t96 = h2 * this->robot_inertial_par(3,0) * this->robot_dim(4,0) *
              this->robot_inertial_par(8,0) * t13 * t14;
        t125 = beta_dot * h1;
        t97_tmp = t125 * this->robot_dim(0,0) * this->robot_inertial_par(3,0) *
                  this->robot_inertial_par(8,0) * phi_r_dot * t14;
        t97 = t97_tmp * 2.0;
        t39 = beta_dot * h2;
        t98_tmp = t39 * this->robot_dim(0,0) * this->robot_inertial_par(3,0) *
                  this->robot_inertial_par(8,0) * phi_r_dot * t14;
        t98 = t98_tmp * 2.0;
        t101_tmp = t125 * this->robot_inertial_par(3,0) * this->robot_dim(4,0) *
                  this->robot_inertial_par(8,0) * phi_r_dot * t14;
        t101 = t101_tmp * 2.0;
        t102_tmp = t39 * this->robot_inertial_par(3,0) * this->robot_dim(4,0) *
                  this->robot_inertial_par(8,0) * phi_r_dot * t14;
        t102 = t102_tmp * 2.0;
        t202 = this->robot_dim(3,0) * this->robot_inertial_par(4,0) *
              this->robot_inertial_par(9,0);
        t106 = t202 * t10 * t37;
        t123 = t125 * this->robot_inertial_par(1,0) * this->robot_dim(4,0) *
              this->robot_inertial_par(7,0) * phi_r_dot * t23 * -2.0;
        t133 = t194 * t13 * t23 / 2.0;
        t39 = this->robot_dim(0,0) * this->robot_inertial_par(3,0) *
              this->robot_inertial_par(8,0) * t10;
        t151 = t39 * t11 * t14;
        t152 = t39 * t12 * t14;
        t39 = this->robot_inertial_par(3,0) * this->robot_dim(4,0) *
              this->robot_inertial_par(8,0) * t10;
        t153 = t39 * t11 * t14;
        t154 = t39 * t12 * t14;
        t172 = t202 * t13 * t37 / 2.0;
        t183 = h1 * this->robot_dim(3,0) * this->robot_dim(4,0) *
              this->robot_inertial_par(9,0) * t13 * t23 * -0.5;
        t70 = h2 * t66;
        t74 = t69 * 2.0;
        t77 = h1 * t69;
        t79 = this->g * this->robot_inertial_par(4,0) * this->robot_inertial_par(9,0) *
              std::cos(phi_r + -this->robot_dim(6,0));
        t87 = this->g * this->robot_dim(0,0) * this->robot_inertial_par(8,0) * t65;
        t88 =
            this->g * this->robot_inertial_par(2,0) * this->robot_inertial_par(6,0) * t65;
        t39 = this->g * this->robot_dim(4,0);
        t89 = t39 * this->robot_inertial_par(8,0) * t65;
        t90 = t39 * this->robot_inertial_par(9,0) * t65;
        t91 = t39 * this->robot_inertial_par(7,0) * t65;
        t92 = h1 * t83;
        t105 = h1 * t73;
        t108 = h1 * t75;
        t111 = t67 / 2.0;
        t155 = t11 * t82;
        t156 = t11 * t83;
        t157 = t82 / 2.0;
        t160 = h1 * t152;
        t161 = h2 * t151;
        t162 = h1 * t154;
        t163 = h2 * t153;
        t164 = t125 * this->robot_inertial_par(4,0) * this->robot_dim(4,0) *
              this->robot_inertial_par(9,0) * phi_r_dot * t64;
        t166 = h1 * t133;
        t167 = t106 / 2.0;
        t182 = h1 * t82 * -0.5;
        t39 = this->robot_inertial_par(4,0) * this->robot_dim(4,0) *
              this->robot_inertial_par(9,0) * t10 * t11 * t64;
        t23 = h1 * this->robot_inertial_par(4,0) * this->robot_dim(4,0) *
              this->robot_inertial_par(9,0) * t13 * t64 / 2.0;
        t65 = h1 * t87;
        t125 = h1 * t88;
        t64 = h1 * t89;
        t10 = h1 * t91;
        t13 = t79 / 2.0;
        t159 = t90 / 2.0;
        t165 = h1 * t157;
        t194 = t11 * -t83;
        t202 = h1 * t90 * -0.5;
        t14 = t155 / 2.0;
        t37 = t39 / 2.0;
        t12 = h1 * t159;
        t39 =
            h1 *
            ((((((((((((((((((((((((((((t67 + t66 * 2.0) + t68 * 2.0) + t74) + t76) +
                                    t82) +
                                  t97_tmp * 4.0) +
                                  t98_tmp * 4.0) +
                                t101_tmp * 4.0) +
                                t102_tmp * 4.0) +
                              h1 * t74) +
                              h1 * t76) +
                            t83 * 2.0) +
                            -(t78 * 2.0)) +
                          -t79) +
                          -(t87 * 2.0)) +
                        -(t88 * 2.0)) +
                        -(t89 * 2.0)) +
                      -t90) +
                      -(t91 * 2.0)) +
                    -t106) +
                    t152 * 2.0) +
                  t154 * 2.0) +
                  -(t151 * 2.0)) +
                -(t153 * 2.0)) +
                -t155) +
              t156 * -2.0) +
              -(t164 * 2.0)) +
            t39) /
            2.0;
        this->b_u_lim_k(0,0) =
            ((((((((((((((((((((((((((((this->u_lim(1,0) + t78) + -t66) + -t68) + -t69) +
                                    -t75) +
                                  t87) +
                                  t88) +
                                t89) +
                                t91) +
                              -t97) +
                              -t98) +
                            -t101) +
                            -t102) +
                          -t77) +
                          t123) +
                        -t111) +
                        -t83) +
                      t13) +
                      t151) +
                    t153) +
                    t156) +
                  t159) +
                  t164) +
                t167) +
                -t157) +
              -t152) +
              -t154) +
            t14) +
            -t37;
        this->b_u_lim_k(1,0) =
            ((((((((((((((((((((((((((((((((((((((((((((((this->u_lim(3,0) + t5) + t6) + t7) +
                                                      t16) +
                                                      t19) +
                                                    t20) +
                                                    t21) +
                                                  t24) +
                                                  t25) +
                                                t26) +
                                                -t4) +
                                              -t8) +
                                              -t9) +
                                            -t15) +
                                            -t17) +
                                          -t18) +
                                          -t22) +
                                        -t27) +
                                        -t28) +
                                      -t29) +
                                      t73) +
                                    -t68) +
                                    -t70) +
                                  t92) +
                                  t93) +
                                t94) +
                                t95) +
                              t96) +
                              t105) +
                            -t111) +
                            t133) +
                          t156) +
                          t160) +
                        t161) +
                        t162) +
                      t163) +
                      t165) +
                    t166) +
                    -t65) +
                  -t125) +
                  -t64) +
                -t10) +
                -t172) +
              t202) +
              t14) +
            -t23) +
            -t39;
        this->b_u_lim_k(2,0) =
            ((((((((((((((((((((((((((((this->u_lim(1,0) + t66) + t68) + t69) + t75) + t77) +
                                  t83) +
                                  t97) +
                                t98) +
                                t101) +
                              t102) +
                              t108) +
                            t111) +
                            -t78) +
                          -t87) +
                          -t88) +
                        -t89) +
                        -t91) +
                      t152) +
                      t154) +
                    t157) +
                    -t13) +
                  -t159) +
                  -t167) +
                -t151) +
                -t153) +
              t194) +
              -t164) +
            -t14) +
            t37;
        this->b_u_lim_k(3,0) =
            ((((((((((((((((((((((((((((((((((((((((((((((this->u_lim(3,0) + t4) + t8) + t9) +
                                                      t15) +
                                                      t17) +
                                                    t18) +
                                                    t22) +
                                                  t27) +
                                                  t28) +
                                                t29) +
                                                -t5) +
                                              -t6) +
                                              -t7) +
                                            -t16) +
                                            -t19) +
                                          -t20) +
                                          -t21) +
                                        -t24) +
                                        -t25) +
                                      -t26) +
                                      t68) +
                                    t70) +
                                    t111) +
                                  -t73) +
                                  t65) +
                                t125) +
                                t64) +
                              t10) +
                              -t92) +
                            -t93) +
                            -t94) +
                          -t95) +
                          -t96) +
                        -t105) +
                        -t133) +
                      t172) +
                      t182) +
                    t183) +
                    t12) +
                  t194) +
                  -t160) +
                -t161) +
                -t162) +
              -t163) +
              -t14) +
            t23) +
            t39;
        this->b_u_lim_k(4,0) =
            ((((((((((((((((((((((((((((-this->u_lim(0,0) + t66) + t68) + t69) + t75) + t77) +
                                  t83) +
                                  t97) +
                                t98) +
                                t101) +
                              t102) +
                              t108) +
                            t111) +
                            -t78) +
                          -t87) +
                          -t88) +
                        -t89) +
                        -t91) +
                      t152) +
                      t154) +
                    t157) +
                    -t13) +
                  -t159) +
                  -t167) +
                -t151) +
                -t153) +
              t194) +
              -t164) +
            -t14) +
            t37;
        this->b_u_lim_k(5,0) =
            ((((((((((((((((((((((((((((((((((((((((((((((t4 + t8) + t9) + t15) +
                                                      t17) +
                                                      t18) +
                                                    t22) +
                                                    t27) +
                                                  t28) +
                                                  t29) +
                                                -this->u_lim(2,0)) +
                                                -t5) +
                                              -t6) +
                                              -t7) +
                                            -t16) +
                                            -t19) +
                                          -t20) +
                                          -t21) +
                                        -t24) +
                                        -t25) +
                                      -t26) +
                                      t68) +
                                    t70) +
                                    t111) +
                                  -t73) +
                                  t65) +
                                t125) +
                                t64) +
                              t10) +
                              -t92) +
                            -t93) +
                            -t94) +
                          -t95) +
                          -t96) +
                        -t105) +
                        -t133) +
                      t172) +
                      t182) +
                    t183) +
                    t12) +
                  t194) +
                  -t160) +
                -t161) +
                -t162) +
              -t163) +
              -t14) +
            t23) +
            t39;
        this->b_u_lim_k(6,0) =
            ((((((((((((((((((((((((((((-this->u_lim(0,0) + t78) + -t66) + -t68) + -t69) +
                                    -t75) +
                                  t87) +
                                  t88) +
                                t89) +
                                t91) +
                              -t97) +
                              -t98) +
                            -t101) +
                            -t102) +
                          -t77) +
                          t123) +
                        -t111) +
                        -t83) +
                      t13) +
                      t151) +
                    t153) +
                    t156) +
                  t159) +
                  t164) +
                t167) +
                -t157) +
              -t152) +
              -t154) +
            t14) +
            -t37;
        this->b_u_lim_k(7,0) =
            ((((((((((((((((((((((((((((((((((((((((((((((t5 + t6) + t7) + t16) +
                                                      t19) +
                                                      t20) +
                                                    t21) +
                                                    t24) +
                                                  t25) +
                                                  t26) +
                                                -this->u_lim(2,0)) +
                                                -t4) +
                                              -t8) +
                                              -t9) +
                                            -t15) +
                                            -t17) +
                                          -t18) +
                                          -t22) +
                                        -t27) +
                                        -t28) +
                                      -t29) +
                                      t73) +
                                    -t68) +
                                    -t70) +
                                  t92) +
                                  t93) +
                                t94) +
                                t95) +
                              t96) +
                              t105) +
                            -t111) +
                            t133) +
                          t156) +
                          t160) +
                        t161) +
                        t162) +
                      t163) +
                      t165) +
                    t166) +
                    -t65) +
                  -t125) +
                  -t64) +
                -t10) +
                -t172) +
              t202) +
              t14) +
            -t23) +
            -t39;
      }

      void b_state_lim_calc(){

        double phi_r=this->q_k(3,0);
        double beta=this->q_k(4,0);
        double phi_r_dot=this->q_k(8,0);
        double beta_dot=this->q_k(9,0);

        double t2;
        double t3;
        t2 = beta_dot * (this->delta_control);
        t3 = (this->delta_control) * phi_r_dot;
        this->b_state_lim_k(0,0) = (-phi_r + this->phi_r_lim(1,0)) - t3;
        this->b_state_lim_k(1,0) = (this->beta_lim(1,0) - beta) - t2;
        this->b_state_lim_k(2,0) = (phi_r + this->phi_r_lim(0,0)) + t3;
        this->b_state_lim_k(3,0) = (this->beta_lim(0,0) + beta) + t2;
      }
       
      void p_cm_tot_rel_wheel_calc(double alpha,double b_gamma){

        double phi_r=this->q_k(3,0);
        double beta=this->q_k(4,0);

        double b_p_cm_tot_rel_wheel_tmp;
        double c_p_cm_tot_rel_wheel_tmp;
        double d_p_cm_tot_rel_wheel_tmp;
        double e_p_cm_tot_rel_wheel_tmp;
        double p_cm_tot_rel_wheel_tmp;
        double t10;
        double t11;
        double t12;
        double t14;
        double t15;
        double t17;
        double t2;
        double t3;
        double t4;
        double t6;
        double t7;
        //     This function was generated by the Symbolic Math Toolbox version 8.7.
        //     11-Aug-2021 12:49:02
        t2 = this->robot_inertial_par(6) * 2.0;
        t3 = this->robot_inertial_par(8) * 2.0;
        t4 = this->robot_inertial_par(7) * 2.0;
        t6 = (beta + phi_r) + this->robot_dim(5);
        t7 = (b_gamma + phi_r) + this->robot_dim(5);
        t10 = std::cos(t6);
        t11 = std::sin(t6);
        t12 = phi_r + -this->robot_dim(6);
        t6 = (phi_r + -alpha) + this->robot_dim(5);
        t14 = std::cos(t6);
        t15 = std::sin(t6);
        t17 = 1.0 / ((((this->robot_inertial_par(9) + t2) + t3) + t4) +
                    this->robot_inertial_par(5) * 2.0);
        p_cm_tot_rel_wheel_tmp = this->robot_dim(3) * this->robot_inertial_par(9);
        b_p_cm_tot_rel_wheel_tmp =
            this->robot_inertial_par(1) * this->robot_inertial_par(7);
        c_p_cm_tot_rel_wheel_tmp = this->robot_dim(4) * this->robot_inertial_par(9);
        d_p_cm_tot_rel_wheel_tmp = this->robot_dim(0) * t3;
        t2 *= this->robot_inertial_par(2);
        e_p_cm_tot_rel_wheel_tmp = this->robot_dim(4) * t3;
        t6 = this->robot_dim(4) * t4;
        this->p_cm_tot_rel_wheel_k(0) =
            -t17 * ((((((((this->robot_inertial_par(3) * this->robot_inertial_par(8) *
                              std::cos(t7) * -2.0 +
                          this->robot_inertial_par(4) * this->robot_inertial_par(9) *
                              std::cos(t12)) -
                          p_cm_tot_rel_wheel_tmp * t11) -
                        b_p_cm_tot_rel_wheel_tmp * t11 * 2.0) +
                        c_p_cm_tot_rel_wheel_tmp * t15) +
                      d_p_cm_tot_rel_wheel_tmp * t15) +
                      t2 * t15) +
                    e_p_cm_tot_rel_wheel_tmp * t15) +
                    t6 * t15);
        this->p_cm_tot_rel_wheel_k(1) =
            t17 * ((((((((-this->robot_inertial_par(4) * this->robot_inertial_par(9) *
                              std::sin(t12) +
                          this->robot_inertial_par(3) * t3 * std::sin(t7)) -
                        p_cm_tot_rel_wheel_tmp * t10) -
                        b_p_cm_tot_rel_wheel_tmp * t10 * 2.0) +
                      c_p_cm_tot_rel_wheel_tmp * t14) +
                      d_p_cm_tot_rel_wheel_tmp * t14) +
                    t2 * t14) +
                    e_p_cm_tot_rel_wheel_tmp * t14) +
                  t6 * t14);

      }
      
      void v_cm_tot_rel_wheel_calc(double alpha,double b_gamma,double h1,double h2){

        double phi_r=this->q_k(3,0);
        double beta=this->q_k(4,0);
        double phi_r_dot=this->q_k(8,0);
        double beta_dot=this->q_k(9,0);

        double b_v_cm_tot_rel_wheel_tmp;
        double c_v_cm_tot_rel_wheel_tmp;
        double d_v_cm_tot_rel_wheel_tmp;
        double e_v_cm_tot_rel_wheel_tmp;
        double f_v_cm_tot_rel_wheel_tmp;
        double g_v_cm_tot_rel_wheel_tmp;
        double h_v_cm_tot_rel_wheel_tmp;
        double i_v_cm_tot_rel_wheel_tmp;
        double j_v_cm_tot_rel_wheel_tmp;
        double k_v_cm_tot_rel_wheel_tmp;
        double l_v_cm_tot_rel_wheel_tmp;
        double m_v_cm_tot_rel_wheel_tmp;
        double n_v_cm_tot_rel_wheel_tmp;
        double o_v_cm_tot_rel_wheel_tmp;
        double t10;
        double t11;
        double t12;
        double t13;
        double t14;
        double t16;
        double t19;
        double t2;
        double t3;
        double t4;
        double t6;
        double t7;
        double v_cm_tot_rel_wheel_tmp;
        //     This function was generated by the Symbolic Math Toolbox version 8.7.
        //     11-Aug-2021 12:49:03
        t2 = this->robot_inertial_par(6) * 2.0;
        t3 = this->robot_inertial_par(8) * 2.0;
        t4 = this->robot_inertial_par(7) * 2.0;
        t6 = (beta + phi_r) + this->robot_dim(5);
        t7 = (b_gamma + phi_r) + this->robot_dim(5);
        t10 = std::cos(t6);
        t11 = std::cos(t7);
        t12 = std::sin(t6);
        t13 = std::sin(t7);
        t14 = phi_r + -this->robot_dim(6);
        t6 = (phi_r + -alpha) + this->robot_dim(5);
        t16 = std::cos(t6);
        t7 = std::sin(t6);
        t19 = 1.0 / ((((this->robot_inertial_par(9) + t2) + t3) + t4) +
                    this->robot_inertial_par(5) * 2.0);
        v_cm_tot_rel_wheel_tmp = beta_dot * h1;
        b_v_cm_tot_rel_wheel_tmp = v_cm_tot_rel_wheel_tmp * this->robot_dim(4);
        c_v_cm_tot_rel_wheel_tmp =
            beta_dot * this->robot_dim(3) * this->robot_inertial_par(9);
        d_v_cm_tot_rel_wheel_tmp = beta_dot * this->robot_inertial_par(1) * t4;
        e_v_cm_tot_rel_wheel_tmp =
            this->robot_dim[0] * this->robot_inertial_par(8) * phi_r_dot;
        f_v_cm_tot_rel_wheel_tmp =
            this->robot_inertial_par(2) * this->robot_inertial_par(6) * phi_r_dot;
        g_v_cm_tot_rel_wheel_tmp =
            this->robot_dim(3) * this->robot_inertial_par(9) * phi_r_dot;
        h_v_cm_tot_rel_wheel_tmp =
            this->robot_dim(4) * this->robot_inertial_par(8) * phi_r_dot;
        i_v_cm_tot_rel_wheel_tmp =
            this->robot_dim(4) * this->robot_inertial_par(9) * phi_r_dot;
        j_v_cm_tot_rel_wheel_tmp =
            this->robot_dim(4) * this->robot_inertial_par(7) * phi_r_dot;
        k_v_cm_tot_rel_wheel_tmp = this->robot_inertial_par(1) * phi_r_dot * t4;
        l_v_cm_tot_rel_wheel_tmp =
            this->robot_inertial_par(4) * this->robot_inertial_par(9) * phi_r_dot;
        m_v_cm_tot_rel_wheel_tmp =
            b_v_cm_tot_rel_wheel_tmp * this->robot_inertial_par(9);
        n_v_cm_tot_rel_wheel_tmp = v_cm_tot_rel_wheel_tmp * this->robot_dim[0] * t3;
        o_v_cm_tot_rel_wheel_tmp = beta_dot * h2 * this->robot_inertial_par(3);
        v_cm_tot_rel_wheel_tmp =
            v_cm_tot_rel_wheel_tmp * this->robot_inertial_par(2) * t2;
        t6 = b_v_cm_tot_rel_wheel_tmp * t3;
        b_v_cm_tot_rel_wheel_tmp *= t4;
        this->v_cm_tot_rel_wheel_k(0) =
            t19 *
            ((((((((((((((((c_v_cm_tot_rel_wheel_tmp * t10 +
                            d_v_cm_tot_rel_wheel_tmp * t10) -
                          e_v_cm_tot_rel_wheel_tmp * t16 * 2.0) -
                          f_v_cm_tot_rel_wheel_tmp * t16 * 2.0) -
                        this->robot_inertial_par(3) * this->robot_inertial_par(8) *
                            phi_r_dot * t13 * 2.0) +
                        g_v_cm_tot_rel_wheel_tmp * t10) -
                      h_v_cm_tot_rel_wheel_tmp * t16 * 2.0) -
                      i_v_cm_tot_rel_wheel_tmp * t16) -
                    j_v_cm_tot_rel_wheel_tmp * t16 * 2.0) +
                    k_v_cm_tot_rel_wheel_tmp * t10) +
                  l_v_cm_tot_rel_wheel_tmp * std::sin(t14)) -
                  o_v_cm_tot_rel_wheel_tmp * this->robot_inertial_par(8) * t13 *
                      2.0) +
                m_v_cm_tot_rel_wheel_tmp * t16) +
                n_v_cm_tot_rel_wheel_tmp * t16) +
              v_cm_tot_rel_wheel_tmp * t16) +
              t6 * t16) +
            b_v_cm_tot_rel_wheel_tmp * t16);
         this->v_cm_tot_rel_wheel_k(1) =
            t19 * ((((((((((((((((c_v_cm_tot_rel_wheel_tmp * t12 +
                                  d_v_cm_tot_rel_wheel_tmp * t12) -
                                e_v_cm_tot_rel_wheel_tmp * t7 * 2.0) -
                                f_v_cm_tot_rel_wheel_tmp * t7 * 2.0) +
                              g_v_cm_tot_rel_wheel_tmp * t12) -
                              h_v_cm_tot_rel_wheel_tmp * t7 * 2.0) -
                            i_v_cm_tot_rel_wheel_tmp * t7) -
                            j_v_cm_tot_rel_wheel_tmp * t7 * 2.0) +
                          this->robot_inertial_par(3) * phi_r_dot * t3 * t11) +
                          k_v_cm_tot_rel_wheel_tmp * t12) -
                        l_v_cm_tot_rel_wheel_tmp * std::cos(t14)) +
                        m_v_cm_tot_rel_wheel_tmp * t7) +
                      n_v_cm_tot_rel_wheel_tmp * t7) +
                      o_v_cm_tot_rel_wheel_tmp * t3 * t11) +
                    v_cm_tot_rel_wheel_tmp * t7) +
                    t6 * t7) +
                  b_v_cm_tot_rel_wheel_tmp * t7);
      }
        
      void J_CoMv_tot_rel_wheel_calc(double alpha,double b_gamma,double h1,double h2){

        double phi_r=this->q_k(3,0);
        double beta=this->q_k(4,0);

        double J_CoMv_tot_rel_wheel_tmp;
        double t11;
        double t12;
        double t14;
        double t15;
        double t17;
        double t18;
        double t19;
        double t2;
        double t20;
        double t21;
        double t25;
        double t26;
        double t28;
        double t3;
        double t4;
        double t5;
        double t7;
        double t8;
        //     This function was generated by the Symbolic Math Toolbox version 8.7.
        //     31-Jul-2021 14:08:50
        t2 = this->robot_dim(0,0) + this->robot_dim(4,0);
        t3 = this->robot_inertial_par(6,0) * 2.0;
        t4 = this->robot_inertial_par(8,0) * 2.0;
        t5 = this->robot_inertial_par(7,0) * 2.0;
        t7 = (beta + phi_r) + this->robot_dim(5,0);
        t8 = (b_gamma + phi_r) + this->robot_dim(5,0);
        t11 = std::cos(t7);
        t12 = std::cos(t8);
        t7 = std::sin(t7);
        t14 = std::sin(t8);
        t15 = phi_r + -this->robot_dim(6,0);
        t8 = (phi_r + -alpha) + this->robot_dim(5,0);
        t17 = this->robot_dim(3,0) * t11;
        t18 = this->robot_inertial_par(1,0) * t11;
        t19 = this->robot_dim(3,0) * t7;
        t20 = this->robot_inertial_par(1,0) * t7;
        t21 = std::cos(t8);
        t11 = std::sin(t8);
        t28 = 1.0 / ((((this->robot_inertial_par(9,0) + t3) + t4) + t5) +
                    this->robot_inertial_par(5,0) * 2.0);
        t7 = this->robot_dim(4,0) * t21;
        t8 = this->robot_dim(4,0) * t11;
        t25 = h1 * t7;
        t26 = h1 * t8;
        this->J_CoMv_tot_rel_wheel_k(0,0) = 0.0;
        this->J_CoMv_tot_rel_wheel_k(1,0) = 0.0;
        this->J_CoMv_tot_rel_wheel_k(0,1) = 0.0;
        this->J_CoMv_tot_rel_wheel_k(1,1) = 0.0;
        this->J_CoMv_tot_rel_wheel_k(0,2) = 0.0;
        this->J_CoMv_tot_rel_wheel_k(1,2) = 0.0;
        J_CoMv_tot_rel_wheel_tmp = this->robot_inertial_par(2,0) * t3;
        this->J_CoMv_tot_rel_wheel_k(0,3) =
            t28 * (((this->robot_inertial_par(9,0) *
                        ((t17 - t7) + this->robot_inertial_par(4,0) * std::sin(t15)) +
                    t5 * (t18 - t7)) -
                    t4 * (this->robot_inertial_par(3,0) * t14 + t2 * t21)) -
                  J_CoMv_tot_rel_wheel_tmp * t21);
        this->J_CoMv_tot_rel_wheel_k(1,3) =
            -t28 *
            (((this->robot_inertial_par(9,0) *
                  ((-t19 + t8) + this->robot_inertial_par(4,0) * std::cos(t15)) -
              t5 * (t20 - t8)) -
              t4 * (this->robot_inertial_par(3,0) * t12 - t2 * t11)) +
            J_CoMv_tot_rel_wheel_tmp * t11);
        J_CoMv_tot_rel_wheel_tmp = h2 * this->robot_inertial_par(3,0);
        t8 = h1 * t2;
        t7 = h1 * this->robot_inertial_par(2,0) * t3;
        this->J_CoMv_tot_rel_wheel_k(0,4) =
            t28 * (((this->robot_inertial_par(9,0) * (t17 + t25) + t5 * (t18 + t25)) -
                    this->robot_inertial_par(8,0) *
                        (J_CoMv_tot_rel_wheel_tmp * t14 - t8 * t21) * 2.0) +
                  t7 * t21);
        this->J_CoMv_tot_rel_wheel_k(1,4) =
            t28 * (((this->robot_inertial_par(9,0) * (t19 + t26) + t5 * (t20 + t26)) +
                    t4 * (J_CoMv_tot_rel_wheel_tmp * t12 + t8 * t11)) +
                  t7 * t11);
      }

      void J_CoMv_dot_tot_rel_wheel_calc(double alpha,double b_gamma,double h1,double h2){
      

        double phi_r=this->q_k(3,0);
        double beta=this->q_k(4,0);
        double phi_r_dot=this->q_k(8,0);
        double beta_dot=this->q_k(9,0);

        double J_CoMv_dot_tot_rel_wheel_tmp;
        double b_J_CoMv_dot_tot_rel_wheel_tmp;
        double t11;
        double t12;
        double t14;
        double t15;
        double t17;
        double t18;
        double t19;
        double t2;
        double t20;
        double t23;
        double t24;
        double t28;
        double t29;
        double t3;
        double t30;
        double t35;
        double t37;
        double t4;
        double t5;
        double t7;
        double t8;
        t2 = this->robot_dim(0,0) + this->robot_dim(4,0);
        t3 = this->robot_inertial_par(6,0) * 2.0;
        t4 = this->robot_inertial_par(8,0) * 2.0;
        t5 = this->robot_inertial_par(7,0) * 2.0;
        t7 = (beta + phi_r) + this->robot_dim(5,0);
        t8 = (b_gamma + phi_r) + this->robot_dim(5,0);
        t11 = std::cos(t7);
        t12 = std::cos(t8);
        t7 = std::sin(t7);
        t14 = std::sin(t8);
        t15 = phi_r + -this->robot_dim(6,0);
        t8 = (phi_r + -alpha) + this->robot_dim(5,0);
        t17 = this->robot_dim(3,0) * t11;
        t18 = this->robot_inertial_par(1,0) * t11;
        t19 = this->robot_dim(3,0) * t7;
        t20 = this->robot_inertial_par(1,0) * t7;
        t23 = std::cos(t8);
        t24 = std::sin(t8);
        t8 = 1.0 / ((((this->robot_inertial_par(9,0) + t3) + t4) + t5) +
                    this->robot_inertial_par(5,0) * 2.0);
        t11 = this->robot_dim(4,0) * t23;
        t28 = this->robot_dim(4,0) * t24;
        t29 = h1 * t11;
        t30 = h1 * t28;
        t7 = beta_dot * t8;
        t35 = t7 * (this->robot_inertial_par(9,0) * t17 + t5 * t18);
        t37 = -(t7 * (this->robot_inertial_par(9,0) * t19 + t5 * t20));
        this->J_CoMv_dot_tot_rel_wheel_k(0,0) = 0.0;
        this->J_CoMv_dot_tot_rel_wheel_k(1,0) = 0.0;
        this->J_CoMv_dot_tot_rel_wheel_k(0,1) = 0.0;
        this->J_CoMv_dot_tot_rel_wheel_k(1,1)= 0.0;
        this->J_CoMv_dot_tot_rel_wheel_k(0,2) = 0.0;
        this->J_CoMv_dot_tot_rel_wheel_k(1,2) = 0.0;
        J_CoMv_dot_tot_rel_wheel_tmp = this->robot_inertial_par(2,0) * t3;
        b_J_CoMv_dot_tot_rel_wheel_tmp = phi_r_dot * t8;
        this->J_CoMv_dot_tot_rel_wheel_k(0,3) =
            t37 +
            b_J_CoMv_dot_tot_rel_wheel_tmp *
                (((this->robot_inertial_par(9,0) *
                      ((-t19 + t28) + this->robot_inertial_par(4,0) * std::cos(t15)) -
                  t5 * (t20 - t28)) -
                  t4 * (this->robot_inertial_par(3,0) * t12 - t2 * t24)) +
                J_CoMv_dot_tot_rel_wheel_tmp * t24);
        this->J_CoMv_dot_tot_rel_wheel_k(1,3) =
            t35 +
            b_J_CoMv_dot_tot_rel_wheel_tmp *
                (((this->robot_inertial_par(9,0) *
                      ((t17 - t11) + this->robot_inertial_par(4,0) * std::sin(t15)) +
                  t5 * (t18 - t11)) -
                  t4 * (this->robot_inertial_par(3,0) * t14 + t2 * t23)) -
                J_CoMv_dot_tot_rel_wheel_tmp * t23);
        J_CoMv_dot_tot_rel_wheel_tmp = h2 * this->robot_inertial_par(3,0);
        t8 = h1 * t2;
        t7 = h1 * this->robot_inertial_par(2,0) * t3;
        this->J_CoMv_dot_tot_rel_wheel_k(0,4) =
            t37 -
            b_J_CoMv_dot_tot_rel_wheel_tmp *
                (((this->robot_inertial_par(9,0) * (t19 + t30) + t5 * (t20 + t30)) +
                  t4 * (J_CoMv_dot_tot_rel_wheel_tmp * t12 + t8 * t24)) +
                t7 * t24);
        this->J_CoMv_dot_tot_rel_wheel_k(1,4) =
            t35 +
            b_J_CoMv_dot_tot_rel_wheel_tmp *
                (((this->robot_inertial_par(9,0) * (t17 + t29) + t5 * (t18 + t29)) -
                  this->robot_inertial_par(8,0) *
                      (J_CoMv_dot_tot_rel_wheel_tmp * t14 - t8 * t23) * 2.0) +
                t7 * t23);
      }

      void alpha_gamma_calc(double beta,double *alpha,double *b_gamma){
        double A;
        double A_tmp;
        double B;
        double B_tmp;
        double cos_gamma;
        double sin_gamma;
        A_tmp = std::cos(beta);
        A = ((((this->robot_dim(0,0) * this->robot_dim(0,0) -
                this->robot_dim(1,0) * this->robot_dim(1,0)) -
              this->robot_dim(2,0) * this->robot_dim(2,0)) -
              this->robot_dim(3,0) * this->robot_dim(3,0)) +
            2.0 * this->robot_dim(1,0) * this->robot_dim(3,0) * A_tmp) /
            (2.0 * this->robot_dim(2,0) *
            (this->robot_dim(3,0) * A_tmp - this->robot_dim(1,0)));
        B_tmp = std::sin(beta);
        B = this->robot_dim(3,0) * B_tmp /
            (this->robot_dim(3,0) * std::cos(beta) - this->robot_dim(1,0));
        sin_gamma = B * B + 1.0;
        cos_gamma = (-A * B + std::sqrt(sin_gamma - A * A)) / sin_gamma;
        sin_gamma = A + B * cos_gamma;
        *b_gamma = atan2(sin_gamma, cos_gamma);
        A = this->robot_dim(3,0) / this->robot_dim(0,0);
        B = this->robot_dim(2,0) / this->robot_dim(0,0);
        *alpha = atan2(A * B_tmp - B * cos_gamma,
                              (this->robot_dim(1,0) / this->robot_dim(0,0) - A * A_tmp) -
                                  B * sin_gamma);
      }

      void knee_jacobians_h1_h2_calc(double beta,double *h1,double *h2){
        double a_tmp;
        double a_tmp_tmp;
        double b_a_tmp;
        double b_a_tmp_tmp;
        double b_knee_jacobians_tmp;
        double c_a_tmp_tmp;
        double c_knee_jacobians_tmp;
        double d_a_tmp_tmp;
        double d_knee_jacobians_tmp;
        double e_a_tmp_tmp;
        double e_knee_jacobians_tmp;
        double f_a_tmp_tmp;
        double f_knee_jacobians_tmp;
        double g_a_tmp_tmp;
        double g_knee_jacobians_tmp;
        double h1_tmp_tmp;
        double h1_tmp_tmp_tmp;
        double h_a_tmp_tmp;
        double h_knee_jacobians_tmp;
        double i_a_tmp_tmp;
        double i_knee_jacobians_tmp;
        double j_a_tmp_tmp;
        double j_knee_jacobians_tmp;
        double k_a_tmp_tmp;
        double k_knee_jacobians_tmp;
        double knee_jacobians_tmp;
        double knee_jacobians_tmp_tmp_tmp;
        double l_a_tmp_tmp;
        h1_tmp_tmp = std::cos(beta);
        h1_tmp_tmp_tmp = std::sin(beta);
        a_tmp_tmp = this->robot_dim(3,0) * h1_tmp_tmp;
        a_tmp = this->robot_dim(1,0) - a_tmp_tmp;
        b_a_tmp_tmp = this->robot_dim(3,0) * this->robot_dim(3,0);
        c_a_tmp_tmp = 2.0 * h1_tmp_tmp;
        d_a_tmp_tmp = this->robot_dim(0,0) * this->robot_dim(0,0);
        e_a_tmp_tmp = this->robot_dim(2,0) * this->robot_dim(2,0);
        f_a_tmp_tmp = c_a_tmp_tmp * this->robot_dim(1,0);
        g_a_tmp_tmp = this->robot_dim(1,0) * this->robot_dim(1,0);
        h_a_tmp_tmp = f_a_tmp_tmp * this->robot_dim(3,0);
        b_a_tmp = (((-d_a_tmp_tmp + g_a_tmp_tmp) - h_a_tmp_tmp) + e_a_tmp_tmp) +
                  b_a_tmp_tmp;
        i_a_tmp_tmp = a_tmp * a_tmp;
        j_a_tmp_tmp = h1_tmp_tmp_tmp * h1_tmp_tmp_tmp;
        k_a_tmp_tmp = b_a_tmp_tmp * j_a_tmp_tmp;
        l_a_tmp_tmp = k_a_tmp_tmp / i_a_tmp_tmp;
        knee_jacobians_tmp = 2.0 * d_a_tmp_tmp;
        b_knee_jacobians_tmp = pow(a_tmp, 3.0);
        c_knee_jacobians_tmp = 4.0 * e_a_tmp_tmp;
        h_a_tmp_tmp = (g_a_tmp_tmp - h_a_tmp_tmp) + b_a_tmp_tmp;
        d_knee_jacobians_tmp = 2.0 * b_a_tmp_tmp;
        e_knee_jacobians_tmp = this->robot_dim(1,0) * b_a_tmp_tmp;
        f_knee_jacobians_tmp = this->robot_dim(2,0) * b_knee_jacobians_tmp;
        g_knee_jacobians_tmp = c_knee_jacobians_tmp * b_knee_jacobians_tmp;
        knee_jacobians_tmp_tmp_tmp = this->robot_dim(3,0) * h1_tmp_tmp_tmp;
        h_knee_jacobians_tmp = this->robot_dim(2,0) * i_a_tmp_tmp;
        i_knee_jacobians_tmp = 2.0 * this->robot_dim(2,0) * i_a_tmp_tmp;
        j_knee_jacobians_tmp = std::sqrt(
            (l_a_tmp_tmp - b_a_tmp * b_a_tmp / (c_knee_jacobians_tmp * i_a_tmp_tmp)) +
            1.0);
        k_knee_jacobians_tmp = j_knee_jacobians_tmp + knee_jacobians_tmp_tmp_tmp *
                                                          b_a_tmp /
                                                          i_knee_jacobians_tmp;
        *h1 =
            -(2.0 * this->robot_dim(0,0) *
              ((h_knee_jacobians_tmp *
                    (((a_tmp_tmp * b_a_tmp / i_knee_jacobians_tmp +
                      e_knee_jacobians_tmp * j_a_tmp_tmp / h_knee_jacobians_tmp) -
                      k_a_tmp_tmp * b_a_tmp / f_knee_jacobians_tmp) +
                    knee_jacobians_tmp_tmp_tmp *
                        ((((((((((pow(this->robot_dim(0,0), 4.0) +
                                  c_a_tmp_tmp * d_a_tmp_tmp * this->robot_dim(1,0) *
                                      this->robot_dim(3,0)) -
                                  knee_jacobians_tmp * e_a_tmp_tmp) -
                                knee_jacobians_tmp * b_a_tmp_tmp) -
                                pow(this->robot_dim(1,0), 4.0)) +
                              c_a_tmp_tmp * pow(this->robot_dim(1,0), 3.0) *
                                  this->robot_dim(3,0)) +
                              f_a_tmp_tmp * e_a_tmp_tmp * this->robot_dim(3,0)) -
                            f_a_tmp_tmp * pow(this->robot_dim(3,0), 3.0)) +
                            pow(this->robot_dim(2,0), 4.0)) -
                          2.0 * e_a_tmp_tmp * b_a_tmp_tmp) +
                          pow(this->robot_dim(3,0), 4.0)) /
                        (g_knee_jacobians_tmp * j_knee_jacobians_tmp)) /
                    (this->robot_dim(0,0) * h_a_tmp_tmp) -
                a_tmp_tmp / this->robot_dim(0,0)) +
              2.0 * this->robot_dim(2,0) * b_a_tmp_tmp * h1_tmp_tmp_tmp *
                  (this->robot_dim(3,0) - this->robot_dim(1,0) * h1_tmp_tmp) *
                  k_knee_jacobians_tmp /
                  (this->robot_dim(0,0) * b_knee_jacobians_tmp *
                    ((l_a_tmp_tmp + 1.0) * (l_a_tmp_tmp + 1.0)))) *
              h_a_tmp_tmp) /
            (a_tmp * (((((((d_knee_jacobians_tmp * (h1_tmp_tmp * h1_tmp_tmp) +
                            d_knee_jacobians_tmp * j_a_tmp_tmp) +
                          d_a_tmp_tmp) +
                          g_a_tmp_tmp) -
                        e_a_tmp_tmp) -
                        b_a_tmp_tmp) -
                      2.0 * this->robot_dim(1,0) * this->robot_dim(3,0) * h1_tmp_tmp) +
                      2.0 * this->robot_dim(2,0) * this->robot_dim(3,0) * h1_tmp_tmp_tmp *
                          j_knee_jacobians_tmp));
          *h2 =
            -(i_a_tmp_tmp *
                  (((this->robot_dim(3,0) * std::cos(beta) *
                        ((((-(this->robot_dim(0,0) * this->robot_dim(0,0)) +
                            this->robot_dim(1,0) * this->robot_dim(1,0)) -
                            2.0 * std::cos(beta) * this->robot_dim(1,0) *
                                this->robot_dim(3,0)) +
                          this->robot_dim(2,0) * this->robot_dim(2,0)) +
                          this->robot_dim(3,0) * this->robot_dim(3,0)) /
                        (2.0 * this->robot_dim(2,0) * (a_tmp * a_tmp)) +
                    e_knee_jacobians_tmp * (h1_tmp_tmp_tmp * h1_tmp_tmp_tmp) /
                        (this->robot_dim(2,0) * (a_tmp * a_tmp))) -
                    b_a_tmp_tmp * (h1_tmp_tmp_tmp * h1_tmp_tmp_tmp) * b_a_tmp /
                        f_knee_jacobians_tmp) +
                  this->robot_dim(3,0) * std::sin(beta) *
                      ((((((((((pow(this->robot_dim(0,0), 4.0) +
                                2.0 * std::cos(beta) *
                                    (this->robot_dim(0,0) * this->robot_dim(0,0)) *
                                    this->robot_dim(1,0) * this->robot_dim(3,0)) -
                                2.0 * (this->robot_dim(0,0) * this->robot_dim(0,0)) *
                                    (this->robot_dim(2,0) * this->robot_dim(2,0))) -
                              2.0 * (this->robot_dim(0,0) * this->robot_dim(0,0)) *
                                  (this->robot_dim(3,0) * this->robot_dim(3,0))) -
                              pow(this->robot_dim(1,0), 4.0)) +
                            2.0 * std::cos(beta) *
                                pow(this->robot_dim(1,0), 3.0) *
                                this->robot_dim(3,0)) +
                            2.0 * std::cos(beta) * this->robot_dim(1,0) *
                                (this->robot_dim(2,0) * this->robot_dim(2,0)) *
                                this->robot_dim(3,0)) -
                          2.0 * std::cos(beta) * this->robot_dim(1,0) *
                              pow(this->robot_dim(3,0), 3.0)) +
                          pow(this->robot_dim(2,0), 4.0)) -
                        2.0 * (this->robot_dim(2,0) * this->robot_dim(2,0)) *
                            (this->robot_dim(3,0) * this->robot_dim(3,0))) +
                        pow(this->robot_dim(3,0), 4.0)) /
                      (g_knee_jacobians_tmp *
                        std::sqrt((b_a_tmp_tmp * (h1_tmp_tmp_tmp * h1_tmp_tmp_tmp) /
                                      (a_tmp * a_tmp) -
                                  b_a_tmp * b_a_tmp /
                                      (c_knee_jacobians_tmp * (a_tmp * a_tmp))) +
                                  1.0))) /
                  h_a_tmp_tmp +
              d_knee_jacobians_tmp * h1_tmp_tmp_tmp *
                  (this->robot_dim(3,0) - this->robot_dim(1,0) * std::cos(beta)) *
                  k_knee_jacobians_tmp /
                  (b_knee_jacobians_tmp *
                  ((l_a_tmp_tmp + 1.0) * (l_a_tmp_tmp + 1.0)))) /
            (b_a_tmp / (2.0 * this->robot_dim(2,0) * a_tmp) -
            knee_jacobians_tmp_tmp_tmp * a_tmp * k_knee_jacobians_tmp / h_a_tmp_tmp);
      }
  
   
    //PRIVATE ATTRIBUTES
    private: 
      physics::WorldPtr world; // Pointer to the world the model is in
      physics::ModelPtr model; // Pointer to the model
      physics::PhysicsEnginePtr physics_engine_ptr; //pointer to the current physics engine
      event::ConnectionPtr updateConnection; // Pointer to the update event connection
      event::ConnectionPtr end_updateConnection; // Pointer to the update event connection
      event::ConnectionPtr before_physics_updateConnection; //Pointer to the before physics update event connection
      double max_step_size;//integration step size
      uint32_t sim_steps_over_control_step; //number of controls per integration step (used to provide the controller with the right state)
      ///////////////////////////////////////////////////
      physics::JointPtr left_wheel_hinge;
      physics::JointPtr left_ah_d_hinge;
      physics::JointPtr left_ah_c_hinge;
      physics::JointPtr left_c_body_hinge;
      physics::JointPtr left_d_body_hinge;

      physics::JointPtr right_wheel_hinge;
      physics::JointPtr right_ah_d_hinge;
      physics::JointPtr right_ah_c_hinge;
      physics::JointPtr right_c_body_hinge;
      physics::JointPtr right_d_body_hinge;

      physics::LinkPtr left_wheel;
      physics::LinkPtr left_ah_link;
      physics::LinkPtr left_d_link;
      physics::LinkPtr left_c_link;
      physics::LinkPtr right_wheel;
      physics::LinkPtr right_ah_link;
      physics::LinkPtr right_d_link;
      physics::LinkPtr right_c_link;
      physics::LinkPtr body;

      physics::JointState left_wheel_hinge_state;
      physics::JointState left_ah_d_hinge_state;
      physics::JointState left_ah_c_hinge_state;
      physics::JointState left_c_body_hinge_state;
      physics::JointState left_d_body_hinge_state;

      physics::JointState right_wheel_hinge_state;
      physics::JointState right_ah_d_hinge_state;
      physics::JointState right_ah_c_hinge_state;
      physics::JointState right_c_body_hinge_state;
      physics::JointState right_d_body_hinge_state;

      physics::LinkState left_wheel_state;
      physics::LinkState left_ah_link_state;
      physics::LinkState left_d_link_state;
      physics::LinkState left_c_link_state;
      physics::LinkState right_wheel_state;
      physics::LinkState right_ah_link_state;
      physics::LinkState right_d_link_state;
      physics::LinkState right_c_link_state;
      physics::LinkState body_state;

      double x_w_left_k;
      double x_w_left_dot_k;
      double z_w_left_k;
      double z_w_left_dot_k;
      double phi_w_left_k;
      double phi_w_left_dot_k;
      double beta_left_k;
      double beta_left_dot_k;

      double x_w_right_k;
      double x_w_right_dot_k;
      double z_w_right_k;
      double z_w_right_dot_k;
      double phi_w_right_k;
      double phi_w_right_dot_k;
      double beta_right_k;
      double beta_right_dot_k;

      double phi_r_k;
      double phi_r_dot_k;

      double alpha_left_k;
      double alpha_right_k;
      double h1_left_k;
      double h1_right_k;
     ///////////////////////////////////////////////////
      double beta_spawn;
      double beta_ref_el;
      double alpha_spawn;
      double gamma_spawn;
      double phi_w_spawn;
      double phi_r_spawn;

      double phi_r0;
      ///////////////////////////////////////////////////
      ignition::math::Vector3d gravity; // world Gravity vector
      double g;
      double mu;// static friction coeff ()
      double delta_control;// controller sampling time 
      MatrixXf Kp=MatrixXf::Zero(2,2); // prop. gain matrix 
      MatrixXf Kd=MatrixXf::Zero(2,2);  // der. gain matrix
      VectorXf u_lim=VectorXf::Zero(4); // torque limits: [C_lb_wheel, C_ub_wheel, C_lb_knee, C_ub_knee]
      VectorXf beta_lim=VectorXf::Zero(2); // knee angle limits: [beta_lb, beta_ub]
      VectorXf phi_r_lim=VectorXf::Zero(2);; // robot pitch angle limits
      double cost_threshold;
      bool forward_prop_state_meas=0;
      bool model_properties_are_manually_set;
      bool use_equiv_central_knee;
      double delta_lambda; //small positive number to ensure normal reaction stays positive
      double current_time=0; // currently not used
      double sim_time;
      uint32_t sim_iteration_count=0;
      uint32_t control_iteration_count=0; //keeps track of the number of control iterations
      uint32_t N=0;
      MatrixXf planar_compensation_Kp=MatrixXf::Zero(2,2);
      MatrixXf planar_compensation_Kd=MatrixXf::Zero(2,2);
      //////////////// model properties //////////////////////
      VectorXf robot_dim=VectorXf::Zero(7);
      VectorXf robot_inertial_par=VectorXf::Zero(15);
      VectorXf robot_knee_el_par=VectorXf::Zero(7);
      ///////////////////////////////////////////////////
      MatrixXf J_CoMv_tot_rel_wheel_k=MatrixXf::Zero(2,5);
      MatrixXf J_CoMv_dot_tot_rel_wheel_k=MatrixXf::Zero(2,5);
      VectorXf p_cm_tot_rel_wheel_k=VectorXf::Zero(2);
      VectorXf v_cm_tot_rel_wheel_k=VectorXf::Zero(2);
      double cost_fun_value_k;
      ///////////////QP matrices ///////////////////////////
      MatrixXf H_k=MatrixXf::Zero(7,7); //quadratic matrix for QP
      VectorXf f_k=VectorXf::Zero(7); //linear vector for QP
      
      MatrixXf A_nonholonomic_k=MatrixXf::Zero(2,7); //matrix of nonholonomic contraint formatted for QP
      MatrixXf A_lambda_dyn_k=MatrixXf::Zero(2,7); //floating base contraint matrix formatted for QP
      MatrixXf A_feas_k=MatrixXf::Zero(1,7);  //feasibility array contraint formatted for QP (1x)-->underactuated system--> joint trajectories cannot be set arbitrarly 
      MatrixXf A_lambda_con_k=MatrixXf::Zero(3,7) ; //constraint matrix on the contact force vector lamdba
      MatrixXf A_u_check_k=MatrixXf::Zero(2,7); //auxiliary matrix for input limits
      MatrixXf B_check_inv_k=MatrixXf::Zero(2,2) ; //auxiliary matrix for input limits
      MatrixXf A_u_lim_k=MatrixXf::Zero(8,7); //constraint matrix for setting input limits
      MatrixXf A_state_lim_k=MatrixXf::Zero(4,7); //constraint matrix for setting state limits (currently not working)
      
      VectorXf b_nonholonomic_k=VectorXf::Zero(2); //right constraint vector for nonholonomic constraints
      VectorXf b_lambda_dyn_k=VectorXf::Zero(2); //right constraint vector
      VectorXf b_feas_k=VectorXf::Zero(1); //right constraint vector
      VectorXf b_lambda_con_k=VectorXf::Zero(3); //right constraint vector
      VectorXf r_u_check_k=VectorXf::Zero(2);//right constraint vector
      VectorXf b_u_lim_k=VectorXf::Zero(8); //right constraint vector
      VectorXf b_state_lim_k=VectorXf::Zero(4); //right constraint vector
      
      MatrixXf A_eq_k=MatrixXf::Zero(A_nonholonomic_k.rows()+A_lambda_dyn_k.rows()+A_feas_k.rows(),A_nonholonomic_k.cols()); //final equality constraint matrix 
      MatrixXf A_ineq_k=MatrixXf::Zero(A_lambda_con_k.rows()+A_u_lim_k.rows()+A_state_lim_k.rows(),A_nonholonomic_k.cols()); //final inequality constraint matrix 
      VectorXf b_eq_k=VectorXf::Zero(b_nonholonomic_k.rows()+b_lambda_dyn_k.rows()+b_feas_k.rows());//final equality constraint right vector
      VectorXf b_ineq_k=VectorXf::Zero(b_lambda_con_k.rows()+b_u_lim_k.rows()+b_state_lim_k.rows()); //final inequality constraint right vector

      //VectorXf b_neg_inf_QP=(-std::numeric_limits<double>::infinity())*(VectorXf::Ones(this->b_ineq_k.rows()));//used to be able to input to the QP solver the inequality constraints mimicking Matlab's solver
      VectorXf b_neg_inf_QP=-10000000000000*(VectorXf::Ones(this->b_ineq_k.rows()));
      // VectorXf b_neg_inf_QP=VectorXf::Constant(this->b_ineq_k.rows(),-std::numeric_limits<double>::infinity());

     /////////////////////////////////////////////////////
      VectorXf Chi_ref_k=VectorXf::Zero(2);//input reference trajectory at the current sampling time (for example CoM cartesian trajectory)
      VectorXf Chi_dot_ref_k=VectorXf::Zero(2); //input derivative reference trajectory at the current sampling time
      VectorXf Chi_ddot_ref_k=VectorXf::Zero(2); //input double derivative reference trajectory at the current sampling time
      VectorXf Chi_k=VectorXf::Zero(2); //current derivate value of the optimization target (for example current CoM position)
      VectorXf Chi_dot_k=VectorXf::Zero(2); //current derivate value of the optimization target 
      
      VectorXf Chi_ddot_command_k=VectorXf::Zero(2); //commanded reference (reference trajectory adjusted with PD)

      /////////////////////////////////////////////////////
      VectorXf x_k=VectorXf::Zero(7);// vector containing solution of QP: [q_p_ddot_k;lambda_k]'
      VectorXf q_k=VectorXf::Zero(10); //current state (from measurements)
      VectorXf u_left_k=VectorXf::Zero(2); //current input  
      VectorXf u_right_k=VectorXf::Zero(2); //current input 
      VectorXf u_middle_k=VectorXf::Zero(2); //current input 



      /////////////////////////////////////////////////////   
      MatrixXf Q_sol; // holds solution for post processing
      MatrixXf U_left; // holds left side inputs throughout simulation for post processing
      MatrixXf U_right; // holds right side inputs throughout simulation for post processing
      MatrixXf U_middle;//hold the control computed using the "equivalent" central knee
      VectorXf dt_control_calc; //holds the control computation time
      VectorXf Time; // holds simulation time
      MatrixXf Chi_ref;// hold the reference Task for the whole sim. time
      MatrixXf Chi_dot_ref;
      MatrixXf Chi_ddot_ref;
      MatrixXf Chi;
      MatrixXf Chi_dot;
      MatrixXf Chi_ddot;
      MatrixXf X_k;//holds all the QP solutions during the simulation-->note that it is really useful to have an estimate of the joint accelerations and ground reactions
      VectorXf QP_cost_value;


  };
  
  // Register this plugin with the simulator
  GZ_REGISTER_MODEL_PLUGIN(ModelController)
}

