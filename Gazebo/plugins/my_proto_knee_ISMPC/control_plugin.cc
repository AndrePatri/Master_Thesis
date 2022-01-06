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

#include <boost/numeric/odeint.hpp> // numerical integration-->use runge_kutta4 class or runge_kutta54_cash_karp with integrate()

using namespace boost::numeric::odeint;

using namespace std; 
using namespace Eigen;

namespace gazebo
{
  class ModelController : public ModelPlugin
  {
    public: ModelController(){//class constructor
    
      this->planar_compensation_Kp<<0,0,0,0;
      this->planar_compensation_Kd<<0,0,0,0;

      this->u_lim<<2.0,20.0;//input limits

      this->max_step_size=0.001;
      this->sim_steps_over_control_step=20;//number of simulation steps per sample step->fixed the control dt as a multiple of the max_step_size
      this->delta_s=max_step_size*this->sim_steps_over_control_step;
      this->sim_time=20; //amount of time to be simulated, in seconds
      this->N=round(this->sim_time/this->delta_s); //total number of control iterations --> necessary to assign solutions for post-processing
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

      //////////////////////////////////////////////////
      this->T_p=1.5;//prediction horizon
      this->N_p=round(this->T_p/this->delta_s);//number of prediction intervals
      this->Time=VectorXf::Zero(this->N+this->N_p+1);//simulation horizon plus last prediction horizon 
      for (int i=1;i<=(this->N+this->N_p);i++){
        this->Time(i)=this->Time(i-1)+this->delta_s;
      }
      this->cost_weights_y1<<15.0,0,0.01,0.01,0;//cost function weights first output (phi_w)
      this->cost_weights_y2<<1000.0,0,0.001,0.001,0;//cost function weights second output (beta)
      // this->cost_weights_y1<<15.0,0,0.01,0.01,0;//cost function weights first output (phi_w)
      // this->cost_weights_y2<<200.0,0,0.0001,0.001,0;//cost function weights second output (beta)

      this->C_pos<<1.0,0;//outputs position matrix (of the single feedback linearized subsystems)
      this->C_vel<<0,1.0;//outputs velocity matrix (of the single feedback linearized subsystems)

      this->Ad<<1.0,this->delta_s,0,1.0;
      this->Bd<<pow(this->delta_s,2)/2.0,this->delta_s;

      this->V_min_k=VectorXf::Zero(2*this->N_p);
      this->V_max_k=VectorXf::Zero(2*this->N_p);
      this->V_k=VectorXf::Zero(2*this->N_p);

      this->Delta_k=MatrixXf::Zero(1,this->N_p);
      this->Y_ref_k=MatrixXf::Zero(this->N_p,2);
      this->dY_ref_k=MatrixXf::Zero(this->N_p,2);
      this->Y_ref=MatrixXf::Zero(this->N+this->N_p+1,2);
      this->dY_ref=MatrixXf::Zero(this->N+this->N_p+1,2);

      this->x_w_freq=0.3;
      this->beta_freq=0.3;
      this->x_w_amplitute=0.5;
      this->beta_amplitute=5.0*M_PI/180.0;
      this->x_w_offset=0;
      this->beta_offset=80.0*M_PI/180.0;
      this->phase_lag_x_w=0.0*M_PI/180.0;
      this->phase_lag_beta=-90.0*M_PI/180.0;

      /////////////////QP///////////////////////
      this->Omega_pos=MatrixXf::Zero(this->N_p,2);
      this->Gamma_pos=MatrixXf::Zero(this->N_p,this->N_p);
      this->Omega_vel=MatrixXf::Zero(this->N_p,2);
      this->Gamma_vel=MatrixXf::Zero(this->N_p,this->N_p);

      MatrixXf Ad_aux=MatrixXf::Zero(2,2);
      MatrixXf Ad_pow=MatrixXf::Identity(2,2);
      double aux=0;

      VectorXf aux_Omega_pos=this->C_pos;
      VectorXf aux_Omega_vel=this->C_vel;

      for (int i=1;i<=(this->N_p);i++){
          for (int j=1;j<=(this->N_p-i+1);j++){
            this->Gamma_pos(i+j-2,j-1)=aux_Omega_pos.dot(this->Bd);
            this->Gamma_vel(i+j-2,j-1)=aux_Omega_vel.dot(this->Bd);
          }    
          aux_Omega_pos=aux_Omega_pos.transpose()*this->Ad;
          this->Omega_pos.row(i-1)=aux_Omega_pos;
          aux_Omega_vel=aux_Omega_vel.transpose()*this->Ad;
          this->Omega_vel.row(i-1)=aux_Omega_vel;
      }

      this->Gamma_pos1=MatrixXf::Zero(this->N_p,(this->N_p)*2);
      this->Gamma_pos2=MatrixXf::Zero(this->N_p,(this->N_p)*2);
      this->Gamma_vel1=MatrixXf::Zero(this->N_p,(this->N_p)*2);
      this->Gamma_vel2=MatrixXf::Zero(this->N_p,(this->N_p)*2);

      this->Gamma_pos1<<this->Gamma_pos,MatrixXf::Zero(this->N_p,this->N_p);
      this->Gamma_pos2<<MatrixXf::Zero(this->N_p,this->N_p),this->Gamma_pos;
      this->Gamma_vel1<<this->Gamma_vel,MatrixXf::Zero(this->N_p,this->N_p);
      this->Gamma_vel2<<MatrixXf::Zero(this->N_p,this->N_p),this->Gamma_vel;

      this->Du=MatrixXf::Identity(this->N_p,this->N_p);
      
      this->D1u=MatrixXf::Zero(this->N_p,(this->N_p)*2);
      this->D2u=MatrixXf::Zero(this->N_p,(this->N_p)*2);
      this->D1u<<this->Du,MatrixXf::Zero(this->N_p,this->N_p);
      this->D2u<<MatrixXf::Zero(this->N_p,this->N_p),this->Du;

      this->Du_diff=MatrixXf::Zero(this->N_p-1,this->N_p);
      for (int i=1;i<=2;i++){
          for (int j=1;j<=(this->N_p-1);j++){
              if (i==1){
                this->Du_diff(j-1,j-1)=1;
              }
              else{
                this->Du_diff(j-1,j)=-1;
              }
            }    
      }

      this->D1u_diff=MatrixXf::Zero(this->N_p-1,(this->N_p)*2);
      this->D2u_diff=MatrixXf::Zero(this->N_p-1,(this->N_p)*2);
      this->D1u_diff<<this->Du_diff,MatrixXf::Zero(this->N_p-1,this->N_p);
      this->D2u_diff<<MatrixXf::Zero(this->N_p-1,this->N_p),this->Du_diff;

      this->Du_diff_prime=MatrixXf::Zero(this->N_p,this->N_p);
      this->Du_diff_prime(0,0)=1;
      this->D1u_diff_prime=MatrixXf::Zero(this->N_p,(this->N_p)*2);
      this->D2u_diff_prime=MatrixXf::Zero(this->N_p,(this->N_p)*2);
      this->D1u_diff_prime<<this->Du_diff_prime,MatrixXf::Zero(this->N_p,this->N_p);
      this->D2u_diff_prime<<MatrixXf::Zero(this->N_p,this->N_p),this->Du_diff_prime;

      this->H=this->cost_weights_y1(0)*this->Gamma_pos1.transpose()*this->Gamma_pos1+this->cost_weights_y1(1)*this->Gamma_vel1.transpose()*this->Gamma_vel1+this->cost_weights_y1(2)*this->D1u.transpose()*this->D1u+
              this->cost_weights_y1(3)*this->D1u_diff.transpose()*this->D1u_diff+this->cost_weights_y1(4)*this->D1u_diff_prime.transpose()*this->D1u_diff_prime+
              this->cost_weights_y2(0)*this->Gamma_pos2.transpose()*this->Gamma_pos2+this->cost_weights_y2(1)*this->Gamma_vel2.transpose()*this->Gamma_vel2+this->cost_weights_y2(2)*this->D2u.transpose()*this->D2u+
              this->cost_weights_y2(3)*this->D2u_diff.transpose()*this->D2u_diff+this->cost_weights_y2(4)*this->D2u_diff_prime.transpose()*this->D2u_diff_prime;
      this->f_k=VectorXf::Zero(2*this->N_p);
      this->f1_diff_k=VectorXf::Zero(2*this->N_p);
      this->f2_diff_k=VectorXf::Zero(2*this->N_p);

      // this->Au_lim_k=MatrixXf::Zero(4,(this->N_p)*2);
      // this->Au_lim_k(0,0)=-1;
      // this->Au_lim_k(1,this->N_p)=-1;
      // this->Au_lim_k(2,0)=1;
      // this->Au_lim_k(3,this->N_p)=1;
      // this->bu_lim_k=VectorXf::Zero(4);

      this->A_bound=MatrixXf::Zero(1,this->N_p);
      this->b_bound=VectorXf::Zero(1);

      this->A_eq_k=MatrixXf::Zero(1,(this->N_p)*2);
      this->b_eq_k=VectorXf::Zero(1);

      // this->Gamma_pos_k=MatrixXf::Zero(4,(this->N_p)*2);
      this->A_ineq_k=MatrixXf::Zero(2,(this->N_p)*2);
      this->A_ineq_k(0,0)=1;
      this->A_ineq_k(1,this->N_p)=1;

      // this->b_ineq_k=VectorXf::Zero(4);
      this->b_ineq_k=VectorXf::Zero(2);
      this->b_ineq_max_k=VectorXf::Zero(2);
      this->b_ineq_max_k=VectorXf::Zero(2);

      this->cost_threshold=1000;
      this->QP_cost_value=VectorXf::Zero(this->N);

      this->VV_k=MatrixXf::Zero(this->N,2*this->N_p);

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
      // Generates reference trajectory
      this->Y_traj_setter();// sets the reference output trajectory
    
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
      }
    
    public:void BeforePhysicsUpdate(){
      // Trying to keep the motion planar kinematically
      ignition::math::Vector3d current_body_angular_vel=this->body->WorldAngularVel();
      ignition::math::Vector3d angular_vel_correction=ignition::math::Vector3d(0, current_body_angular_vel.Y(),0);
      this->body->SetAngularVel(angular_vel_correction);	//sets to zero vertical angular velocity of body link
      // ignition::math::Vector3d current_body_linear_vel=this->body->WorldCoGLinearVel();
      // ignition::math::Vector3d linear_vel_correction=ignition::math::Vector3d(current_body_linear_vel.X(), 0,current_body_linear_vel.Z());
      // this->body->SetLinearVel(linear_vel_correction);	//sets planar motion

    }


    public: void OnUpdate(){// Called by the world update START event (@ each sim update --> control sample)
          
     
      /////////////////////////////////////////////////////////////////////////////////////////
      //Reading current time (referred to the end of simulation iteration)
      this->current_time=(this->world->SimTime()).Double();//gets current sim_time as double
      this->sim_iteration_count=this->world->Iterations();//gets current simulation iteration count
      //(it is always the end time of the control step)
      if (current_time>this->sim_time) {
        this->world->SetPaused(true);
        this->WriteToTxT();//write simulation results to .txt file	
        //this->world->Reset();
        //this->world->Stop();

      }//if current time is greater than the required simulation horizon, stop/pause the simulation 
      // (the simulation will never exceed sim_time, due to how sim_time is output)
     
      else{ //simulate
        
        if (((this->sim_iteration_count-1)%this->sim_steps_over_control_step)==0){//Corresponds to a control interval-->compute and apply a new control
         
          this->control_iteration_count=((this->sim_iteration_count-1)/this->sim_steps_over_control_step)+1;
          // this->Time(this->control_iteration_count,0)=this->control_iteration_count*this->delta_s;

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

            this->Y_traj_sampler();//samples the current reference trajectory
            this->lambda_b_T_zeta_u_calc();
            this->f_k_calc();
            this->M_C_bar_k_calc();
            this->V_min_max_k_calc();
            this->QPMatrixCalc();// computes all necessary QP matrices
            this->V_k_calc();
            // // this->cost_fnctn_val_calc();// computes the cost function-->necessary to check if the QP goes crazy; if it does, do not apply any input or apply the previous one
            // // this->QP_cost_value(this->control_iteration_count-1,0)=this->cost_fun_value_k;
            this->ISMPC_input_calc_and_application();

            end= chrono::high_resolution_clock::now();
            std::chrono::duration<double> dt_control = end - start;

            this->dt_control_calc(this->control_iteration_count-1,0)=dt_control.count();

            this->PrintStuff();
          }

        }
        else{//not in a control interval-->use previous input

          this->left_wheel_hinge->SetForce( 0,  this->U_left(this->control_iteration_count-1,0));
          this->left_ah_d_hinge->SetForce( 0,   this->U_left(this->control_iteration_count-1,1));
          this->right_wheel_hinge->SetForce( 0, this->U_right(this->control_iteration_count-1,0)); 
          this->right_ah_d_hinge->SetForce( 0,  this->U_right(this->control_iteration_count-1,1));

        }

      }
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
       
      printf("q_k:\n");
      for(int i=0;i<(this->q_k).rows();i++)  // loop 3 times for three lines
        {
          for(int j=0;j<(this->q_k).cols();j++)  // loop for the three elements on the line
          {
              cout<<this->q_k(i,j);  // display the current element out of the array
              printf("\t");
          }
          cout<<endl;  // when the inner loop is done, go to a new line
        }
        printf("\n");

      // printf("alpha_spawn:\t%f",this->alpha_spawn);
      // printf("\n");

  
      // printf("gamma_spawn:\t%f",this->gamma_spawn);
      // printf("\n");

      // printf("N_p:\t%d",this->N_p);
      // printf("\n");


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

      // printf("Omega_pos:\n");
      // for(int i=0;i<(this->Omega_pos).rows();i++)  // loop 3 times for three lines
      //   {
      //     for(int j=0;j<(this->Omega_pos).cols();j++)  // loop for the three elements on the line
      //     {
      //         cout<<this->Omega_pos(i,j);  // display the current element out of the array
      //         printf("\t");
      //     }
      //     cout<<endl;  // when the inner loop is done, go to a new line
      //   }
      //   printf("\n");

      // printf("Gamma_pos:\n");
      // for(int i=0;i<(this->Gamma_pos).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->Gamma_pos).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->Gamma_pos(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");

      // printf("Ad:\n");
      // for(int i=0;i<(this->Ad).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->Ad).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->Ad(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");
      // printf("Bd:\n");
      // for(int i=0;i<(this->Bd).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->Bd).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->Bd(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");
      // printf("C:\n");
      // for(int i=0;i<(this->C_pos).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->C_pos).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->C_pos(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");
      // printf("Gamma_pos1:\n");
      // for(int i=0;i<(this->Gamma_pos1).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->Gamma_pos1).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->Gamma_pos1(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");

      // printf("Gamma_pos2:\n");
      // for(int i=0;i<(this->Gamma_pos2).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->Gamma_pos2).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->Gamma_pos2(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");

      // printf("D1u:\n");
      // for(int i=0;i<(this->D1u).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->D1u).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->D1u(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");

      // printf("D2u:\n");
      // for(int i=0;i<(this->D2u).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->D2u).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->D2u(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");

      // printf("D2u:\n");
      // for(int i=0;i<(this->D2u).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->D2u).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->D2u(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");

      // printf("D2u:\n");
      // for(int i=0;i<(this->D2u).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->D2u).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->D2u(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");


      
      // printf("Du_diff:\n");
      // for(int i=0;i<(this->Du_diff).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->Du_diff).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->Du_diff(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");

      // printf("D1u_diff:\n");
      // for(int i=0;i<(this->D1u_diff).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->D1u_diff).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->D1u_diff(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");

      // printf("D2u_diff:\n");
      // for(int i=0;i<(this->D2u_diff).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->D2u_diff).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->D2u_diff(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");

      // printf("Du_diff_prime:\n");
      // for(int i=0;i<(this->Du_diff_prime).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->Du_diff_prime).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->Du_diff_prime(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");
      // printf("D1u_diff_prime:\n");
      // for(int i=0;i<(this->D1u_diff_prime).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->D1u_diff_prime).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->D1u_diff_prime(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");
      // printf("D2u_diff_prime:\n");
      // for(int i=0;i<(this->D2u_diff_prime).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->D2u_diff_prime).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->D2u_diff_prime(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");

      // printf("H:\n");
    
      // for(int i=0;i<(H).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(H).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<H(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");

      // printf("f_k:\n");
    
      // for(int i=0;i<(this->f_k).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->f_k).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->f_k(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");

      // printf("Y_ref:\n");
      // for(int i=0;i<(this->Y_ref).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->Y_ref).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->Y_ref(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");

      printf("Y_ref_k:\n");
      for(int i=0;i<(this->Y_ref_k).rows();i++)  // loop 3 times for three lines
      {
        for(int j=0;j<(this->Y_ref_k).cols();j++)  // loop for the three elements on the line
        {
            cout<<this->Y_ref_k(i,j);  // display the current element out of the array
            printf("\t");
        }
        cout<<endl;  // when the inner loop is done, go to a new line
      }
      printf("\n");

      // printf("dY_ref:\n");
      // for(int i=0;i<(this->dY_ref).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->dY_ref).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->dY_ref(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");

      // printf("dY_ref_k:\n");
      // for(int i=0;i<(this->dY_ref_k).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->dY_ref_k).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->dY_ref_k(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");

      // printf("Time\n");
      // for(int i=0;i<(this->Time).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->Time).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->Time(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");

      // printf("M_bar_k:\n");
      // for(int i=0;i<(this->M_bar_k).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->M_bar_k).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->M_bar_k(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");

      // printf("C_bar_k:\n");
      // for(int i=0;i<(this->C_bar_k).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->C_bar_k).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->C_bar_k(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");

      // printf("V_min_k:\n");
      // for(int i=0;i<(this->V_min_k).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->V_min_k).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->V_min_k(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");

      // printf("V_max_k:\n");
      // for(int i=0;i<(this->V_max_k).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->V_max_k).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->V_max_k(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");

      // printf("V_k:\n");
      // for(int i=0;i<(this->V_k).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->V_k).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->V_k(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");

      // printf("lambda_u_k:\t%f",this->lambda_u_k);
      // printf("\n");

      // printf("b_u_k:\t%f",this->b_u_k);
      // printf("\n");

      // printf("phi_r_st_eq_k:\t%f",this->phi_r_st_eq_k);
      // printf("\n");

      // printf("dphi_r_st_eq_k_dbeta:\t%f",this->dphi_r_st_eq_k_dbeta);
      // printf("\n");

      // printf("ddphi_r_st_eq_k_ddbeta:\t%f",this->ddphi_r_st_eq_k_ddbeta);
      // printf("\n");

      // printf("l_p_eq_k:\t%f",this->l_p_eq_k);
      // printf("\n");

      // printf("I_p_eq_k:\t%f",this->I_p_eq_k);
      // printf("\n");

      // printf("m_p_eq:\t%f",this->m_p_eq);
      // printf("\n");

      // printf("T_k:\n");
      // for(int i=0;i<(this->T_k).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->T_k).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->T_k(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");

      // printf("zeta_k:\n");
      // for(int i=0;i<(this->zeta_k).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->zeta_k).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->zeta_k(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");

      // printf("A_eq:\n");
      // for(int i=0;i<(this->A_eq_k).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->A_eq_k).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->A_eq_k(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");
      
      // printf("b_eq:\n");
      // for(int i=0;i<(this->b_eq_k).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->b_eq_k).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->b_eq_k(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");

      // printf("A_ineq:\n");
      // for(int i=0;i<(this->A_ineq_k).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->A_ineq_k).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->A_ineq_k(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");

      // printf("b_ineq:\n");
      // for(int i=0;i<(this->b_ineq_k).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->b_ineq_k).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->b_ineq_k(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");

      printf("u_k:\n");
      for(int i=0;i<(this->u_middle_k).rows();i++)  // loop 3 times for three lines
      {
        for(int j=0;j<(this->u_middle_k).cols();j++)  // loop for the three elements on the line
        {
            cout<<this->u_middle_k(i,j);  // display the current element out of the array
            printf("\t");
        }
        cout<<endl;  // when the inner loop is done, go to a new line
      }
      printf("\n");
      
      // printf("Delta:\n");
      // for(int i=0;i<(this->Delta_k).rows();i++)  // loop 3 times for three lines
      // {
      //   for(int j=0;j<(this->Delta_k).cols();j++)  // loop for the three elements on the line
      //   {
      //       cout<<this->Delta_k(i,j);  // display the current element out of the array
      //       printf("\t");
      //   }
      //   cout<<endl;  // when the inner loop is done, go to a new line
      // }
      // printf("\n");
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //private method 

    private:void WriteToTxT(){

      ofstream myfile;

      myfile.open ("~/IIT_CRIS_HHCM/Gazebo/sim_results/ISMPC/Time.txt");
      myfile << this->Time.head(this->N+1); 
      myfile.close();
      myfile.open ("~/IIT_CRIS_HHCM/Gazebo/sim_results/ISMPC/Q_sol.txt");
      myfile << this->Q_sol; 
      myfile.close();
      myfile.open ("~/IIT_CRIS_HHCM/Gazebo/sim_results/ISMPC/U_left.txt");
      myfile << this->U_left; 
      myfile.close();
      myfile.open ("~/IIT_CRIS_HHCM/Gazebo/sim_results/ISMPC/U_right.txt");
      myfile << this->U_right; 
      myfile.close();

      MatrixXf Y_ref_aux=MatrixXf::Zero(this->N+1,2);
      Y_ref_aux<<(this->Y_ref.col(0)).head(this->N+1),(this->Y_ref.col(1)).head(this->N+1);
      MatrixXf dY_ref_aux=MatrixXf::Zero(this->N+1,2);
      dY_ref_aux<<(this->dY_ref.col(0)).head(this->N+1),(this->dY_ref.col(1)).head(this->N+1);

      myfile.open ("~/IIT_CRIS_HHCM/Gazebo/sim_results/ISMPC/Y_ref.txt");
      myfile << Y_ref_aux; 
      myfile.close();
      myfile.open ("~/IIT_CRIS_HHCM/Gazebo/sim_results/ISMPC/dY_ref.txt");
      myfile <<dY_ref_aux; 
      myfile.close();
      myfile.open ("~/IIT_CRIS_HHCM/Gazebo/sim_results/ISMPC/dt_control.txt");
      myfile << this->dt_control_calc; 
      myfile.close();
      myfile.open ("~/IIT_CRIS_HHCM/Gazebo/sim_results/ISMPC/cost.txt");
      myfile << this->QP_cost_value; 
      myfile.close();
      myfile.open ("~/IIT_CRIS_HHCM/Gazebo/sim_results/ISMPC/U_middle.txt");
      myfile << this->U_middle; 
      myfile.close();
      myfile.open ("~/IIT_CRIS_HHCM/Gazebo/sim_results/ISMPC/VV_k.txt");
      myfile << this->VV_k; 
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

        double k1t=30.0;
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
    //PRIVATE method 

    private: void initial_q_setter(){
        
        ignition::math::Pose3d left_wheel_pose;
        left_wheel_pose.Set(0, 0, this->robot_inertial_par(0,0), 0, -phi_r0, 0);
        this->model->SetLinkWorldPose(left_wheel_pose,this->left_wheel); 
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

    //  // OBTAINING beta and phi_r derivatives numerically to correct for strange velocity errors
    //   if (this->control_iteration_count==0){//before first sample time
    //   }
      
    //   else{
    //     this->x_w_left_dot_k=(this->Q_sol((this->control_iteration_count),0)-this->Q_sol((this->control_iteration_count-1),0))/(this->delta_s);
    //     this->z_w_left_dot_k=(this->Q_sol((this->control_iteration_count),1)-this->Q_sol((this->control_iteration_count-1),1))/(this->delta_s);
    //     this->phi_w_left_dot_k=(this->Q_sol((this->control_iteration_count),2)-this->Q_sol((this->control_iteration_count-1),2))/(this->delta_s);
    //     this->phi_r_dot_k=(this->Q_sol((this->control_iteration_count),3)-this->Q_sol((this->control_iteration_count-1),3))/(this->delta_s);
    //     this->beta_left_dot_k=(this->Q_sol((this->control_iteration_count),4)-this->Q_sol((this->control_iteration_count-1),4))/(this->delta_s);
    //     // this->beta_right_dot_k=(this->Q_sol((this->control_iteration_count),4)-this->Q_sol((this->control_iteration_count-1),4))/(this->delta_s);
    //     this->Q_sol.row(this->control_iteration_count)<<this->x_w_left_k,
    //                                                     this->z_w_left_k,
    //                                                     this->phi_w_left_k,
    //                                                     this->phi_r_k,
    //                                                     this->beta_left_k,
    //                                                     this->x_w_left_dot_k,
    //                                                     this->z_w_left_dot_k,
    //                                                     this->phi_w_left_dot_k,
    //                                                     this->phi_r_dot_k,
    //                                                     this->beta_left_dot_k;
      

    //   }

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
    ///////////////////////////////////////////////////////////////////////////////////////////
    //private method for the actual controller-->to be called inside Load and Onupdate (similar to Matlab)
    private: void ISMPC_input_calc_and_application(){

      VectorXf Delta_left=VectorXf::Zero(2);
      VectorXf Delta_dot_left=VectorXf::Zero(2);
      VectorXf Delta_right=VectorXf::Zero(2);
      VectorXf Delta_dot_right=VectorXf::Zero(2);

      this->u_middle_k=(this->M_bar_k*this->v_k+this->C_bar_k)/2.0;

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

      if ( this->cost_fun_value_k>this->cost_threshold || isnan(this->u_middle_k(0,0)) || isnan(this->u_middle_k(1,0)) || abs(this->q_k(3,0))>=M_PI/180.0*70.0){ 
        //if the cost threshold is exceeded or if the QP cannot compute a solution or the robot fell, do not apply any input
        this->u_left_k<<0,0;
        this->u_right_k<<0,0;
        this->u_middle_k<<0,0;
      }

      this->U_middle.row(this->control_iteration_count-1)=this->u_middle_k.transpose();//assigning inputs
      this->U_left.row(this->control_iteration_count-1)=(this->u_left_k).transpose();
      this->U_right.row(this->control_iteration_count-1)=(this->u_right_k).transpose();

      this->left_wheel_hinge->SetForce(0,  this->u_left_k(0,0));
      this->left_ah_d_hinge->SetForce(0,   this->u_left_k(1,0));
      this->right_wheel_hinge->SetForce(0, this->u_right_k(0,0)); 
      this->right_ah_d_hinge->SetForce(0,  this->u_right_k(1,0));
  
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    //PRIVATE method for computing all necessary QP matrices

    private: void QPMatrixCalc(){
      
      for (int i=1;i<=(this->N_p);i++){
        this->Delta_k(0,i-1)=exp(-this->lambda_u_k*this->delta_s*(i-1));
      }
      this->Delta_k=this->Delta_k*(exp(-this->lambda_u_k*this->delta_s)-1)/this->lambda_u_k;
      
      // this->bu_lim_k(0)=-this->v1_min_k;
      // this->bu_lim_k(1)=-this->v2_min_k;
      // this->bu_lim_k(2)=this->v1_max_k;
      // this->bu_lim_k(3)=this->v2_max_k;

      this->A_bound=this->Delta_k*this->b_u_k;
      this->b_bound=this->zeta_u_k;

      this->A_eq_k<<this->A_bound,MatrixXf::Zero(1,this->N_p);
      // this->A_ineq_k<<this->Au_lim_k;
      this->b_eq_k=this->b_bound;
      // this->b_ineq_k=this->bu_lim_k;
      this->b_ineq_max_k=this->v_max_k;
      this->b_ineq_min_k=this->v_min_k;
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////
    //PRIVATE method for getting the current value of the reference output over the prediction horizon (and its derivative) 

    private: void Y_traj_setter(){

      for (int j=1;j<=(this->N+this->N_p+1);j++){
        // Uncomment for sinusoidal ref. trajectory

        this->Y_ref.row(j-1)<<this->x_w_offset/this->robot_inertial_par(0)+this->x_w_amplitute/this->robot_inertial_par(0)*sin(2*M_PI*this->x_w_freq*this->Time(j-1)+this->phase_lag_x_w),this->beta_offset+this->beta_amplitute*sin(2*M_PI*this->beta_freq*this->Time(j-1)+this->phase_lag_beta);
        this->dY_ref.row(j-1)<<(2*M_PI*this->x_w_freq*this->x_w_amplitute/this->robot_inertial_par(0))*cos(2*M_PI*this->x_w_freq*this->Time(j-1)+this->phase_lag_x_w),2*M_PI*this->beta_freq*this->beta_amplitute*cos(2*M_PI*this->beta_freq*this->Time(j-1)+this->phase_lag_beta); 
       
        // // uncomment for hybrid re. trajectory
        // this->Y_ref.row(j-1)<<((atan(this->Time(j-1)-this->sim_time/2.0))-(atan(-this->sim_time/2.0)))/this->robot_inertial_par(0),
        //                       this->beta_offset+this->beta_amplitute*sin(2*M_PI*this->beta_freq*this->Time(j-1)+this->phase_lag_beta);
        // this->dY_ref.row(j-1)<< (1.0/(1.0+pow(this->Time(j-1)-this->sim_time/2.0,2)))/this->robot_inertial_par(0),
        //                        2*M_PI*this->beta_freq*this->beta_amplitute*cos(2*M_PI*this->beta_freq*this->Time(j-1)+this->phase_lag_beta);
      }    

       
    }
    //PRIVATE method for getting the current value of the reference output over the prediction horizon (and its derivative) 

    private: void Y_traj_sampler(){
      VectorXf Y_ref_K_aux1=this->Y_ref.col(0);
      VectorXf Y_ref_K_aux2=this->Y_ref.col(1);
      VectorXf dY_ref_K_aux1=this->dY_ref.col(0);
      VectorXf dY_ref_K_aux2=this->dY_ref.col(1);

      this->Y_ref_k<<Y_ref_K_aux1.segment(this->control_iteration_count,this->N_p),
                     Y_ref_K_aux2.segment(this->control_iteration_count,this->N_p);
      this->dY_ref_k<<dY_ref_K_aux1.segment(this->control_iteration_count,this->N_p),
                     dY_ref_K_aux2.segment(this->control_iteration_count,this->N_p);               
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
      }
    ///////////////////////////////////////////////////////////////////////////////////////////
    //PRIVATE method for computing the QP solution V_k ( [\ddot{q}_p;lambda] )

    private: void V_k_calc(){

      // this->V_k=solveQP_hpipm(this->H, this->f_k, this->A_eq_k, this->b_eq_k, this->A_ineq_k,this->b_neg_inf_QP, this->b_ineq_k); 
      this->V_k=solveQP_hpipm(this->H, this->f_k, this->A_eq_k, this->b_eq_k, this->A_ineq_k,this->b_ineq_min_k, this->b_ineq_max_k); 
      (this->VV_k).row(this->control_iteration_count-1)<<(this->V_k).transpose();//assigning inputs

      this->v_k<<this->V_k(0),this->V_k(this->N_p);
      this->v1_k_m1=this->v_k(0);
      this->v2_k_m1=this->v_k(1);

    }
    //////////////////////////////////////////////////////////////////////////////////////////
   
    private: void lambda_b_T_zeta_u_calc(){

      this->m_p_eq=this->robot_inertial_par(9)+2*(this->robot_inertial_par(6)+this->robot_inertial_par(7)+this->robot_inertial_par(8)); // m_body+2*(m_ah+m_d+m_c);
      this->l_p_eq_k_calc();
      this->I_p_eq_k_calc();

      this->lambda_u_k=sqrt(this->m_p_eq*this->l_p_eq_k*this->g/(this->I_p_eq_k+this->m_p_eq*this->l_p_eq_k*(this->l_p_eq_k+this->robot_inertial_par(0))));
      double b_u_k_aux=0;
      b_u_k_aux=this->robot_inertial_par(0)*(this->m_p_eq*(this->l_p_eq_k+this->robot_inertial_par(0))+this->robot_inertial_par(0)*(this->robot_inertial_par(5)+this->robot_inertial_par(10)/pow(this->robot_inertial_par(0),2)))/(2.0*this->lambda_u_k*(this->I_p_eq_k+this->m_p_eq*this->l_p_eq_k*(this->l_p_eq_k+this->robot_inertial_par(0))));
      this->b_u_k=b_u_k_aux;
      this->T_k<<1.0/2.0,-1.0/(2.0*this->lambda_u_k),1.0/2.0,1.0/(2.0*this->lambda_u_k);

      this->phi_r_st_eq_k_calc();

      this->theta_r_theta_r_dot_k<<this->q_k(3)-this->phi_r_st_eq_k,this->q_k(8)-this->dphi_r_st_eq_k_dbeta*this->q_k(9);
      this->zeta_k=this->T_k*this->theta_r_theta_r_dot_k;
      this->zeta_u_k<<this->zeta_k(1);
    }
    //////////////////////////////////////////////////////////////////////////////////////////
    private: void f_k_calc(){
      VectorXf q1_bar_k=VectorXf::Zero(2);
      VectorXf q2_bar_k=VectorXf::Zero(2);
      q1_bar_k(0)=this->q_k(2);
      q1_bar_k(1)=this->q_k(7);
      q2_bar_k(0)=this->q_k(4);
      q2_bar_k(1)=this->q_k(9);

      this->f1_diff_k(0)=-this->v1_k_m1;
      this->f2_diff_k(this->N_p)=-this->v2_k_m1;

      this->f_k=(this->cost_weights_y1(0)*(this->Omega_pos*q1_bar_k-this->Y_ref_k.col(0)).transpose()*this->Gamma_pos1).transpose()+
                 (this->cost_weights_y1(1)*(this->Omega_vel*q1_bar_k-this->dY_ref_k.col(0)).transpose()*this->Gamma_vel1).transpose()+ this->cost_weights_y1(4)*this->f1_diff_k+
                (this->cost_weights_y2(0)*(this->Omega_pos*q2_bar_k-this->Y_ref_k.col(1)).transpose()*this->Gamma_pos2).transpose()+
                (this->cost_weights_y2(1)*(this->Omega_vel*q2_bar_k-this->dY_ref_k.col(1)).transpose()*this->Gamma_vel2).transpose()+ this->cost_weights_y2(4)*this->f2_diff_k;
    }
    //////////////////////////////////////////////////////////////////////////////////////////
    private: void M_C_bar_k_calc(){
          this->M_hat_k_calc();
          this->C_hat_k_calc();

          double m11=M_hat_k(0,0);
          double m12=M_hat_k(0,1);
          double m13=M_hat_k(0,2);
          double m22=M_hat_k(1,1);
          double m23=M_hat_k(1,2);
          double m33=M_hat_k(2,2);
          double c1=C_hat_k(0);
          double c2=C_hat_k(1);
          double c3=C_hat_k(2);

          double h1;
          double h2;
          this->knee_jacobians_h1_h2_calc(this->q_k(4), &h1, &h2);

          this->M_bar_k<<(m11+(m12*(m11-m12))/(m22-m12)),(m13+(m12*(m13-m23))/(m22-m12)),
                       (m13+m11*h1+((m23+h1*m12)*(m11-m12))/(m22-m12)),(m33+m13*h1+((m23+h1*m12)*(m13-m23))/(m22-m12));
          this->C_bar_k<<c1+m12*(c1-c2)/(m22-m12),
                       c3+h1*c1+(m23+h1*m12)*(c1-c2)/(m22-m12);
    }
    //////////////////////////////////////////////////////////////////////////////////////////
    private: void V_min_max_k_calc(){

      VectorXf u_bar_lim=2*this->u_lim;
      MatrixXf M_bar_inv=(this->M_bar_k).inverse();
      this->v_min_k=-M_bar_inv*(u_bar_lim+this->C_bar_k);
      this->v_max_k=M_bar_inv*(u_bar_lim-this->C_bar_k);

      // this->v1_min_k=v_min(0);
      // this->v1_max_k=v_max(0);
      // this->v2_min_k=v_min(1);
      // this->v2_max_k=v_max(1);
      
      // this->V_min_k.segment(0,this->N_p)= this->v1_min_k*MatrixXf::Constant(this->N_p,1,1.0);
      // this->V_min_k.segment(this->N_p,this->N_p)= this->v2_min_k*MatrixXf::Constant(this->N_p,1,1.0);
      
      // this->V_max_k.segment(0,this->N_p)= this->v1_max_k*MatrixXf::Constant(this->N_p,1,1.0);
      // this->V_max_k.segment(this->N_p,this->N_p)= this->v2_max_k*MatrixXf::Constant(this->N_p,1,1.0);
    }
    
    //Private "low-level" methods
    private: 

      void M_hat_k_calc() {

        double phi_r=this->q_k(3);
        double beta=this->q_k(4);

        double A;
        double B;
        double a_tmp_tmp;
        double b_a_tmp_tmp;
        double c_a_tmp_tmp;
        double cos_gamma;
        double h1;
        double h2;
        double sin_gamma;
        double t11;
        double t12;
        double t15;
        double t16;
        double t25;
        double t26;
        double t29;
        double t4;
        double t45;
        double t46;
        double t48;
        double t5;
        double t50;
        double t52;
        double t53;
        double t6;
        double t65;
        double t72;
        double t73;
        double t74;
        double t75;
        double t9;
        //     This function was generated by the Symbolic Math Toolbox version 8.7.
        //     30-Sep-2021 13:00:33
        t15 = cos(beta);
        t16 = sin(beta);
        t9 = this->robot_dim(3) * t15;
        t48 = this->robot_dim(1) - t9;
        a_tmp_tmp = this->robot_dim(3) * this->robot_dim(3);
        sin_gamma = 2.0 * t15;
        b_a_tmp_tmp = this->robot_dim(0) * this->robot_dim(0);
        t11 = this->robot_dim(2) * this->robot_dim(2);
        A = sin_gamma * this->robot_dim(1);
        t12 = this->robot_dim(1) * this->robot_dim(1);
        B = A * this->robot_dim(3);
        c_a_tmp_tmp = -b_a_tmp_tmp + t12;
        t45 = ((c_a_tmp_tmp - B) + t11) + a_tmp_tmp;
        t52 = t48 * t48;
        t65 = t16 * t16;
        t72 = a_tmp_tmp * t65;
        t26 = t72 / t52;
        t29 = 2.0 * b_a_tmp_tmp;
        t4 = pow(t48, 3.0);
        t5 = 4.0 * t11;
        t73 = (t12 - B) + a_tmp_tmp;
        t74 = 2.0 * a_tmp_tmp;
        t46 = this->robot_dim(1) * a_tmp_tmp;
        t53 = this->robot_dim(2) * t4;
        t29 = (((((((((pow(this->robot_dim(0), 4.0) +
                      sin_gamma * b_a_tmp_tmp * this->robot_dim(1) *
                          this->robot_dim(3)) -
                      t29 * t11) -
                    t29 * a_tmp_tmp) -
                    pow(this->robot_dim(1), 4.0)) +
                  sin_gamma * pow(this->robot_dim(1), 3.0) *
                      this->robot_dim(3)) +
                  A * t11 * this->robot_dim(3)) -
                A * pow(this->robot_dim(3), 3.0)) +
                pow(this->robot_dim(2), 4.0)) -
              2.0 * t11 * a_tmp_tmp) +
              pow(this->robot_dim(3), 4.0);
        cos_gamma = t5 * t4;
        t6 = this->robot_dim(3) * t16;
        B = this->robot_dim(2) * t52;
        sin_gamma = 2.0 * this->robot_dim(2) * t52;
        A = sqrt((t26 - t45 * t45 / (t5 * t52)) + 1.0);
        t75 = this->robot_dim(3) - this->robot_dim(1) * t15;
        t50 = A + t6 * t45 / sin_gamma;
        t25 = 2.0 * this->robot_dim(1) * this->robot_dim(3) * t15;
        h1 = -(2.0 * this->robot_dim(0) *
              ((B *
                    (((t9 * t45 / sin_gamma + t46 * t65 / B) - t72 * t45 / t53) +
                      t6 * t29 / (cos_gamma * A)) /
                    (this->robot_dim(0) * t73) -
                t9 / this->robot_dim(0)) +
                2.0 * this->robot_dim(2) * a_tmp_tmp * t16 * t75 * t50 /
                    (this->robot_dim(0) * t4 * ((t26 + 1.0) * (t26 + 1.0)))) *
              t73) /
            (t48 *
              (((((((t74 * (t15 * t15) + t74 * t65) + b_a_tmp_tmp) + t12) - t11) -
                a_tmp_tmp) -
                t25) +
              2.0 * this->robot_dim(2) * this->robot_dim(3) * t16 * A));
        h2 =
            -(t52 *
                  (((this->robot_dim(3) * t15 *
                        (((c_a_tmp_tmp -
                            2.0 * t15 * this->robot_dim(1) * this->robot_dim(3)) +
                          t11) +
                          a_tmp_tmp) /
                        (2.0 * this->robot_dim(2) * (t48 * t48)) +
                    t46 * (t16 * t16) / (this->robot_dim(2) * (t48 * t48))) -
                    a_tmp_tmp * (t16 * t16) * t45 / t53) +
                  this->robot_dim(3) * t16 * t29 /
                      (cos_gamma * sqrt((a_tmp_tmp * (t16 * t16) / (t48 * t48) -
                                              t45 * t45 / (t5 * (t48 * t48))) +
                                              1.0))) /
                  t73 +
              t74 * t16 * t75 * t50 / (t4 * ((t26 + 1.0) * (t26 + 1.0)))) /
            (t45 / (2.0 * this->robot_dim(2) * t48) - t6 * t48 * t50 / t73);
        //  DESCRIPTION:
        //  computes, given a linkage beta, the corresponding alpha and gamma,
        //  employing proper kinematic relationships and controls
        //  PARAMETERS:
        //  beta-->linkage angle (see the system's drawing for clarifications)
        B = t9 - this->robot_dim(1);
        A = ((((b_a_tmp_tmp - t12) - t11) - a_tmp_tmp) + t25) /
            (2.0 * this->robot_dim(2) * B);
        B = t6 / B;
        sin_gamma = B * B + 1.0;
        cos_gamma = (-A * B + sqrt(sin_gamma - A * A)) / sin_gamma;
        sin_gamma = A + B * cos_gamma;
        t75 = atan2(sin_gamma, cos_gamma);
        B = this->robot_dim(3) / this->robot_dim(0);
        A = this->robot_dim(2) / this->robot_dim(0);
        B = atan2(B * t16 - A * cos_gamma,
                          (this->robot_dim(1) / this->robot_dim(0) - B * t15) -
                              A * sin_gamma);
        t4 = this->robot_inertial_par(12) * 2.0;
        t5 = h1 * h1;
        t6 = h2 * h2;
        t9 = this->robot_inertial_par(3) * this->robot_inertial_par(3);
        t11 = this->robot_inertial_par(2) * this->robot_inertial_par(2);
        t12 = this->robot_dim(4) * this->robot_dim(4);
        t15 = cos(B + beta);
        t16 = sin(B + t75);
        sin_gamma = cos((beta + phi_r) + this->robot_dim(5));
        t25 = sin((t75 + phi_r) + this->robot_dim(5));
        t26 = this->robot_inertial_par(9) * a_tmp_tmp;
        t29 = this->robot_inertial_par(7) *
              (this->robot_inertial_par(1) * this->robot_inertial_par(1)) * 2.0;
        t45 = this->robot_dim(3) * this->robot_dim(4) * this->robot_inertial_par(9) *
              t15;
        t46 = this->robot_dim(3) * this->robot_inertial_par(9) * sin_gamma;
        A = cos((phi_r + -B) + this->robot_dim(5));
        t48 = sin((-B + this->robot_dim(5)) + this->robot_dim(6));
        c_a_tmp_tmp = this->robot_inertial_par(1) * this->robot_dim(4) *
                      this->robot_inertial_par(7) * t15;
        cos_gamma = c_a_tmp_tmp * 2.0;
        t50 = this->robot_inertial_par(1) * this->robot_inertial_par(7) * sin_gamma *
              2.0;
        t53 = this->robot_dim(3) * this->robot_inertial_par(4) *
              this->robot_inertial_par(9) *
              sin((beta + this->robot_dim(5)) + this->robot_dim(6));
        t52 = h1 * t45;
        t65 = this->robot_dim(4) * this->robot_inertial_par(9) * A;
        t72 = this->robot_dim(0) * this->robot_inertial_par(8) * A * 2.0;
        t73 = this->robot_inertial_par(2) * this->robot_inertial_par(6) * A * 2.0;
        t74 = this->robot_dim(4) * this->robot_inertial_par(8) * A * 2.0;
        t75 = this->robot_dim(4) * this->robot_inertial_par(7) * A * 2.0;
        B = h1 * this->robot_inertial_par(8);
        sin_gamma = h1 * this->robot_dim(0);
        A = h2 * this->robot_inertial_par(3);
        cos_gamma =
            ((((((((((((((((((((t4 + this->robot_inertial_par(13) * h2 * 2.0) +
                              -(this->robot_inertial_par(11) * h1 * 2.0)) +
                              t26) +
                            t29) +
                            h2 * this->robot_inertial_par(8) * t9 * 2.0) +
                          -(sin_gamma * this->robot_dim(4) *
                            this->robot_inertial_par(8) * 4.0)) +
                          -(B * b_a_tmp_tmp * 2.0)) +
                        -(h1 * this->robot_inertial_par(6) * t11 * 2.0)) +
                        -(B * t12 * 2.0)) +
                      -(h1 * this->robot_inertial_par(9) * t12)) +
                      -(h1 * this->robot_inertial_par(7) * t12 * 2.0)) +
                    t52) +
                    -t45) +
                  -cos_gamma) +
                  h1 * cos_gamma) +
                h2 * this->robot_dim(0) * this->robot_inertial_par(3) *
                    this->robot_inertial_par(8) * t16 * 2.0) +
                A * this->robot_dim(4) * this->robot_inertial_par(8) * t16 * 2.0) +
              -(sin_gamma * this->robot_inertial_par(3) *
                this->robot_inertial_par(8) * t16 * 2.0)) +
              -(h1 * this->robot_inertial_par(3) * this->robot_dim(4) *
                this->robot_inertial_par(8) * t16 * 2.0)) +
            -t53) +
            -(h1 * this->robot_inertial_par(4) * this->robot_dim(4) *
              this->robot_inertial_par(9) * t48);
        sin_gamma =
            this->robot_inertial_par(0) *
            (((((((t46 + t50) + A * this->robot_inertial_par(8) * t25 * -2.0) +
                h1 * t65) +
                h1 * t72) +
              h1 * t73) +
              h1 * t74) +
            h1 * t75);
        B = -(this->robot_inertial_par(0) *
              ((((((((this->robot_inertial_par(3) * this->robot_inertial_par(8) *
                          t25 * 2.0 +
                      -t46) +
                    -t50) +
                    t65) +
                  -(this->robot_inertial_par(4) * this->robot_inertial_par(9) *
                    sin(phi_r + -this->robot_dim(6)))) +
                  t72) +
                t73) +
                t74) +
              t75));
        this->M_hat_k(0,0) = this->robot_inertial_par(10) * 2.0 +
                  this->robot_inertial_par(0) * this->robot_inertial_par(0) *
                      ((((this->robot_inertial_par(6) * 2.0 +
                          this->robot_inertial_par(8) * 2.0) +
                          this->robot_inertial_par(9)) +
                        this->robot_inertial_par(7) * 2.0) +
                        this->robot_inertial_par(5) * 2.0);
        this->M_hat_k(1,0) = B;
        this->M_hat_k(2,0) = sin_gamma;
        this->M_hat_k(0,1) = B;
        A = this->robot_dim(0) * this->robot_dim(4) * this->robot_inertial_par(8);
        this->M_hat_k(1,1) = ((((((((((((((((((this->robot_inertial_par(11) * 2.0 +
                                    this->robot_inertial_par(13) * 2.0) +
                                    this->robot_inertial_par(14)) +
                                  t4) +
                                  t26) +
                                t29) -
                                t45 * 2.0) -
                              t53 * 2.0) +
                              this->robot_inertial_par(6) * t11 * 2.0) +
                            this->robot_inertial_par(8) * b_a_tmp_tmp * 2.0) +
                            this->robot_inertial_par(8) * t9 * 2.0) +
                          this->robot_inertial_par(8) * t12 * 2.0) +
                          this->robot_inertial_par(9) * t12) +
                        this->robot_inertial_par(7) * t12 * 2.0) +
                        this->robot_inertial_par(4) * this->robot_inertial_par(4) *
                            this->robot_inertial_par(9)) +
                      A * 4.0) +
                      this->robot_dim(0) * this->robot_inertial_par(3) *
                          this->robot_inertial_par(8) * t16 * 4.0) +
                    this->robot_inertial_par(3) * this->robot_dim(4) *
                        this->robot_inertial_par(8) * t16 * 4.0) -
                    c_a_tmp_tmp * 4.0) +
                  this->robot_inertial_par(4) * this->robot_dim(4) *
                      this->robot_inertial_par(9) * t48 * 2.0;
        this->M_hat_k(2,1) = cos_gamma;
        this->M_hat_k(0,2) = sin_gamma;
        this->M_hat_k(1,2) = cos_gamma;
        sin_gamma = this->robot_inertial_par(8) * t5;
        B = h1 * h2;
        this->M_hat_k(2,2) = ((((((((((((((t4 + t26) + t29) + t52 * 2.0) +
                              this->robot_inertial_par(11) * t5 * 2.0) +
                            this->robot_inertial_par(13) * t6 * 2.0) +
                            this->robot_inertial_par(6) * t5 * t11 * 2.0) +
                          sin_gamma * b_a_tmp_tmp * 2.0) +
                          this->robot_inertial_par(8) * t6 * t9 * 2.0) +
                        sin_gamma * t12 * 2.0) +
                        this->robot_inertial_par(9) * t5 * t12) +
                      this->robot_inertial_par(7) * t5 * t12 * 2.0) +
                      A * t5 * 4.0) +
                    h1 * this->robot_inertial_par(1) * this->robot_dim(4) *
                        this->robot_inertial_par(7) * t15 * 4.0) -
                    B * this->robot_dim(0) * this->robot_inertial_par(3) *
                        this->robot_inertial_par(8) * t16 * 4.0) -
                  B * this->robot_inertial_par(3) * this->robot_dim(4) *
                      this->robot_inertial_par(8) * t16 * 4.0;
      }

      void C_hat_k_calc(){

          double phi_r=this->q_k(3);
          double beta=this->q_k(4);
          double phi_r_dot=this->q_k(8);
          double beta_dot=this->q_k(9);
                                  
          double A;
          double B;
          double a_tmp;
          double a_tmp_tmp;
          double b_a_tmp;
          double b_gamma;
          double b_h1_tmp;
          double c_a_tmp;
          double c_h1_tmp;
          double cos_gamma;
          double d_a_tmp;
          double d_h1_tmp;
          double e_a_tmp;
          double e_h1_tmp;
          double f_a_tmp;
          double f_h1_tmp;
          double h1;
          double h1_tmp;
          double h2;
          double sin_gamma;
          double t15;
          double t16;
          double t20;
          double t22;
          double t23;
          double t24;
          double t4;
          double t5;
          double t6;
          double t7;
          double t8;
          double t9;
          //     This function was generated by the Symbolic Math Toolbox version 8.7.
          //     30-Sep-2021 13:00:34
          t23 = cos(beta);
          t24 = sin(beta);
          t7 = this->robot_dim(3) * t23;
          a_tmp = this->robot_dim(1) - t7;
          t8 = this->robot_dim(3) * this->robot_dim(3);
          sin_gamma = 2.0 * t23;
          t9 = this->robot_dim(0) * this->robot_dim(0);
          t15 = this->robot_dim(2) * this->robot_dim(2);
          A = sin_gamma * this->robot_dim(1);
          t16 = this->robot_dim(1) * this->robot_dim(1);
          B = A * this->robot_dim(3);
          a_tmp_tmp = -t9 + t16;
          b_a_tmp = ((a_tmp_tmp - B) + t15) + t8;
          c_a_tmp = a_tmp * a_tmp;
          d_a_tmp = t24 * t24;
          e_a_tmp = t8 * d_a_tmp;
          f_a_tmp = e_a_tmp / c_a_tmp;
          h1_tmp = 2.0 * t9;
          t4 = pow(a_tmp, 3.0);
          t5 = 4.0 * t15;
          b_h1_tmp = (t16 - B) + t8;
          c_h1_tmp = 2.0 * t8;
          d_h1_tmp = this->robot_dim(1) * t8;
          e_h1_tmp = this->robot_dim(2) * t4;
          h1_tmp = (((((((((pow(this->robot_dim(0), 4.0) +
                            sin_gamma * t9 * this->robot_dim(1) * this->robot_dim(3)) -
                          h1_tmp * t15) -
                          h1_tmp * t8) -
                        pow(this->robot_dim(1), 4.0)) +
                        sin_gamma * pow(this->robot_dim(1), 3.0) *
                            this->robot_dim(3)) +
                      A * t15 * this->robot_dim(3)) -
                      A * pow(this->robot_dim(3), 3.0)) +
                    pow(this->robot_dim(2), 4.0)) -
                    2.0 * t15 * t8) +
                  pow(this->robot_dim(3), 4.0);
          cos_gamma = t5 * t4;
          t6 = this->robot_dim(3) * t24;
          B = this->robot_dim(2) * c_a_tmp;
          sin_gamma = 2.0 * this->robot_dim(2) * c_a_tmp;
          A = sqrt((f_a_tmp - b_a_tmp * b_a_tmp / (t5 * c_a_tmp)) + 1.0);
          t20 = this->robot_dim(3) - this->robot_dim(1) * t23;
          t22 = A + t6 * b_a_tmp / sin_gamma;
          f_h1_tmp = 2.0 * this->robot_dim(1) * this->robot_dim(3) * t23;
          h1 =
              -(2.0 * this->robot_dim(0) *
                ((B *
                      (((t7 * b_a_tmp / sin_gamma + d_h1_tmp * d_a_tmp / B) -
                        e_a_tmp * b_a_tmp / e_h1_tmp) +
                      t6 * h1_tmp / (cos_gamma * A)) /
                      (this->robot_dim(0) * b_h1_tmp) -
                  t7 / this->robot_dim(0)) +
                2.0 * this->robot_dim(2) * t8 * t24 * t20 * t22 /
                    (this->robot_dim(0) * t4 * ((f_a_tmp + 1.0) * (f_a_tmp + 1.0)))) *
                b_h1_tmp) /
              (a_tmp *
              (((((((c_h1_tmp * (t23 * t23) + c_h1_tmp * d_a_tmp) + t9) + t16) - t15) -
                  t8) -
                f_h1_tmp) +
                2.0 * this->robot_dim(2) * this->robot_dim(3) * t24 * A));
          h2 = -(c_a_tmp *
                    (((this->robot_dim(3) * t23 *
                            (((a_tmp_tmp -
                              2.0 * t23 * this->robot_dim(1) * this->robot_dim(3)) +
                              t15) +
                            t8) /
                            (2.0 * this->robot_dim(2) * (a_tmp * a_tmp)) +
                        d_h1_tmp * (t24 * t24) /
                            (this->robot_dim(2) * (a_tmp * a_tmp))) -
                      t8 * (t24 * t24) * b_a_tmp / e_h1_tmp) +
                      this->robot_dim(3) * t24 * h1_tmp /
                          (cos_gamma *
                          sqrt((t8 * (t24 * t24) / (a_tmp * a_tmp) -
                                      b_a_tmp * b_a_tmp / (t5 * (a_tmp * a_tmp))) +
                                    1.0))) /
                    b_h1_tmp +
                c_h1_tmp * t24 * t20 * t22 /
                    (t4 * ((f_a_tmp + 1.0) * (f_a_tmp + 1.0)))) /
              (b_a_tmp / (2.0 * this->robot_dim(2) * a_tmp) -
                t6 * a_tmp * t22 / b_h1_tmp);
          //  DESCRIPTION:
          //  computes, given a linkage beta, the corresponding alpha and gamma,
          //  employing proper kinematic relationships and controls
          //  PARAMETERS:
          //  beta-->linkage angle (see the system's drawing for clarifications)
          B = t7 - this->robot_dim(1);
          A = ((((t9 - t16) - t15) - t8) + f_h1_tmp) / (2.0 * this->robot_dim(2) * B);
          B = t6 / B;
          sin_gamma = B * B + 1.0;
          cos_gamma = (-A * B + sqrt(sin_gamma - A * A)) / sin_gamma;
          sin_gamma = A + B * cos_gamma;
          b_gamma = atan2(sin_gamma, cos_gamma);
          B = this->robot_dim(3) / this->robot_dim(0);
          A = this->robot_dim(2) / this->robot_dim(0);
          h1_tmp = atan2(B * t24 - A * cos_gamma,
                                (this->robot_dim(1) / this->robot_dim(0) - B * t23) -
                                    A * sin_gamma);
          t4 = beta_dot * beta_dot;
          t5 = h1 * h1;
          t6 = h2 * h2;
          t7 = phi_r_dot * phi_r_dot;
          t8 = cos(h1_tmp + b_gamma);
          t9 = sin(h1_tmp + beta);
          t15 = cos((beta + this->robot_dim(5)) + this->robot_dim(6));
          t16 = cos((b_gamma + phi_r) + this->robot_dim(5));
          B = sin((beta + phi_r) + this->robot_dim(5));
          t20 = cos(phi_r + -this->robot_dim(6));
          t22 = this->robot_dim(3) * this->robot_inertial_par(9) * B;
          t23 = cos((-h1_tmp + this->robot_dim(5)) + this->robot_dim(6));
          t24 = sin((phi_r + -h1_tmp) + this->robot_dim(5));
          B *= this->robot_inertial_par(1) * this->robot_inertial_par(7);
          sin_gamma = B * 2.0;
          d_a_tmp =
              this->robot_dim(3) * this->robot_dim(4) * this->robot_inertial_par(9);
          A = d_a_tmp * t4;
          e_a_tmp = this->robot_inertial_par(1) * this->robot_dim(4) *
                    this->robot_inertial_par(7);
          cos_gamma = e_a_tmp * t4;
          b_h1_tmp = this->robot_dim(0) * this->robot_inertial_par(8);
          c_h1_tmp = -t22 + -sin_gamma;
          d_h1_tmp = this->robot_inertial_par(2) * this->robot_inertial_par(6);
          e_h1_tmp = this->robot_inertial_par(3) * this->robot_inertial_par(8);
          f_h1_tmp = this->robot_dim(4) * this->robot_inertial_par(8);
          a_tmp = this->robot_dim(4) * this->robot_inertial_par(9);
          a_tmp_tmp = this->robot_dim(4) * this->robot_inertial_par(7);
          b_a_tmp = h1 * this->robot_dim(4);
          c_a_tmp = h1 * this->robot_dim(0);
          f_a_tmp = h2 * this->robot_inertial_par(3);
          this->C_hat_k(0) =
              this->robot_inertial_par(0) *
              ((t7 * (((((((c_h1_tmp + b_h1_tmp * t24 * 2.0) + d_h1_tmp * t24 * 2.0) -
                          e_h1_tmp * t16 * 2.0) +
                        this->robot_inertial_par(4) * this->robot_inertial_par(9) *
                            t20) +
                        f_h1_tmp * t24 * 2.0) +
                      a_tmp * t24) +
                      a_tmp_tmp * t24 * 2.0) +
                t4 * ((((((c_h1_tmp + b_h1_tmp * t5 * t24 * 2.0) +
                          d_h1_tmp * t5 * t24 * 2.0) -
                        e_h1_tmp * t6 * t16 * 2.0) +
                        f_h1_tmp * t5 * t24 * 2.0) +
                      a_tmp * t5 * t24) +
                      a_tmp_tmp * t5 * t24 * 2.0)) -
              beta_dot * phi_r_dot *
                  (((((((t22 * 2.0 + B * 4.0) +
                        c_a_tmp * this->robot_inertial_par(8) * t24 * 4.0) +
                        h1 * this->robot_inertial_par(2) * this->robot_inertial_par(6) *
                            t24 * 4.0) +
                      f_a_tmp * this->robot_inertial_par(8) * t16 * 4.0) +
                      b_a_tmp * this->robot_inertial_par(8) * t24 * 4.0) +
                    b_a_tmp * this->robot_inertial_par(9) * t24 * 2.0) +
                    b_a_tmp * this->robot_inertial_par(7) * t24 * 4.0));
          b_h1_tmp = this->g * this->robot_dim(4);
          c_h1_tmp = this->robot_dim(0) * this->robot_inertial_par(3) *
                    this->robot_inertial_par(8) * t4;
          d_h1_tmp = this->robot_inertial_par(3) * this->robot_dim(4) *
                    this->robot_inertial_par(8) * t4;
          e_h1_tmp = beta_dot * this->robot_dim(3);
          f_h1_tmp = beta_dot * h1;
          a_tmp = beta_dot * h2;
          a_tmp_tmp = ((this->g * t22 + this->g * sin_gamma) + -(A * t5 * t9)) +
                      -(cos_gamma * t5 * t9 * 2.0);
          b_a_tmp = this->robot_dim(3) * this->robot_inertial_par(4) *
                    this->robot_inertial_par(9);
          this->C_hat_k(1) =
              ((((((((((((((((((((((((a_tmp_tmp - this->g * this->robot_dim(0) *
                                                      this->robot_inertial_par(8) *
                                                      t24 * 2.0) -
                                    this->g * this->robot_inertial_par(2) *
                                        this->robot_inertial_par(6) * t24 * 2.0) +
                                    this->g * this->robot_inertial_par(3) *
                                        this->robot_inertial_par(8) * t16 * 2.0) -
                                  this->g * this->robot_inertial_par(4) *
                                      this->robot_inertial_par(9) * t20) -
                                  b_h1_tmp * this->robot_inertial_par(8) * t24 * 2.0) -
                                b_h1_tmp * this->robot_inertial_par(9) * t24) -
                                b_h1_tmp * this->robot_inertial_par(7) * t24 * 2.0) -
                              b_a_tmp * t4 * t15) +
                              A * t9) +
                            cos_gamma * t9 * 2.0) -
                            c_h1_tmp * t5 * t8 * 2.0) +
                          c_h1_tmp * t6 * t8 * 2.0) -
                          d_h1_tmp * t5 * t8 * 2.0) +
                        d_h1_tmp * t6 * t8 * 2.0) +
                        this->robot_inertial_par(4) * this->robot_dim(4) *
                            this->robot_inertial_par(9) * t4 * t5 * t23) -
                      e_h1_tmp * this->robot_inertial_par(4) *
                          this->robot_inertial_par(9) * phi_r_dot * t15 * 2.0) +
                      e_h1_tmp * this->robot_dim(4) * this->robot_inertial_par(9) *
                          phi_r_dot * t9 * 2.0) +
                    beta_dot * this->robot_inertial_par(1) * this->robot_dim(4) *
                        this->robot_inertial_par(7) * phi_r_dot * t9 * 4.0) +
                    f_h1_tmp * this->robot_dim(0) * this->robot_inertial_par(3) *
                        this->robot_inertial_par(8) * phi_r_dot * t8 * 4.0) +
                  a_tmp * this->robot_dim(0) * this->robot_inertial_par(3) *
                      this->robot_inertial_par(8) * phi_r_dot * t8 * 4.0) +
                  f_h1_tmp * this->robot_inertial_par(3) * this->robot_dim(4) *
                      this->robot_inertial_par(8) * phi_r_dot * t8 * 4.0) +
                a_tmp * this->robot_inertial_par(3) * this->robot_dim(4) *
                    this->robot_inertial_par(8) * phi_r_dot * t8 * 4.0) +
                f_h1_tmp * this->robot_dim(3) * this->robot_dim(4) *
                    this->robot_inertial_par(9) * phi_r_dot * t9 * 2.0) +
              f_h1_tmp * this->robot_inertial_par(1) * this->robot_dim(4) *
                  this->robot_inertial_par(7) * phi_r_dot * t9 * 4.0) +
              f_h1_tmp * this->robot_inertial_par(4) * this->robot_dim(4) *
                  this->robot_inertial_par(9) * phi_r_dot * t23 * -2.0;
          b_h1_tmp = h1_tmp * h1;
          c_h1_tmp = this->robot_knee_el_par[5] * h1;
          d_h1_tmp = this->robot_knee_el_par[6] * h2;
          e_h1_tmp = b_gamma * h2;
          f_h1_tmp = this->g * h1;
          a_tmp = f_h1_tmp * this->robot_dim(4);
          B = h1 * this->robot_dim(3) * this->robot_dim(4) *
              this->robot_inertial_par(9);
          sin_gamma = h1 * this->robot_inertial_par(1) * this->robot_dim(4) *
                      this->robot_inertial_par(7);
          c_a_tmp = c_a_tmp * this->robot_inertial_par(3) * this->robot_inertial_par(8);
          A = h2 * this->robot_dim(0) * this->robot_inertial_par(3) *
              this->robot_inertial_par(8);
          cos_gamma = h1 * this->robot_inertial_par(3) * this->robot_dim(4) *
                      this->robot_inertial_par(8);
          f_a_tmp = f_a_tmp * this->robot_dim(4) * this->robot_inertial_par(8);
          this->C_hat_k(2) =
              ((((((((((((((((((((((((((((((((a_tmp_tmp +
                                              h1_tmp * this->robot_knee_el_par[1] *
                                                  2.0) -
                                            this->robot_knee_el_par[1] *
                                                this->robot_knee_el_par[5] * 2.0) -
                                            this->robot_knee_el_par[1] *
                                                this->robot_knee_el_par[4] * 2.0) -
                                          this->robot_knee_el_par[3] *
                                              this->robot_knee_el_par[4] * 2.0) +
                                          beta * this->robot_knee_el_par[1] * 2.0) +
                                        beta * this->robot_knee_el_par[3] * 2.0) +
                                        b_h1_tmp * this->robot_knee_el_par[0] * 2.0) -
                                      c_h1_tmp * this->robot_knee_el_par[0] * 2.0) +
                                      b_h1_tmp * this->robot_knee_el_par[1] * 2.0) +
                                    h1_tmp * h2 * this->robot_knee_el_par[0] * 2.0) -
                                    c_h1_tmp * this->robot_knee_el_par[1] * 2.0) -
                                  this->robot_knee_el_par[5] * h2 *
                                      this->robot_knee_el_par[0] * 2.0) -
                                  this->robot_knee_el_par[4] * h1 *
                                      this->robot_knee_el_par[1] * 2.0) +
                                beta * h1 * this->robot_knee_el_par[1] * 2.0) -
                                this->robot_knee_el_par[6] * h1 *
                                    this->robot_knee_el_par[0] * 2.0) -
                              d_h1_tmp * this->robot_knee_el_par[0] * 2.0) -
                              d_h1_tmp * this->robot_knee_el_par[2] * 2.0) +
                            b_gamma * h1 * this->robot_knee_el_par[0] * 2.0) +
                            e_h1_tmp * this->robot_knee_el_par[0] * 2.0) +
                          e_h1_tmp * this->robot_knee_el_par[2] * 2.0) +
                          f_h1_tmp * this->robot_dim(0) * this->robot_inertial_par(8) *
                              t24 * 2.0) +
                        f_h1_tmp * this->robot_inertial_par(2) *
                            this->robot_inertial_par(6) * t24 * 2.0) +
                        this->g * h2 * this->robot_inertial_par(3) *
                            this->robot_inertial_par(8) * t16 * 2.0) +
                      a_tmp * this->robot_inertial_par(8) * t24 * 2.0) +
                      a_tmp * this->robot_inertial_par(9) * t24) +
                    a_tmp * this->robot_inertial_par(7) * t24 * 2.0) +
                    b_a_tmp * t7 * t15) -
                  d_a_tmp * t7 * t9) -
                  e_a_tmp * t7 * t9 * 2.0) -
                c_a_tmp * t7 * t8 * 2.0) -
                A * t7 * t8 * 2.0) -
              cos_gamma * t7 * t8 * 2.0) +
              (((((((((f_a_tmp * t7 * t8 * -2.0 - B * t4 * t9) - B * t7 * t9) -
                    sin_gamma * t4 * t9 * 2.0) -
                    sin_gamma * t7 * t9 * 2.0) +
                  h1 * this->robot_inertial_par(4) * this->robot_dim(4) *
                      this->robot_inertial_par(9) * t7 * t23) -
                  c_a_tmp * t4 * t6 * t8 * 2.0) -
                A * t4 * t5 * t8 * 2.0) -
                cos_gamma * t4 * t6 * t8 * 2.0) -
              f_a_tmp * t4 * t5 * t8 * 2.0);
      }

      void l_p_eq_k_calc(){//equivalent pendulum length
        
        double beta=this->q_k(4);

        double A;
        double A_tmp;
        double B;
        double B_tmp;
        double alpha;
        double b_A_tmp;
        double b_d_CoM_partial_tmp;
        double b_gamma;
        double c_A_tmp;
        double c_d_CoM_partial_tmp;
        double cos_gamma;
        double d_CoM_partial_tmp;
        double sin_gamma;
        double t2;
        double t6;
        double t7;
        double t8;
        //  DESCRIPTION:
        //  computes, given a linkage beta, the corresponding alpha and gamma,
        //  employing proper kinematic relationships and controls
        //  PARAMETERS:
        //  beta-->linkage angle (see the system's drawing for clarifications)
        A_tmp = cos(beta);
        b_A_tmp = this->robot_dim(0) * this->robot_dim(0);
        c_A_tmp = this->robot_dim(3) * this->robot_dim(3);
        A = ((((b_A_tmp - this->robot_dim(1) * this->robot_dim(1)) -
              this->robot_dim(2) * this->robot_dim(2)) -
              c_A_tmp) +
            2.0 * this->robot_dim(1) * this->robot_dim(3) * A_tmp) /
            (2.0 * this->robot_dim(2) *
            (this->robot_dim(3) * A_tmp - this->robot_dim(1)));
        B_tmp = sin(beta);
        B = this->robot_dim(3) * B_tmp /
            (this->robot_dim(3) * cos(beta) - this->robot_dim(1));
        sin_gamma = B * B + 1.0;
        cos_gamma = (-A * B + sqrt(sin_gamma - A * A)) / sin_gamma;
        sin_gamma = A + B * cos_gamma;
        b_gamma = atan2(sin_gamma, cos_gamma);
        A = this->robot_dim(3) / this->robot_dim(0);
        B = this->robot_dim(2) / this->robot_dim(0);
        alpha = atan2(A * B_tmp - B * cos_gamma,
                              (this->robot_dim(1) / this->robot_dim(0) - A * A_tmp) -
                                  B * sin_gamma);
        t2 = this->robot_dim(4) * this->robot_dim(4);
        t6 = this->robot_inertial_par(8) * this->robot_inertial_par(8);
        t7 = this->robot_inertial_par(9) * this->robot_inertial_par(9);
        t8 = this->robot_inertial_par(7) * this->robot_inertial_par(7);
        A = ((this->robot_inertial_par(9) + this->robot_inertial_par(6) * 2.0) +
            this->robot_inertial_par(8) * 2.0) +
            this->robot_inertial_par(7) * 2.0;
        A = 1.0 / (A * A);
        B = this->robot_dim(0) * this->robot_dim(4);
        sin_gamma = B * this->robot_inertial_par(8);
        cos_gamma = this->robot_inertial_par(2) * this->robot_dim(4) *
                    this->robot_inertial_par(6);
        A_tmp = this->robot_inertial_par(3) * this->robot_dim(4);
        B_tmp = A_tmp * this->robot_inertial_par(8);
        d_CoM_partial_tmp = this->robot_inertial_par(4) * this->robot_dim(4);
        b_d_CoM_partial_tmp = this->robot_dim(3) * this->robot_dim(4);
        c_d_CoM_partial_tmp = this->robot_inertial_par(1) * this->robot_dim(4);
        this->l_p_eq_k=sqrt(
            ((A * (((((((((((((((((((t2 * t6 * 4.0 + t2 * t7) + t2 * t8 * 4.0) +
                                  b_A_tmp * t6 * 4.0) +
                                  c_A_tmp * t7) +
                                this->robot_inertial_par(3) *
                                    this->robot_inertial_par(3) * t6 * 4.0) +
                                this->robot_inertial_par(1) *
                                    this->robot_inertial_par(1) * t8 * 4.0) +
                              this->robot_inertial_par(4) *
                                  this->robot_inertial_par(4) * t7) +
                              this->robot_inertial_par(2) *
                                  this->robot_inertial_par(2) *
                                  (this->robot_inertial_par(6) *
                                  this->robot_inertial_par(6)) *
                                  4.0) +
                            B * t6 * 8.0) +
                            this->robot_inertial_par(8) *
                                this->robot_inertial_par(9) * t2 * 4.0) +
                          this->robot_inertial_par(7) * this->robot_inertial_par(8) *
                              t2 * 8.0) +
                          this->robot_inertial_par(7) * this->robot_inertial_par(9) *
                              t2 * 4.0) +
                        this->robot_dim(0) * this->robot_inertial_par(2) *
                            this->robot_inertial_par(6) *
                            this->robot_inertial_par(8) * 8.0) +
                        sin_gamma * this->robot_inertial_par(9) * 4.0) +
                      this->robot_inertial_par(1) * this->robot_dim(3) *
                          this->robot_inertial_par(9) * this->robot_inertial_par(7) *
                          4.0) +
                      sin_gamma * this->robot_inertial_par(7) * 8.0) +
                    cos_gamma * this->robot_inertial_par(8) * 8.0) +
                    cos_gamma * this->robot_inertial_par(9) * 4.0) +
                  cos_gamma * this->robot_inertial_par(7) * 8.0) -
              A * sin((beta + this->robot_dim(5)) + this->robot_dim(6)) *
                  (this->robot_dim(3) * this->robot_inertial_par(4) * t7 * 2.0 +
                  4.0 * this->robot_inertial_par(1) * this->robot_inertial_par(4) *
                      this->robot_inertial_par(7) * this->robot_inertial_par(9))) +
            A * sin(beta - b_gamma) *
                (this->robot_dim(3) * this->robot_inertial_par(3) *
                      this->robot_inertial_par(8) * this->robot_inertial_par(9) *
                      4.0 +
                  this->robot_inertial_par(1) * this->robot_inertial_par(3) *
                      this->robot_inertial_par(8) * this->robot_inertial_par(7) *
                      8.0)) +
            (((A * sin(alpha + b_gamma) *
                  ((((this->robot_dim(0) * this->robot_inertial_par(3) * t6 * 8.0 +
                      A_tmp * t6 * 8.0) +
                      this->robot_inertial_par(2) * this->robot_inertial_par(3) *
                          this->robot_inertial_par(6) * this->robot_inertial_par(8) *
                          8.0) +
                    B_tmp * this->robot_inertial_par(9) * 4.0) +
                    B_tmp * this->robot_inertial_par(7) * 8.0) +
              A * sin((-alpha + this->robot_dim(5)) + this->robot_dim(6)) *
                  ((((d_CoM_partial_tmp * t7 * 2.0 +
                      this->robot_dim(0) * this->robot_inertial_par(4) *
                          this->robot_inertial_par(8) * this->robot_inertial_par(9) *
                          4.0) +
                      this->robot_inertial_par(2) * this->robot_inertial_par(4) *
                          this->robot_inertial_par(6) * this->robot_inertial_par(9) *
                          4.0) +
                    d_CoM_partial_tmp * this->robot_inertial_par(8) *
                        this->robot_inertial_par(9) * 4.0) +
                    d_CoM_partial_tmp * this->robot_inertial_par(9) *
                        this->robot_inertial_par(7) * 4.0)) -
              A * cos(alpha + beta) *
                  (((((((((b_d_CoM_partial_tmp * t7 * 2.0 +
                          c_d_CoM_partial_tmp * t8 * 8.0) +
                          this->robot_dim(0) * this->robot_dim(3) *
                              this->robot_inertial_par(8) *
                              this->robot_inertial_par(9) * 4.0) +
                        this->robot_dim(0) * this->robot_inertial_par(1) *
                            this->robot_inertial_par(8) *
                            this->robot_inertial_par(7) * 8.0) +
                        this->robot_inertial_par(2) * this->robot_dim(3) *
                            this->robot_inertial_par(6) *
                            this->robot_inertial_par(9) * 4.0) +
                      this->robot_inertial_par(1) * this->robot_inertial_par(2) *
                          this->robot_inertial_par(6) * this->robot_inertial_par(7) *
                          8.0) +
                      b_d_CoM_partial_tmp * this->robot_inertial_par(8) *
                          this->robot_inertial_par(9) * 4.0) +
                    c_d_CoM_partial_tmp * this->robot_inertial_par(8) *
                        this->robot_inertial_par(7) * 8.0) +
                    b_d_CoM_partial_tmp * this->robot_inertial_par(9) *
                        this->robot_inertial_par(7) * 4.0) +
                  c_d_CoM_partial_tmp * this->robot_inertial_par(9) *
                      this->robot_inertial_par(7) * 4.0)) -
            this->robot_inertial_par(3) * this->robot_inertial_par(4) *
                this->robot_inertial_par(8) * this->robot_inertial_par(9) * A *
                cos((b_gamma + this->robot_dim(5)) + this->robot_dim(6)) *
                4.0));
      }

      void I_p_eq_k_calc(){//equivalent pendulum Y-axis rotational moment of inertia
        
        double beta=this->q_k(4);

        double A;
        double A_tmp;
        double B;
        double B_tmp;
        double J_p_eq_tmp;
        double b_A_tmp;
        double b_J_p_eq_tmp;
        double b_gamma;
        double c_A_tmp;
        double c_J_p_eq_tmp;
        double cos_gamma;
        double d_J_p_eq_tmp;
        double e_J_p_eq_tmp;
        double sin_gamma;
        double t10;
        double t11;
        double t12;
        double t16;
        double t19;
        double t6;
        double t7;
        double t8;
        double t9;
        //     This function was generated by the Symbolic Math Toolbox version 8.7.
        //     09-May-2021 10:13:26
        //  DESCRIPTION:
        //  computes, given a linkage beta, the corresponding alpha and gamma,
        //  employing proper kinematic relationships and controls
        //  PARAMETERS:
        //  beta-->linkage angle (see the system's drawing for clarifications)
        A_tmp = cos(beta);
        b_A_tmp = this->robot_dim(0) * this->robot_dim(0);
        c_A_tmp = this->robot_dim(3) * this->robot_dim(3);
        A = ((((b_A_tmp - this->robot_dim(1) * this->robot_dim(1)) -
              this->robot_dim(2) * this->robot_dim(2)) -
              c_A_tmp) +
            2.0 * this->robot_dim(1) * this->robot_dim(3) * A_tmp) /
            (2.0 * this->robot_dim(2) *
            (this->robot_dim(3) * A_tmp - this->robot_dim(1)));
        B_tmp = sin(beta);
        B = this->robot_dim(3) * B_tmp /
            (this->robot_dim(3) * cos(beta) - this->robot_dim(1));
        sin_gamma = B * B + 1.0;
        cos_gamma = (-A * B + sqrt(sin_gamma - A * A)) / sin_gamma;
        sin_gamma = A + B * cos_gamma;
        b_gamma = atan2(sin_gamma, cos_gamma);
        A = this->robot_dim(3) / this->robot_dim(0);
        B = this->robot_dim(2) / this->robot_dim(0);
        A = atan2(A * B_tmp - B * cos_gamma,
                          (this->robot_dim(1) / this->robot_dim(0) - A * A_tmp) -
                              B * sin_gamma);
        t6 = this->robot_inertial_par(3) * this->robot_inertial_par(3);
        t7 = this->robot_inertial_par(1) * this->robot_inertial_par(1);
        t8 = this->robot_inertial_par(2) * this->robot_inertial_par(2);
        t9 = this->robot_inertial_par(4) * this->robot_inertial_par(4);
        t10 = this->robot_dim(4) * this->robot_dim(4);
        t11 = cos(A + beta);
        t12 = sin(A + b_gamma);
        t16 = sin((beta + this->robot_dim(5)) + this->robot_dim(6));
        t19 = sin(beta + -b_gamma);
        A = sin((-A + this->robot_dim(5)) + this->robot_dim(6));
        B = this->robot_inertial_par(6) * this->robot_inertial_par(8);
        sin_gamma = this->robot_inertial_par(6) * this->robot_inertial_par(9);
        cos_gamma = this->robot_inertial_par(6) * this->robot_inertial_par(7);
        A_tmp = this->robot_inertial_par(8) * this->robot_inertial_par(9);
        B_tmp = this->robot_inertial_par(7) * this->robot_inertial_par(8);
        J_p_eq_tmp = this->robot_inertial_par(7) * this->robot_inertial_par(9);
        b_J_p_eq_tmp = this->robot_inertial_par(2) * this->robot_dim(4) *
                      this->robot_inertial_par(6);
        c_J_p_eq_tmp = this->robot_dim(0) * this->robot_inertial_par(3);
        d_J_p_eq_tmp = c_J_p_eq_tmp * this->robot_inertial_par(8);
        e_J_p_eq_tmp = this->robot_dim(3) * this->robot_inertial_par(4);
        this->I_p_eq_k=(((((((((((((((((((((((((((((((((((((((((this->robot_inertial_par(6) *
                                                            this->robot_inertial_par
                                                                [11] *
                                                            4.0 +
                                                        this->robot_inertial_par(8) *
                                                            this->robot_inertial_par
                                                                [11] *
                                                            4.0) +
                                                      this->robot_inertial_par(6) *
                                                          this->robot_inertial_par
                                                              [13] *
                                                          4.0) +
                                                      this->robot_inertial_par(9) *
                                                          this->robot_inertial_par
                                                              [11] *
                                                          2.0) +
                                                    this->robot_inertial_par(6) *
                                                        this->robot_inertial_par
                                                            [14] *
                                                        2.0) +
                                                    this->robot_inertial_par(7) *
                                                        this->robot_inertial_par(11) *
                                                        4.0) +
                                                  this->robot_inertial_par(6) *
                                                      this->robot_inertial_par(12) *
                                                      4.0) +
                                                  this->robot_inertial_par(8) *
                                                      this->robot_inertial_par(13) *
                                                      4.0) +
                                                this->robot_inertial_par(9) *
                                                    this->robot_inertial_par(13) *
                                                    2.0) +
                                                this->robot_inertial_par(8) *
                                                    this->robot_inertial_par(14) *
                                                    2.0) +
                                              this->robot_inertial_par(9) *
                                                  this->robot_inertial_par(14)) +
                                              this->robot_inertial_par(7) *
                                                  this->robot_inertial_par(13) *
                                                  4.0) +
                                            this->robot_inertial_par(8) *
                                                this->robot_inertial_par(12) * 4.0) +
                                            this->robot_inertial_par(7) *
                                                this->robot_inertial_par(14) * 2.0) +
                                          this->robot_inertial_par(9) *
                                              this->robot_inertial_par(12) * 2.0) +
                                          this->robot_inertial_par(7) *
                                              this->robot_inertial_par(12) * 4.0) +
                                        B * b_A_tmp * 4.0) +
                                        B * t6 * 4.0) +
                                      B * t8 * 4.0) +
                                      B * t10 * 4.0) +
                                    sin_gamma * c_A_tmp * 2.0) +
                                    sin_gamma * t8 * 2.0) +
                                  sin_gamma * t9 * 2.0) +
                                  sin_gamma * t10 * 2.0) +
                                cos_gamma * t7 * 4.0) +
                                cos_gamma * t8 * 4.0) +
                              cos_gamma * t10 * 4.0) +
                              A_tmp * b_A_tmp * 2.0) +
                            A_tmp * c_A_tmp * 2.0) +
                            A_tmp * t6 * 2.0) +
                          A_tmp * t9 * 2.0) +
                          B_tmp * b_A_tmp * 4.0) +
                        B_tmp * t6 * 4.0) +
                        B_tmp * t7 * 4.0) +
                      J_p_eq_tmp * c_A_tmp * 2.0) +
                      J_p_eq_tmp * t7 * 2.0) +
                    J_p_eq_tmp * t9 * 2.0) -
                    this->robot_dim(0) * this->robot_inertial_par(2) *
                        this->robot_inertial_par(6) * this->robot_inertial_par(8) *
                        8.0) +
                  this->robot_dim(0) * this->robot_dim(4) *
                      this->robot_inertial_par(6) * this->robot_inertial_par(8) *
                      8.0) -
                  this->robot_inertial_par(1) * this->robot_dim(3) *
                      this->robot_inertial_par(9) * this->robot_inertial_par(7) *
                      4.0) -
                b_J_p_eq_tmp * this->robot_inertial_par(8) * 8.0) +
                ((((((((((((((((((((((b_J_p_eq_tmp * this->robot_inertial_par(9) *
                                          -4.0 -
                                      b_J_p_eq_tmp * this->robot_inertial_par(7) *
                                          8.0) +
                                    this->robot_inertial_par(3) *
                                        this->robot_inertial_par(4) *
                                        this->robot_inertial_par(8) *
                                        this->robot_inertial_par(9) *
                                        cos((b_gamma + this->robot_dim(5)) +
                                                  this->robot_dim(6)) *
                                        4.0) +
                                    c_J_p_eq_tmp * this->robot_inertial_par(6) *
                                        this->robot_inertial_par(8) * t12 * 8.0) +
                                  this->robot_dim(0) * this->robot_dim(3) *
                                      this->robot_inertial_par(8) *
                                      this->robot_inertial_par(9) * t11 * 4.0) +
                                  d_J_p_eq_tmp * this->robot_inertial_par(9) * t12 *
                                      4.0) +
                                d_J_p_eq_tmp * this->robot_inertial_par(7) * t12 *
                                    8.0) +
                                this->robot_dim(0) * this->robot_inertial_par(1) *
                                    this->robot_inertial_par(8) *
                                    this->robot_inertial_par(7) * t11 * 8.0) -
                              this->robot_dim(0) * this->robot_inertial_par(4) *
                                  this->robot_inertial_par(8) *
                                  this->robot_inertial_par(9) * A * 4.0) -
                              this->robot_inertial_par(2) *
                                  this->robot_inertial_par(3) *
                                  this->robot_inertial_par(6) *
                                  this->robot_inertial_par(8) * t12 * 8.0) +
                            this->robot_inertial_par(2) * this->robot_dim(3) *
                                this->robot_inertial_par(6) *
                                this->robot_inertial_par(9) * t11 * 4.0) +
                            this->robot_inertial_par(1) *
                                this->robot_inertial_par(2) *
                                this->robot_inertial_par(6) *
                                this->robot_inertial_par(7) * t11 * 8.0) -
                          e_J_p_eq_tmp * this->robot_inertial_par(6) *
                              this->robot_inertial_par(9) * t16 * 4.0) -
                          this->robot_dim(3) * this->robot_inertial_par(3) *
                              this->robot_inertial_par(8) *
                              this->robot_inertial_par(9) * t19 * 4.0) -
                        this->robot_inertial_par(2) * this->robot_inertial_par(4) *
                            this->robot_inertial_par(6) *
                            this->robot_inertial_par(9) * A * 4.0) -
                        this->robot_inertial_par(1) * this->robot_inertial_par(3) *
                            this->robot_inertial_par(8) *
                            this->robot_inertial_par(7) * t19 * 8.0) -
                      e_J_p_eq_tmp * this->robot_inertial_par(8) *
                          this->robot_inertial_par(9) * t16 * 4.0) +
                      this->robot_inertial_par(3) * this->robot_dim(4) *
                          this->robot_inertial_par(6) * this->robot_inertial_par(8) *
                          t12 * 8.0) -
                    this->robot_dim(3) * this->robot_dim(4) *
                        this->robot_inertial_par(6) * this->robot_inertial_par(9) *
                        t11 * 4.0) -
                    e_J_p_eq_tmp * this->robot_inertial_par(9) *
                        this->robot_inertial_par(7) * t16 * 4.0) +
                  this->robot_inertial_par(1) * this->robot_inertial_par(4) *
                      this->robot_inertial_par(9) * this->robot_inertial_par(7) *
                      t16 * 4.0) -
                  this->robot_inertial_par(1) * this->robot_dim(4) *
                      this->robot_inertial_par(6) * this->robot_inertial_par(7) *
                      t11 * 8.0) +
                this->robot_inertial_par(4) * this->robot_dim(4) *
                    this->robot_inertial_par(6) * this->robot_inertial_par(9) * A *
                    4.0)) /
              (((this->robot_inertial_par(6) * 2.0 +
                  this->robot_inertial_par(8) * 2.0) +
                this->robot_inertial_par(9)) +
                this->robot_inertial_par(7) * 2.0);
      }

      void phi_r_st_eq_k_calc(){//equivalent pendulum Y-axis rotational moment of inertia

        double beta=this->q_k(4);

        double A;
        double B;
        double G1;
        double G2;
        double a_tmp;
        double a_tmp_tmp;
        double a_tmp_tmp_tmp;
        double b_a_tmp;
        double b_a_tmp_tmp;
        double b_h1_tmp;
        double b_h1_tmp_tmp_tmp;
        double c_a_tmp;
        double c_a_tmp_tmp;
        double c_h1_tmp;
        double cos_gamma;
        double d_G1_beta;
        double d_a_tmp;
        double d_a_tmp_tmp;
        double d_h1_tmp;
        double e_a_tmp_tmp;
        double e_h1_tmp;
        double f_h1_tmp;
        double g_h1_tmp;
        double h1;
        double h1_tmp;
        double h1_tmp_tmp;
        double h1_tmp_tmp_tmp;
        double h2;
        double h_h1_tmp;
        double i_h1_tmp;
        double j_h1_tmp;
        h1_tmp_tmp = cos(beta);
        h1_tmp_tmp_tmp = sin(beta);
        a_tmp_tmp_tmp = this->robot_dim(3) * h1_tmp_tmp;
        a_tmp = this->robot_dim(1) - a_tmp_tmp_tmp;
        a_tmp_tmp = this->robot_dim(3) * this->robot_dim(3);
        A = 2.0 * h1_tmp_tmp;
        b_a_tmp_tmp = this->robot_dim(0) * this->robot_dim(0);
        c_a_tmp_tmp = this->robot_dim(2) * this->robot_dim(2);
        d_a_tmp_tmp = A * this->robot_dim(1);
        e_a_tmp_tmp = this->robot_dim(1) * this->robot_dim(1);
        B = d_a_tmp_tmp * this->robot_dim(3);
        G1 = -b_a_tmp_tmp + e_a_tmp_tmp;
        b_a_tmp = ((G1 - B) + c_a_tmp_tmp) + a_tmp_tmp;
        h2 = a_tmp * a_tmp;
        c_a_tmp = h1_tmp_tmp_tmp * h1_tmp_tmp_tmp;
        d_G1_beta = a_tmp_tmp * c_a_tmp;
        d_a_tmp = d_G1_beta / h2;
        h1_tmp = 2.0 * b_a_tmp_tmp;
        b_h1_tmp = pow(a_tmp, 3.0);
        c_h1_tmp = 4.0 * c_a_tmp_tmp;
        d_h1_tmp = (e_a_tmp_tmp - B) + a_tmp_tmp;
        e_h1_tmp = 2.0 * a_tmp_tmp;
        f_h1_tmp = this->robot_dim(1) * a_tmp_tmp;
        g_h1_tmp = this->robot_dim(2) * b_h1_tmp;
        h1_tmp = (((((((((pow(this->robot_dim(0), 4.0) +
                          A * b_a_tmp_tmp * this->robot_dim(1) * this->robot_dim(3)) -
                        h1_tmp * c_a_tmp_tmp) -
                        h1_tmp * a_tmp_tmp) -
                      pow(this->robot_dim(1), 4.0)) +
                      A * pow(this->robot_dim(1), 3.0) * this->robot_dim(3)) +
                    d_a_tmp_tmp * c_a_tmp_tmp * this->robot_dim(3)) -
                    d_a_tmp_tmp * pow(this->robot_dim(3), 3.0)) +
                  pow(this->robot_dim(2), 4.0)) -
                  2.0 * c_a_tmp_tmp * a_tmp_tmp) +
                pow(this->robot_dim(3), 4.0);
        h_h1_tmp = c_h1_tmp * b_h1_tmp;
        b_h1_tmp_tmp_tmp = this->robot_dim(3) * h1_tmp_tmp_tmp;
        B = this->robot_dim(2) * h2;
        A = 2.0 * this->robot_dim(2) * h2;
        cos_gamma = sqrt((d_a_tmp - b_a_tmp * b_a_tmp / (c_h1_tmp * h2)) + 1.0);
        i_h1_tmp = this->robot_dim(3) - this->robot_dim(1) * h1_tmp_tmp;
        G2 = cos_gamma + b_h1_tmp_tmp_tmp * b_a_tmp / A;
        j_h1_tmp = 2.0 * this->robot_dim(1) * this->robot_dim(3) * h1_tmp_tmp;
        h1 =
            -(2.0 * this->robot_dim(0) *
              ((B *
                    (((a_tmp_tmp_tmp * b_a_tmp / A + f_h1_tmp * c_a_tmp / B) -
                      d_G1_beta * b_a_tmp / g_h1_tmp) +
                    b_h1_tmp_tmp_tmp * h1_tmp / (h_h1_tmp * cos_gamma)) /
                    (this->robot_dim(0) * d_h1_tmp) -
                a_tmp_tmp_tmp / this->robot_dim(0)) +
              2.0 * this->robot_dim(2) * a_tmp_tmp * h1_tmp_tmp_tmp * i_h1_tmp * G2 /
                  (this->robot_dim(0) * b_h1_tmp *
                    ((d_a_tmp + 1.0) * (d_a_tmp + 1.0)))) *
              d_h1_tmp) /
            (a_tmp *
            (((((((e_h1_tmp * (h1_tmp_tmp * h1_tmp_tmp) + e_h1_tmp * c_a_tmp) +
                  b_a_tmp_tmp) +
                  e_a_tmp_tmp) -
                c_a_tmp_tmp) -
                a_tmp_tmp) -
              j_h1_tmp) +
              2.0 * this->robot_dim(2) * this->robot_dim(3) * h1_tmp_tmp_tmp *
                  cos_gamma));
        h2 =
            -(h2 *
                  (((this->robot_dim(3) * h1_tmp_tmp *
                        (((G1 - 2.0 * h1_tmp_tmp * this->robot_dim(1) *
                                    this->robot_dim(3)) +
                          c_a_tmp_tmp) +
                          a_tmp_tmp) /
                        (2.0 * this->robot_dim(2) * (a_tmp * a_tmp)) +
                    f_h1_tmp * (h1_tmp_tmp_tmp * h1_tmp_tmp_tmp) /
                        (this->robot_dim(2) * (a_tmp * a_tmp))) -
                    a_tmp_tmp * (h1_tmp_tmp_tmp * h1_tmp_tmp_tmp) * b_a_tmp /
                        g_h1_tmp) +
                  this->robot_dim(3) * h1_tmp_tmp_tmp * h1_tmp /
                      (h_h1_tmp *
                        sqrt((a_tmp_tmp * (h1_tmp_tmp_tmp * h1_tmp_tmp_tmp) /
                                      (a_tmp * a_tmp) -
                                  b_a_tmp * b_a_tmp / (c_h1_tmp * (a_tmp * a_tmp))) +
                                  1.0))) /
                  d_h1_tmp +
              e_h1_tmp * h1_tmp_tmp_tmp * i_h1_tmp * G2 /
                  (b_h1_tmp * ((d_a_tmp + 1.0) * (d_a_tmp + 1.0)))) /
            (b_a_tmp / (2.0 * this->robot_dim(2) * a_tmp) -
            b_h1_tmp_tmp_tmp * a_tmp * G2 / d_h1_tmp);
        
        d_a_tmp_tmp = a_tmp_tmp_tmp - this->robot_dim(1);
        A = ((((b_a_tmp_tmp - e_a_tmp_tmp) - c_a_tmp_tmp) - a_tmp_tmp) + j_h1_tmp) /
            (2.0 * this->robot_dim(2) * d_a_tmp_tmp);
        B = b_h1_tmp_tmp_tmp / d_a_tmp_tmp;
        d_a_tmp_tmp = B * B + 1.0;
        cos_gamma = (-A * B + sqrt(d_a_tmp_tmp - A * A)) / d_a_tmp_tmp;
        A += B * cos_gamma;
        d_a_tmp_tmp = this->robot_dim(3) / this->robot_dim(0);
        B = this->robot_dim(2) / this->robot_dim(0);
        d_a_tmp_tmp = atan2(d_a_tmp_tmp * h1_tmp_tmp_tmp - B * cos_gamma,
                                    (this->robot_dim(1) / this->robot_dim(0) -
                                    d_a_tmp_tmp * h1_tmp_tmp) -
                                        B * A) -
                      this->robot_dim(5);
        b_a_tmp = sin(d_a_tmp_tmp);
        B = atan2(A, cos_gamma) + this->robot_dim(5);
        f_h1_tmp = cos(B);
        A = beta + this->robot_dim(5);
        g_h1_tmp = sin(A);
        j_h1_tmp = this->robot_dim(3) * this->robot_inertial_par(9) +
                  2.0 * this->robot_inertial_par(1) * this->robot_inertial_par(7);
        a_tmp = (((2.0 * this->robot_dim(0) * this->robot_inertial_par(8) +
                  2.0 * this->robot_inertial_par(2) * this->robot_inertial_par(6)) +
                  2.0 * this->robot_dim(4) * this->robot_inertial_par(8)) +
                this->robot_dim(4) * this->robot_inertial_par(9)) +
                2.0 * this->robot_dim(4) * this->robot_inertial_par(7);
        cos_gamma = this->robot_inertial_par(4) * this->robot_inertial_par(9);
        h_h1_tmp = 2.0 * this->robot_inertial_par(3) * this->robot_inertial_par(8);
        G1 = this->g * (((g_h1_tmp * j_h1_tmp + b_a_tmp * a_tmp) -
                        cos_gamma * cos(this->robot_dim(6))) +
                        h_h1_tmp * f_h1_tmp);
        e_h1_tmp = cos(d_a_tmp_tmp);
        d_h1_tmp = sin(B);
        c_a_tmp = cos(A);
        B = 2.0 * h2 * this->robot_inertial_par(3) * this->robot_inertial_par(8);
        d_G1_beta =
            this->g * ((c_a_tmp * j_h1_tmp + h1 * e_h1_tmp * a_tmp) - B * d_h1_tmp);
        G2 =
            this->g *
            (((cos(beta + this->robot_dim(5)) *
                  (this->robot_dim(3) * this->robot_inertial_par(9) +
                    2.0 * this->robot_inertial_par(1) * this->robot_inertial_par(7)) -
              e_h1_tmp * a_tmp) -
              cos_gamma * sin(this->robot_dim(6))) -
            h_h1_tmp * d_h1_tmp);
        cos_gamma =
            this->g *
            ((h1 * b_a_tmp * a_tmp -
              sin(beta + this->robot_dim(5)) *
                  (this->robot_dim(3) * this->robot_inertial_par(9) +
                  2.0 * this->robot_inertial_par(1) * this->robot_inertial_par(7))) -
            B * f_h1_tmp);
        this->phi_r_st_eq_k = -atan(G1 / G2);
        h_h1_tmp = G1 * G1 + G2 * G2;
        i_h1_tmp = d_G1_beta * G2 - G1 * cos_gamma;
        this->dphi_r_st_eq_k_dbeta = -1.0 / h_h1_tmp * i_h1_tmp;
        A = h1 * h1;
        d_a_tmp_tmp =
            2.0 * this->g * this->robot_inertial_par(3) * this->robot_inertial_par(8);
        B = h2 * h2;
        this->ddphi_r_st_eq_k_ddbeta = -(
            -(2.0 * (G1 + G2)) / (h_h1_tmp * h_h1_tmp) * (d_G1_beta + cos_gamma) *
                i_h1_tmp +
            1.0 / h_h1_tmp *
                (((-this->g * b_a_tmp * a_tmp * A - d_a_tmp_tmp * f_h1_tmp * B) -
                  this->g * g_h1_tmp * j_h1_tmp) *
                    G2 -
                G1 * ((this->g * e_h1_tmp * a_tmp * A + d_a_tmp_tmp * d_h1_tmp * B) -
                      this->g * c_a_tmp * j_h1_tmp)));
      }

      void alpha_gamma_calc(double beta,double *alpha,double *b_gamma){
        double A;
        double A_tmp;
        double B;
        double B_tmp;
        double cos_gamma;
        double sin_gamma;
        A_tmp = cos(beta);
        A = ((((this->robot_dim(0,0) * this->robot_dim(0,0) -
                this->robot_dim(1,0) * this->robot_dim(1,0)) -
              this->robot_dim(2,0) * this->robot_dim(2,0)) -
              this->robot_dim(3,0) * this->robot_dim(3,0)) +
            2.0 * this->robot_dim(1,0) * this->robot_dim(3,0) * A_tmp) /
            (2.0 * this->robot_dim(2,0) *
            (this->robot_dim(3,0) * A_tmp - this->robot_dim(1,0)));
        B_tmp = sin(beta);
        B = this->robot_dim(3,0) * B_tmp /
            (this->robot_dim(3,0) * cos(beta) - this->robot_dim(1,0));
        sin_gamma = B * B + 1.0;
        cos_gamma = (-A * B + sqrt(sin_gamma - A * A)) / sin_gamma;
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
        h1_tmp_tmp = cos(beta);
        h1_tmp_tmp_tmp = sin(beta);
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
        j_knee_jacobians_tmp = sqrt(
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
                  (((this->robot_dim(3,0) * cos(beta) *
                        ((((-(this->robot_dim(0,0) * this->robot_dim(0,0)) +
                            this->robot_dim(1,0) * this->robot_dim(1,0)) -
                            2.0 * cos(beta) * this->robot_dim(1,0) *
                                this->robot_dim(3,0)) +
                          this->robot_dim(2,0) * this->robot_dim(2,0)) +
                          this->robot_dim(3,0) * this->robot_dim(3,0)) /
                        (2.0 * this->robot_dim(2,0) * (a_tmp * a_tmp)) +
                    e_knee_jacobians_tmp * (h1_tmp_tmp_tmp * h1_tmp_tmp_tmp) /
                        (this->robot_dim(2,0) * (a_tmp * a_tmp))) -
                    b_a_tmp_tmp * (h1_tmp_tmp_tmp * h1_tmp_tmp_tmp) * b_a_tmp /
                        f_knee_jacobians_tmp) +
                  this->robot_dim(3,0) * sin(beta) *
                      ((((((((((pow(this->robot_dim(0,0), 4.0) +
                                2.0 * cos(beta) *
                                    (this->robot_dim(0,0) * this->robot_dim(0,0)) *
                                    this->robot_dim(1,0) * this->robot_dim(3,0)) -
                                2.0 * (this->robot_dim(0,0) * this->robot_dim(0,0)) *
                                    (this->robot_dim(2,0) * this->robot_dim(2,0))) -
                              2.0 * (this->robot_dim(0,0) * this->robot_dim(0,0)) *
                                  (this->robot_dim(3,0) * this->robot_dim(3,0))) -
                              pow(this->robot_dim(1,0), 4.0)) +
                            2.0 * cos(beta) *
                                pow(this->robot_dim(1,0), 3.0) *
                                this->robot_dim(3,0)) +
                            2.0 * cos(beta) * this->robot_dim(1,0) *
                                (this->robot_dim(2,0) * this->robot_dim(2,0)) *
                                this->robot_dim(3,0)) -
                          2.0 * cos(beta) * this->robot_dim(1,0) *
                              pow(this->robot_dim(3,0), 3.0)) +
                          pow(this->robot_dim(2,0), 4.0)) -
                        2.0 * (this->robot_dim(2,0) * this->robot_dim(2,0)) *
                            (this->robot_dim(3,0) * this->robot_dim(3,0))) +
                        pow(this->robot_dim(3,0), 4.0)) /
                      (g_knee_jacobians_tmp *
                        sqrt((b_a_tmp_tmp * (h1_tmp_tmp_tmp * h1_tmp_tmp_tmp) /
                                      (a_tmp * a_tmp) -
                                  b_a_tmp * b_a_tmp /
                                      (c_knee_jacobians_tmp * (a_tmp * a_tmp))) +
                                  1.0))) /
                  h_a_tmp_tmp +
              d_knee_jacobians_tmp * h1_tmp_tmp_tmp *
                  (this->robot_dim(3,0) - this->robot_dim(1,0) * cos(beta)) *
                  k_knee_jacobians_tmp /
                  (b_knee_jacobians_tmp *
                  ((l_a_tmp_tmp + 1.0) * (l_a_tmp_tmp + 1.0)))) /
            (b_a_tmp / (2.0 * this->robot_dim(2,0) * a_tmp) -
            knee_jacobians_tmp_tmp_tmp * a_tmp * k_knee_jacobians_tmp / h_a_tmp_tmp);
      }
  
   
    //PRIVATE ATTRIBUTES
    private: 
      ////////////////////////////////////////////////////
      physics::WorldPtr world; // Pointer to the world the model is in
      physics::ModelPtr model; // Pointer to the model
      physics::PhysicsEnginePtr physics_engine_ptr; //pointer to the current physics engine
      event::ConnectionPtr updateConnection; // Pointer to the update event connection
      event::ConnectionPtr end_updateConnection; // Pointer to the update event connection
      event::ConnectionPtr before_physics_updateConnection; //Pointer to the before physics update event connection

      double max_step_size;//integration step size
      uint32_t sim_steps_over_control_step; //number of controls per integration step (used to provide the controller with the right state)

      ignition::math::Vector3d gravity; // world Gravity vector
      double g;
      ///////////////////////////////////////////////////
      double delta_s;// controller sampling time 
      double T_p;//prediction horizon
      uint32_t N_p;//number of prediction intervals
      VectorXf cost_weights_y1=VectorXf::Zero(5);//cost function weights first output (phi_w)
      VectorXf cost_weights_y2=VectorXf::Zero(5);//cost function weights second output (beta)
     
      VectorXf C_pos=VectorXf::Zero(2);//outputs position matrix (of the single feedback linearized subsystems)ù
      VectorXf C_vel=VectorXf::Zero(2);//outputs velocity matrix (of the single feedback linearized subsystems)

      MatrixXf Ad=MatrixXf::Zero(2,2);
      VectorXf Bd=VectorXf::Zero(2);

      VectorXf v_k=VectorXf::Zero(2); //current linearized input
      double v1_k_m1=0;//linearized inputs at t_(k-1)
      double v2_k_m1=0;

      double v1_min_k=0;
      double v1_max_k=0;
      double v2_min_k=0;
      double v2_max_k=0;
      VectorXf v_min_k=VectorXf::Zero(2); 
      VectorXf v_max_k=VectorXf::Zero(2);
      VectorXf V_min_k;
      VectorXf V_max_k;
      VectorXf V_k;

      MatrixXf M_bar_k=MatrixXf::Zero(2,2);
      VectorXf C_bar_k=VectorXf::Zero(2);
      MatrixXf M_hat_k=MatrixXf::Zero(3,3);
      VectorXf C_hat_k=VectorXf::Zero(3);

      double lambda_u_k=0;
      double phi_r_st_eq_k=0;
      double dphi_r_st_eq_k_dbeta=0;
      double ddphi_r_st_eq_k_ddbeta=0;

      VectorXf theta_r_theta_r_dot_k=VectorXf::Zero(2);
      double b_u_k=0;
      VectorXf zeta_u_k=VectorXf::Zero(1);
      VectorXf zeta_k=VectorXf::Zero(2);
      MatrixXf Delta_k;
      MatrixXf T_k=MatrixXf::Zero(2,2);

      MatrixXf Y_ref_k;
      MatrixXf dY_ref_k;
      MatrixXf Y_ref;
      MatrixXf dY_ref;

      MatrixXf VV_k;//holds all the QP solutions during the simulation-->note that it is really useful to have an estimate of the joint accelerations and ground reactions
      ////////////////////////////////////////////////////
      bool forward_prop_state_meas=0;
      bool model_properties_are_manually_set;
      bool use_equiv_central_knee;
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
      double m_p_eq=0;
      double l_p_eq_k=0;
      double I_p_eq_k=0;
      //////////////////////////////////////////////////
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
      //////////////////////////////////////////////////
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
      /////////////////QP///////////////////////
      MatrixXf Omega_pos;
      MatrixXf Omega_vel;
      MatrixXf Gamma_pos;
      MatrixXf Gamma_vel;
      MatrixXf Gamma_pos1;
      MatrixXf Gamma_vel1;
      MatrixXf Gamma_pos2;
      MatrixXf Gamma_vel2;
      MatrixXf Du;
      MatrixXf D1u;
      MatrixXf D2u;
      MatrixXf Du_diff;
      MatrixXf D1u_diff;
      MatrixXf D2u_diff;
      MatrixXf Du_diff_prime;
      MatrixXf D1u_diff_prime;
      MatrixXf D2u_diff_prime;
      MatrixXf H;
      VectorXf f_k;
      VectorXf f1_diff_k;
      VectorXf f2_diff_k;

      MatrixXf Au_lim_k;
      VectorXf bu_lim_k;
      MatrixXf A_bound;
      VectorXf b_bound;

      MatrixXf A_eq_k;
      MatrixXf A_ineq_k;
      VectorXf b_eq_k;
      VectorXf b_ineq_k;
      VectorXf b_ineq_max_k;
      VectorXf b_ineq_min_k;

      double cost_fun_value_k;
      double cost_threshold;

      VectorXf b_neg_inf_QP=-10000000*(VectorXf::Ones(this->b_ineq_k.rows()));

      ////////////////////////////////////////////////////
      VectorXf u_lim=VectorXf::Zero(2); // torque limits: [Cm1_lim,Cm2_lim]
      VectorXf beta_lim=VectorXf::Zero(2); // knee angle limits: [beta_lb, beta_ub]
      VectorXf phi_r_lim=VectorXf::Zero(2);; // robot pitch angle limits
      /////////////////////////////////////////////////////
      VectorXf q_k=VectorXf::Zero(10); //current state (from measurements)
      VectorXf q_k_eq=VectorXf::Zero(10); //current equivalent state (from measurements)

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
      VectorXf QP_cost_value;
      /////////////////////////////////////////////////////
      double x_w_freq;
      double beta_freq;
      double x_w_amplitute;
      double beta_amplitute;
      double x_w_offset;
      double beta_offset;
      double phase_lag_x_w;
      double phase_lag_beta;
  };
  
  // Register this plugin with the simulator
  GZ_REGISTER_MODEL_PLUGIN(ModelController)
}

