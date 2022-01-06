function [figures_handles,graphs_handles]=graph_plot_Gazebo_task_tracker()
    
    [Q_sol_gazebo,time_gazebo,X_k,U_middle,~,~,Chi_gazebo,Chi_dot_gazebo,Chi_ref_gazebo,Chi_dot_ref_gazebo,~,~,~]=LoadDataFromGazeboSim_task_tracker();
    
    f1=figure('Name','State evolution');
    f1.Visible = 'on';
%     f1.WindowState ='maximized';
    tiledlayout(5,1);

    plt1=nexttile;
    ax1=gca;
    hold(ax1)
    graph1_1=plot(ax1,time_gazebo,Q_sol_gazebo(:,1),'Color','b','LineStyle','-','LineWidth',1);
    grid on
    xlabel('t\,[s]','Interpreter','latex')
    ylabel('[m]','Interpreter','latex')
    title(plt1,'$x_w$','Interpreter','latex')

    plt3=nexttile;
    ax3=gca;
    hold(ax3)
    graph3_1=plot(ax3,time_gazebo,Q_sol_gazebo(:,2),'Color','b','LineStyle','-','LineWidth',1);
    grid on
    xlabel('t\,[s]','Interpreter','latex')
    ylabel('[m]','Interpreter','latex')
    title(plt3,'$z_w$','Interpreter','latex')

    plt5=nexttile;
    ax5=gca;
    hold(ax5)
    graph5_1=plot(ax5,time_gazebo,Q_sol_gazebo(:,3)*180/pi,'Color','b','LineStyle','-','LineWidth',1);
    grid on
    xlabel('t\,[s]','Interpreter','latex')
    ylabel('[deg]','Interpreter','latex')
    title(plt5,'$\phi_w$','Interpreter','latex')

    plt7=nexttile;
    ax7=gca;
    hold(ax7)
    graph7_1=plot(ax7,time_gazebo,Q_sol_gazebo(:,4)*180/pi,'Color','b','LineStyle','-','LineWidth',1);
    grid on
    xlabel('t\,[s]','Interpreter','latex')
    ylabel('[deg]','Interpreter','latex')
    title(plt7,'$\phi_r$','Interpreter','latex')

    plt9=nexttile;
    ax9=gca;
    hold(ax9)
    graph9_1=plot(ax9,time_gazebo,Q_sol_gazebo(:,5)*180/pi,'Color','b','LineStyle','-','LineWidth',1);
    grid on
    xlabel('t\,[s]','Interpreter','latex')
    ylabel('[deg]','Interpreter','latex')
    title(plt9,'$\beta$','Interpreter','latex')

    f2=figure('Name','State velocity evolution');
    f2.Visible = 'on';
%     f2.WindowState ='maximized';
    tiledlayout(5,1);

    plt2=nexttile;
    ax2=gca;
    hold(ax2)
    graph2_1=plot(ax2,time_gazebo,Q_sol_gazebo(:,6),'Color','b','LineStyle','-','LineWidth',1);
    grid on
    xlabel('t\,[s]','Interpreter','latex')
    ylabel('[m/s]','Interpreter','latex')
    title(plt2,'$\dot{x}_w$','Interpreter','latex')

    plt4=nexttile;
    ax4=gca;
    hold(ax4)
    graph4_1=plot(ax4,time_gazebo,Q_sol_gazebo(:,7),'Color','b','LineStyle','-','LineWidth',1);
    grid on
    xlabel('t\,[s]','Interpreter','latex')
    ylabel('[m/s]','Interpreter','latex')
    title(plt4,'$\dot{z}_w$','Interpreter','latex')

    plt6=nexttile;
    ax6=gca;
    hold(ax6)
    graph6_1=plot(ax6,time_gazebo,Q_sol_gazebo(:,8)*180/pi,'Color','b','LineStyle','-','LineWidth',1);
    grid on
    xlabel('t\,[s]','Interpreter','latex')
    ylabel('[deg/s]','Interpreter','latex')
    title(plt6,'$\dot{\phi}_w$','Interpreter','latex')

    plt8=nexttile;
    ax8=gca;
    hold(ax8)
    graph8_1=plot(ax8,time_gazebo,Q_sol_gazebo(:,9)*180/pi,'Color','b','LineStyle','-','LineWidth',1);
    grid on
    xlabel('t\,[s]','Interpreter','latex')
    ylabel('[deg/s]','Interpreter','latex')
    title(plt8,'$\dot{\phi}_r$','Interpreter','latex')

    plt10=nexttile;
    ax10=gca;
    hold(ax10)
    graph10_1=plot(ax10,time_gazebo,Q_sol_gazebo(:,10)*180/pi,'Color','b','LineStyle','-','LineWidth',1);
    grid on
    xlabel('t\,[s]','Interpreter','latex')
    ylabel('[deg/s]','Interpreter','latex')
    title(plt10,'$\dot{\beta}$','Interpreter','latex')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ACTUATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    f4=figure('Name','Actuation input');
    f4.Visible = 'on';
%     f4.WindowState ='maximized';
    % Top two plots
    tiledlayout(2,1);
    
    plt14=nexttile;
    ax11=gca;
    hold(ax11)
    graph14_1=stairs(ax11,time_gazebo(1:(end-1)),U_middle(:,1),'Color','b','LineStyle','-','LineWidth',1);
    grid on
    xlabel('t\,[s]','Interpreter','latex')
    ylabel('[N\,m]','Interpreter','latex')
    title(plt14,'$C_{m1}$','Interpreter','latex')

    plt15=nexttile;
    ax12=gca;
    hold(ax12)
    graph15_1=stairs(ax12,time_gazebo(1:(end-1)),U_middle(:,2),'Color','b','LineStyle','-','LineWidth',1);
    grid on
    xlabel('t\,[s]','Interpreter','latex')
    ylabel('[N\,m]','Interpreter','latex')
    title(plt15,'$C_{m2}$','Interpreter','latex')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%GROUND REACTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    f6=figure('Name','Ground reactions');
    f6.Visible = 'on';
%     f6.WindowState ='maximized';
    tiledlayout(2,1);

    nexttile
    graph20=plot(time_gazebo(1:(end-1)),X_k(:,6),'Color','b','LineStyle','-','LineWidth',1);
    grid on
    xlabel('t\,[s]','Interpreter','latex')
    ylabel('[N]','Interpreter','latex')
    title('Horizontal reaction at ground','Interpreter','latex')

    nexttile
    graph21=plot(time_gazebo(1:(end-1)),X_k(:,7),'Color','b','LineStyle','-','LineWidth',1);
    grid on
    xlabel('t\,[s]','Interpreter','latex')
    ylabel('[N]','Interpreter','latex')
    title('Vertical reaction at ground','Interpreter','latex')

  
    f8=figure('Name','Relative CoM trajectory');
    f8.Visible = 'on';
    f8.WindowState ='maximized';
    tiledlayout(4,1);

    nexttile
    ax12=gca;
    hold(ax12)
    graph26=plot(ax12,time_gazebo(2:(end-1)),Chi_gazebo(:,1),'Color','b','LineStyle','-','LineWidth',1);
    graph30=plot(ax12,time_gazebo(2:(end-1)),Chi_ref_gazebo(:,1),'LineStyle','--');
    

    grid on
    xlabel('t\,[s]','Interpreter','latex')
    ylabel('[m]','Interpreter','latex')
    title('Relative CoM x position','Interpreter','latex')
    legend('Measured','Reference','Interpreter','latex');

    nexttile
    ax12=gca;
    hold(ax12)
    graph27=plot(ax12,time_gazebo(2:(end-1)),Chi_gazebo(:,2),'Color','b','LineStyle','-','LineWidth',1);
    graph31=plot(ax12,time_gazebo(2:(end-1)),Chi_ref_gazebo(:,2),'LineStyle','--');


    grid on
    xlabel('t\,[s]','Interpreter','latex')
    ylabel('[m]','Interpreter','latex')
    title('Relative CoM z position','Interpreter','latex')
    legend('Measured','Reference','Interpreter','latex');

    nexttile
    ax12=gca;
    hold(ax12)
    graph28=plot(ax12,time_gazebo(2:(end-1)),Chi_dot_gazebo(:,1),'Color','b','LineStyle','-','LineWidth',1);
    graph32=plot(ax12,time_gazebo(2:(end-1)),Chi_dot_ref_gazebo(:,1),'LineStyle','--');


    grid on
    xlabel('t\,[s]','Interpreter','latex')
    ylabel('[m/s]','Interpreter','latex')
    title('Relative CoM x velocity','Interpreter','latex')
    legend('Measured','Reference','Interpreter','latex');

    nexttile
    ax12=gca;
    hold(ax12)
    graph29=plot(ax12,time_gazebo(2:(end-1)),Chi_dot_gazebo(:,2),'Color','b','LineStyle','-','LineWidth',1);
    graph33=plot(ax12,time_gazebo(2:(end-1)),Chi_dot_ref_gazebo(:,2),'LineStyle','--');


    grid on
    xlabel('t\,[s]','Interpreter','latex')
    ylabel('[m/s]','Interpreter','latex')
    title('Relative CoM z velocity','Interpreter','latex')
    legend('Measured','Reference','Interpreter','latex');

    graphs_handles=[f1,f2,graph1_1,graph2_1,graph3_1,graph4_1,graph5_1,graph6_1,graph7_1,graph8_1,graph9_1,graph10_1,...
        f4,graph14_1,graph15_1,...
        f6,graph20,graph21...
        f8,graph26,graph27,graph28,graph29,...
        graph30,graph31,graph32,graph33]';
    figures_handles=[f1,f2,f4,f6,f8]';

end