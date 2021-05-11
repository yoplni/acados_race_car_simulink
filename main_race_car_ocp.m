%
% Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
% Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
% Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
% Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
%
% This file is part of acados.
%
% The 2-Clause BSD License
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.;
%
% Author: Daniel Kloeser
% Ported by Thomas Jespersen (thomasj@tkjelectronics.dk), TKJ Electronics
%

%% Example of the frc_racecars in simulation without obstacle avoidance:
%% This example is for the optimal racing of the frc race cars. The model is a simple bicycle model and the lateral acceleration is constraint in order to validate the model assumptions.
%% The simulation starts at s=-2m until one round is completed(s=8.71m). The beginning is cut in the final plots to simulate a 'warm start'. 

clear all

% model_path = fullfile(pwd,'..','pendulum_on_cart_model');
% addpath(model_path)

check_acados_requirements()


import casadi.*
addpath('helper_functions', 'tracks');

% track_file = 'LMS_Track.txt';
% [Sref, ~, ~, ~, ~] = getTrack(track_file);

Sref = [0;0.0625000000000000;0.125000000000000;0.187500000000000;0.250000000000000;0.312500000000000;0.375000000000000;0.437500000000000;0.500000000000000;0.562500000000000;0.625000000000000;0.687500000000000;0.750000000000000;0.812500000000000;0.875000000000000;0.937500000000000;1;1.02453380000000;1.04906770000000;1.07360150000000;1.09813530000000;1.12266920000000;1.14720300000000;1.17173690000000;1.19627070000000;1.22080450000000;1.24533840000000;1.26987220000000;1.29440600000000;1.31893990000000;1.34347370000000;1.36800760000000;1.39254140000000;1.41707520000000;1.44160910000000;1.46614290000000;1.49067670000000;1.51521060000000;1.53974440000000;1.56427830000000;1.58881210000000;1.61334590000000;1.63787980000000;1.66241360000000;1.68694740000000;1.71148130000000;1.73601510000000;1.76054900000000;1.78508280000000;1.80961660000000;1.83415050000000;1.85868430000000;1.88321810000000;1.90775200000000;1.93228580000000;1.95681960000000;1.98135350000000;2.00588730000000;2.03042120000000;2.05495500000000;2.07948880000000;2.10402270000000;2.12855650000000;2.15309030000000;2.17762420000000;2.20215800000000;2.22669190000000;2.25122570000000;2.27575950000000;2.30029340000000;2.32482720000000;2.34936100000000;2.37389490000000;2.39842870000000;2.42296260000000;2.44749640000000;2.47203020000000;2.49656410000000;2.52109790000000;2.54563170000000;2.57016560000000;2.59469940000000;2.61923330000000;2.64376710000000;2.66830090000000;2.69283480000000;2.71736860000000;2.74190240000000;2.76643630000000;2.79097010000000;2.81550400000000;2.84003780000000;2.86457160000000;2.88910550000000;2.91363930000000;2.93817310000000;2.96270700000000;2.98724080000000;3.01177460000000;3.03630850000000;3.06084230000000;3.08537620000000;3.10991000000000;3.13444380000000;3.15897770000000;3.18351150000000;3.20804530000000;3.23257920000000;3.25711300000000;3.28164690000000;3.30618070000000;3.33071450000000;3.35524840000000;3.41774840000000;3.48024840000000;3.54274840000000;3.60524840000000;3.66774840000000;3.73024840000000;3.79274840000000;3.85524840000000;3.87978220000000;3.90431600000000;3.92884990000000;3.95338370000000;3.97791760000000;4.00245140000000;4.02698520000000;4.05151910000000;4.07605290000000;4.10058670000000;4.12512060000000;4.14965440000000;4.17418830000000;4.19872210000000;4.22325590000000;4.24778980000000;4.31028980000000;4.37278980000000;4.43528980000000;4.49778980000000;4.56028980000000;4.62278980000000;4.68528980000000;4.74778980000000;4.77232360000000;4.79685740000000;4.82139130000000;4.84592510000000;4.87045890000000;4.89499280000000;4.91952660000000;4.94406050000000;4.96859430000000;4.99312810000000;5.01766200000000;5.04219580000000;5.06672960000000;5.09126350000000;5.11579730000000;5.14033120000000;5.16486500000000;5.18939880000000;5.21393270000000;5.23846650000000;5.26300030000000;5.28753420000000;5.31206800000000;5.33660190000000;5.36113570000000;5.38566950000000;5.41020340000000;5.43473720000000;5.45927100000000;5.48380490000000;5.50833870000000;5.53287260000000;5.59537260000000;5.65787260000000;5.72037260000000;5.78287260000000;5.84537260000000;5.90787260000000;5.97037260000000;6.03287260000000;6.05740640000000;6.08194020000000;6.10647410000000;6.13100790000000;6.15554170000000;6.18007560000000;6.20460940000000;6.22914320000000;6.25367710000000;6.27821090000000;6.30274480000000;6.32727860000000;6.35181240000000;6.37634630000000;6.40088010000000;6.42541390000000;6.44994780000000;6.47448160000000;6.49901550000000;6.52354930000000;6.54808310000000;6.57261700000000;6.59715080000000;6.62168460000000;6.64621850000000;6.67075230000000;6.69528620000000;6.71982000000000;6.74435380000000;6.76888770000000;6.79342150000000;6.81795530000000;6.88045530000000;6.94295530000000;7.00545530000000;7.06795530000000;7.13045530000000;7.19295530000000;7.25545530000000;7.31795530000000;7.38045530000000;7.44295530000000;7.50545530000000;7.56795530000000;7.63045530000000;7.69295530000000;7.75545530000000;7.81795530000000;7.84248920000000;7.86702300000000;7.89155690000000;7.91609070000000;7.94062450000000;7.96515840000000;7.98969220000000;8.01422600000000;8.03875990000000;8.06329370000000;8.08782750000000;8.11236140000000;8.13689520000000;8.16142910000000;8.18596290000000;8.21049670000000;8.28192530000000;8.35335390000000;8.42478240000000;8.49621100000000;8.56763960000000;8.63906820000000;8.71049670000000];

% load track parameters
% [s00, ~, ~, ~, kapparef] = getTrack(track_file);

s00 = Sref;

kapparef = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;0;0;0;0;0;0;0;0;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;0;0;0;0;0;0;0;0;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;0;0;0;0;0;0;0;0;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;-4;0;0;0;0;0;0;0;0];

%% Solver parameters
compile_interface = 'auto';
codgen_model = 'true';
nlp_solver = 'sqp_rti'; % sqp, sqp_rti
qp_solver = 'partial_condensing_hpipm';
    % full_condensing_hpipm, partial_condensing_hpipm, full_condensing_qpoases
nlp_solver_exact_hessian = 'false'; % false=gauss_newton, true=exact
qp_solver_cond_N = 5; % for partial condensing
regularize_method = 'no_regularize';
%regularize_method = 'project';
%regularize_method = 'mirror';
%regularize_method = 'convexify';
% integrator type
sim_method = 'erk'; % erk, irk, irk_gnsf

%% horizon parameters
N = 50;
T = 1.0; % time horizon length

%% model dynamics
[model, constraint] = bicycle_model(s00, kapparef);

% Define initial conditions
model.x0 = [-2, 0, 0, 0, 0, 0];

nx = length(model.x);
nu = length(model.u);

%% model to create the solver
ocp_model = acados_ocp_model();

%% acados ocp model
ocp_model.set('name', model.name);
ocp_model.set('T', T);

% symbolics
ocp_model.set('sym_x', model.x);
ocp_model.set('sym_u', model.u);
ocp_model.set('sym_xdot', model.xdot);
%ocp_model.set('sym_z', model.z);
%ocp_model.set('sym_p', model.p);

% dynamics
if (strcmp(sim_method, 'erk'))
    ocp_model.set('dyn_type', 'explicit');
    ocp_model.set('dyn_expr_f', model.f_expl_expr);
else % irk irk_gnsf
    ocp_model.set('dyn_type', 'implicit');
    ocp_model.set('dyn_expr_f', model.f_impl_expr);
end

% constraintsJbx = zeros(1,nx);
nbx = 1;
Jbx = zeros(nbx,nx);
Jbx(1,2) = 1;
ocp_model.set('constr_Jbx', Jbx);
ocp_model.set('constr_lbx', -12);
ocp_model.set('constr_ubx', 12);

nbu = 2;
Jbu = zeros(nbu,nu);
Jbu(1,1) = 1;
Jbu(2,2) = 1;
ocp_model.set('constr_Jbu', Jbu);
ocp_model.set('constr_lbu', [model.dthrottle_min, model.ddelta_min]);
ocp_model.set('constr_ubu', [model.dthrottle_max, model.ddelta_max]);

%ocp_model.set('constr_type', 'bgh');
nh = 5;
ocp_model.set('constr_expr_h', constraint.expr);
ocp_model.set('constr_lh', [...
                                constraint.along_min,...
                                constraint.alat_min,...
                                model.n_min,...
                                model.throttle_min,...
                                model.delta_min,...
                            ]);
ocp_model.set('constr_uh', [...
                                constraint.along_max,...
                                constraint.alat_max,...
                                model.n_max,...
                                model.throttle_max,...
                                model.delta_max,...
                            ]);
%ocp_model.set('constr_expr_h_e', constraint.expr);
%ocp_model.set('constr_lh_e', 0);
%ocp_model.set('constr_uh_e', 0);

% Configure constraint slack variables
nsh = 2;
Jsh = zeros(nh, nsh);
Jsh(1,1) = 1;
Jsh(3,2) = 1;
ocp_model.set('constr_Jsh', Jsh);
% Set cost on slack
% L1 slack (linear term)
ocp_model.set('cost_zl', 100 * ones(nsh,1));
ocp_model.set('cost_zu', 100 * ones(nsh,1));
% L2 slack (squared term)
ocp_model.set('cost_Zl', 0 * ones(nsh,nsh));
ocp_model.set('cost_Zu', 0 * ones(nsh,nsh));

% set intial condition
ocp_model.set('constr_x0', model.x0);

% cost = define linear cost on x and u
%ocp_model.set('cost_expr_ext_cost', model.expr_ext_cost);
%ocp_model.set('cost_expr_ext_cost_e', model.expr_ext_cost_e);

ocp_model.set('cost_type', 'linear_ls');
ocp_model.set('cost_type_e', 'linear_ls');

% number of outputs is the concatenation of x and u
ny = nx + nu;
ny_e = nx;

% The linear cost contributions is defined through Vx, Vu and Vz
Vx = zeros(ny, nx);
Vx_e = zeros(ny_e, nx);
Vu = zeros(ny, nu);

Vx(1:nx,:) = eye(nx);
Vx_e(1:nx,:) = eye(nx);
Vu(nx+1:end,:) = eye(nu);
ocp_model.set('cost_Vx', Vx);
ocp_model.set('cost_Vx_e', Vx_e);
ocp_model.set('cost_Vu', Vu);

% Define cost on states and input
Q = diag([ 1e-1, 1e-8, 1e-8, 1e-8, 1e-3, 5e-3 ]);
R = eye(nu);
R(1, 1) = 1e-3;
R(2, 2) = 5e-3;
Qe = diag([ 5e0, 1e1, 1e-8, 1e-8, 5e-3, 2e-3 ]);

unscale = N / T;
W = unscale * blkdiag(Q, R);
W_e = Qe / unscale;
ocp_model.set('cost_W', W);
ocp_model.set('cost_W_e', W_e);

% set intial references
y_ref = zeros(ny,1);
y_ref_e = zeros(ny_e,1);
y_ref(1) = 1; % set reference on 's' to 1 to push the car forward (progress)
ocp_model.set('cost_y_ref', y_ref);
ocp_model.set('cost_y_ref_e', y_ref_e);

% ... see ocp_model.model_struct to see what other fields can be set

%% acados ocp set opts
ocp_opts = acados_ocp_opts();
%ocp_opts.set('compile_interface', compile_interface);
%ocp_opts.set('codgen_model', codgen_model);
ocp_opts.set('param_scheme_N', N);
ocp_opts.set('nlp_solver', nlp_solver);
ocp_opts.set('nlp_solver_exact_hessian', nlp_solver_exact_hessian); 
ocp_opts.set('sim_method', sim_method);
ocp_opts.set('sim_method_num_stages', 4);
ocp_opts.set('sim_method_num_steps', 3);
ocp_opts.set('qp_solver', qp_solver);
%ocp_opts.set('regularize_method', regularize_method);
ocp_opts.set('qp_solver_cond_N', qp_solver_cond_N);
ocp_opts.set('nlp_solver_tol_stat', 1e-4);
ocp_opts.set('nlp_solver_tol_eq', 1e-4);
ocp_opts.set('nlp_solver_tol_ineq', 1e-4);
ocp_opts.set('nlp_solver_tol_comp', 1e-4);
% ... see ocp_opts.opts_struct to see what other fields can be set

%% create ocp solver
ocp = acados_ocp(ocp_model, ocp_opts);

%% Simulate - call ocp solver
dt = T / N;
Tf = 10.00;  % maximum simulation time[s]
Nsim = round(Tf / dt);
sref_N = 3;  % reference for final reference progress

% initialize data structs
% simX = zeros(Nsim, nx);
% simU = zeros(Nsim, nu);
% s0 = model.x0(1);
tcomp_sum = 0;
tcomp_max = 0;

%% call ocp solver
% update initial state
% ocp.set('constr_x0', model.x0);
% ocp.set('constr_lbx', model.x0, 0)
% ocp.set('constr_ubx', model.x0, 0)

% set trajectory initialization
ocp.set('init_x', model.x0' * ones(1,N+1));
ocp.set('init_u', zeros(nu, N));
ocp.set('init_pi', zeros(nx, N));

% simulate
% for i = 1:Nsim
%     % update reference
%     sref = s0 + sref_N;
%     for j = 0:(N-1)
%         yref = [s0 + (sref - s0) * j / N, 0, 0, 0, 0, 0, 0, 0];
% %         yref=[1,0,0,0,0,0,0,0]
%         ocp.set('cost_y_ref', yref, j);   
%     end
%     % yref_N = [sref, 0, 0, 0, 0, 0];
%     yref_N = [0, 0, 0, 0, 0, 0];
% %     yref_N = np.array([0,0,0,0,0,0]);
%     ocp.set('cost_y_ref_e', yref_N);

    % solve ocp
    t = tic();

    ocp.solve();
    status = ocp.get('status'); % 0 - success
    if status ~= 0
        % borrowed from acados/utils/types.h
        %statuses = {
        %    0: 'ACADOS_SUCCESS',
        %    1: 'ACADOS_FAILURE',
        %    2: 'ACADOS_MAXITER',
        %    3: 'ACADOS_MINSTEP',
        %    4: 'ACADOS_QP_FAILURE',
        %    5: 'ACADOS_READY'
        error(sprintf('acados returned status %d in closed loop iteration %d. Exiting.', status, i));
    end
    %ocp.print('stat')

    elapsed = toc(t);

    % manage timings
    tcomp_sum = tcomp_sum + elapsed;
    if elapsed > tcomp_max
        tcomp_max = elapsed;
    end

    % get solution
    utraj = ocp.get('u', 0);
    xtraj = ocp.get('x', 0);
%     for j = 1:nx
%         simX(i, j) = x0(j);
%     end
%     for j = 1:nu
%         simU(i, j) = u0(j);
%     end
%% Commented out - might need to revert
    % update initial condition
    x0 = ocp.get('x', 1);
    % update initial state
    ocp.set('constr_x0', x0);    
    ocp.set('constr_lbx', x0, 0);
    ocp.set('constr_ubx', x0, 0);
    % s0 = x0(1);

%     % check if one lap is done and break and remove entries beyond
%     if s0 > Sref(end) + 0.1
%         % find where vehicle first crosses start line
%         N0 = find(diff(sign(simX(:, 1))));
%         N0 = N0(1);
%         Nsim = i - N0 + 1;  % correct to final number of simulation steps for plotting
%         simX = simX(N0:i, :);
%         simU = simU(N0:i, :);
% %         break
%     end
% end

%% Plots
% t = linspace(0.0, Nsim * dt, Nsim);
% 
% % Plot results
% figure(1);
% subplot(2,1,1);
% plot(t, simU(:,1), 'r');
% hold on;
% plot(t, simU(:,2), 'g');
% hold off;
% title('closed-loop simulation');
% legend('dD','ddelta');
% xlabel('t');
% ylabel('u');
% grid;
% xlim([t(1), t(end)]);
% 
% subplot(2,1,2);
% plot(t, simX);
% xlabel('t');
% ylabel('x');
% legend('s','n','alpha','v','D','delta');
% grid;
% xlim([t(1), t(end)]);
% 
% % Plot track
% figure(2);
% plotTrackProj(simX, track_file);
% 
% % Plot alat
% alat = zeros(Nsim,1);
% for i = 1:Nsim
%     alat(i) = full(constraint.alat(simX(i,:),simU(i,:)));
% end
% figure(3);
% plot(t, alat);
% line([t(1), t(end)], [constraint.alat_min, constraint.alat_min], 'LineStyle', '--', 'Color', 'k');
% line([t(1), t(end)], [constraint.alat_max, constraint.alat_max], 'LineStyle', '--', 'Color', 'k');
% xlabel('t');
% ylabel('alat');
% xlim([t(1), t(end)]);


%% go embedded
% to generate templated C code
% ocp.generate_c_code;
