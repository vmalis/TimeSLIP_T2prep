function TimeSlip_t2prep(tissues, TEeffective_values, TI_values)
%--------------------------------------------------------------------------
%
%   Function to simulate and plot Time-SLIP MRI sequences with T2prep module
%
%   INPUT:  
%       tissues (cell array)                - Cell array with tissue names and properties
%                                             Each row contains {'Name', T1, T2}
%                                             where T1 and T2 are in milliseconds.
%       TEeffective_values (1D array)       - Array of TE effective values for T2 preparation (double)
%       TI_values (1D array)                - Array of inversion times (double)
%
%   OUTPUT: 
%       This function plots the results for each tissue and saves figures
%       as PDFs for each sequence order.
%
%   EXAMPLE:
%       tissues = {
%           'Water', 2500, 1300;  
%           'Milk', 1400, 80; 
%           'Cream', 850, 50;
%       };
%       TEeffective_values = [50, 100, 150, 200, 300, 400];
%       TI_values = 20:100:4000;
%       TimeSlip_t2prep(tissues, TEeffective_values, TI_values);
%
%   SUBROUTINES:
%       - simulate_case
%       - apply_TimeSlip_first
%       - apply_T2prep_first
%       - apply_TimeSlip_only
%       - plot_universal_results
%       - plot_SIR
%       - store_results
%       - initialize_arrays
%       - relaxation
%       - T2prep
%       - rotation
%--------------------------------------------------------------------------
% Vadim Malis (vmalis@ucsd.edu) UC San Diego
%--------------------------------------------------------------------------

    % Simulation parameters
    M0 = 1;  % Equilibrium magnetization
    Rad = [1, 0, 0; 0, 1, 0; 0, 0, -1];  % Adiabatic inversion rotation matrix (inverts Mz)
    Rno = eye(3);  % No pulse (identity matrix)
    sequence_orders = {'TimeSlipFirst', 'T2prepFirst','NoPrep'};
    
    % Preallocate results structure
    results = struct();
    
    % Loop over tissues
    for t = 1:size(tissues, 1)
        tissue_name = tissues{t, 1};
        T1 = tissues{t, 2};
        T2 = tissues{t, 3};
        
        fprintf('Simulating for %s (T1=%d ms, T2=%d ms)\n', tissue_name, T1, T2);
        
        % Loop over sequence orders
        for order_idx = 1:length(sequence_orders)
            sequence_order = sequence_orders{order_idx};
            fprintf('Running %s sequence\n', sequence_order);
            
            % Loop over TEeffective values
            for te_idx = 1:length(TEeffective_values)
                TEeffective = TEeffective_values(te_idx);
                
                % Initialize arrays to store results
                [Mx_case1, My_case1, Mz_case1, Signal_case1, ...
                 Mx_case2, My_case2, Mz_case2, Signal_case2,] = initialize_arrays(length(TI_values));
                
                % Loop over TI values
                for idx = 1:length(TI_values)
                    TI = TI_values(idx);
                    
                    %% Case 1: Two adiabatic inversion pulses
                    [Mx_case1(idx), My_case1(idx), Mz_case1(idx), Signal_case1(idx)] = ...
                        simulate_case(T1, T2, TEeffective, TI, M0, Rad, sequence_order, true);
                    
                    %% Case 2: One adiabatic inversion pulse
                    [Mx_case2(idx), My_case2(idx), Mz_case2(idx), Signal_case2(idx)] = ...
                        simulate_case(T1, T2, TEeffective, TI, M0, Rad, sequence_order, false);
            

                end
                
                % Store results for both cases
                results = store_results(results, tissue_name, TEeffective, sequence_order, ...
                                        TI_values, Mx_case1, My_case1, Mz_case1, Signal_case1, ...
                                        Mx_case2, My_case2, Mz_case2, Signal_case2);
            end
        end
    end
    
    % Plot the results for both sequence orders
    plot_universal_results(results, TEeffective_values, tissues, 'TimeSlipFirst');
    plot_universal_results(results, TEeffective_values, tissues, 'T2prepFirst');
    plot_SIR(results, TEeffective_values, tissues, 'TimeSlipFirst')
    plot_SIR(results, TEeffective_values, tissues, 'T2prepFirst')

end

function [Mx, My, Mz, Signal] = simulate_case(T1, T2, TEeffective, TI, M0, Rad, sequence_order, double_inversion)
    % Initialize magnetization vector
    M = [0; 0; M0];
    
    % Determine sequence order
    switch sequence_order
        case 'TimeSlipFirst'
            M = apply_TimeSlip_first(M, T1, T2, TEeffective, TI, M0, Rad, double_inversion);
        case 'T2prepFirst'
            M = apply_T2prep_first(M, T1, T2, TEeffective, TI, M0, Rad, double_inversion);
        case 'NoPrep'
            M = apply_TimeSlip_only(M, T1, T2, TEeffective, TI, M0, Rad, double_inversion);
        otherwise
            error('Invalid sequence order.');
    end
    
    % Readout (signal is proportional to transverse magnetization)
    Signal = sqrt(M(1)^2 + M(2)^2);
    Mx = M(1);
    My = M(2);
    Mz = M(3);
end

function M = apply_TimeSlip_first(M, T1, T2, TEeffective, TI, M0, Rad, double_inversion)
    % Apply the sequence: TimeSlip -> TI -> T2prep
    M = Rad * M;  % First inversion pulse
    M = relaxation(M, T1, T2, 20, M0);
    
    if double_inversion
        M = Rad * M;  % Second inversion pulse if applicable
    end
    
    M = relaxation(M, T1, T2, TI, M0);  % Relaxation over TI
    M = T2prep(M, T2, TEeffective, M0);  % T2prep module
    M = rotation(M, 90, 'y');  % 90-degree excitation pulse
end


function M = apply_T2prep_first(M, T1, T2, TEeffective, TI, M0, Rad, double_inversion)
    % Apply the sequence: T2prep -> TimeSlip -> TI
    M = T2prep(M, T2, TEeffective, M0);  % T2prep module
    M = Rad * M;  % First inversion pulse
    M = relaxation(M, T1, T2, 20, M0);
    
    if double_inversion
        M = Rad * M;  % Second inversion pulse if applicable
    end
    
    M = relaxation(M, T1, T2, TI, M0);  % Relaxation over TI
    M = rotation(M, 90, 'y');  % 90-degree excitation pulse
end

function M = apply_TimeSlip_only(M, T1, T2, TEeffective, TI, M0, Rad, double_inversion)
    % Apply the sequence: TimeSlip -> TI -> T2prep
    M = Rad * M;  % First inversion pulse
    M = relaxation(M, T1, T2, 20, M0);
    
    if double_inversion
        M = Rad * M;  % Second inversion pulse if applicable
    end
    M = relaxation(M, T1, T2, TI, M0);  % Relaxation over TI

    M = rotation(M, 90, 'y');  % 90-degree excitation pulse
end

function plot_universal_results(results, TEeffective_values, tissues, sequence_name)
    for t = 1:size(tissues, 1)
        tissue_name = tissues{t, 1};
        figure('Name', [tissue_name, ' - ', sequence_name, ' Signal Comparison'],'Color', 'white', 'Units', 'normalized', 'OuterPosition', [0 0 0.3 0.5]);
        for te_idx = 1:length(TEeffective_values)
            TEeffective = TEeffective_values(te_idx);
            data = results.(tissue_name).(sequence_name).(['TE' num2str(TEeffective)]);
            dataNoPrep = results.(tissue_name).('NoPrep').(['TE' num2str(TEeffective)]);
            subplot(2, ceil(length(TEeffective_values)/2), te_idx);
            plot(data.TI, abs(data.Mx_case1), 'r', 'LineWidth', 1.5); hold on;
            plot(data.TI, abs(data.Mx_case2), 'b', 'LineWidth', 1.5);
            legend({'NS+S IRs', 'NS IR'}, 'Location', 'southeast', 'Interpreter', 'latex');
            ylim([0, 1]);
            axis square;
            title(sprintf('$T_{2}$ Prep TE = %d ms', TEeffective), 'Interpreter', 'latex');
            xlabel('Inversion Time TI (ms)', 'Interpreter', 'latex');
            ylabel('Signal Magnitude', 'Interpreter', 'latex');
            box on
            grid on;
            % Set tick labels to use LaTeX interpreter
            set(gca, 'TickLabelInterpreter', 'latex');
        end

        switch sequence_name
            case 'TimeSlipFirst'
                sequence_description = 'Time-SLIP $\rightarrow T_{2}\mathrm{Prep}$';
            case 'T2prepFirst'
                sequence_description = '$T_{2}\mathrm{Prep} \rightarrow$ Time-SLIP';
            case 'NoPrep'
                sequence_description = 'Time-SLIP';
            otherwise
                error('Invalid sequence order.');
        end
        sgtitle([strrep(tissue_name, '_', ' '), ' : ', sequence_description], 'Interpreter', 'latex');
        export_fig(strcat(sequence_name,tissue_name,'_SIR.pdf'))

    end
end



function plot_SIR(results, TEeffective_values, tissues, sequence_name)
    % Create a figure for the plots with specified properties
    figure('Name', [sequence_name, ' - Signal Comparison'], 'Color', 'white', 'Units', 'normalized', 'OuterPosition', [0 0 0.5 0.5]);

    % Get the MATLAB default color order
    colors = lines(size(tissues, 1));  % Standard MATLAB colors for different tissues

    % Loop over TEeffective values
    for te_idx = 1:length(TEeffective_values)
        TEeffective = TEeffective_values(te_idx);
        
        subplot(2, ceil(length(TEeffective_values) / 2), te_idx);
        hold on;  % Hold on to plot multiple tissues on the same subplot
        
        % Loop over tissues
        for t = 1:size(tissues, 1)
            tissue_name = tissues{t, 1};
            data = results.(tissue_name).(sequence_name).(['TE' num2str(TEeffective)]);
            dataNoPrep = results.(tissue_name).('NoPrep').(['TE' num2str(TEeffective)]);
            
            % Plot each tissue with a different color from the standard color order
            plot(data.TI, abs(abs(data.Mx_case1) - abs(data.Mx_case2)) / dataNoPrep.Mx_case1(end), ...
                'Color', colors(t, :), 'LineWidth', 1.5);
        end
        
        % Set legend for tissues
        legend(tissues(:, 1), 'Location', 'northeast', 'Interpreter', 'latex');
        
        ylim([0, 1]);
        axis square;
        title(sprintf('$T_{2}$ Prep TE = %d ms', TEeffective), 'Interpreter', 'latex');
        xlabel('Inversion Time TI (ms)', 'Interpreter', 'latex');
        ylabel('Signal Magnitude', 'Interpreter', 'latex');
        grid on;
        box on
        % Set tick labels to use LaTeX interpreter
        set(gca, 'TickLabelInterpreter', 'latex');
    end

    % Sequence description for the plot title
    switch sequence_name
        case 'TimeSlipFirst'
            sequence_description = 'Time-SLIP $\rightarrow T_{2}\mathrm{Prep}$';
        case 'T2prepFirst'
            sequence_description = '$T_{2}\mathrm{Prep} \rightarrow$ Time-SLIP';
        case 'NoPrep'
            sequence_description = 'Time-SLIP';
        otherwise
            error('Invalid sequence order.');
    end

    % Set the overall title for the figure
    sgtitle(sequence_description, 'Interpreter', 'latex');
    export_fig(strcat(sequence_name,'_SIR.pdf'))
end



function results = store_results(results, tissue_name, TEeffective, sequence_order, TI_values, Mx_case1, My_case1, Mz_case1, Signal_case1, Mx_case2, My_case2, Mz_case2, Signal_case2)
    % Store results for both cases
    results.(tissue_name).(sequence_order).(['TE' num2str(TEeffective)]).TI = TI_values;
    results.(tissue_name).(sequence_order).(['TE' num2str(TEeffective)]).Signal_case1 = Signal_case1;
    results.(tissue_name).(sequence_order).(['TE' num2str(TEeffective)]).Mx_case1 = Mx_case1;
    results.(tissue_name).(sequence_order).(['TE' num2str(TEeffective)]).My_case1 = My_case1;
    results.(tissue_name).(sequence_order).(['TE' num2str(TEeffective)]).Mz_case1 = Mz_case1;
    results.(tissue_name).(sequence_order).(['TE' num2str(TEeffective)]).Signal_case2 = Signal_case2;
    results.(tissue_name).(sequence_order).(['TE' num2str(TEeffective)]).Mx_case2 = Mx_case2;
    results.(tissue_name).(sequence_order).(['TE' num2str(TEeffective)]).My_case2 = My_case2;
    results.(tissue_name).(sequence_order).(['TE' num2str(TEeffective)]).Mz_case2 = Mz_case2;
end

function [Mx_case1, My_case1, Mz_case1, Signal_case1, Mx_case2, ...
    My_case2, Mz_case2, Signal_case2] = initialize_arrays(len)

    % Initialize arrays to store results
    Mx_case1 = zeros(len, 1);
    My_case1 = zeros(len, 1);
    Mz_case1 = zeros(len, 1);
    Signal_case1 = zeros(len, 1);
    Mx_case2 = zeros(len, 1);
    My_case2 = zeros(len, 1);
    Mz_case2 = zeros(len, 1);
    Signal_case2 = zeros(len, 1);
end

function M = relaxation(M, T1, T2, dt, M0)
    % Applies relaxation over time dt
    E1 = exp(-dt / T1);
    E2 = exp(-dt / T2);
    Mx = M(1) * E2;
    My = M(2) * E2;
    Mz = M(3) * E1 + M0 * (1 - E1);
    M = [Mx; My; Mz];
end

function M = T2prep(M, T2, TEeffective, M0)
    % Simulates the T2prep module (MLEV4)
    M = rotation(M, 90, 'y');
    M = relaxation(M, Inf, T2, TEeffective / 2, M0);
    M = rotation(M, 180, 'x');
    M = relaxation(M, Inf, T2, TEeffective / 2, M0);
    M = rotation(M, -90, 'y');
end

function M = rotation(M, angle_deg, axis)
    % Rotates magnetization vector M by angle_deg around specified axis
    angle_rad = deg2rad(angle_deg);
    switch axis
        case 'x'
            R = [1, 0, 0; 0, cos(angle_rad), -sin(angle_rad); 0, sin(angle_rad), cos(angle_rad)];
        case 'y'
            R = [cos(angle_rad), 0, sin(angle_rad); 0, 1, 0; -sin(angle_rad), 0, cos(angle_rad)];
        case 'z'
            R = [cos(angle_rad), -sin(angle_rad), 0; sin(angle_rad), cos(angle_rad), 0; 0, 0, 1];
        otherwise
            error('Invalid rotation axis');
    end
    M = R * M;
end