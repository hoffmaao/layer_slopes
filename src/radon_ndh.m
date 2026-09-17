function [optimum_angle rad_result angles x_rad radon_snr] = radon_ndh(xaxis,yaxis,data,angle_thresh,plotter,method);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Method Descriptions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) This method determines the optimal angle by searching for the angle
% stack with the maximum amplitude peak (corrects for number of data
% samples in stack and edge effects using a gaussian).
%
% 2) This method determines the optimal angle by maximizing the return
% power for a given stack (only would work if looking at amplitudes that
% oscillate around 0)
%
% 3) This method looks at the variance for stacks at different angles, and
% selects the angle that produces the stack with the maximum variance
% (Normalized values, and corrected for number of data samples in a stack,
% edge effects using a boxcar).
%
% 4) This method determines the optimal angle by comparing the normalized
% slant stack with the normalized original trace. The angle that minimizes
% the difference between them is selected as the optimal trace.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


pause_duration = 0.1;

if exist('angle_thresh') == 0
    angle_thresh = 0;
end

if exist('plotter') == 0
    plotter = 0;
end

if exist('method') == 0
    method = 4;
end

if angle_thresh == 0
    angle_thresh = 90;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Horizontal Interpolation to set Aspect Ratio to 1:1 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% FIX: griddedInterpolant below requires ascending sample points, so a
%%% descending y axis (elevation rather than depth) threw "Sample points
%%% must be sorted in ascending order" before ever reaching the
%%% time_increases_downward == 0 branch that was meant to handle it. Solve
%%% in the ascending frame and mirror the answer back: flipping the y axis
%%% is a pure reordering: it flips the axis and the data together, so the
%%% (y, value) association is preserved and the dip is unchanged.
y_was_descending = 0;
if numel(yaxis) > 1 && yaxis(2) < yaxis(1)
    y_was_descending = 1;
    yaxis = fliplr(yaxis(:).');
    data = flipud(data);
end

xstep = xaxis(2)-xaxis(1);
ystep = yaxis(2)-yaxis(1);

if xstep > ystep;
    new_x = xaxis(1):ystep:xaxis(end);
    new_y = yaxis;
    [y_0 x_0] = ndgrid(yaxis,xaxis);
    [interp_y interp_x] = ndgrid(new_y,new_x);
    A = griddedInterpolant(y_0,x_0,data,'spline');
    new_im = A(interp_y,interp_x);
elseif ystep > xstep
    new_x = xaxis;
    new_y = yaxis(1):xstep:yaxis(end);
    [y_0 x_0] = ndgrid(yaxis,xaxis);
    [interp_y interp_x] = ndgrid(new_y,new_x);
    A = griddedInterpolant(y_0,x_0,data,'spline');
    new_im = A(interp_y,interp_x);
end

if ystep ~= xstep
    xaxis = new_x;
    yaxis = new_y;
    data = new_im;
    xstep = xaxis(2)-xaxis(1);
    ystep = yaxis(2)-yaxis(1);
end



%%

if yaxis(end) > yaxis(1)
    time_increases_downward = 1;
else
    time_increases_downward = 0;
end


d_theta = 0.1;

if time_increases_downward == 1
    %angles = 180-angle_thresh:1:180+angle_thresh; % Left looking, to down looking, to right looking
    angles = 90-angle_thresh:d_theta:90+angle_thresh; % Top looking, to left looking, to down looking;
    dip_angles = angles-90;
    dip_angles(find(dip_angles > 90)) = -1*(180-dip_angles(find(dip_angles > 90))); 
    
else
    %angles = [360-angle_thresh:1:360 1:1:angle_thresh]; % Left looking, to down looking, to right looking
    angles = 270-angle_thresh:d_theta:270+angle_thresh; % Up looking, to left looking, to down looking
    dip_angles = angles-270;
    dip_angles(find(dip_angles < -90)) = dip_angles(find(dip_angles < -90))+180; 
    
end

%%% FIX: sign convention. As computed above, dip_angles is the NEGATIVE of
%%% the geologic dip: a layer that deepens with increasing x came back with
%%% a negative angle, verified against synthetics in tests/test_sign.m.
%%% Nick's full NDH_MatlabTools version of RollingRadon compensates for
%%% this with `slopegrid = slopegrid*-1` at the very end, but the public
%%% SlopeExtraction_Radar release dropped that line, so every slope it
%%% returns is sign-flipped. Fix it here instead, so the function's own
%%% contract is right and anyone calling radon_ndh directly gets the
%%% documented convention:
%%%
%%%     POSITIVE dip = the layer gets DEEPER with increasing x.
%%%
dip_angles = -dip_angles;

real_angles = angles;


% Compute the correction function for the final radon transform (removes
% the contribution of the number of cells in the stack), plus the window
% filters.
%
%%% SPEED: none of these depend on the DATA, only on its size and the angle
%%% list, yet the released code rebuilt all of them - including two extra
%%% Radon transforms and a per-column normpdf loop - on every single call.
%%% In a rolling-window solver that is the same work tens of thousands of
%%% times over. Cache them against (size, angles); the arithmetic is
%%% unchanged.
persistent c_key c_rad_corr c_gauss0 c_circle c_radgauss0 c_gauss c_boxcar
key = [size(data) real_angles(1) real_angles(end) numel(real_angles)];

if isempty(c_key) || ~isequal(key, c_key)
    correction_image = ones(size(data));
    rad_correction = radon(correction_image,real_angles);
    rad_correction = value2value(rad_correction,0,1);

    % Generate the gauss filter for method 0;
    filt_sd = 2.5;
    gauss_a = normpdf(-filt_sd:filt_sd*2/(length(data(1,:))-1):filt_sd);
    gauss_b = normpdf(-filt_sd:filt_sd*2/(length(data(:,1))-1):filt_sd);
    gauss_filter_0 = gauss_b'*gauss_a;
    gauss_filter_0(find(gauss_filter_0 < 0.01)) = 0;
    circle_filter = gauss_filter_0;
    circle_filter(find(circle_filter > 0)) = 1;

    rad_correction_gaussian0 = radon(gauss_filter_0,real_angles);
    rad_correction_gaussian0 = value2value(rad_correction_gaussian0,0,1);

    % Generate the gauss filter for method 1;
    gauss_filter = zeros(size(rad_correction));
    filt_sd = 1.6;
    for i = 1:length(rad_correction(1,:))
        indecies = find(rad_correction(:,i) ~= 1);
        filt_length = indecies(end)-indecies(1);
        gauss_filter(indecies,i) = normpdf(-filt_sd:filt_sd*2/filt_length:filt_sd);
    end
    boxcar_filter = gauss_filter;
    filt_val = 0.1;
    boxcar_filter(find(gauss_filter > filt_val)) = 1;
    boxcar_filter(find(gauss_filter <= filt_val)) = 0;

    c_key = key;
    c_rad_corr = rad_correction;
    c_gauss0 = gauss_filter_0;
    c_circle = circle_filter;
    c_radgauss0 = rad_correction_gaussian0;
    c_gauss = gauss_filter;
    c_boxcar = boxcar_filter;
else
    rad_correction = c_rad_corr;
    gauss_filter_0 = c_gauss0;
    circle_filter = c_circle;
    rad_correction_gaussian0 = c_radgauss0;
    gauss_filter = c_gauss;
    boxcar_filter = c_boxcar;
end
column_filter = ones(length(data(:,1)),1)*sum(data,1);

% FIX: this divided the shifted data by the max of the UNSHIFTED data, so
% for log-power input (every value negative) the divisor was negative and
% the whole window came back sign-flipped. Methods 3 and 4 normalise each
% Radon column by its own max, which is not invariant to that flip, so a
% dB radargram could select a different angle than the same data in linear
% units. Normalise onto [0,1] properly instead.
data_min = min(min(data));
data_max = max(max(data));
if data_max > data_min
    data = (data - data_min)/(data_max - data_min);
else
    data = zeros(size(data));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SNR Evaluation                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimization Criteria for the Radon Transform %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if method == 0
    %[rad_result x_rad] = radon(data.*gauss_filter_0,real_angles); %Pre_filter data with gaussian
    [rad_result x_rad] = radon(data.*circle_filter,real_angles); %Pre_filter data with circle window
    
    %[rad_result x_rad] = radon(data./column_filter.*gauss_filter_0,real_angles); %Pre_filter data with gaussian and column_filter
    %rad_result = rad_result./rad_correction_gaussian0;
    [xcol_final ycol_final] = find(max(max(abs(rad_result))) == abs(rad_result));
    
    optimum_angle = dip_angles(ycol_final);
    

%     %%%%% Restrict based on certainty
%     [temp_rad] = radon(mean(mean(data))*ones(size(data)).*gauss_filter_0,real_angles);
%     temp_rad = temp_rad./rad_correction;
% 
%     pick_val = max(max(abs(rad_result-temp_rad)));
%     group_mean = mean(max(abs(rad_result-temp_rad),[],1));
%     group_std = std(max(abs(rad_result-temp_rad),[],1));
%     
%     error_thresh = 3;
%     if pick_val < group_mean + error_thresh*group_std
%         optimum_angle = NaN;
%     end
    
elseif method == 1
[rad_result x_rad] = radon(data,real_angles);
    rad_result = rad_result./rad_correction.*gauss_filter;
    [xcol_final ycol_final] = find(max(max(abs(rad_result))) == abs(rad_result));
    optimum_angle = dip_angles(ycol_final);
    
elseif method == 2
[rad_result x_rad] = radon(data,real_angles);
    optimum_angle = dip_angles(find(sum(rad_result,1) == max(sum(rad_result,1))));
    
elseif method == 3
    [rad_result x_rad] = radon(data-mean(mean(data)),real_angles);
    rad_result = rad_result./rad_correction.*boxcar_filter;
    for i = 1:length(rad_result(1,:))
        rad_result(:,i) = rad_result(:,i)/max(rad_result(:,i));
    end
    optimum_angle = dip_angles(find(var(rad_result,0,1) == max(var(rad_result,0,1))));
    
elseif method == 4
    [rad_result x_rad] = radon(data-mean(mean(data)),real_angles);
    rad_result = rad_result;%./rad_correction;
    for i = 1:length(rad_result(1,:))
        rad_result(:,i) = rad_result(:,i)/max(rad_result(:,i));
    end
    sample_length = length(data(:,1));
    length_dif = length(rad_result(:,1))-sample_length;
    comp_trace = data(:,round(length(data(1,:))/2));
    comp_trace = comp_trace/max(comp_trace);
    trace_range = round(length_dif/2):round(length_dif/2)+length(comp_trace)-1;
    rad_result = rad_result(trace_range,:);
    %rad_result = rad_result - comp_trace*ones(1,length(rad_result(1,:)));
    optimum_angle = dip_angles(find(sum(rad_result,1) == max(sum(rad_result,1))));
end

%%% FIX: the criteria above use find(x == max(x)), which returns SEVERAL
%%% indices when angles tie and NONE at all when the window is featureless
%%% (max of an all-NaN vector is NaN, and NaN == NaN is false). The callers
%%% assign into a scalar element, so a tie raised a size error and an empty
%%% result raised "Index exceeds array bounds". Take the first candidate,
%%% and abstain with NaN when there is no candidate - RollingRadon already
%%% treats NaN as "no slope here".
if isempty(optimum_angle)
    optimum_angle = NaN;
else
    optimum_angle = optimum_angle(1);
end

radon_snr = max(max(abs(rad_result))) - min(max(abs(rad_result)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Plot Portion - Ignore if plotter == 0  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if plotter == 1
    % Find the points that define the margins
    if length(xaxis(:,1)) > 1
        xaxis = xaxis';
    end
    if length(yaxis(:,1)) > 1
        yaxis = yaxis';
    end
    margins = combvec(xaxis,yaxis)';
    margins1 = [];
    margin_angles = [];
    for i = 1:length(margins)
        if margins(i,1) ~= min(xaxis) & margins(i,1) ~= max(xaxis) & margins(i,2) ~= min(yaxis) & margins(i,2) ~= max(yaxis)
        else
            margins1 = [margins1; margins(i,:)];
            margin_angles = [margin_angles; rad2deg(atan2(margins1(end,2)-mean(yaxis),margins1(end,1)-mean(xaxis)))];
        end
    end
    
    if time_increases_downward == 1;
        margin_angles = -margin_angles;
    end
    
    for i = 1:length(margin_angles)
        if margin_angles(i) < 0
            margin_angles(i) = 360-margin_angles(i)*-1;
        end
    end
    
    rad_result_sums = sum(rad_result,1);
    rad_result_var = var(rad_result,0,1);
    rad_result_compare = sum(rad_result,1);
    %Initialize the plot
    subplot(4,1,1:2)
    imagesc(xaxis,yaxis,data)
    set(gca,'YDir','reverse')
    axis equal
    colormap(gray)
    hold all
    
    if exist('trace_range') == 0
        trace_range = 1:length(x_rad);
    end
    
    for i = 1:length(angles)
        % Find the points that edge points that define a line with the
        % slope of interest.
        edge_angles = [angles(i) angles(i)+90 angles(i)+180 angles(i)+270];
        % Converts angles to a 0 to 360 scale
        edge_angles(find(edge_angles < 0)) = edge_angles(find(edge_angles < 0))*-1 + 360;
        edge_angles(find(edge_angles > 360)) = edge_angles(find(edge_angles > 360)) - 360;
        
        edge_ind(1) = find_nearest(margin_angles,edge_angles(1));
        edge_ind(2) = find_nearest(margin_angles,edge_angles(2));
        edge_ind(3) = find_nearest(margin_angles,edge_angles(3));
        edge_ind(4) = find_nearest(margin_angles,edge_angles(4));
        
        
        subplot(4,1,1:2)
        line1 = plot(margins1([edge_ind(2) edge_ind(4)],1),margins1([edge_ind(2) edge_ind(4)],2),'Color','blue','LineWidth',1.5);
        line2 = plot(margins1(edge_ind(2),1),margins1(edge_ind(2),2),'p','Color','blue','MarkerFaceColor','blue','MarkerSize',15);
        xlim([min(xaxis) max(xaxis)])
        ylim([min(yaxis) max(yaxis)])
        
        subplot(4,1,3)
        line3 = plot(x_rad(trace_range),rad_result(:,i),'Color','blue','LineWidth',2);
        if min(min(rad_result)) ~= max(max(rad_result))
            ylim([min(min(rad_result)) max(max(rad_result))])
        end
        
        subplot(4,2,7)
        hold off
        name1 = text(0,0,num2str(dip_angles(i)),'FontSize',24);
        xlim([-1 2])
        ylim([-1 1])
        
        subplot(4,2,8)
        if method == 2
            name2 = plot(angles(1:i),rad_result_sums(1:i),'Color','blue','LineWidth',2);
            if min(min(rad_result)) ~= max(max(rad_result))
                ylim([min(min(rad_result)) max(max(rad_result))])
            end
        elseif method == 3
            name2 = plot(angles(1:i),rad_result_var(1:i),'Color','blue','LineWidth',2);
            if min(min(rad_result)) ~= max(max(rad_result))
                ylim([0 max(max(rad_result))])
            end
        elseif method == 4
            name2 = plot(angles(1:i),rad_result_compare(1:i),'Color','blue','LineWidth',2);
            if min(min(rad_result)) ~= max(max(rad_result))
                ylim([min(min(rad_result_compare)) max(max(rad_result_compare))])
            end
        end
        
        xlim([min(angles) max(angles)])
        
        
        
        pause(pause_duration)
        
        if i < length(angles)
            subplot(4,1,1:2)
            delete(line1)
            delete(line2)
            subplot(4,1,3)
            if angles(i) ~= 0
                delete(line3)
            end
            subplot(4,2,7)
            delete(name1)
            if method > 1
                subplot(4,2,8)
                delete(name2)
            end
        end
        
    end
end

end