function [nx ny output outindsx outindsy] = regrid(xaxis,yaxis,data,nx,ny,subsample_largergrid,interp_type);
% (C) Nick Holschuh - Penn State University - 2016 (Nick.Holschuh@gmail.com)
% This does a 2d interpolation for gridded data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The inputs are as follows:
%
% xaxis - The xaxis for the input data grid
% yaxis - The yaxis for the input data grid
% data - The input data grid
%
% There are several ways to use the following inputs:
% 1) nx - the desired xaxis for the output grid
%    ny - the desired yaxis for the output grid
%
% 2) nx - set to 0 - the output data will have equal grid spacing in the x
%         and y direction, set to be the finer grid spacing of the original data
%    ny - [];
%
% 2) nx - set to 1 - the data is assumed to be ice penetrating radar data,
%           and the source yaxis is assumed to be time in seconds. This is
%           immediately converted to depth, and the data are interpolated to an even
%           spacing in the x and y direction that is defined by spacing that is 4x
%           the nyquist frequency of the target data (sample rate provided in ny)
%    ny - center frequency of the radar data.
%
% subsample_largergrid - flag, 0 or 1, to indicate if you want the
%   regridded product to only take the samples of the larger grid that
%   overlap.
% interp_type - accepts 'spline', 'linear','nearest','next','cubic'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
skipinterp = 0;

if min(size(xaxis)) > 1
    matflag = 1;
else
    matflag = 0;
end

if exist('subsample_largergrid') == 0
    subsample_largergrid = 0;
end

if exist('interp_type') == 0
    interp_type = 'linear';
end



%%%%%%%%%%%%%%%%%%%%%%%%%%% This only maters in interpolating radar data
%%%%%%%%%%%%%%%%%%%%%%%%%%% specifically
if length(nx) == 1
    % This automatically regrids the data to the finer spacing of the two
    % source spacings
    xstep = xaxis(2)-xaxis(1);
    ystep = yaxis(2)-yaxis(1);
    
    %%%%%%%%%%%%%%%%%%%%%%%% This was for interpolating radargrams,
    %%%%%%%%%%%%%%%%%%%%%%%% specifically.
    %%% NOTE: only the STEP is converted to a length here, not yaxis
    %%% itself, so the y axis is still two-way travel time downstream.
    %%% Everything below has to keep that distinction straight.
    y_is_twtt = (nx == 1 & ystep < 1e-5);
    if y_is_twtt
        cice_import
        ystep = ystep*cice/2;
    end

    if xstep > ystep;
        nx_temp = xaxis(1):ystep:xaxis(end);
        ny_temp = yaxis;
    elseif ystep > xstep
        nx_temp = xaxis;
        ny_temp = yaxis(1):xstep:yaxis(end);
    else
        skipinterp = 1;
        nx_temp = xaxis;
        ny_temp = yaxis;
    end

    % This regrids the data to ~20 samples per wavelength in ice.
    %
    % FIX: the target spacing is a LENGTH, (c_ice/f)/20. The original used
    % (1/f)/20, which is a TIME, and then applied it to both axes. On a
    % distance axis that asks for (track length)/(8e-11 s) samples -- for a
    % 2.5 km line, 3e13 elements -- so the call died in the allocator
    % before the guard below could reject the grid. The x axis therefore
    % gets the length step, and the y axis gets it converted back into
    % seconds whenever y is still two-way travel time.
    if nx == 1
        f = ny;
        cice_import
        dl_target = (cice/f)/20;             % target spacing, metres
        if y_is_twtt
            dy_target = dl_target*2/cice;    % metres -> seconds
        else
            dy_target = dl_target;
        end

        % Size the candidate grids arithmetically, so an unreasonable
        % request is rejected by the guard below rather than attempted.
        n_x2 = floor((xaxis(end)-xaxis(1))/dl_target)+1;
        n_y2 = floor((yaxis(end)-yaxis(1))/dy_target)+1;

        if n_x2 <= length(nx_temp) | n_y2 <= length(ny_temp)
            nx_temp = xaxis(1):dl_target:xaxis(end);
            ny_temp = yaxis(1):dy_target:yaxis(end);
        end
    end

    nx = nx_temp;
    ny = ny_temp;    
end

if skipinterp == 0
    
    if matflag == 0
        [y_0 x_0] = ndgrid(yaxis,xaxis);
        
        if subsample_largergrid == 1
            outindsx = find_nearest(nx,min(xaxis)):find_nearest(nx,max(xaxis));
            outindsy = find_nearest(ny,min(yaxis)):find_nearest(ny,max(yaxis));
            nx = nx(outindsx);
            ny = ny(outindsy);
        else
           outindsx = 1:length(nx);
           outindsy = 1:length(ny);
        end
        
        [interp_y interp_x] = ndgrid(ny,nx);
    else
        y_0 = yaxis;
        x_0 = xaxis;
        interp_y = ny;
        interp_x = nx;
        
        outindsx = NaN;
        outindsy = NaN;
    end
    
    A = griddedInterpolant(y_0,x_0,data,interp_type);
    output = A(interp_y,interp_x);
    if min(min(isnan(output))) == 1
        A = griddedInterpolant(y_0,x_0,data,'linear');
        output = A(interp_y,interp_x);
    end
else
    output = data;
end
    
end