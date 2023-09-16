function acf=get_ACF(xys,dt,param)               
% Compute the velocity auto correlation function of trajectories
% 
% Syntax:
%   get_ACF;
%   acf=get_ACF;
%   acf=get_ACF(xys);
%   acf=get_ACF(xys,dt);
%   acf=get_ACF(xys,dt,param);
% 
% input:
%   xys (optional): a N*1 cell array, Each cell contains individual cell trajectory data matrice. 
%   The format of trajectory matrice must be in the the sequence of cell id, time frame, x, y, z ;      
%   Each trajectories  must have the same time length, Nt.
%   if xys is not provided, the "get_trajfile.m" will be called to
%   select and input trajectories file. (see get_trajfile for more
%   details). 
%   dt  (optional): time step size of trajectory data,
%   if not included, input the desired value through pop-up window. 
%   param (optional) : parameters setting , if not included a default value
%   of settting will be used.
%    param.dim :dimensionality of trajectories (either 2 or 3) (default: 2);
%    param.tlag : frame lag for calculate displacements for ACF (default: 1); 
%    param.showfig : if true, ACF result will be plotted into a figure. 
%    default: true).
%    param.markertype : Specify marker types in generated ACF graph
%    (default: 'bo').
%    param.outfigurenum : Specify figure number (postive integer) for ACF
%    graph(default : 303).
%    param.saveres : if true, ACFs of will be output to user-specified
%    excel file. Both individual ACFs and ensemble-averaged ACF will be saved in two different sheet where 
%    first column is time lags. (default: true).  
% 
% output: 
%   acf (optional): velocity autocorrelation function of trajectories, a tt * N matrix  
%   where tt is depdent on param.tlag. Maximum of tt is equal to Nt-2 only when param.tlag is set to 1;
%   Each column is acf profile of individual trajectories over time lag (1:tt-1)*dt*param.tlag 
% 
% developed  by-
%  Pei-Hsun Wu, Ph.D. 
%  Johns Hopkins University
%
%%%%%%%%%%%%%%%%%%%

% main program
% check input variable
     if nargin==0
        xys=get_trajfile;
     end
    if nargin<=1;
         answer=inputdlg('Trajectory time step size','input time step size');
         dt=str2double(answer);
        
    end
    if nargin<=2;
        param.dim=2;
        param.tlag=1;    
        param.saveres=1;
        param.showfig=1;       
        param.markertype='bo';
        param.outfigurenum=303;           
    end
    
    acf=[];
    tlag=param.tlag;
    for k=1:length(xys);
        xy=xys{k};
        xy=xy(1:tlag:end,:);
        dxy=xy(2:end,1:param.dim)-xy(1:end-1,1:param.dim) ;           
        M=length(dxy(:,1));       
        cef=[]; tll=[];
        for k=0:M-1
            cef(k+1)=sum(sum(dxy(1:M-k,:).*dxy(1+k:M,:),2))/(M-k);  % using inner dot (taking vector)           
            tll(k+1)=k; % unit of time lag
        end 
        acf=[acf, cef(:)] ;  
    end    
    
    if param.showfig
       figure(param.outfigurenum);             
           m33=mean(acf,2);  % average profile of ACF      
           m33=m33/m33(3);  % normalization           
           tll=tll(:); % unit of time lag
           cc0=tll>=2 & tll< 50  ; % set time lag range, get rid of first two points       
           m33=m33(cc0); tll0=tll(cc0);
           if param.aviv
%                 m33 = smooth(m33, 33);
                m33 = abs(m33);
% Fit a fourth-order polynomial to the data
p = polyfit(tll0*dt*tlag, log(m33), 5);

% Generate fitted curve points
fitted_curve_x = linspace(min(tll0*dt*tlag), max(tll0*dt*tlag), 100);
fitted_curve_y = exp(polyval(p, fitted_curve_x));

% Normalize the fitted curve so that it starts at 1
% fitted_curve_y = fitted_curve_y / fitted_curve_y(1);

% Plot the fitted curve
if param.aviv3
    semilogy(fitted_curve_x, fitted_curve_y, 'color',param.markertype,'Marker',param.markertype2, 'LineWidth', 1);
else
    semilogy(fitted_curve_x, fitted_curve_y, param.markertype, 'LineWidth', 1);
end
           end
%            semilogy(tll0*dt*tlag,(m33),param.markertype,'linewidth',2); hold on;

% % Fit a linear function to the log-transformed y values
% fit_linear = fit(tll0*dt*tlag, log(m33), 'poly1');
% 
% % Generate fitted line points
% fitted_line_x = linspace(min(tll0*dt*tlag), max(tll0*dt*tlag), 100);
% fitted_line_y = exp(feval(fit_linear, fitted_line_x));
% 
% % Plot the fitted line
% semilogy(fitted_line_x, fitted_line_y, '-', 'LineWidth', 1);
%            xlim([0 750])
%            ylim([1e-3 1e0])
           
            bjff3;  
           hold on;    
    end
    
    if param.saveres
          [filename, pathname] = uiputfile( ...       
         {'*.xlsx',  'excel files (*.xlsx)'; ...
           '*.xls','excel file (*.xls)'}, ...             
           'output ACF to a excel file','ACF.xlsx');               
        xlswrite([pathname,filename],[tll(:)*dt*tlag acf],'individual ACF');
        xlswrite([pathname,filename],[tll(:)*dt*tlag ,mean(acf,2)],'average ACF');
        delete_extra_sheet(pathname,filename);
    end
    
    if nargout==0
        clear acf
    end
end
    
