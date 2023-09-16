function [feq,bin]=get_dR_PDF(xys,tlag,param)
% Calculate the probability density funciton (PDF) of displacements    
% 
% Syntax:
%   get_dR_PDF;
%   [feq,bin]=get_dR_PDF;
%   [feq,bin]=get_dR_PDF(xys);
%   [feq,bin]=get_dR_PDF(xys,tlag);
%   [feq,bin]=get_dR_PDF(xys,tlag,param);
% 
% input:
%   xys (optional): a N*1 cell array, Each cell contains individual cell trajectory data matrice. 
%   The format of trajectory matrice must be in the the sequence of cell id, time frame, x, y, z ;      
%   Each trajectories  must have the same time length, Nt.
%   if xys is not provided, the "get_trajfile.m" will be called to
%   select and input trajectories file. (see get_trajfile for more
%   details). 
%   tlag  (optional): frame lag for calculate displacements; 
%   if not included, input the desired value through pop-up window. 
%   param (optional) : parameters setting. if not included, default values of setting will be used.
%    param.dim : dimensionality of trajectories (either 2 or 3) (default: 2);
%    param.showfig : if true, displacement-PDF result will be plotted into a figure. 
%    default: true).
%    param.markertype : Specify marker types in generated displacement-PDF graph
%    (default: 'bo').
%    param.outfigurenum : Specify figure number (postive integer) for displacement-PDF
%    graph(default : 302).
%    param.saveres : if true, displacement-PDF and displacements of trajectories 
%    will be output to user-specified excel file. (default: true).  
%    param.dxmax : maximum displacment value for calculation of PDF)
%    (default: 50)
%    param.binn  : number of displacement bins between 0 and param.dxmax for calculation of PDF)
%    (default: 70)
%
% output: 
%   bin (optional): Displacement bins ( a monotonic vector)  
%   feq (optional): Probability of displacement occurrence at difffernet displacement bins   
%   
% 
% developed  by-
%  Pei-Hsun Wu, Ph.D. 
%  Johns Hopkins University
%
%%%%%%%%%%%%%%%%%%%

% main program
% check input variable% 
        if nargin==0
            xys=get_trajfile;    
        end
        if isempty(xys)
            xys=get_trajfile;    
        end
       if nargin <=1
           answer=inputdlg('input frame lag to compute dR','input frame lag');
           tlag=str2double(answer);
       end
        if  nargin <=2 % set the bin size
            param.dxmax=50;
            param.binn=70;            
            param.showfig=1;
            param.saveres=1;
            param.dim=2;            
            param.outfigurenum=302;
            param.markertype='bo';
        end
        dxyout=[];     
         for k=1:length(xys);  
                xy=xys{k};         
                dxy=xy(1+tlag:end,1:param.dim)-xy(1:end-tlag,1:param.dim) ;  
                dxyout=[dxyout,dxy];             
         end     
            binw=param.dxmax/param.binn;
            bin=linspace(binw/2,param.dxmax,param.binn);
            [feq,bin]=hist(abs(dxyout(:)),bin);
            cc=feq>=0; % highlight all.
            feq=feq/sum(feq)/binw;          
            bin=bin(cc);
            feq=feq(cc);
            bin=bin(:);
            feq=feq(:);    
            
             %%% plot result
            if param.showfig
                if param.aviv
                    feq = smooth(feq, 15);
                    if param.aviv2 == 1
                        figure(param.outfigurenum); 
                        semilogy(bin,feq,'--s','color',param.markertype,'linewidth',2); 
                        xlim([0 60])
                        bjff3;      
                        hold on;
                    elseif param.aviv2 == 2
                        figure(param.outfigurenum); 
                        semilogy(bin,feq,'--x','color',param.markertype,'linewidth',2); 
                        xlim([0 60])
                        bjff3;      
                        hold on;
                    else
                        figure(param.outfigurenum); 
                        semilogy(bin,feq,'--','color',param.markertype,'linewidth',2); 
                        xlim([0 60])
                        bjff3;      
                        hold on;
                    end
                else
                    
                    figure(param.outfigurenum); 
                    if param.aviv3 == 1
                    semilogy(bin,feq,param.markertype,'linewidth',2);
                    else
                    semilogy(bin,feq,'color',param.markertype,'linewidth',2);
                    end
                    xlim([0 60])
                    bjff3;      
                    hold on;
                end
            end
            %%% output results
            if param.saveres
             [filename, pathname] = uiputfile( ...       
                 {'*.xlsx',  'excel files (*.xlsx)'; ...
                   '*.xls','excel file (*.xls)'}, ...             
                   'output pdf of displacement','PDF-dR.xlsx');            
               
                xlswrite([pathname,filename],[bin feq],'pdf data');
                xlswrite([pathname,filename],[dxyout],'displacement data');
                
                delete_extra_sheet(pathname,filename)
             end
    
            if nargout==0
                clear
            end
    
end