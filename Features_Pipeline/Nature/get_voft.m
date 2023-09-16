function [feq,bin]=get_voft(xys,tlag0,param)
%  calculate the velocity over time
        if nargin==0
            xys=get_trajfile;    
        end
        if isempty(xys)
            xys=get_trajfile;    
        end
       if nargin <=1
           tlag0=1;
       end
        if  nargin <=2 % set the bin size             
            param.showfig=1;
            param.saveres=1;
            param.dim=2;
            param.outfigurenum=300;
            param.markertype='b-';
        end
        vout=[];
         for k=1:length(xys);  
                xy=xys{k};         
                dxy=xy(1+tlag0:end,1:param.dim)-xy(1:end-tlag0,1:param.dim) ;  
                [~,dr]=cart2pol(dxy(:,1),dxy(:,2)) ;
                dr=dr/tlag0;
                vout=[vout,dr];             
         end     
                        
            if param.saveres
             [filename, pathname] = uiputfile( ...       
                 {'*.xlsx',  'excel files (*.xlsx)'; ...
                   '*.xls','excel file (*.xls)'}, ...             
                   'save average velocity profile','velocity over time.xlsx');                
                
                xlswrite([pathname,filename],[mean(vout,2)],'mean velocity over time');
                delete_extra_sheet(pathname,filename)
             end
            if param.showfig
                figure(param.outfigurenum); 
                plot(mean(vout,2),param.markertype);    
                bjff3;      
                hold on;
            end
            if nargout==0
                clear
            end
    
end