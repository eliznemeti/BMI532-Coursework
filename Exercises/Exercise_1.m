clear
close all
clc

% Question 1.1

 DATA = LogisticMap(4, 0.5, 10000);
 find_bifurcation(DATA)
 
 function graph_values = LogisticMap(r, x0, N)
     graph_values=[];
     for r = linspace(0,r,N)        
         xold = x0;                     
        
         % checking for the steady state
         for i=1:2000 
             xnew=((xold-xold^2)*r);
             
             xold=xnew;
         end
         
         xss=xnew;
         for i=1:1000
             xnew=((xold-xold^2)*r);
             xold=xnew;
             
             % forming the matrix to create the plot with
             graph_values(1,length(graph_values)+1)=r;
             graph_values(2,length(graph_values))=xnew;
             
             if(abs(xnew-xss)<.0001)
                 break
             end
         end
     end
 
     % plot the values for a bifurcation plot and add details
     plot(graph_values(1,:), graph_values(2,:), '.', 'LineWidth', .1, 'MarkerSize',1.2,...
     'Color',[0 0 0])                   %the 1 1 1 vector represents white color
     set(gca, 'color', 'w', 'xcolor', 'k', 'ycolor', 'k')
     set(gcf, 'color', 'w')
     xlabel('r')
     ylabel('Steady-State (x*)')
 
     savefig('bifurcation.fig')
 
 end
 
 %% Question 1.2
 
 function bifurcation_list = find_bifurcation(data)
 bifurcation_list = [];
 i=1;
     for j=data(1,3:end) % start at 3 to avoid the first 1 columnds of 0 counting as a bifurcation 
         if isempty(bifurcation_list)
             if abs(j - data(1,i+1)) == 0 % this if statement identifies the bifurcation point, and the values where the branches emerge
                 bifurcation_list = j;
             end
         elseif length(bifurcation_list) == 1
             if (abs(j - data(1,i+4)) == 0) % this if statement identifies the 2nd bifurcation point, and the values where 4 branches now emerge
                 bifurcation_list = [bifurcation_list j];
             end
         else
             break
         end
         i = i+1; 
     end
 
 
 end


