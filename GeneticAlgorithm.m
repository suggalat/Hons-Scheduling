
function [bestConfig]= GeneticAlgorithm(g,p,N)
  
    had_matrix = hadamard(N);
    init_pop = had_matrix(:,1:p);
    bestConfig=zeros(N,1);
    bestFitness = 0;
    for g_i =1:g 
        for b=1:p
            currentFitness = fitness(init_pop(:,b));
            if(currentFitness>bestFitness)
                bestFitness = currentFitness;                
                bestConfig = init_pop(:,b)               
            end            
        end
        sel_pop=selection(init_pop);    
        cross_pop=sel_pop;
        prob_cross=0.5;
        for i=1:2:p
           [cross_pop(:,i),cross_pop(:,i+1)] = crossover(sel_pop(:,i),sel_pop(i+1),prob_cross);
        end
        mut_pop =cross_pop;
        prob_mut=0.5;
        for i=1:p
            mut_pop(:,i) = mutation(cross_pop(:,i),prob_mut);
        end
        init_pop=mut_pop
    end 

end

 