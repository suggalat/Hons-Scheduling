function [cross_c1,cross_c2] =crossover(c1,c2,prob_cross)
    load('ga.mat');
    rand_val=rand(1,1);
    cross_c1=c1;
    cross_c2=c2;
    if(rand_val<prob_cross)
        cross_pos = randi(N-1);
        cross_c1 = [c1(1:cross_pos); c2(cross_pos+1:N)];
        cross_c2 = [c2(1:cross_pos); c1(cross_pos+1:N)];        
    end
end
