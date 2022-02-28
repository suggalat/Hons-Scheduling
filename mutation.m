function [c2] = mutation(c1,prob_mut)
    c2=c1;
    load('ga.mat');
    if(prob_mut>rand(1,1))
        rand_index=randi(N);
        if(c2(rand_index)==-1)
            c2(rand_index)=1;
        else
            c2(rand_index)=-1;
        end
    end
end