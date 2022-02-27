function [c2] = mutation(c1,prob_mut)
    c2=c1;
    if(prob_mut>rand(0,1))
        rand_index=randi(N);
        if(c2(rand_index)==-1)
            c2(rand_index)=1;
        else
            c2(rand_index)=-1;
        end
    end
end