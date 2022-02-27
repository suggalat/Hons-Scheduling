function [sel_pop] =selection(init_pop)
    fit_sum=0;
    p=size(init_pop,2);
    for i=1:p
        fit_sum = fit_sum+fitness(init_pop(:,i));
    end
    prob =zeros(1,p);
    cum_prob=zeros(1,p);
    for i = 1:p
        prob(i)=fitness(init_pop(:,i))/fit_sum;
        if(i>1)
            cum_prob(i)=prob(i)+cum_prob(i-1);
        else
            cum_prob(i) = prob(i);
        end        
    end
    sel_pop=zeros(size(init_pop));
    for i=1:p
       rand_val = rand(1,1);
       for j=1:p
           if(cum_prob(j)>rand_val)
               sel_pop(:,i)=init_pop(:,j);
           end
       end
    end
end