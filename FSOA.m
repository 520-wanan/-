  function [Best_score, Best_pos, curve] = FSOA(N, T, lb, ub, D, objective)
    habitat.Position = [];
    habitat.Cost = [];
    VarSize = [1, D]; 
    if numel(lb) == 1
    lb = repmat(lb, 1, D);
    end
    if numel(ub) == 1
    ub = repmat(ub, 1, D);
    end

    pop = repmat(habitat, N, 1);

    for i = 1:N
        pop(i).Position = unifrnd(lb, ub, VarSize); 
        pop(i).Cost = objective(pop(i).Position); 
    end

    [~, SortOrder] = sort([pop.Cost]);
    pop = pop(SortOrder);

    BestCost = zeros(T, 1);
    BestCost(1) = pop(1).Cost;
    Best_pos = pop(1).Position;

    % main loop
    for it = 2:T
        newpop = pop; 
        X_best = pop(1).Position;
        X_SecondBest = pop(2).Position;
        
        for i = 1:N
            r = rand();
            r1 = rand();
            r2 = rand();
            r3 = rand();
            r4 = rand();
            r5 = rand();
            r6 = rand(); 
            op1 = rand();
            op2 = trnd(10);
            op3 = levyFlight(1);
            op4 = 0.1*gamma(r);
            S=3 / (1 + exp(-2 * (it - T/2) / T)) - 1;

             if r1 <=1/4
                l = op1;
            elseif r2 >1/4&&r2<2/4
                l = op2;
            elseif r3 >2/4&&r3<3/4
                l = op3;
            else
                l = op4;
             end

               A=round(unifrnd(1,N));
            while A==i
                A=round(unifrnd(1,N));
            end
            B=round(unifrnd(1,N));
            while B==i || B==A
                B=round(unifrnd(1,N));
            end

            if i <= 10*N / 30
                     newpop(i).Position = newpop(i).Position + l * (pop(A).Position - pop(B).Position) ;         
            elseif i > 10*N / 30 && i <=  20*N / 30
                if r4 < 0.5
                      newpop(i).Position =pop(i).Position + l * ((pop(A).Position - pop(B).Position)) +  0.1*r5 * (pop(i).Position- newpop(i).Position);
                else
                     newpop(i).Position = X_best *S;
                end
            else 
                 if r6<0.5
                    newpop(i).Position=X_best-cos(2 * pi * rand(1, D)).*(X_best- pop(i).Position);
                else
                     newpop(i).Position=X_SecondBest-cos(2 * pi * rand(1, D)).*(X_SecondBest -pop(i).Position); 
                 end
            end

            
for d = 1:D
if newpop(i).Position(d) < lb(d)
    newpop(i).Position(d) = lb(d) + (lb(d) - newpop(i).Position(d)); 
    newpop(i).Position(d) = min(newpop(i).Position(d), ub(d));
elseif newpop(i).Position(d) > ub(d)
    newpop(i).Position(d) = ub(d) - (newpop(i).Position(d) - ub(d)); 
    newpop(i).Position(d) = max(newpop(i).Position(d), lb(d));
end
end
            
            newpop(i).Cost = objective(newpop(i).Position);

            if newpop(i).Cost < pop(i).Cost
                pop(i) = newpop(i);
                if pop(i).Cost < pop(1).Cost
                    X_SecondBest = X_best; 
                    X_best = pop(i).Position; 
                end
            end
        end
         
        [~, SortOrder] = sort([pop.Cost]);
        pop = pop(SortOrder);
        BestCost(it) = pop(1).Cost;
        Best_pos = pop(1).Position;
    end
    Best_score = BestCost(end); 
    curve = BestCost; 
end

function o = levyFlight(D)
    beta = 3 / 2;
    sigma = (gamma(1 + beta) * sin(pi * beta / 2) / (gamma((1 + beta) / 2) * beta * 2^((beta - 1) / 2)))^(1 / beta);
    u = randn(1, D) * sigma;
    v = randn(1, D);
    step = u ./ abs(v).^(1 / beta);
    o = step;
end   