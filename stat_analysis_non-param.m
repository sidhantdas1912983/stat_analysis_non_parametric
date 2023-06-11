[p,tbl] = friedman(ans(:,2:end),1)
%%
col25 = PS(:,1);
col50 = PS(:,2);
col75 = PS(:,3);
col100 = PS(:,4);
col125 = PS(:,5);
%%
%Only if u need to find pvalues for a pairwise experiment using wilcoxon signed rank test
[p1,v1] = signrank(col25,col50);
[p2,v2] = signrank(col25,col75);
[p3,v3] = signrank(col25,col100);
[p4,v4] = signrank(col25,col125);
[p5,v5] = signrank(col50,col75);
[p6,v6] = signrank(col50,col100);
[p7,v7] = signrank(col50,col125);
[p8,v8] = signrank(col75,col100);
[p9,v9] = signrank(col75,col125);
[p10,v10] = signrank(col100,col125);

%%
p = [p1; p2; p3; p4;p5;p6;p7;p8;p9;p10 ];%can be modified as per your requirement
q = sort(p);
old_result =zeros(length(p),1);

for i = 1:length(p)
    if p(i)>0.05
        old_result(i) = 0; 
    else
        old_result(i) = 1;
    end
end

%%
fp = 0.05;
q = sort(p);
if length(unique(q)) == length(q)
% Finds adjusted pvalue and test hypothesis when unique p values are present    
    x = zeros(length(q),1);
    h = zeros(length(q),1);
    for i = 1:length(q)
        x(i) = (fp*(i/length(p)));
        if x(i)>q(i)
            h(i) = 1;%%reject null hypothesis : statistically different
        else
            h(i) =0;%%%%accept null hypothesis : statistically similar
        end
    end
    for i=1:length(q)
        ap(i) = (q(i)*length(q))/Q(i);
    end
else
% Finds adjusted pvalue and test hypothesis when identical p values are present
    Q = benjamin_hochberg_equalp(q);
    x = zeros(length(q),1);
    ap = zeros(length(q),1);
    h = zeros(length(q),1);
    for i = 1:length(q)
        x(i) = (fp*(Q(i)/length(p)));
        if x(i)>q(i)
            h(i) = 1;%%reject null hypothesis : statistically different
        else
            h(i) =0;%%%%accept null hypothesis : statistically similar
        end
    end
    for i=1:length(q)
        ap(i) = (q(i)*length(q))/Q(i);
    end
    
end
us = [];
for i = 1:length(p)
    us(i)=find(q==p(i),1);
end
h = h(us);
x = x(us);
ap = ap(us);
a = ["25-100";"25-125";"25-75";"50-125";"50-100";"25-50";"50-75";"75-120";"75-100";"100-125"];
Summary = ["pairwise p-vals", "old p-values", "new critical vals","adjusted p vals","result";a p x ap h];
    

%just unsorting the sorted array

 
function [rankwise_matrix] = benjamin_hochberg_equalp(s)
    [v,w,rak]= unique(s);
        turtoise =1;
        rankwise_matrix=zeros(length(rak),1);
        hare = turtoise+1;
        for i = 1:length(rak)-1
            if(rak(turtoise)==rak(hare)&&turtoise == length(rak)-1)
                rankwise_matrix(turtoise) = hare;
                rankwise_matrix(hare) = hare;
            elseif(rak(turtoise)==rak(hare)&&turtoise ~= length(rak)-1)
                hare = hare+1;    
            else
                rankwise_matrix(turtoise) = turtoise;
                dif = hare-turtoise-1;
                if(dif~=0)
                    rankwise_matrix(turtoise:hare-1,1)=hare-1;    
                end
                   turtoise = turtoise+dif+1;
                   rankwise_matrix(turtoise) = hare;
                   hare = turtoise+1; 

            end
        end
end        
