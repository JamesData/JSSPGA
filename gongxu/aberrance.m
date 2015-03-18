function ChromNew=aberrance(Chrom,MUTR)

% Chrom=[1 3 2 3 1 2 1 3 2; 
%     1 1 2 3 3 1 2 3 2;
%     1 3 2 3 2 2 1 3 1;
%     1 3 3 3 1 2 1 2 2;
%     2 3 2 3 1 2 1 3 1;
%     1 3 2 1 1 2 3 3 2 ;
%     3 3 2 1 2 2 1 3 1;
%     1 3 2 2 1 2 1 3 3;
%     1 3 3 2 2 2 1 3 1 ;
%     1 3 2 1 3 2 1 3 2 ]; 
%  
%   MUTR=0.2;
%  

%初始化
[NIND,WNumber]=size(Chrom);
ChromNew=Chrom;

for i=1:NIND    
    %是否变异
    if MUTR>rand;
        
        %变异位置
        Pos1=unidrnd(WNumber);
        Pos2=unidrnd(WNumber);
        
        %变异位置不相同
        while Pos1==Pos2      
            Pos2=unidrnd(WNumber);
        end 
        
        %取一个个体
        S=Chrom(i,:); 
        
        %交换
        temp=S(Pos1);
        S(Pos1)=S(Pos2);
        S(Pos2)=temp; 
        
        %数据放入新群
        ChromNew(i,:)=S;
    end
end

 