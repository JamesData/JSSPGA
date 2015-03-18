function ChromNew=across(Chrom,XOVR,PNumber)

%  Chrom=[1 3 2 3 1 2 1 3 2 4 4 4; 
%      1 1 2 4 3 3 4 1 2 4 3 2;
%  ]; 
%    XOVR=0.7;
%  PNumber=4;

%新种群
[NIND,WNumber]=size(Chrom);
ChromNew=Chrom;
S11=zeros(1,WNumber);
S22=S11;

%随机选择交叉个体
SelNum=randperm(NIND);

SelPNumber=randperm(PNumber);

 Pos1=unidrnd(PNumber);
 Pos2=unidrnd(PNumber);
        while Pos1==Pos2      
            Pos2=unidrnd(PNumber);
        end
        
        if Pos1>Pos2
            temp=Pos1;
            Pos1=Pos2;
            Pos2=temp;
        end
  
SelPNumber=SelPNumber(Pos1:Pos2);
        
%交叉个体组个数
Num=NIND/2;
Num=2*fix(Num);
for i=1:2:Num
    if XOVR>rand; 
       
        %取交换的个体
        S1=Chrom(SelNum(i),:);
        S2=Chrom(SelNum(i+1),:);
                
        for j=1:WNumber
            if ismember(S1(j),SelPNumber)==1
               S11(j)=S1(j);
            else
                S11(j)=0;
            end
            if ismember(S2(j),SelPNumber)==1
               S22(j)=S2(j);
            else
                S22(j)=0;
            end
        end
        
        k=1;t=1;
        for j=1:WNumber            
           if ismember(S2(j),SelPNumber)==0
             while(S11(k)~=0)
               k=k+1;
             end
               S11(k)=S2(j);
           end
           
           if ismember(S1(j),SelPNumber)==0
             while(S22(t)~=0)
               t=t+1;
             end
               S22(t)=S1(j);
           end
        end

   
        
        %数据放入新群
        ChromNew(SelNum(i),:)=S11;
        ChromNew(SelNum(i+1),:)=S22;
    end
 end