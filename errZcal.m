mydata=xlsread('Static Test.xlsx');
z=zeros(4,length(mydata)-1);
for i =1:length(mydata)-1
   %% Assuming using the Axis angle
   z(:,i)=[mydata(i,5),mydata(i,6),mydata(i,7),sqrt(mydata(i,5)^2+mydata(i,6)^2+mydata(i,7)^2)]';
end
err=var(z(4,:));