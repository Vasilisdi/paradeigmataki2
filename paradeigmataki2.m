clear all;clc;close all;
%%
%vasilis dimitriou 5380
%example 
%project in FEM
%part1 - initial experiment every support is pinned
%%
%------------part of meshing-----------------------------------------------
%nodes
j=1; c0=10; f=10000;
y3_temp=linspace(-50,50,c0);scale2=1/10;scale1=1/100;
y3_temp=[y3_temp']; y3=[]; x3=[]; c=10;
for i=1:length(linspace(-85,85,c))
y3=[y3 y3_temp];
end
for i=1:length(y3_temp)
x3=[x3;linspace(-85,85,c)];
end

[c r]=size(x3);numnode=c*r; Elem3=[];
%elements
ie=1; Elem3=[];  NodeID3=[(1:length(reshape(x3,1,numnode))) ; reshape(x3',1,numnode) ; reshape(y3',1,numnode)]';
for j=1:c-1           %c-1
     for i=1:r-1
Elem3=[Elem3; i+(r)*(j-1)  i+1+(j-1)*(r)  (i)+(j)*(r)  ;...
     i+1+(j-1)*(r)  i+1+(j)*(r)  i+(j)*(r) ];  
     end
end

XYtotf=[NodeID3(Elem3(:,1),[2 3]) NodeID3(Elem3(:,2),[2 3]) NodeID3(Elem3(:,3),[2 3])];
figure(500)
for i=1:length(XYtotf)
    plot([XYtotf(i,1) XYtotf(i,3) XYtotf(i,5)],[XYtotf(i,2),XYtotf(i,4),XYtotf(i,6)],'g');hold on
end


Nodes=NodeID3;Elements=Elem3;
v=0.3; E=210000; t=4;Ktot=zeros(length(Nodes)*2);
D=E/(1-v^2)*[1 v 0;v 1 0;0 0 (1-v)/2]; As=[];
for i=1:length(Elements)
    b= [ Nodes(find(Nodes(:,1)==Elements(i,2)),3)-Nodes(find(Nodes(:,1)==Elements(i,3)),3);...
        Nodes(find(Nodes(:,1)==Elements(i,3)),3)-Nodes(find(Nodes(:,1)==Elements(i,1)),3);...
        Nodes(find(Nodes(:,1)==Elements(i,1)),3)-Nodes(find(Nodes(:,1)==Elements(i,2)),3)];
    c=[Nodes(find(Nodes(:,1)==Elements(i,3)),2)-Nodes(find(Nodes(:,1)==Elements(i,2)),2);...
        Nodes(find(Nodes(:,1)==Elements(i,1)),2)-Nodes(find(Nodes(:,1)==Elements(i,3)),2);...
        Nodes(find(Nodes(:,1)==Elements(i,2)),2)-Nodes(find(Nodes(:,1)==Elements(i,1)),2)];
   
    A=1/2*det([1 Nodes(find(Nodes(:,1)==Elements(i,1)),2)  Nodes(find(Nodes(:,1)==Elements(i,1)),3);...
        1 Nodes(find(Nodes(:,1)==Elements(i,2)),2) Nodes(find(Nodes(:,1)==Elements(i,2)),3);...
        1 Nodes(find(Nodes(:,1)==Elements(i,3)),2) Nodes(find(Nodes(:,1)==Elements(i,3)),3)]);
    B=(1/(2*A))*[b(1,1) 0 b(2,1) 0 b(3,1) 0;0 c(1,1) 0 c(2,1) 0 c(3,1);c(1,1) b(1,1) c(2,1) b(2,1) c(3,1) b(3,1)];
    ke=(B'*D*B)*t*A;
    DOF_temp=[2*find(Nodes(:,1)==Elements(i,1))-1,2*find(Nodes(:,1)==Elements(i,1)),...
        2*find(Nodes(:,1)==Elements(i,2))-1,2*find(Nodes(:,1)==Elements(i,2)),...
        2*find(Nodes(:,1)==Elements(i,3))-1,2*find(Nodes(:,1)==Elements(i,3))];
    Ktot(DOF_temp,DOF_temp)=Ktot(DOF_temp,DOF_temp)+ke;
end

%%
%--------------------Static Load Modeling----------------------------------

%determine the boundary condition - stable-movable
HoldBoundaryNodes=find(Nodes(:,2)==-85); HoldRestiDOFs=zeros(2*length(HoldBoundaryNodes),1);
for i=1:length(HoldBoundaryNodes)
HoldRestiDOFs(2*i-1)=2*HoldBoundaryNodes(i)-1;
HoldRestiDOFs(2*i)=2*HoldBoundaryNodes(i);
end
actDOFs=(1:2*length(Nodes))'; actDOFs(HoldRestiDOFs)=[];
K=Ktot; 



%determine the Loaded condition
My=10000000; %Nmm
HoldLoadedNodes=find(Nodes(:,2)==85);HoldMeanDistperElem=zeros(1,length(HoldLoadedNodes)-1);
for i=1:length(HoldLoadedNodes)-1
    HoldMeanDistperElem(i)=(Nodes(HoldLoadedNodes(i),3)+Nodes(HoldLoadedNodes(i+1),3))/2; 
end
alph=(My/2)*(1/(sum(HoldMeanDistperElem.^2))); Fbtw=zeros(length(HoldMeanDistperElem),1);Fh=zeros(length(HoldLoadedNodes),1);
Fbtw=(alph*HoldMeanDistperElem)';
LNodes=2*HoldLoadedNodes;
for i=1:length(HoldLoadedNodes)
  if i==1
      Fh(i)=Fbtw(i)/2;
  end
  if i==length(HoldLoadedNodes)
      Fh(i)=Fbtw(i-1)/2;
  end
  if i~=1 && i~=length(HoldLoadedNodes)
      Fh(i)=(Fbtw(i)/2)+(Fbtw(i-1)/2);
  end
end
F=zeros(2*length(Nodes),1);F(2*HoldLoadedNodes-1,1)=Fh;
ua=K(actDOFs,actDOFs)\F(actDOFs);
u=zeros(2*length(Nodes),1);u(actDOFs,1)=ua;   %eliminatiion method
XYtotfnew=zeros(length(Elements(:,1)),6);
Nodesst=zeros(length(Nodes),3);Nodesst(:,[2 3])=Nodes(:,[2 3])+reshape(u,2,length(Nodes))';Nodesst(:,1)=Nodes(:,1);
for i=1:length(Elements(:,1))
XYtotfnew(i,[1 2])=[Nodesst(find(Nodes(:,1)==Elements(i,1)),[2 3])];
XYtotfnew(i,[3 4])=[Nodesst(find(Nodes(:,1)==Elements(i,2)),[2 3])];
XYtotfnew(i,[5 6])=[Nodesst(find(Nodes(:,1)==Elements(i,3)),[2 3])];
end
figure(500)
for i=1:length(XYtotfnew)
    plot([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],'b');hold on
end
title('kampsh katharh-prasino prin mple meta')
for i=1:length(HoldBoundaryNodes)
plot(Nodesst(HoldBoundaryNodes(i),2),Nodesst(HoldBoundaryNodes(i),3),'ro')
end
for i=1:length(XYtotf)
    plot([XYtotf(i,1) XYtotf(i,3) XYtotf(i,5)],[XYtotf(i,2),XYtotf(i,4),XYtotf(i,6)],'g');hold on
end
eps=zeros(3,size(Elements,1)); tension=zeros(3,size(Elements,1)); u2=u;
for i=1:length(Elements)
    A=1/2*det([1 Nodes(find(Nodes(:,1)==Elements(i,1)),2)  Nodes(find(Nodes(:,1)==Elements(i,1)),3);...
        1 Nodes(find(Nodes(:,1)==Elements(i,2)),2) Nodes(find(Nodes(:,1)==Elements(i,2)),3); 1 Nodes(find(Nodes(:,1)==Elements(i,3)),2) Nodes(find(Nodes(:,1)==Elements(i,3)),3)]);
    b= [ Nodes(find(Nodes(:,1)==Elements(i,2)),3)-Nodes(find(Nodes(:,1)==Elements(i,3)),3);...
        Nodes(find(Nodes(:,1)==Elements(i,3)),3)-Nodes(find(Nodes(:,1)==Elements(i,1)),3);...
        Nodes(find(Nodes(:,1)==Elements(i,1)),3)-Nodes(find(Nodes(:,1)==Elements(i,2)),3)];
    c=[Nodes(find(Nodes(:,1)==Elements(i,3)),2)-Nodes(find(Nodes(:,1)==Elements(i,2)),2);...
        Nodes(find(Nodes(:,1)==Elements(i,1)),2)-Nodes(find(Nodes(:,1)==Elements(i,3)),2);...
        Nodes(find(Nodes(:,1)==Elements(i,2)),2)-Nodes(find(Nodes(:,1)==Elements(i,1)),2)]; 
    B=(1/(2*A))*[b(1,1) 0 b(2,1) 0 b(3,1) 0;0 c(1,1) 0 c(2,1) 0 c(3,1);c(1,1) b(1,1) c(2,1) b(2,1) c(3,1) b(3,1)];ke=(B'*D*B)*t*A;
    eps(:,i)=B*[u2(2*(find(Nodes(:,1)==Elements(i,1)))-1);...
        u2(2*(find(Nodes(:,1)==Elements(i,1))));...
        u2(2*(find(Nodes(:,1)==Elements(i,2)))-1);...
        u2(2*(find(Nodes(:,1)==Elements(i,2))));...
        u2(2*(find(Nodes(:,1)==Elements(i,3)))-1);...
        u2(2*(find(Nodes(:,1)==Elements(i,3))))]; tension(:,i)=D*eps(:,i);   %per element
end
sigx=zeros(length(Nodes),1);sigy=zeros(length(Nodes),1);dutot=zeros(length(Nodes),1);tau=zeros(length(Nodes),1);
ex=zeros(length(Nodes),1);ey=zeros(length(Nodes),1);gama=zeros(length(Nodes),1);vonmises=zeros(length(Nodes),1);
u1=reshape(u,2,length(Nodes))';ihold=zeros(1,length(Nodes(:,1)));
for i=1:length(Nodes)
    ihold(i)=0;
  for j=1:size(Elements,1)
   if (Elements(j,1)==Nodes(i,1) || Elements(j,2)==Nodes(i,1) || Elements(j,3)==Nodes(i,1))
       sigx(i)=sigx(i)+tension(1,j);sigy(i)=sigy(i)+tension(2,j);tau(i)=tau(i)+tension(3,j);
       ex(i)=ex(i)+eps(1,j);ey(i)=ey(i)+eps(2,j);gama(i)=gama(i)+eps(3,j);
       ihold(i)=ihold(i)+1;      %#of elems. sourrounding one node 
       dutot(i)=sqrt((u1(i,1)^2+u1(i,2)^2));
       vonmises(i)=sqrt(sigx(i)^2-(sigx(i)*sigy(i))+(sigy(i)^2)+(3*(gama(i)^2)));
   end
  end
   sigx(i)=sigx(i)/ihold(i);sigy(i)=sigy(i)/ihold(i);tau(i)=tau(i)/ihold(i);  %tension mean val
   ex(i)=ex(i)/ihold(i);ey(i)=ey(i)/ihold(i);gama(i)=gama(i)/ihold(i);        %deformation mean val
   vonmises(i)=vonmises(i)/ihold(i);
end
sigma1=zeros(length(Nodes),1);sigma2=zeros(length(Nodes),1);
for i=1:size(Nodes,1)
sigma1(i,1)=((sigx(i)+sigy(i))/2)+sqrt(((sigx(i)-sigy(i))/2)^2+tau(i)^2);
sigma2(i,1)=((sigx(i)+sigy(i))/2)-sqrt(((sigx(i)-sigy(i))/2)^2+tau(i)^2);
end
colormap;Sigx=zeros(length(Elements(:,1)),3);
figure(2);title('normal stress {\sigma}_x[N/{mm}^2]');
for i=1:size(Elements,1)
Sigx(i,:)=[sigx(find(Nodes(:,1)==Elements(i,1))) sigx(find(Nodes(:,1)==Elements(i,2)))...
     sigx(find(Nodes(:,1)==Elements(i,3)))];
end
for i=1:size(Elements,1)
  patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Sigx(i,:)');
  hold on
end
axis equal;colorbar;
colormap;Sigy=zeros(length(Elements(:,1)),3);
figure(3);title('normal stress {\sigma}_y[N/{mm}^2]');
for i=1:size(Elements,1)
Sigy(i,:)=[sigy(find(Nodes(:,1)==Elements(i,1))) sigy(find(Nodes(:,1)==Elements(i,2)))...
     sigy(find(Nodes(:,1)==Elements(i,3))) ];
end
for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Sigy(i,:)');
   hold on
end
axis equal;colorbar;
colormap;Tau=zeros(length(Elements(:,1)),3);
figure(4);title('shear stress {\tau}_{xy}[N/{mm}^2]');
for i=1:size(Elements,1)
Tau(i,:)=[tau(find(Nodes(:,1)==Elements(i,1))) tau(find(Nodes(:,1)==Elements(i,2)))...
    tau(find(Nodes(:,1)==Elements(i,3))) ];
end
 for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Tau(i,:)');
   hold on
end
axis equal;colorbar;
colormap;Du=zeros(length(Elements(:,1)),3); 
 figure(5);title('total du[mm]')
for i=1:size(Elements,1)
Du(i,:)=[dutot(find(Nodes(:,1)==Elements(i,1))) dutot(find(Nodes(:,1)==Elements(i,2)))...
     dutot(find(Nodes(:,1)==Elements(i,3))) ];
end
 for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Du(i,:)');
   hold on
 end
axis equal;colorbar;
 colormap;Ex=zeros(length(Elements(:,1)),3);
figure(6);title('strain displacement ex'); 
for i=1:size(Elements,1)
Ex(i,:)=[ex(find(Nodes(:,1)==Elements(i,1))) ex(find(Nodes(:,1)==Elements(i,2)))...
     ex(find(Nodes(:,1)==Elements(i,3))) ];
end
for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Ex(i,:)');
   hold on
end
axis equal;colorbar;
colormap;Ey=zeros(length(Elements(:,1)),3);
figure(7);title('strain displacement ey');
for i=1:size(Elements,1)
Ey(i,:)=[ey(find(Nodes(:,1)==Elements(i,1))) ey(find(Nodes(:,1)==Elements(i,2)))...
     ey(find(Nodes(:,1)==Elements(i,3)))];
end
for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Ey(i,:)');
   hold on
end
axis equal;colorbar;
colormap;Gam=zeros(length(Elements(:,1)),3);
figure(8);title('strain displacement {\gamma}_{xy}');
for i=1:size(Elements,1)
Gam(i,:)=[gama(find(Nodes(:,1)==Elements(i,1))) gama(find(Nodes(:,1)==Elements(i,2)))...
    gama(find(Nodes(:,1)==Elements(i,3)))];
end
for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Gam(i,:)');
   hold on
end
axis equal;colorbar;
colormap;Sig1=zeros(length(Elements(:,1)),3);
figure(9);title(' {\sigma}_{1}[N/mm^2]');
for i=1:size(Elements,1)
Sig1(i,:)=[sigma1(find(Nodes(:,1)==Elements(i,1))) sigma1(find(Nodes(:,1)==Elements(i,2)))...
     sigma1(find(Nodes(:,1)==Elements(i,3))) ];
end
for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Sig1(i,:)');
   hold on
end
axis equal;colorbar;
colormap;Sig2=zeros(length(Elements(:,1)),3);
figure(10);title(' {\sigma}_{2}[N/mm^2]');
for i=1:size(Elements,1)
Sig2(i,:)=[sigma2(find(Nodes(:,1)==Elements(i,1))) sigma2(find(Nodes(:,1)==Elements(i,2)))...
     sigma2(find(Nodes(:,1)==Elements(i,3))) ];
end
for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Sig2(i,:)');
   hold on
end
axis equal;colorbar;
colormap;VM=zeros(length(Elements(:,1)),3);
figure(25);title('vonmises[N/mm^2]');
for i=1:size(Elements,1)
VM(i,:)=[vonmises(find(Nodes(:,1)==Elements(i,1))) vonmises(find(Nodes(:,1)==Elements(i,2)))...
     vonmises(find(Nodes(:,1)==Elements(i,3)))];
end
for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],VM(i,:)');
   hold on
end
axis equal;colorbar;
figure(22);
for i=1:length(Nodes(:,1))
    dx=sigma1(i).*cos((1/2)*atan((2*tau(i))./(sigx(i)-sigy(i))));
    dy=sigma1(i).*sin((1/2)*atan((2*tau(i))./(sigx(i)-sigy(i))));
   quiver(Nodes(i,2),Nodes(i,3),dx*scale1,dy*scale1,'b')
   hold on;
end
title(['print arrows for {\sigma}_{1} scale for arrows:' num2str(scale1)])
figure(23);
for i=1:length(Nodes(:,1))
    dx=sigma2(i).*cos(pi/2-((1/2)*atan((2*tau(i))./(sigx(i)-sigy(i)))));
    dy=sigma2(i).*sin(pi/2-((1/2)*atan((2*tau(i))./(sigx(i)-sigy(i)))));
   quiver(Nodes(i,2),Nodes(i,3),dx*scale1,dy*scale1,'b')
   hold on;
end
title(['print arrows for {\sigma}_{2} scale for arrows:' num2str(scale1)])
%%
%in txt file
[Ysx,Isx]=max(sigx);[Ysy,Isy]=max(sigy);[Yt,It]=max(tau);[Ys1,Is1]=max(sigma1);[Yvm,Ivm]=max(vonmises);
[Ydu,Idu]=max(dutot);[Yex,Iex]=max(ex);[Yey,Iey]=max(ey);[Yg,Ig]=max(gama);
fpr=[Ysx,Nodes(Isx,:);Ysy,Nodes(Isy,:);Yt,Nodes(It,:);Ys1,Nodes(Is1,:);Yvm,Nodes(Ivm,:);...
    Ydu,Nodes(Idu,:);Yex,Nodes(Iex,:);Yey,Nodes(Iey,:);Yg,Nodes(Ig,:)];


nod=fopen('AEM5380_2nd_A_ResultsBending_MeshingCir.txt','w+t');
fprintf(nod,  'Maximum_Value #Code_of_NODE      X[mm]    Y[mm]  \n');fprintf('\n');
fprintf(nod,'  sigma_{x}[N/mm2] \n');fprintf('\n');
fprintf(nod,'  %f      %.0f        %f     %f  \n',[ fpr(1,:) ] );fprintf('\n');
fprintf(nod,'  sigma_{y}[N/mm2] \n');fprintf('\n');
fprintf(nod,'  %f      %.0f        %f     %f  \n',[ fpr(2,:) ] );fprintf('\n');
fprintf(nod,'  tau_{xy}[N/mm2] \n');fprintf('\n');
fprintf(nod,'  %f      %.0f        %f     %f  \n',[ fpr(3,:) ] );fprintf('\n');
fprintf(nod,'  principal stress sigma_{1}[N/mm2] \n');fprintf('\n');
fprintf(nod,'  %f      %.0f        %f     %f  \n',[ fpr(4,:) ] );fprintf('\n');
fprintf(nod,'  Max_Von_Mises_stress[N/mm2] \n');fprintf('\n');
fprintf(nod,'  %f      %.0f        %f     %f  \n',[ fpr(5,:) ] );fprintf('\n');
fprintf(nod,'  du \n');fprintf('\n');
fprintf(nod,'  %f         %.0f        %f     %f  \n',[ fpr(6,:) ] );fprintf('\n');
fprintf(nod,'  epsilon_{x} \n');fprintf('\n');
fprintf(nod,'  %f         %.0f        %f     %f  \n',[ fpr(7,:) ] );fprintf('\n');
fprintf(nod,'  epsilon_{y} \n');fprintf('\n');
fprintf(nod,'  %f         %.0f        %f     %f  \n',[ fpr(8,:) ] );fprintf('\n');
fprintf(nod,'  gamma_{xy} \n');fprintf('\n');
fprintf(nod,'  %f         %.0f        %f     %f  \n',[ fpr(9,:) ] );fprintf('\n');
alpha_kappa=Ysx/mean(abs(sigx(find(Nodes(:,2)==Nodes(Isx,2)))));
fprintf(nod,'  alpha_{kappa} \n');fprintf('\n');
fprintf(nod,'      %f  \n',[ alpha_kappa ] );fprintf('\n');
fclose(nod);

%figure(1)
%for i=1:length(Elements)
%  patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],tension(1,i));
%  hold on
%end
%%



F=zeros(2*length(Nodes),1);F(2*HoldLoadedNodes-1,1)=f;
F(2*HoldLoadedNodes(1)-1,1)=f/2;F(2*HoldLoadedNodes(length(HoldLoadedNodes))-1,1)=f/2;
%u=K\F;                                       %penalty method

ua=K(actDOFs,actDOFs)\F(actDOFs);
u=zeros(2*length(Nodes),1);u(actDOFs,1)=ua;   %eliminatiion method
XYtotfnew=zeros(length(Elements(:,1)),6);
Nodesst=zeros(length(Nodes),3);Nodesst(:,[2 3])=Nodes(:,[2 3])+reshape(u,2,length(Nodes))';Nodesst(:,1)=Nodes(:,1);
for i=1:length(Elements(:,1))
XYtotfnew(i,[1 2])=[Nodesst(find(Nodes(:,1)==Elements(i,1)),[2 3])];
XYtotfnew(i,[3 4])=[Nodesst(find(Nodes(:,1)==Elements(i,2)),[2 3])];
XYtotfnew(i,[5 6])=[Nodesst(find(Nodes(:,1)==Elements(i,3)),[2 3])];
end
figure(600)
for i=1:length(XYtotfnew)
    plot([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],'b');hold on
end
title('efelkusmos katharos-prasino prin mple meta')
for i=1:length(HoldBoundaryNodes)
plot(Nodesst(HoldBoundaryNodes(i),2),Nodesst(HoldBoundaryNodes(i),3),'ro')
end
for i=1:length(XYtotf)
    plot([XYtotf(i,1) XYtotf(i,3) XYtotf(i,5)],[XYtotf(i,2),XYtotf(i,4),XYtotf(i,6)],'g');hold on
end

eps=zeros(3,size(Elements,1)); tension=zeros(3,size(Elements,1)); u2=u;
for i=1:length(Elements)
    A=1/2*det([1 Nodes(find(Nodes(:,1)==Elements(i,1)),2)  Nodes(find(Nodes(:,1)==Elements(i,1)),3);...
        1 Nodes(find(Nodes(:,1)==Elements(i,2)),2) Nodes(find(Nodes(:,1)==Elements(i,2)),3); 1 Nodes(find(Nodes(:,1)==Elements(i,3)),2) Nodes(find(Nodes(:,1)==Elements(i,3)),3)]);
    b= [ Nodes(find(Nodes(:,1)==Elements(i,2)),3)-Nodes(find(Nodes(:,1)==Elements(i,3)),3);...
        Nodes(find(Nodes(:,1)==Elements(i,3)),3)-Nodes(find(Nodes(:,1)==Elements(i,1)),3);...
        Nodes(find(Nodes(:,1)==Elements(i,1)),3)-Nodes(find(Nodes(:,1)==Elements(i,2)),3)];
    c=[Nodes(find(Nodes(:,1)==Elements(i,3)),2)-Nodes(find(Nodes(:,1)==Elements(i,2)),2);...
        Nodes(find(Nodes(:,1)==Elements(i,1)),2)-Nodes(find(Nodes(:,1)==Elements(i,3)),2);...
        Nodes(find(Nodes(:,1)==Elements(i,2)),2)-Nodes(find(Nodes(:,1)==Elements(i,1)),2)]; 
    B=(1/(2*A))*[b(1,1) 0 b(2,1) 0 b(3,1) 0;0 c(1,1) 0 c(2,1) 0 c(3,1);c(1,1) b(1,1) c(2,1) b(2,1) c(3,1) b(3,1)];ke=(B'*D*B)*t*A;
   eps(:,i)=B*[u2(2*(find(Nodes(:,1)==Elements(i,1)))-1);...
        u2(2*(find(Nodes(:,1)==Elements(i,1))));...
        u2(2*(find(Nodes(:,1)==Elements(i,2)))-1);...
        u2(2*(find(Nodes(:,1)==Elements(i,2))));...
        u2(2*(find(Nodes(:,1)==Elements(i,3)))-1);...
        u2(2*(find(Nodes(:,1)==Elements(i,3))))]; tension(:,i)=D*eps(:,i);   %per element
end
sigx=zeros(length(Nodes),1);sigy=zeros(length(Nodes),1);dutot=zeros(length(Nodes),1);tau=zeros(length(Nodes),1);
ex=zeros(length(Nodes),1);ey=zeros(length(Nodes),1);gama=zeros(length(Nodes),1);vonmises=zeros(length(Nodes),1);
u1=reshape(u,2,length(Nodes))';ihold=zeros(1,length(Nodes(:,1)));
for i=1:length(Nodes)
    ihold(i)=0;
  for j=1:size(Elements,1)
   if (Elements(j,1)==Nodes(i,1) || Elements(j,2)==Nodes(i,1) || Elements(j,3)==Nodes(i,1))
       sigx(i)=sigx(i)+tension(1,j);sigy(i)=sigy(i)+tension(2,j);tau(i)=tau(i)+tension(3,j);
       ex(i)=ex(i)+eps(1,j);ey(i)=ey(i)+eps(2,j);gama(i)=gama(i)+eps(3,j);
       ihold(i)=ihold(i)+1;      %#of elems. sourrounding one node 
       dutot(i)=sqrt((u1(i,1)^2+u1(i,2)^2));
       vonmises(i)=sqrt(sigx(i)^2-(sigx(i)*sigy(i))+(sigy(i)^2)+(3*(gama(i)^2)));
   end
  end
   sigx(i)=sigx(i)/ihold(i);sigy(i)=sigy(i)/ihold(i);tau(i)=tau(i)/ihold(i);  %tension mean val
   ex(i)=ex(i)/ihold(i);ey(i)=ey(i)/ihold(i);gama(i)=gama(i)/ihold(i);        %deformation mean val
   vonmises(i)=vonmises(i)/ihold(i);
end
sigma1=zeros(length(Nodes),1);sigma2=zeros(length(Nodes),1);
for i=1:size(Nodes,1)
sigma1(i,1)=((sigx(i)+sigy(i))/2)+sqrt(((sigx(i)-sigy(i))/2)^2+tau(i)^2);
sigma2(i,1)=((sigx(i)+sigy(i))/2)-sqrt(((sigx(i)-sigy(i))/2)^2+tau(i)^2);
end
colormap;Sigx=zeros(length(Elements(:,1)),3);
figure(11);title('normal stress {\sigma}_x[N/{mm}^2]');
for i=1:size(Elements,1)
Sigx(i,:)=[sigx(find(Nodes(:,1)==Elements(i,1))) sigx(find(Nodes(:,1)==Elements(i,2)))...
     sigx(find(Nodes(:,1)==Elements(i,3))) ];
end
for i=1:size(Elements,1)
  patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Sigx(i,:)');
  hold on
end
axis equal;colorbar;
colormap;Sigy=zeros(length(Elements(:,1)),3);
figure(12);title('normal stress {\sigma}_y[N/{mm}^2]');
for i=1:size(Elements,1)
Sigy(i,:)=[sigy(find(Nodes(:,1)==Elements(i,1))) sigy(find(Nodes(:,1)==Elements(i,2)))...
    sigy(find(Nodes(:,1)==Elements(i,3))) ];
end
for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Sigy(i,:)');
   hold on
end
axis equal;colorbar;
colormap;Tau=zeros(length(Elements(:,1)),3);
figure(13);title('shear stress {\tau}_{xy}[N/{mm}^2]');
for i=1:size(Elements,1)
Tau(i,:)=[tau(find(Nodes(:,1)==Elements(i,1))) tau(find(Nodes(:,1)==Elements(i,2)))...
     tau(find(Nodes(:,1)==Elements(i,3))) ];
end
 for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Tau(i,:)');
   hold on
end
axis equal;colorbar;
colormap;Du=zeros(length(Elements(:,1)),3); 
figure(14);title('total du[mm]')
for i=1:size(Elements,1)
Du(i,:)=[dutot(find(Nodes(:,1)==Elements(i,1))) dutot(find(Nodes(:,1)==Elements(i,2)))...
    dutot(find(Nodes(:,1)==Elements(i,3))) ];
end
 for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Du(i,:)');
   hold on
 end
axis equal;colorbar;
 colormap;Ex=zeros(length(Elements(:,1)),3);
figure(15);title('strain displacement ex'); 
for i=1:size(Elements,1)
Ex(i,:)=[ex(find(Nodes(:,1)==Elements(i,1))) ex(find(Nodes(:,1)==Elements(i,2)))...
     ex(find(Nodes(:,1)==Elements(i,3))) ];
end
for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Ex(i,:)');
   hold on
end
axis equal;colorbar;
colormap;Ey=zeros(length(Elements(:,1)),3);
figure(16);title('strain displacement ey');
for i=1:size(Elements,1)
Ey(i,:)=[ey(find(Nodes(:,1)==Elements(i,1))) ey(find(Nodes(:,1)==Elements(i,2)))...
     ey(find(Nodes(:,1)==Elements(i,3))) ];
end
for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Ey(i,:)');
   hold on
end
axis equal;colorbar;
colormap;Gam=zeros(length(Elements(:,1)),3);
figure(17);title('strain displacement {\gamma}_{xy}');
for i=1:size(Elements,1)
Gam(i,:)=[gama(find(Nodes(:,1)==Elements(i,1))) gama(find(Nodes(:,1)==Elements(i,2)))...
     gama(find(Nodes(:,1)==Elements(i,3)))];
end
for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Gam(i,:)');
   hold on
end
axis equal;colorbar;
colormap;Sig1=zeros(length(Elements(:,1)),3);
figure(18);title(' {\sigma}_{1}[N/mm^2]');
for i=1:size(Elements,1)
Sig1(i,:)=[sigma1(find(Nodes(:,1)==Elements(i,1))) sigma1(find(Nodes(:,1)==Elements(i,2)))...
     sigma1(find(Nodes(:,1)==Elements(i,3))) ];
end
for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Sig1(i,:)');
   hold on
end
axis equal;colorbar;
colormap;Sig2=zeros(length(Elements(:,1)),3);
figure(19);title(' {\sigma}_{2}[N/mm^2]');
for i=1:size(Elements,1)
Sig2(i,:)=[sigma2(find(Nodes(:,1)==Elements(i,1))) sigma2(find(Nodes(:,1)==Elements(i,2)))...
     sigma2(find(Nodes(:,1)==Elements(i,3))) ];
end
for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],Sig2(i,:)');
   hold on
end
axis equal;colorbar;
colormap;VM=zeros(length(Elements(:,1)),3);
figure(24);title('vonmises[N/mm^2]');
for i=1:size(Elements,1)
VM(i,:)=[vonmises(find(Nodes(:,1)==Elements(i,1))) vonmises(find(Nodes(:,1)==Elements(i,2)))...
    vonmises(find(Nodes(:,1)==Elements(i,3)))];
end
for i=1:size(Elements,1)
   patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],VM(i,:)');
   hold on
end
axis equal;colorbar;
figure(20);
for i=1:length(Nodes(:,1))
    dx=sigma1(i).*cos((1/2)*atan((2*tau(i))./(sigx(i)-sigy(i))));
    dy=sigma1(i).*sin((1/2)*atan((2*tau(i))./(sigx(i)-sigy(i))));
   quiver(Nodes(i,2),Nodes(i,3),dx*scale2,dy*scale2,'b')
   hold on;
end
title(['print arrows for {\sigma}_{1} scale for arrows:' num2str(scale2)])
figure(21);
for i=1:length(Nodes(:,1))
    dx=sigma2(i).*cos(pi/2-((1/2)*atan((2*tau(i))./(sigx(i)-sigy(i)))));
    dy=sigma2(i).*sin(pi/2-((1/2)*atan((2*tau(i))./(sigx(i)-sigy(i)))));
   quiver(Nodes(i,2),Nodes(i,3),dx*scale2,dy*scale2,'b')
   hold on;
end
title(['print arrows for {\sigma}_{2} scale for arrows:' num2str(scale2)])
%%
%in txt file
[Ysx,Isx]=max(sigx);[Ysy,Isy]=max(sigy);[Yt,It]=max(tau);[Ys1,Is1]=max(sigma1);[Yvm,Ivm]=max(vonmises);
[Ydu,Idu]=max(dutot);[Yex,Iex]=max(ex);[Yey,Iey]=max(ey);[Yg,Ig]=max(gama);
fpr=[Ysx,Nodes(Isx,:);Ysy,Nodes(Isy,:);Yt,Nodes(It,:);Ys1,Nodes(Is1,:);Yvm,Nodes(Ivm,:);...
    Ydu,Nodes(Idu,:);Yex,Nodes(Iex,:);Yey,Nodes(Iey,:);Yg,Nodes(Ig,:)];


nod=fopen('AEM5380_2nd_A_ResultsTensile_MeshingCir.txt','w+t');
fprintf(nod,  'Maximum_Value #Code_of_NODE      X[mm]    Y[mm]  \n');fprintf('\n');
fprintf(nod,'  sigma_{x}[N/mm2] \n');fprintf('\n');
fprintf(nod,'  %f         %.0f        %f     %f  \n',[ fpr(1,:) ] );fprintf('\n');
fprintf(nod,'  sigma_{y}[N/mm2] \n');fprintf('\n');
fprintf(nod,'  %f         %.0f        %f     %f  \n',[ fpr(2,:) ] );fprintf('\n');
fprintf(nod,'  tau_{xy}[N/mm2] \n');fprintf('\n');
fprintf(nod,'  %f         %.0f        %f     %f  \n',[ fpr(3,:) ] );fprintf('\n');
fprintf(nod,'  principal stress sigma_{1}[N/mm2] \n');fprintf('\n');
fprintf(nod,'  %f         %.0f        %f     %f  \n',[ fpr(4,:) ] );fprintf('\n');
fprintf(nod,'  Max_Von_Mises_stress[N/mm2] \n');fprintf('\n');
fprintf(nod,'  %f         %.0f        %f     %f  \n',[ fpr(5,:) ] );fprintf('\n');
fprintf(nod,'  du \n');fprintf('\n');
fprintf(nod,'  %f         %.0f        %f     %f  \n',[ fpr(6,:) ] );fprintf('\n');
fprintf(nod,'  epsilon_{x} \n');fprintf('\n');
fprintf(nod,'  %f         %.0f        %f     %f  \n',[ fpr(7,:) ] );fprintf('\n');
fprintf(nod,'  epsilon_{y} \n');fprintf('\n');
fprintf(nod,'  %f         %.0f        %f     %f  \n',[ fpr(8,:) ] );fprintf('\n');
fprintf(nod,'  gamma_{xy} \n');fprintf('\n');
fprintf(nod,'  %f         %.0f        %f     %f  \n',[ fpr(9,:) ] );fprintf('\n');
alpha_kappa=Ysx/mean(abs(sigx(find(Nodes(:,2)==Nodes(Isx,2)))));
fprintf(nod,'  alpha_{kappa} \n');fprintf('\n');
fprintf(nod,'      %f  \n',[ alpha_kappa ] );fprintf('\n');
fclose(nod);clc;


%figure(1000)
%for i=1:length(Elements)
%  patch([XYtotfnew(i,1) XYtotfnew(i,3) XYtotfnew(i,5)],[XYtotfnew(i,2),XYtotfnew(i,4),XYtotfnew(i,6)],tension(1,i));
%  hold on
%end