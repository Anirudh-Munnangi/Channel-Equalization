tic
clear all
close all
clc
%% initialisation of program variables
mu=0.001;
ch=[0.2 0.3 1 0.3];
sn=20;
outer=50;
w=rand(1,4);
%% initialisation of random input
N=2000;
ip=sign(randn(1,N));
% created 1 X N input random sequence
% created decisions
%% multiplication with channel coefficients
r(1,:)=ip.*ch(1,1);
r(2,:)=ip.*ch(1,2);
r(3,:)=ip.*ch(1,3);
r(4,:)=ip.*ch(1,4);
% created 4 X N input samples
%% noise addition
rn=awgn(r,sn);
%% Adaptation part 1 first iteration
% multiplying with the filter coefficients in bulk
fm1=rn(1,:).*w(1,1);
fm2=rn(2,:).*w(1,2);
fm3=rn(3,:).*w(1,3);
fm4=rn(4,:).*w(1,4);
% calculating the total in bulk
fmt=fm1+fm2+fm3+fm4;        
% actual adaptive algorithm implementation in bulk
for i=1:length(ip)
    ferr(i)=abs(fmt(i)-ip(i));
    ferrsq(i)=(ferr(i))^2;
    wtest(i,:)=w;
    %w=w+mu.*fmt(i).*abs(ferr(i));
    w=w+(mu/(fmt(i)*(fmt(i))')).*fmt(i).*abs(ferr(i));
end

%% Adaptation part with outer iteration and input as the averaged vector
for out=1:outer
    % multiplying with the filter coefficients in bulk
    fm1=rn(1,:).*w(1,1);
    fm2=rn(2,:).*w(1,2);
    fm3=rn(3,:).*w(1,3);
    fm4=rn(4,:).*w(1,4);
    % calculating the total in bulk
    fmt=fm1+fm2+fm3+fm4;
    % actual adaptive algorithm implementation in bulk
    for i=1:length(ip)
        err(out,i)=abs(fmt(i)-ip(i));
        ferrsq(i)=(err(out,i))^2;
        wtest2(i,1:4,out)=w;
        w=w+(mu/(fmt(i)*(fmt(i))')).*fmt(i).*abs(err(out,i));
    end
    ep(out)=mean((err(out,:)).^2);
end
w
figure
plot(ep);
title(' Mean Square Error Curve ');
xlabel('No of Iterations ');
ylabel('Error Magnitude');
%% ber plot part
N=10000;
new=sign(randn(1,N));
ctr=1;
 for sn=-10:1:15
     % multiplying with the channel gains
     rte(1,:)=new.*ch(1,1);
     rte(2,:)=new.*ch(1,2);  
     rte(3,:)=new.*ch(1,3);
     rte(4,:)=new.*ch(1,4);
     % addition of noise
     rnte=awgn(rte,sn);
     % multiplying with the filter gains
     fm1te=rnte(1,:).*w(1,1);
     fm2te=rnte(2,:).*w(1,2);
     fm3te=rnte(3,:).*w(1,3);
     fm4te=rnte(4,:).*w(1,4);
     % summing the individual delay elements
     fmte=fm1te+fm2te+fm3te+fm4te;
     tr=sign(fmte);
     ber=0;
     for i=1:length(new)
%         if(fmte(i) ~= new(i))
%             ber=ber+1;
%         end
          if(tr(i) ~= new(i))
              ber=ber+1;
          end
     end
     ber_t(ctr)=mean(ber);
     ctr=ctr+1;
 end
 xax=linspace(-10,15,26);
 figure
 plot(xax,log10(ber_t));
  title('Ber Vs Snr Plot');
 xlabel(' SNR value ');
 ylabel(' Ber Value ');
 figure
 semilogy(xax,ber_t);
 toc
