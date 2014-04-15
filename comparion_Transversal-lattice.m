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
wn=w;
wr=w;
p       = 4;                % filter order
lambda  = 1.0;              % forgetting factor
laminv  = 1/lambda;
delta   = 1.0;              % initialization parameter

% Filter Initialization
% w       = zeros(p,1);       % filter coefficients
P       = delta*eye(p);     % inverse correlation matrix
% e       = x*0;              % error signal
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
figure
plot(rn(2,:),rn(1,:),'r*');
title(' Scatter Plot before Equalisation ');
xlabel(' In phase ');
ylabel(' Quadrature phase ');
%% Adaptation part with outer iteration and input as the averaged vector
for out=1:outer
    % TYPE NEW
    % multiplying with the filter coefficients in bulk
        fm1=rn(1,:).*w(1,1);
        fm2=rn(2,:).*w(1,2);
        fm3=rn(3,:).*w(1,3);
        fm4=rn(4,:).*w(1,4);
        
        nfm1=rn(1,:).*wn(1,1);
        nfm2=rn(2,:).*wn(1,2);
        nfm3=rn(3,:).*wn(1,3);
        nfm4=rn(4,:).*wn(1,4);
        
        rfm1=rn(1,:);
        rfm2=rn(2,:);
        rfm3=rn(3,:);
        rfm4=rn(4,:);
        
        % calculating the total in bulk
        fmt=fm1+fm2+fm3+fm4;
        nfmt=nfm1+nfm2+nfm3+nfm4;
        rfmt=rfm1+rfm2+rfm3+rfm4;
       % tr=sign(fmt);
        
        % actual adaptive algorithm implementation in bulk
        for i=1:length(ip)
            err(out,i)=abs(fmt(i)-ip(i));
            nerr(out,i)=abs(nfmt(i)-ip(i));
            ferrsq(i)=(err(out,i))^2;
            nferrsq(i)=(nerr(out,i))^2;
            wtest2(i,1:4,out)=w;
            wntest2(i,1:4,out)=wn;
            w=w+mu.*fmt(i).*abs(err(out,i));
            wn=wn+(mu/(nfmt(i)*(nfmt(i))')).*nfmt(i).*abs(nerr(out,i));
         end
   ep(out)=mean((err(out,:)).^2);
   nep(out)=mean((nerr(out,:)).^2);
   for m = p:length(rfmt)
            % Acquire chunk of data
               y = rfmt(m:-1:m-p+1);

            % Error signal equation
               er(out,m) = ip(m)-wr*(y');
    
            % Parameters for efficiency
               Pi = P*(y');
    
            % Filter gain vector update
                k = (Pi)/(lambda+y*Pi);

            % Inverse correlation matrix update
                P = (P - k*y*P)*laminv;

            % Filter coefficients adaption
                        wrtest2(m,1:4,out)=wr;
                wr = wr + ((er(out,m)*k)');

    end
        epr(out)=mean((er(out,:)).^2);
end
w
wn
wr
figure
plot(fm2,fm1,'r*');
title(' Scatter Plot after Equalisation: LMS ');
xlabel(' In phase ');
ylabel(' Quadrature phase ');

figure
plot(nfm2,nfm1,'r*');
title(' Scatter Plot after Equalisation : NLMS ');
xlabel(' In phase ');
ylabel(' Quadrature phase ');
figure
hold on
plot(ep,'r');
plot(nep,'b');
plot(epr,'k');
title('Mean Error comparision of algorithms');
legend('LMS','NLMS','RLS');
xlabel('No of training iterations');
ylabel('Magnitude');
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
     
     nfm1te=rnte(1,:).*wn(1,1);
     nfm2te=rnte(2,:).*wn(1,2);
     nfm3te=rnte(3,:).*wn(1,3);
     nfm4te=rnte(4,:).*wn(1,4);
     
     rfm1te=rnte(1,:);
     rfm2te=rnte(2,:);
     rfm3te=rnte(3,:);
     rfm4te=rnte(4,:);
     
     % summing the individual delay elements
     fmte=fm1te+fm2te+fm3te+fm4te;
     nfmte=nfm1te+nfm2te+nfm3te+nfm4te;
     rfmte=rfm1te+rfm2te+rfm3te+rfm4te;
     tr=sign(fmte);
     ntr=sign(nfmte);
     rtr=sign(rfmte);
     ber=0;
     nber=0;
     rber=0;
     for i=1:length(new)
%         if(fmte(i) ~= new(i))
%             ber=ber+1;
%         end
          if(tr(i) ~= new(i))
              ber=ber+1;
          end
          if(ntr(i) ~= new(i))
              nber=nber+1;
          end

     end
     for m = p:length(rfmte)
            % Acquire chunk of data
               y = rfmte(m:-1:m-p+1);
               op(m)=wr*(y');
               opr(m)=sign(op(m));
               if(new(m) ~= opr(m))
                  rber=rber+1;
               end
    end
     ber_t(ctr)=mean(ber);
     nber_t(ctr)=mean(nber);
     rber_t(ctr)=mean(rber);
     ctr=ctr+1;
 end
 xax=linspace(-10,15,26);
 figure
 hold on
 plot(xax,log10(ber_t),'r');
 plot(xax,log10(nber_t),'b');
 plot(xax,log10(rber_t),'k');
 title(' BER plot comparision ');
 xlabel(' S N R  VALUE ');
 ylabel(' BER in logarithmic scale ');
 legend('LMS','NLMS','RLS');
%  figure
%  semilogy(xax,ber_t,nber_t);
 toc
