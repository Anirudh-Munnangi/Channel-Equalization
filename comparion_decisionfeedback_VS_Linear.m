 tic
clear all
close all
N=2000;
ii=sign(randn(1,N));
wi=rand(1,4);
fei=rand(1,4);
clc
%% initialisation of program variables
mu=0.001;
ch=[0.2 0.3 1 0.3];
sn=20;
outer=50;
w=wi;%rand(1,4);
w1=wi;w2=wi;
fe=fei;%rand(1,4);
fe1=fe;fe2=fe;

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
ip=ii;%sign(randn(1,N));
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
% TYPE NEW (NAN with probability of 0.3  >> working sometimes( mu=0.001) with mean error getting constant at 1 :(  )

% multiplying with the filter coefficients in bulk
fm1(1,:)=rn(1,:).*w(1,1);
fm2(1,:)=rn(2,:).*w(1,2);
fm3(1,:)=rn(3,:).*w(1,3);
fm4(1,:)=rn(4,:).*w(1,4);

em1(1,:)=rn(1,:).*w(1,1);
em2(1,:)=rn(2,:).*w(1,2);
em3(1,:)=rn(3,:).*w(1,3);
em4(1,:)=rn(4,:).*w(1,4);

% calculating the total in bulk

fmt(1,:)=fm1(1,:)+fm2(1,:)+fm3(1,:)+fm4(1,:);

emt(1,:)=em1(1,:)+em2(1,:)+em3(1,:)+em4(1,:);
% tr=sign(fmt);

% actual adaptive algorithm implementation in bulk
for i=1:length(ip)
    err(1,i)=abs(fmt(1,i)-ip(i));
     err1(1,i)=abs(emt(1,i)-ip(i));
    wtest(i,:)=w;
    w=w+mu.*fmt(i).*abs(err(1,i));
     w1=w1+mu.*emt(i).*abs(err1(1,i));
end

dec(1,:)=sign(fmt);

dec1(1,:)=sign(fmt);
%% Adaptation part with outer iteration and input as the averaged vector
for out=2:outer
    
    % multiplying with the filter coefficients in bulk
    fm1(out,:)=rn(1,:).*w(1,1);
    fm2(out,:)=rn(2,:).*w(1,2);
    fm3(out,:)=rn(3,:).*w(1,3);
    fm4(out,:)=rn(4,:).*w(1,4);
    
    t(1,:)=(dec(out-1,:)).*fe(1,1);
    t(2,:)=(dec(out-1,:)).*fe(1,2);
    t(3,:)=(dec(out-1,:)).*fe(1,3);
    t(4,:)=(dec(out-1,:)).*fe(1,4);
    
    
    fmt(out,:)=fm1(out,:)+fm2(out,:)+fm3(out,:)+fm4(out,:)+t(1,:)+t(2,:)+t(3,:)+t(4,:);
    %     % tr=sign(fmt);
    fmc=fmt(out,:);
    
    em1(out,:)=rn(1,:).*w1(1,1);
    em2(out,:)=rn(2,:).*w1(1,2);
    em3(out,:)=rn(3,:).*w1(1,3);
    em4(out,:)=rn(4,:).*w1(1,4);
    
    t1(1,:)=(dec1(out-1,:)).*fe1(1,1);
    t1(2,:)=(dec1(out-1,:)).*fe1(1,2);
    t1(3,:)=(dec1(out-1,:)).*fe1(1,3);
    t1(4,:)=(dec1(out-1,:)).*fe1(1,4);
    
    
    emt(out,:)=em1(out,:)+em2(out,:)+em3(out,:)+em4(out,:)+t1(1,:)+t1(2,:)+t1(3,:)+t1(4,:);
    %     % tr=sign(fmt);
    emc=emt(out,:);
    % actual adaptive algorithm implementation in bulk
    for i=1:length(ip)
        %
        
        err(out,i)=abs(fmc(i)-ip(i));
        e=err(out,i);
        w=w+(mu/(fmc(i)*(fmc(i))')).*fmc(i).*abs(err(out,i));
        fe=fe+(mu/(fmc(i)*(fmc(i))')).*fmc(i).*abs(err(out,i));
        
         err1(out,i)=abs(emc(i)-ip(i));
        e1=err1(out,i);
        w1=w1+(mu).*emc(i).*abs(err1(out,i));
        fe1=fe1+(mu).*emc(i).*abs(err1(out,i));
        
    end
    dec(out,:)=sign(fmc);
     ep(out)=mse(err(out,:));
     
     dec1(out,:)=sign(emc);
     epLMS(out)=mse(err1(out,:));
     
end
    
for out=1:outer
    % TYPE NEW
    
    % multiplying with the filter coefficients in bulk
        rm1=rn(1,:);
        rm2=rn(2,:);
        rm3=rn(3,:);
        rm4=rn(4,:);
        
        % calculating the total in bulk
        rmt=rm1+rm2+rm3+rm4;
       % tr=sign(fmt);
        
        % actual adaptive algorithm implementation in bulk
        for m = p:length(rmt)
            % Acquire chunk of data
               y = rmt(m:-1:m-p+1);

            % Error signal equation
               e(out,m) = ip(m)-w2*(y');
    
            % Parameters for efficiency
               Pi = P*(y');
    
            % Filter gain vector update
                k = (Pi)/(lambda+y*Pi);

            % Inverse correlation matrix update
                P = (P - k*y*P)*laminv;

            % Filter coefficients adaption
                       
                w2= w2 + ((e(out,m)*k)');
           
        end
   epRLS(out)=mean((e(out,:)).^2);
end

%% ber plot part
N=10000;
new=sign(randn(1,N));
ctr=1;
k=0;
fe;
c=1;
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
    fmq1(c,:)=rnte(1,:).*w(1,1);
    
      em1te=rnte(1,:).*w1(1,1);
    em2te=rnte(2,:).*w1(1,2);
    em3te=rnte(3,:).*w1(1,3);
    em4te=rnte(4,:).*w1(1,4);
    emq1(c,:)=rnte(1,:).*w1(1,1);
    
    rm1te=rnte(1,:).*w2(1,1);
    rm2te=rnte(2,:).*w2(1,2);
    rm3te=rnte(3,:).*w2(1,3);
    rm4te=rnte(4,:).*w2(1,4);
    rmq1(c,:)=rnte(1,:).*w2(1,1);
    
    if(k==1)
        s1(1,:)=sign(fm1te).*fe(1,1);
        s1(2,:)=sign(fm2te).*fe(1,2);
        s1(3,:)=sign(fm3te).*fe(1,3);
        s1(4,:)=sign(fm4te).*fe(1,4);
        fmte=fm1te+fm2te+fm3te+fm4te+s1(1,:)+s1(2,:)+s1(3,:)+s1(4,:);
        
        s11(1,:)=sign(em1te).*fe1(1,1);
        s11(2,:)=sign(em2te).*fe1(1,2);
        s11(3,:)=sign(em3te).*fe1(1,3);
        s11(4,:)=sign(em4te).*fe1(1,4);
        emte=em1te+em2te+em3te+em4te+s11(1,:)+s11(2,:)+s11(3,:)+s11(4,:);
        
        s111(1,:)=sign(rm1te).*fe2(1,1);
        s111(2,:)=sign(rm2te).*fe2(1,2);
        s111(3,:)=sign(rm3te).*fe2(1,3);
        s111(4,:)=sign(rm4te).*fe2(1,4);
        rmte=rm1te+rm2te+rm3te+rm4te+s111(1,:)+s111(2,:)+s111(3,:)+s111(4,:);
    else
        fmte=fm1te+fm2te+fm3te+fm4te;
        
        emte=em1te+em2te+em3te+em4te;
        
        rmte=rm1te+rm2te+rm3te+rm4te;
        k=1;
    end
    tr=sign(fmte);
    ber=0;
    
    tr1=sign(emte);
    ber1=0;
    
    tr2=sign(rmte);
    ber2=0;
    for i=1:length(new)
        
        if(tr(i) ~= new(i))
            ber=ber+1;
        end
    end
    ber_t(ctr)=mean(ber);
        
    for i=1:length(new)
        
        if(tr1(i) ~= new(i))
            ber1=ber1+1;
        end
    end
    ber_t1(ctr)=mean(ber1);
    
    for i=1:length(new)
        
        if(tr2(i) ~= new(i))
            ber2=ber2+1;
        end
    end
    ber_t2(ctr)=mean(ber2);
    ctr=ctr+1;
    c=c+1;
end


%LINEAR MODEL
%% initialisation of program variables
mu=0.001;
ch=[0.2 0.3 1 0.3];
sn=20;
outer=50;
w=wi;%rand(1,4);
wn=wi;%w;
wr=wi;%w;
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
ip=ii;%sign(randn(1,N));
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
%             ferrsq(i)=(err(out,i))^2;
%             nferrsq(i)=(nerr(out,i))^2;
%             wtest2(i,1:4,out)=w;
%             wntest2(i,1:4,out)=wn;
            w=w+mu.*fmt(i).*abs(err(out,i));
            wn=wn+(mu/(nfmt(i)*(nfmt(i))')).*nfmt(i).*abs(nerr(out,i));
         end
   epL(out)=mean((err(out,:)).^2);
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
%           if(rtr(i) ~= new(i))
%               rber=rber+1;
%           end
          % done till here
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
     ber_tL(ctr)=mean(ber);
     nber_tL(ctr)=mean(nber);
     rber_t(ctr)=mean(rber);
     ctr=ctr+1;
 end
 

 
figure
hold on
plot(ep,'r');
plot(epL,'b');
legend('Descision Feedback NLMS', 'Linear Model NLMS');
hold off
title('Mean Square Error Versus Iterations');
xlabel('No of Iterations');
ylabel('Mean Square Error');

figure
hold on
plot(epLMS,'r');
plot(epL,'b');
legend('Descision Feedback LMS', 'Linear Model LMS');
hold off
title('Mean Square Error Versus Iterations');
xlabel('No of Iterations');
ylabel('Mean Square Error');

figure
hold on
plot(epRLS,'r');
plot(epr,'b');
legend('Descision Feedback RLS', 'Linear Model RLS');
hold off
title('Mean Square Error Versus Iterations');
xlabel('No of Iterations');
ylabel('Mean Square Error');

 xax=linspace(-10,15,26);
 figure
 hold on
 plot(xax,log10(ber_t),'r');
 plot(xax,log10(nber_tL),'b');
  title(' BER plot comparision NLMS ');
 xlabel(' S N R  VALUE ');
 ylabel(' BER in logarithmic scale ');
%  legend('Decision Feedback','Linear  Model');
 legend('Linear  Model','Decision Feedback');
 
 xax=linspace(-10,15,26);
 figure
 hold on
 plot(xax,log10(ber_t1),'r');
 plot(xax,log10(ber_tL),'b');
  title(' BER plot comparision LMS ');
 xlabel(' S N R  VALUE ');
 ylabel(' BER in logarithmic scale ');
 %legend('Decision Feedback','Linear  Model');
 legend('Linear  Model','Decision Feedback');
 
 xax=linspace(-10,15,26);
 figure
 hold on
 plot(xax,log10(ber_t2),'r');
 plot(xax,log10(rber_t),'b');
  title(' BER plot comparision RLS ');
 xlabel(' S N R  VALUE ');
 ylabel(' BER in logarithmic scale ');
% legend('Decision Feedback','Linear  Model');
 legend('Linear  Model','Decision Feedback');
% xax=linspace(-10,15,26);
% figure
% plot(xax,log10(ber_t));
% figure
% semilogy(xax,ber_t);
 toc
