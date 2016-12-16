function Lab6_7()




        %%%% MOHIT KUMAR AHUJA %%%%
        %%%%     MSCV_2016     %%%%



% Exercise 1: 
         
% Question 1.1  

% ButterWorth Filter (Low-pass butterworth filter)

[x,y] = butter(3,0.5,'low');
[H, W]= freqz(x,y);
figure();
subplot(2,2,1);
plot(W/pi, abs(H)); 
title('Low-Pass Butterworth Filter');

% High-pass butterworth filter

[x,y] = butter(3,0.5, 'high');
[H, W]= freqz(x,y);
subplot(2,2,2);
plot(W/pi, abs(H)); 
title('High-pass Butterworth Filter');

% Band-pass butterworth filter

[x,y] = butter(3,[0.4 0.6]); %pass band between 0.4-0.6
[H, W]= freqz(x,y);
subplot(2,2,3);
plot(W/pi, abs(H));
title('Band-pass Butterworth Filter');

% Band-stop butterworth filter

[x,y] = butter(3,[0.4 0.6] , 'stop'); %stop band between 0.4-0.6
[H, W]= freqz(x,y);
subplot(2,2,4);
plot(W/pi, abs(H));
title('Band-Stop - Butterworth');

% Chebychev-I Filter

% Low-pass Chebychev-I filter

[x,y] = cheby1(3,1,0.5,'low');  
[H, W]= freqz(x,y);
figure();
subplot(2,2,1);
plot(W/pi, abs(H)); 
title('Low-pass Chebychev-I filter');

% High-pass Chebychev-I filter

[x,y] = cheby1(3,1, 0.5, 'high');
[H, W]= freqz(x,y);
subplot(2,2,2);
plot(W/pi, abs(H));
title('High-pass Chebychev-I filter');

% Band-pass Chebychev-I filter

[x,y] = cheby1(3,1, [0.4 0.6]); %pass band between 0.4-0.6
[H, W]= freqz(x,y);
subplot(2,2,3);
plot(W/pi, abs(H));
title('Band-pass Chebychev-I filter');

% Band-stop Chebychev-I filter

[x,y] = cheby1(3,1, [0.4 0.6] , 'stop'); %stop band between 0.4-0.6
[H, W]= freqz(x,y);
subplot(2,2,4);
plot(W/pi, abs(H)); 
title('Band-stop Chebychev-I filter');



% Question 1.2

% With order 3
 
figure();
[x,y] = cheby1(3,1,0.5,'low');
[H, W]= freqz(x,y);
subplot(2,2,1);
plot(W/pi, abs(H));
title('Low-pass Chebyshev-I filter with order 3');

%With order 5

[x,y] = cheby1(5,1,0.5,'low');
[H, W]= freqz(x,y);
subplot(2,2,2);
plot(W/pi, abs(H)); 
title('Low-pass Chebyshev-I filter with order 5');

%with order 10

[x,y] = cheby1(10,1,0.5,'low');
[H, W]= freqz(x,y);
subplot(2,2,3);
plot(W/pi, abs(H));
title('Low-pass Chebyshev-I filter with order 10');

%with order 20

[x,y] = cheby1(20,1,0.5,'low');
[H, W]= freqz(x,y);
subplot(2,2,4);
plot(W/pi, abs(H));
title('Low-pass Chebyshev-I filter with order 20');


% By increasing the order, the gradient of low pass becomes more steep, 
% which means it gets close to look like an ideal low pass filter, but 
% at the same time the ripple caused by the filter increasing with the 
% order number.







% Exercise 2: 

% Question 2.1

% Calling dirac function previously defined in previous labs

drc=dirac(40, 20); %dirac with n=20 and N=40
figure();   
plot(drc); 
title('Dirac Function'); xlabel('k'); ylabel('x(k)');


%Setting value of scaling

scl = 0.5; 
Ts = 1; 
alp = scl*Ts;  
a_s = exp(-alp) ;


%Anti-Causal part of smoothing filter

d_a_s = zeros(length(drc),1);
d_l=length(drc)-2:-1:1;
for i =  d_l 
 d_a_s(i) = scl*alp*a_s*dirac(i+1)+(2*a_s)*d_a_s(i+1)-(a_s^2)*d_a_s(i+2);
end
figure()
stem (d_a_s) ;
title('Anti-causal Smoothing'); 


%Causal part of smoothing

d_c_s = zeros(length(drc),1);
for i = 3:length(drc);
 d_c_s(i) = -scl*alp*a_s*dirac(i-1)+(2*a_s)*d_a_s(i-1)-(a_s^2)*d_a_s(i-2);
end    
figure()
stem (d_c_s) ;
title('Causal Smmothing'); 


% Step input Signal

step10= step(40,10)
step30= step(40,30)
stepp=step10-step30;
figure()
stem(stepp)

% Causal Derivative function

s_c = zeros(length (stepp),1) ;
for i = 3 : length(stepp)
 s_c(i) = stepp(i)+a_s*(alp-1)*stepp(i-1)+(2*a_s)*s_c(i-1)-(a_s^2)*s_c(i-2) ;
end
figure()
stem (s_c) ;
title('Causal Deravative');

% Anti-causal derivative function

s_a = zeros(length (stepp),1);
s_l = length(stepp)-2 : -1 : 1 ; 
for i = s_l  
 s_a(i) = a_s*(alp+1)*step(i+1)-(a_s^2)*step(i+2)+(2*a_s)*s_a(i+1)-(a_s^2)*s_a(i+2) ;  
end
figure()
stem (s_a);
title('Anti-Causal Deravative');






% Exercise 3: 

% Question 3

brbra = imread('C:\Users\mohitkumarahuja\Documents\GitHub\DSP-TP-1617\barbara.gif');
figure();
imshow(brbra);  
case1 = zeros(size(brbra)); 
case2 = zeros(size(brbra)); 


for i = 1:size(brbra, 2)
    img1 = brbra(:,i);
    
    r_c = zeros(length (img1),1) ;
    for i = 3 : length(img1)
     r_c(i) = img1(i)+a_s*(alp-1)*img1(i - 1)+(2*a_s)*r_c(i-1)-(a_s^2)*r_c(i-2) ;
    end
    
    r_a = zeros(length (img1 ),1) ;
    b_l = length(img1)-2 : -1 : 1 ;
    for i =  b_l
     r_a(i) = a_s*(alp+1)*img1(i+1)-(a_s^2)*img1(i+2)+(2*a_s)*r_a(i+1)-(a_s^2)*r_a(i+2) ;
    end
    
    response = r_c + r_a;
    
    case1(:,i) = response;            
end
figure();
imshow (case1, []);

for i = 1:size(brbra, 2)

    img2 = brbra(i,:);
    
    r_c = zeros(length (img2),1) ;
    for i = 3 : length(img2)
     r_c(i) = img2(i)+a_s*(alp-1)*img2(i - 1)+(2*a_s)*r_c(i-1)-(a_s^2)*r_c(i-2) ;
    end
    
    r_a = zeros(length (img2 ),1) ;
    b_l = length(img2)-2 : -1 : 1 ;
    for i =  b_l
     r_a(i) = a_s*(alp+1)*img2(i+1)-(a_s^2)*img2(i+2)+(2*a_s)*r_a(i+1)-(a_s^2)*r_a(i+2) ;
    end
    response2 = r_c + r_a;
    
    case2(i,:) = response2;   
end
figure();
imshow (case2, []);






function S1 =  dirac(n,N)                % Function Defination 
 
    if ((n<1)||(n>N))
        
        disp('n is greater than N-1');   % Display error
        S1= 0;
    
    else
        
        s = zeros(1,N);                   %Dirac f/n
        s(n) = 1 ;
        S1 = s;
           
    end
  
end




function S3 =  step(n,N)                 % Function Defination 
    
    if ((n<1)||(n>N))
        
        disp('n is greater than N-1');   % Display error
        S3= 0;
        
    else
        
        s = zeros(N,1);
        
        for i = n+1:N
            
            s(i) = 1 ;
            
        end
        
        S3 = s;

        figure(2)
        subplot(3,1,1)
        stem(S3) ;                       % Unit Step
        xlabel('X'); ylabel('Y')
           
    end
  
end 




end