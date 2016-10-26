function DSP_LAB_2()


      


        %MOHIT KUMAR AHUJA%
        %MSCV_2016%




% Reminder :

fs = 20 ;                    % Sampling freq
f = 1;                       % Frequency
t = (0:1/fs:10); 
S = sin(2*pi*f*t);
figure(1)
plot(t,S);
title('Sin Function')
xlabel('T'); ylabel('Sin(2*pi*f*t1)')
         
t1 = (0:1/fs:20);
S1 = sin(2*pi*(f/fs)*t1);
figure(2)
plot(t1,S1);
title('Sin Function')   
xlabel('T'); ylabel('Sin(2*pi*(f/fs)*t2)')
         
         



% Exercise 1: 
         
% Question 1.1   

x = step(4,20);                % plotting Signal with Step

for k= 1:19
   
    y(k) = x(k)/2+ (x(k+1))/2; % it is non causal because of the future 
                               % values(k+1)

end


subplot(3,1,2); 
stem(S)
title('Non Causal System')  
xlabel('x(k)'); ylabel('y(k)')
  



% Question 1.2

for k= 2:20
    
    y(k) = x(k)/2+ (x(k-1))/2;  % Causal Signal by shifting its value to 
                                % past by putting k-1
end
        
subplot(3,1,3); 
stem(S)
title('Causal System')  
xlabel('x(k) '); ylabel(' y(k) ') 










% Exercise 2 :

% Question 2.1

S2 = x;                          %x is a step function

for i = 2:1:20
    
    S2(i) = S2(i-1)+x(i);
    
end

figure(3)        
subplot(2,1,1);
stem(S2)  
title('Accumulation')  
xlabel('x(k)'); ylabel('y(k)')
        

% Unstable system because it keep adding previous values till our range and
% the signal will keep increasing with increment of 1.
        
   



% Question 2.2

D = Dirac(4,20);                 % D is a step function
S3 = D;

for i = 2:1:20
    
    S3(i) = S3(i-1)+D(i); 
    
end
    
subplot(2,1,2); 
stem(S3)  
title('Stable System')  
xlabel('D(k)'); ylabel('y(k)')
      

% Because of the value of Dirac function i.e 1 on specific  index, but it
% is 0 previously and at the next step also which ressults in constant for 
% all set of values.




         
% Question 2.3

D = Dirac(4,20);
S4 = D;

for i = 2:1:20
    
    S4(i) = D(i)+ 2*(S4(i-1));    % Exponential increasing(Non-Zero)
    
end
       
figure(4)
subplot(2,1,1); 
stem(S4)  
title('Unstable System')  
xlabel('D(k)'); ylabel('y(k)')

 
        


% Question 2.4

D = Dirac(4,20);
S5 = D;
for i = 2:1:20
    S5(i) = D(i)+ (S5(i-1)/3);
end
subplot(2,1,2); 
stem(S5)  
title('Stable System')  
xlabel('D(k)'); ylabel('y(k)')
        
% Exponentially decaying(Value approaching to zero) : stable sytem.
        
   








% Exercise 3 :

% Question 3.1 
        
x1 = [0 0 0 0 1 2 3 4 5 0 0 0 0 0 0 0 0 0 0];
x2 = [0 0 0 0 0 0 0 0 0 4 3 2 1 0 0 0 0 0 0];
y1(1) = 0;
y2(1) = 0;


for i = 2:1:19-1
    
    y1(i) = 3*x1(i-1)-2*x1(i)+x1(i+1);
    
end
        
figure(5);
subplot(2,1,1); 
stem(y1)  
title('System X1')  
xlabel('x(1) '); ylabel('y(1)')

%System with input X1





%  Question 3.2

for i = 2:1:19-1
    
    y2(i) = 3*x2(i-1)-2*x2(i)+x2(i+1);
    
end

subplot(2,1,2); 
stem(y2)  
title('System X1')  
xlabel('x(2)'); ylabel('y(2)')

%System with input X2





% Question 3.3 

h = [1,-2,3];
x = x1+x2; 
y1 = conv(x,h);
figure(6)
stem(y1);
title('Convolution of two signal(Adding first and then convolving)')  
xlabel('X = conv(x1 + x2)'); ylabel('y')

% Adding first and then convolving


% Question 3.4 

y2 = conv(x1,h)+conv(x2,h);   
figure(7)
stem(y1);
title('Convolution of two signal(Convolving first and then adding)')  
xlabel('X = conv(x1) + conv(x2)'); ylabel('y)')
        
% Convolving first and then adding
% System is linear and invariant
 
       
        
        
 
 
 





% Defining all Functions %

  
function S1 =  Dirac(n,N)                % Function Defination 
 
    if ((n<1)||(n>N))
        
        disp('n is greater than N-1');   % Display error
        S1= 0;
    
    else
        
        s = zeros(1,N);
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