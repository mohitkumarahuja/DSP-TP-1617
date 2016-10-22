function Assignment_Exr_1()




      %MOHIT KUMAR AHUJA%
      %MSCV_2016%






% Exercise 1 :  

% Question 1.1

d = Dirac(0,20);                    %Dirac Function 
figure(1); 
stem(d);  
title('Dirac Function/Ques 1.1');   %Title for plot  
xlabel('n'); ylabel('delta(k-n)');


% Question 1.2

Stp =Step(10,20);                   %Step Function
figure(2); 
stem(Stp); 
title ('Step Function n=10/Ques 1.2');%Title for plot 
xlabel('n'); ylabel('H(k-n)');      % H is a Step Function

% Question 1.3

Rmp = Ramp(2,10,20);                %Ramp Function
figure(3); 
stem(Rmp); 
title ('Ramp Function a=2, n=10/Ques 1.3');%Title for plot 
xlabel('n'); ylabel('P(k-n)');     %P is a Ramp Function

% Question 1.4

Gmt = Geo(2, 10, 20);
figure(4); 
stem(Gmt); 
title ('Geometric Function a=2, n=10/Ques 1.4');%Title for plot 
xlabel('n'); ylabel('G(k-n)');      %G is a Geometric Function

% Question 1.5

Bx = box(3, 10, 20);
figure(5); 
stem(Bx); 
title ('Box Function a=3, n=10/Ques 1.5');%Title for plot 
xlabel('n'); ylabel('B(k-n)');      %B is a Box Function

% Question 1.6.1

Snfn_1 = Sinfn(1,10,100);
figure(6); 
stem(Snfn_1); 
title ('Sin Function fs=100/Ques 1.6.1');%Title for plot 
xlabel('n'); ylabel('sin(2*pi*f*t)');

% Question 1.6.2

Snfn_2 = Sinfn(2, 10, 1000);
figure(7); 
stem(Snfn_2); 
title ('Sin Function fs=1000/Ques 1.6.2');%Title for plot 
xlabel('n'); ylabel('sin(2*pi*f*t)');

% Question 1.6.3

Snfn_3 = Sinfn(2, 10, 30);        
figure(8); 
stem(Snfn_3); 
title ('Sin Function fs=30/Ques 1.6.3');%Title for plot 
xlabel('n'); ylabel('sin(2*pi*f*t)');

end







% Defining all Functions %

% Question 1.1

function d = Dirac(n, N)              %Function Creation

if  n >(N-1)
    disp('n is greater than N-1 ');   %Displayig error
else
    M = zeros(1,N);
    M(n+1) = 1; 
    d = M;

end
end


% Question 1.2

function Stp = Step(n, N)             %Function Creation

if  n >(N-1)
    disp('n is greater than N-1 ');   %Displayig error
else
    S = zeros(1,N);

for a = n+1:N 
    S(a) = 1;   
end
 
Stp = S;

end
end


% Question 1.3

function Rmp = Ramp(a,n, N)          %Function Creation


t= n+1:1:N;
size(t);


if  n > (N-1)
    Disp('n is greater than N-1 ');  %Displayig error
    
else
    R = zeros(1,N);
    
for j = n:N
    R(j) = a*(j-n);
    
end

Rmp = R;

end

end


% Question 1.4

function Gmt = Geo(a, n, N)          %Function Creation

if  n > (N-1)
    disp('n is greater than N-1 ');  %Displayig error
    
else
    G = zeros(1,N);
    for i = n+1:N; 
    G(i) = 1;
    
end

k=0;

for j = n+1:N 
    G(j) = (a^k).*G(j);
    k=k+1;
    
end
 
Gmt = G;

end
end


% Question 1.5

function Bx =  box(a,n,N)             %Function Creation

    if n >(N-1)
            Disp('n is greater than N-1'); %Displayig error 
    else
            
            s = zeros(1,N);  
            for j = (n-a):(n+a)
                s(j) = 1 ;
            end
            Bx = s;
                  
    end
  
end

% Question 1.6

function Snfn = Sinfn(n,f,Ts)         %Function Creation

t = (0:1/Ts:n/f);
Snfn = sin(2*pi*f*t);

end

