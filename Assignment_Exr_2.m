function Assignment_Exr_2 ()

     




      %MOHIT KUMAR AHUJA%
      %MSCV_2016%







%  Exercise 2 : Random Signals

    %Question 2.1 :

G1 = Gaussian(100);                                        % Function
figure(1); 
hist (G1);
histfit(G1);
xlabel ('Xn')
ylabel ('N')
title('Histogram of gaussian function')                    % Title for plot 
    
G2 = Gaussi(1000);                                         % Function
figure(2); 
hist (G2);
histfit(G2);
xlabel ('Xn')
ylabel ('N')
title('Histogram of gaussian function with 1000 samples')  % Title for plot

G3 = Gauss(10000);                                         % Function
figure(3); 
hist (G3);
histfit(G3);                                              
xlabel ('Xn')                                              
ylabel ('N')
title('Histogram of gaussian function with 10000 samples') % Title for plot

 % As it is clear from the output diagram, that the increase in sampling
 % frequency will lead to the histogram looks like Gaussian Distribution


function G1 = Gaussian(N)                                %Function Creation

i = mean(randn(1,N));
j = std(randn(1,N));
G1= i.*randn(1,N) + j;

end

function G2 = Gaussi(N)                                  %Function Creation

i = mean(randn(1,N));
j = std(randn(1,N));
G2= i.*randn(1,N) + j;

end
 
function G3 = Gauss(N)                                   %Function Creation
 
i = mean(randn(1,N));
j = std(randn(1,N));
G3= i.*randn(1,N) + j;

end








   % Question 2.2 :

 G4 = rand(1,1000);                    % Value of Xu=1000
 SD = 1; 
 mn = 0;
 figure(4); 
 hist(G4,100)                          % Plotting Histogram
 title('Histogram with xu=1000')       % Title for plot
 xlabel('Xu')
 ylabel('U')
 
 Thrtcl_dist = unifpdf(G4,mn,SD);
 figure(5); 
 stem(G4,Thrtcl_dist);
 title('Thr Dstrbn with xu=1000')
 xlabel('Xu')
 ylabel('U')
 
 G5 = rand(1,10000);                    % Value of Xu=10000
 figure(6);  
 hist(G5,100)                         % Plotting Histogram
 title('Histogram with xu=10000')     % Title for plot
 xlabel('Xu')
 ylabel('U')
 
 Thrtcl_dist = unifpdf(G5,mn,SD);
 figure(7); 
 stem(G5,Thrtcl_dist);
 title('Thr dstrbn with xu=10000')  % Title for plot
 xlabel('Xu')
 ylabel('U')
 
 % The change is in prospective of Intensity with uniform distribution
 
 
 
 
 
 
 
 
    % Question 2.3 :  

 [r1,lags1] = xcorr(G1);                 %Cross-correlation
 figure(8);
 plot(lags1,r1)
 title('Autocorrelation xn=1000 ')       % Normal/Gaussian random process
 xlabel('Array of Random numbers')
 ylabel('Intensity at each index)')
 

 % As seen in the output figure that it has a peak in the middle and 
 % variation on the sides.

 [r2,lags2] = xcorr(G2);
 figure(9);
 plot(lags2,r2)
 title('Autocorrelation xn=1000')       % Normal/Gaussian random process
 xlabel('Xn')
 ylabel('N')
 
 % As the number of samples are increasing(1000), there is very small
 % changes in the result(very less variation in the sides)
 % Changes cn be better seen on zooming 
  
 [r3,lags3] = xcorr(G3);
 figure(10);
 plot(lags3,r3)
 title('Autocorrelation xn=10000')       % Normal/Gaussian random process
 xlabel('Xn')
 ylabel('N')
 
 % As the number of samples are increasing(10000), there is very small
 % changes in the result(very less variation in the sides)
 % Changes cn be better seen on zooming 
 
 [r4,lags4] = xcorr(G4);
 figure(11);
 plot(lags4,r4)
 title('Autocorrelation xu=1000')      % For Random process U
 xlabel('Xu')
 ylabel('U')
 
 % The result we get now is just the oppposite as we get earlier(there is a
 % huge variation in the centre but not in the sides)
 
 
 [r5,lags5] = xcorr(G5);
 figure(12);
 plot(lags5,r5)
 title('Autocorrelation xu=10000')      % For Random process U
 xlabel('Xu')
 ylabel('U')
 
 % As the number of samples are increasing(10000), there is very small
 % changes in the result(very less variation in the sides)
  
  
 
  
  
  
 
 
    % Question 2.4 :
    
    
  s1 = round(rand(1,50));              % Generating binary random signal 1
  s2 = round(rand(1,50));              % Generating binary random signal 2
  s3 = round(rand(1,50));              % Generating binary random signal 3
  s=[zeros(1,15) s1 zeros(1,15) s2 zeros(1,15) s3 zeros(1,15)]; % Generating Whole signal with shift of 15
 
 figure(13);
 plot(s);                              %Plotting signal
 title('Whole Signal with shift of 15')
 xlabel('X')
 ylabel('Intensity')
 
 figure(14);
 [r1,lags1] = xcorr(s,s1);        %cross correlation of whole signal and s1
 plot(lags1,r1)                   %Plotting cross correlation of s1
 title('Cross-correlation of random signal s1 with the Whole Signal')
 xlabel('Random signal 1')
 ylabel('Intensity')
 
 figure(15);
 [r1,lags1] = xcorr(s,s2);        %cross correlation of whole signal and s2
 plot(lags1,r1)                   %Plotting cross correlation of s2
 title('Cross-correlation of random signal s2 with the Whole Signal')
 xlabel('Random signal 2')
 ylabel('Intensity')
 
 
 figure(16);
 [r1,lags1] = xcorr(s,s3);        %cross correlation of whole signal and s3
 plot(lags1,r1)                   %Plotting cross correlation of s3
 title('Cross-correlation of random signal s3 with the Whole Signal')
 xlabel('Random signal 3')
 ylabel('Intensity')
 
 % The peaks for each Random signal while cross-correlation of that signal 
 % provides us a method to find signal patterns.
 
 end









