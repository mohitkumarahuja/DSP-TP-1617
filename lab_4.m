function lab_4()







            %%%%% Mohit Kumar Ahuja%%%%
            %%%%%     MSCV 2016    %%%%
            


            
            
% Exercise 1: DFT

% Question 1.1: 
 
f = 5;
fs = 50;
t = 0: 1/fs : 1;
xn = sin(2 * pi * f * t);
N = length(xn);
fr = (-N/2 : N/2-1)* fs/N;
xf = fftshift(fft(xn));
figure(1);
subplot(221); plot(t, xn); 
title('Sin Signal'); 
xlabel('Time(sec)'); ylabel('Amplitude');
subplot(222); plot(fr, abs(xf)); 
title('Magnitude'); 
xlabel('Frequency'); ylabel('|X(f)|');
subplot(223); plot(fr, real(xf)); 
title('Real'); 
xlabel('Frequency'); ylabel('Re(X(f))');
subplot(224); plot(fr, imag(xf)); 
title('Imaginary'); 
xlabel('Frequency'); ylabel('Im(X(f))');


% Question 1.2: 

xn2 = cos(2 * pi * f * t);
N = length(xn2);
fr = (-N/2 : N/2-1)* fs/N;
xf = fftshift(fft(xn2));
figure(2);
subplot(221); plot(t, xn2); 
title('Cos Signal'); 
xlabel('Time(sec)'); ylabel('Amplitude');
subplot(222); plot(fr, abs(xf)); 
title('Magnitude'); 
xlabel('Frequency'); ylabel('|X(f)|');
subplot(223); plot(fr, real(xf)); 
title('Real'); 
xlabel('Frequency'); ylabel('Re(X(f))');
subplot(224); plot(fr, imag(xf)); 
title('Imaginary'); 
xlabel('Frequency'); ylabel('Im(X(f))');



% Question 1.3:

t2 = 0: 1/fs : 2;
xn3 = square(2 * pi * f * t2);
N = length(xn3);
fr = (-N/2 : N/2-1)* fs/N;
xf = fftshift(fft(xn3));
figure(3);
subplot(221); plot(t2, xn3); 
title('Cos Signal'); 
xlabel('Time(sec)'); ylabel('Amplitude');
subplot(222); plot(fr, abs(xf)); 
title('Magnitude'); 
xlabel('Frequency'); ylabel('|X(f)|');
subplot(223); plot(fr, real(xf)); 
title('Real'); 
xlabel('Frequency'); ylabel('Re(X(f))');
subplot(224); plot(fr, imag(xf)); 
title('Imaginary'); 
xlabel('Frequency'); ylabel('Im(X(f))');


% Question 1.4:

xn4 = randn(1,10000);
N = length(xn4);
fr = (-N/2 : N/2-1)* fs/N;
xf = fftshift(fft(xn4));
figure(4);
subplot(221); plot(xn4); 
title('Cos Signal'); 
xlabel('Time(sec)'); ylabel('Amplitude');
subplot(222); plot(fr, abs(xf)); 
title('Magnitude'); 
xlabel('Frequency'); ylabel('|X(f)|');
subplot(223); plot(fr, real(xf)); 
title('Real'); 
xlabel('Frequency'); ylabel('Re(X(f))');
subplot(224); plot(fr, imag(xf)); 
title('Imaginary'); 
xlabel('Frequency'); ylabel('Im(X(f))');





% Exercise 2: Sampling

% Question 2.1: 

f1 = 5;
f2 = 20;
fs = [10 20 25 40 50 100 150];
for i = 1: length(fs)
    n = 0: 1/fs(i) : 1;
    xn5 = (3*cos(2*pi*f1*n)+(4*sin(2*pi*f2*n)));
end

% Question 2.4:

N = 1000;
fr5 = (-N/2 : N/2-1)* fs(i)/N;
xf5 = fftshift(fft(xn5,1000));

% Question 2.2:
figure(5);
title('Signal in Time Domain'); 
subplot(221); plot( n , xn5); 
title('Signal'); 
xlabel('Time(sec)'); ylabel('Amplitude');
subplot(222); plot(fr5, abs(xf5)); 
title('Magnitude'); 
xlabel('Frequency'); ylabel('|X(f)|');
subplot(223); plot(fr5, real(xf5)); 
title('Real'); 
xlabel('Frequency'); ylabel('Re(X(f))');
subplot(224); plot(fr5, imag(xf5)); 
title('Imaginary'); 
xlabel('Frequency'); ylabel('Im(X(f))');



% Question 2.3:

% Aliasing effect: 
% When the sampling freqency(fs) is less than the twise of the highest 
% frequency there will be aliasing effect.
% And when our sampling frequency(fs) is equal to twise highest frequency,
% there will be no aliasing. so to avoid aliasing we should put the value
% ofour sampling frequency either twise or more than twise of the signal
% frequency. i.e when our fs >= 40hz, there will be no aliasing.






% Exercise 3: 1D DFT for image classification


% Question 3.1: 

image_folder='C:\Users\mohitkumarahuja\Downloads\DSP-TP-1617-master\images\1D-DFT'; 
cd (image_folder)
files=dir;
amnt_of_img = length(files)-2; 
 for i = 1:amnt_of_img
     img{i} = imread(files(i+2).name);  
 end
 
[y,x] = size(img{1});            %storing the size of the first index image 
 for idx = 2:amnt_of_img
     [a,b]=size(img{idx});       %Comparing the size of all the images in 
     if a*b<x*y                  %the directory and saving the lowest size 
     y=a;                        %in x and y
     x=b;
     end
 end
 
 p=1;q=1;
 for i=1:amnt_of_img
    I = imread(files(i+2).name);
    I = double(I);
    [a,b,c]= size(I);         % loading the size of the image in a,b and c
        if c==4               % checking if the tif is 4 dimensional or not
            I=I(:,:,1);       % if four dimension then making it 1.
        end
    I=I/(max(I(:)));          %normalizing the image
    I_resized = imresize(I, [y x]); %resizing the image to the lowest size 
    z=round(y/2/2);                 %taking the middle row from the centre of image
    img_1D =I_resized(z,:);     %storing the z row into a new matrix

    N=x;
    fr = (-N/2 : N/2-1);
    x1 = ifftshift(fft(img_1D)); 
    
    % Using DFt transfrom on the 1d Profile
    % We know that in Bar code images only vertical straight white and 
    % black lines will be there so other rows should have the same spectrum
    % so taking that into consideration, i will take two more rows to 
    % double cross my results. one from the centre and the other one from 
    % the middle of centre and bottom of the screen.
    
    
    %%%% We can increase the lines to crosscheck our results %%%%
    %%%%     The more the lines = the better the accuracy    %%%%
    
    
    % Hence calculating DFT of the two other 1D profiles
    img_1D2 = I_resized(z+10,:);
    img_1D3 = I_resized(z+25,:); 
    x2 = ifftshift(fft(img_1D2));
    x3 = ifftshift(fft(img_1D3));
    
    % Finding the absolute values of each one 1D profile  
    v1=abs(x1);
    v2=abs(x2);
    v3=abs(x3);
    
    % Taking maximum values of each 1D profile
    threshold1=max(v1);
    threshold2=abs(2*max(v1)-max(v2)-max(v3)); % Storing the difference of the three 1D profiles maximum value 
    threshold(i) = threshold1*threshold2; % Multiplication of the two thresholds and creating an offset

    if (threshold(i)<48) % According to the threshold seen at the output, from obeservation a   
                         % threshold less than 48 gives the highest percentage of Barcode Image
        Brcd_image(p)=i ;% If image is barcode, the index of the image is stored 
        p=p+1; 
    else
        NonBrcd_image(q)=i; % If image is barcode, the index of the image is stored
        q=q+1;
    end
 end
 figure;
 stem(threshold);         %showing the threshold of each image in one graph
 disp('Barcode');
 disp(Brcd_image); 
 disp('Non barcode');
 disp(NonBrcd_image); %displaying the images that are barcode and the ones nonbarcode
 Barcode=[1,2,6,44:54]; 
 NonBarcode = [3,4,5,7:43];
 
 %Looking for the images that are wrongly classified
 x=length(Barcode);
 u=0;v=0;
 for i=1:length(Barcode)
   for j=1:length(Brcd_image)
     if Brcd_image(j)==Barcode(i)
         u=u+1;
     end
   end
 end
 
 for i=1:length(NonBarcode)
   for j=1:length(Brcd_image)
     if Brcd_image(j)==NonBarcode(i)
         v=v+1;
     end
   end
 end
 %storing the wrong images in v
 v=v+x-u;
 
 %Calculating the percentage of accuracy
 prcnt_of_acc= ((amnt_of_img-v)/amnt_of_img)*100;
 sprintf('Percentage of correct distinction between the Barcode and NonBarcode images is: %d percent', round(prcnt_of_acc))
 
 
 
        % Sorry for the late Submission %
        
        
end





