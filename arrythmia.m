
%Load data

        close all;
        clear all;
        Data = load('207.txt')

        x = Data(:,1);
        y = Data(:,2);
        figure;
        plot (x, y);
        xlabel('Time (s)');
        ylabel('Amplitude (V)');
        title('RAW DATA')

%FFT_RAW_SIGNAL

       xdft = fft(y);    
       
      % sampling frequency 
       t_diff_vect = diff(x);
       DT = mean(t_diff_vect);
      
    
       Fs = round(1/DT)
       DF = Fs/length(x);
       freq = 0:DF:Fs/2;
       xdft = xdft(1:length(xdft)/2+1);
       figure;
       subplot(2,1,1);
       plot(freq,20*log10(abs(xdft)));
       xlabel('Frequency (Hz)');
       ylabel('Amplitude (db)');
       title('RAW DATA FFT')
       ylim([-40 80]);
 
%FILTER_SECTION

       [b,a] = butter(2,[59.9,60.1]/(Fs/2),'stop');
       filtered55 = filter(b,a,y);
    
       [d,c] = butter(10,100/(Fs/2),'low');
       filtered105 = filter(d,c,filtered55);
       
       [f,e] = butter(3,0.5/(Fs/2),'high');
       dataOut = filter(f,e,filtered105);

       
 %FFT_FILTERED_SIGNAL
 
       xdft2 = fft(dataOut);
       xdft2 = xdft2(1:length(xdft2)/2+1);
       subplot(2,1,2);
       plot(freq,20*log10(abs(xdft2)));
       xlabel('Frequency (Hz)');
       ylabel('Amplitude (db)');    
       title('FILTERED SIGNAL FFT')
       ylim([-40 80]);      
       
 %Peak detector     
 
       max=0;
       for i = 1:size(dataOut)-(3*Fs)
 
           for k = 1:(3*Fs)
                temp = dataOut(i+k);
                
               if(temp > max)
                max = temp;
               end                
           end
         tresholds(i) = max;   
         max=0;
         
       end

       peak_treshold = 0.4*mean(tresholds)

       figure;
       plot(x(1:end-(3*Fs)),tresholds);
       xlabel('Time (s)');
       ylabel('Amplitude (V)');
       title('R-PEAK ENVELOPE')
       
       figure;
       MAXW = 0.01;
       MPD = 0.3;
       [pks,locs] = findpeaks(dataOut,x,'MinPeakHeight', peak_treshold,'MinPeakWidth', MAXW,'MinPeakDistance',MPD);
       plot(x,dataOut,locs,pks,'o');
       xlabel('Time (s)');
       ylabel('Amplitude (v)');
       title('PEAKS FOUND')

       
%Peak distances

       for i=2:size(locs)
          peak_diffs(i-1) = locs(i)-locs(i-1);
       end
       
        locs = locs(2:end);
    
       figure;
       plot(locs, peak_diffs*1000);
       xlabel('Time (s)');
       ylabel('R-R (ms)');
       title('RR-TIME INTERVALS')
       
%Simple HR plot

       HR = 60./peak_diffs;
       
       figure;
       plot(locs, HR);
       xlabel('Time (s)');
       ylabel('HR (bpm)');
       title('HEART RATE')
       ylim([40 220]);       
       
       %Statistics

       RR_mean = mean(peak_diffs)
       HR_mean = mean(HR)
       
       RR_STD = std(peak_diffs)
       HR_STD = std(HR)

%PoinCare

       figure;

           plot(1000*peak_diffs(1:end-1), 1000*peak_diffs(2:end),'.');
       %end
       xlabel('R-R(n) (ms)');
       ylabel('R-R(n+1) (ms)');
       title('POINCARE HRV')
       
%Arrhythmia

RR_category = ones(1,(length(peak_diffs)));


      for i=2:(length(peak_diffs)-1)
          
          
          rr1 = peak_diffs(i-1);
          rr2 = peak_diffs(i);
          rr3 = peak_diffs(i+1);

          %Ventricular flutter/fibrillation
          if(rr2 < 0.6 && (1.8*rr2 < rr1));
             RR_category(i) = 3;
             
             for k = i+1:(size(peak_diffs)-1)
                 
                 rr1k = peak_diffs(k-1);
                 rr2k = peak_diffs(k);
                 rr3k = peak_diffs(k+1);
                 
                if((rr1k < 0.7 && rr2k < 0.7 && rr3k < 0.7) || (rr1k+rr2k+rr3k < 1.7))
                    RR_category(i) = 3;
                else %Less than 4 cycles = normal
                    if(RR_category(i-4)== 1 && RR_category(i-3)== 3 && RR_category(i-2)== 3 && RR_category(i-1)== 3)
                        RR_category(i-3) = 1;
                        RR_category(i-2) = 1;
                        RR_category(i-1) = 1;
                    end    
                end
             end
          end 
         
          %Premature ventricular contractions
          if(   ((1.15*rr2 < rr1)&&(1.15*rr2 < rr3))...
             ||((abs(rr1-rr2)< 0.3) && (rr1 < 0.8 && rr2 < 0.8)...
             && (rr3 > 1.2*((rr1+rr2)/2)))...
             ||(abs(rr2-rr3)< 0.3)&&(rr2 < 0.8 && rr3 < 0.8)...
             && (rr1 > 1.2 * ((rr2+rr3)/2))    )
                RR_category(i) = 2;
          end
          
         %2° heart block 
          if((2.2 < rr2 && rr2 < 3.0 )&&((abs(rr1-rr2) < 0.2) || (abs(rr2-rr3) < 0.2)))
              RR_category(i) = 4;
          end
      end
       
    
       figure;
       plot(locs, RR_category, 'o')
       
       xlabel('Time (s)');
       ylabel('RR category');
       title('HEART BEAT CATEGORY')
       ylim([0 5]);  

       yticks([1 2 3 4])
       yticklabels({'Normal','Premature ventricular contractions','Ventricular flutter / fibrillation','2° heart block'})
       
  
       cat1_tot = 0;
       cat2_tot = 0;
       cat3_tot = 0;
       cat4_tot = 0;
       %Sum of categories
      for i =1:length(peak_diffs)
         if (RR_category(i) == 1) 
            cat1_tot = cat1_tot + 1;
         end
         if (RR_category(i) == 2) 
            cat2_tot = cat2_tot + 1;
         end
         if (RR_category(i) == 3) 
            cat3_tot = cat3_tot + 1;
         end
         if (RR_category(i) == 4) 
            cat4_tot = cat4_tot + 1;          
         end
      end
      
       cat1_tot
       cat2_tot
       cat3_tot
       cat4_tot
       
       