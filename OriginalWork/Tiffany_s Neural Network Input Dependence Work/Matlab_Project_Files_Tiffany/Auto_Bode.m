function [Amplitude_save,phase_save,Amplitude_ratio_stdev] = Auto_Bode(X,Y,N_test_cycles,N_points_per_cycle,N_mean_percent,N_initial_discard_percent)

    % David Turcic - ME Department - Portland State University
    
    Xin    = X;
    Xout   = Y;

    N_discard_cycles     = ceil(N_initial_discard_percent /100*N_test_cycles);
    N_test_cycles_AD     = N_test_cycles - N_discard_cycles;

    N_mean_discard     = ceil(N_mean_percent/100*N_test_cycles_AD);

    istart              = N_discard_cycles*N_points_per_cycle + 1 ;
    iend                = N_test_cycles*N_points_per_cycle;

    Ydata_in  = Xin(istart:iend);
    Ydata_out = Xout(istart:iend);

%          figure
%          plot(Xin)
%          title(['Xin    frequency = ',num2str(jj)],'Interpreter','none')
%          hold on;

    istart2 = 1;

    for ii=1:N_test_cycles_AD
        iend2 = istart2 + N_points_per_cycle - 1;
        [ymax,iymax] = max(Ydata_in(istart2:iend2));
        [ymin,iymin] = min(Ydata_in(istart2:iend2));
        iymax = iymax + istart2 + istart - 2;
        iymin = iymin + istart2 + istart - 2;

%           plot(iymax,ymax,'or')
%           plot(iymin,ymin,'og')

        iymax_in_save(ii) = iymax;
        iymin_in_save(ii) = iymin;
        ymax_in_save(ii)  = ymax;
        ymin_in_save(ii)  = ymin;

        istart2 = istart2 + N_points_per_cycle;
    end

%          figure
%          plot(Xout)
%          title(['Xout    frequency = ',num2str(jj)],'Interpreter','none')
%          hold on;

    istart2 = 1;

    for ii=1:N_test_cycles_AD
        iend2 = istart2 + N_points_per_cycle - 1;
        [ymax,iymax] = max(Ydata_out(istart2:iend2));
        [ymin,iymin] = min(Ydata_out(istart2:iend2));
        iymax = iymax + istart2 + istart - 2;
        iymin = iymin + istart2 + istart - 2;

%           plot(iymax,ymax,'or')
%           plot(iymin,ymin,'og')

        iymax_out_save(ii) = iymax;
        iymin_out_save(ii) = iymin;
        ymax_out_save(ii)  = ymax;
        ymin_out_save(ii)  = ymin;

        istart2 = istart2 + N_points_per_cycle;
    end

    mean_start = N_mean_discard+1;
    mean_end   = N_test_cycles_AD - N_mean_discard;

    phase          = -(iymax_out_save - iymax_in_save)*360/N_points_per_cycle;
    phase          = sort(phase);
    phase_m        = mean(phase((mean_start):(mean_end)));
    phase_save     = phase_m;

    Amplitude_in   = abs(ymax_in_save - ymin_in_save);
    Amplitude_in_vec = Amplitude_in((mean_start):(mean_end));   % do this pre-sort to get accurate stdev (arbitrarily cuts off first and last instead of largest and smallest in order to have same number of data points)
    Amplitude_in   = sort(Amplitude_in);
    Amplitude_in_m = mean(Amplitude_in((mean_start):(mean_end)));

    Amplitude_out   = abs(ymax_out_save - ymin_out_save);
    Amplitude_out_vec = Amplitude_out((mean_start):(mean_end));
    Amplitude_out   = sort(Amplitude_out);
    Amplitude_out_m = mean(Amplitude_out((mean_start):(mean_end))); 
    
    Amplitude_ratio_vec     = Amplitude_out_vec./Amplitude_in_vec;
    Amplitude_ratio_stdev   = std(Amplitude_ratio_vec);

    Amplitude_save  = Amplitude_out_m/Amplitude_in_m;

end
