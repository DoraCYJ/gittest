clear ;
close all;
%% 
% Path = 'C:\Users\13704\Desktop\jesd_param_matlab'; 
Path = pwd; % 存放数据的文件夹路径(绝对路径)
Mode = 'NAmp'; % 输入是否接放大器, Amp(有放大器),NAmp(无放大器)
% Signal_Type = 'Single';
Signal_Type = 'QPSK';
f = 600; % MHz % 载波频率
file_num = 0 ;
FileType = 'csv';
FileName = strcat(Path,'\',Mode,'\',Signal_Type,'\f',num2str(f),'m',num2str(file_num),'.',FileType);
FileName_noise = strcat(Path,'\',Mode,'\',Signal_Type,'\f','_noise',num2str(f),'m',num2str(file_num),'.',FileType);
fs = 4800; % MHz 采样频率
N  = fs/f;
% T  = 1e-3; % 采样持续时间, s
% Nt = fs*1e6*T;
% freq =(0:Nt/2)*fs/Nt;

M = readmatrix(FileName);
Trig = M(:,3);
index = find(Trig==1)+1;
% rd_en = M(:,36);
% index_rd_en = find(rd_en==1)+1;
% Trig_Data1 = M(index_rd_en,4:35);
Trig_Data1 = M(index:end,4:35);
row = 1:32;
number_l = row(1:16);
number_h = row(17:32);
number   = [fliplr(number_h),fliplr(number_l)];

sl = Trig_Data1(:,1:16);
sh = Trig_Data1(:,17:32);

s = Trig_Data1;

% s = [sl,sh];
% s = [fliplr(sh),fliplr(sl)];
% s = fliplr(Trig_Data1); % before
% s    = [fliplr(sl(:,1:16)),fliplr(sh(:,1:16))];
%% 
Num = numel(s);
signal = reshape(s.',1,Num);

%% 
% y = signal(1,((N*0+1):(N*16))+5000);
y = signal(1,((N*0+1):(N*8))+5000);
figure
% subplot(1,2,2);
plot(y, 'b');
title('时域(一小段)');

F_range = 40; % MHz
Nt = length(signal);
freq =(0:Nt/2)*fs/Nt;
index_f1 = max(find(freq<( f-(F_range/2) )));
index_f2 = min(find(freq>( f+(F_range/2) )));
figure
S_FFT = abs(fft(signal,Nt));
S_FFT = S_FFT(1:Nt/2+1);
% [Value,Position] = sort(S_FFT(1:Nt/2+1),'descend');
subplot(1,2,1);
plot( freq(index_f1:index_f2),log10(S_FFT(index_f1:index_f2)), 'b' );
% plot( freq(Rang_Lower:Rang_Upper),S_FFT((Rang_Lower:Rang_Upper)), 'b' );
xlabel('MHz');
legend('无放大器','location','northwest'); 
title('频域');
grid on;
y = signal(1,((N*0+1):(N*20))+5000);
% figure
subplot(1,2,2);
plot(y, 'b');
title('时域(一小段)');



 %% 
% 
% clear ;
% close all;
% %% 
% % Path = 'C:\Users\13704\Desktop\jesd_param_matlab'; 
% Path = pwd; % 存放数据的文件夹路径(绝对路径)
% Mode = 'NAmp'; % 输入是否接放大器, Amp(有放大器),NAmp(无放大器)
% Signal_Type = 'Single';
% f = 600; % MHz % 载波频率
% file_num = 0 ;
% FileType = 'csv';
% FileName = strcat(Path,'\',Mode,'\',Signal_Type,'\f',num2str(f),'m',num2str(file_num),'.',FileType);
% FileName_noise = strcat(Path,'\',Mode,'\',Signal_Type,'\f','_noise',num2str(f),'m',num2str(file_num),'.',FileType);
% fs = 4800; % MHz 采样频率
% N  = fs/f;
% % T  = 1e-3; % 采样持续时间, s
% % Nt = fs*1e6*T;
% % freq =(0:Nt/2)*fs/Nt;
% 
% M = readmatrix(FileName);
% Trig = M(:,3);
% index = find(Trig==1)+1;
% rd_en = M(:,36);
% index_rd_en = find(rd_en==1)+1;
% Trig_Data1 = M(index_rd_en,4:35);
% sl = Trig_Data1(:,1:16);
% sh = Trig_Data1(:,17:32);
% 
% s = [sl,sh];
% 
% %% 
% Num = numel(s);
% signal = reshape(s.',1,Num);
% 
% %% 
% % y = signal(1,((N*0+1):(N*16))+5000);
% y = signal(1,((N*0+1):(N*16))+5000);
% figure
% % subplot(1,2,2);
% plot(y, 'b');
% title('时域(一小段)');
% 
% Nt = length(signal);
% freq =(0:Nt/2)*fs/Nt;
% 
% figure
% S_FFT = abs(fft(signal,Nt));
% [Value,Position] = sort(S_FFT(1:Nt/2+1),'descend');
% subplot(1,2,1);
% plot( freq,log10(S_FFT(1:Nt/2+1)), 'b' );
% % plot( freq(Rang_Lower:Rang_Upper),S_FFT((Rang_Lower:Rang_Upper)), 'b' );
% xlabel('MHz');
% legend('无放大器','location','northwest'); 
% title('频域');
% grid on;
% y = signal(1,((N*0+1):(N*20))+5000);
% % figure
% subplot(1,2,2);
% plot(y, 'b');
% title('时域(一小段)');


