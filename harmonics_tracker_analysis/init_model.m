%[text] ## General Settings
%[text] ### Simulink model initialization
close all
clear all
clc
beep off
pm_addunit('percent', 0.01, '1');
options = bodeoptions;
options.FreqUnits = 'Hz';
simlength = 2;
s=tf('s');

model_name = 'harmonic_tracker';
%[text] ### PWM and sampling time and data length storage
fPWM = 10e3;
fPWM_AFE = fPWM; % PWM frequency 
tPWM_AFE = 1/fPWM_AFE;
fPWM_INV = fPWM; % PWM frequency 
tPWM_INV = 1/fPWM_INV;
fPWM_DAB = fPWM; % PWM frequency 
tPWM_DAB = 1/fPWM_DAB;

double_sampling = 0;

if double_sampling 
    ts_afe = 1/fPWM_AFE/2;
    ts_inv = 1/fPWM_INV/2;
    ts_dab = 1/fPWM_DAB/2;
else
    ts_afe = 1/fPWM_AFE;
    ts_inv = 1/fPWM_INV;
    ts_dab = 1/fPWM_DAB;
end

ts_battery = ts_dab;
tc = ts_inv/200;
decimation_tc = 10;

z_dab = tf('z',ts_dab);
z_afe = tf('z',ts_afe);
z_inv = tf('z',ts_inv);

t_misura = simlength;
% t_misura = 2;
Nc = ceil(t_misura/tc)/decimation_tc;
Ns_battery = ceil(t_misura/ts_battery);
Ns_dab = ceil(t_misura/ts_dab);
Ns_afe = ceil(t_misura/ts_afe);
Ns_inv = ceil(t_misura/ts_inv);

dead_time_INV = 0;
delay_pwm = 0;
%[text] #### Phase shift filter for Q component derivation at 50Hz and 80Hz
frequency_set = 50;
omega_set = 2*pi*frequency_set;

% gain for active sogi
kepsilon = 2;

a = 1 + 2*pi*frequency_set*ts_inv;
b = 1 - 2*pi*frequency_set*ts_inv;
phase_shift_filter_gain = 1;
phase_shit_filter_d = phase_shift_filter_gain * (1-a*z_inv^-1)/(1-b*z_inv^-1);
flt_dq = 2/(s/omega_set + 1)^2;
flt_dq_d = c2d(flt_dq,ts_inv);
%[text] #### Phase shift from SOGI
delta_sogi = 0.05;
sogi_flt_alpha = delta_sogi*omega_set*s/(s^2 + delta_sogi*omega_set*s + omega_set^2);
sogi_flt_beta = delta_sogi*omega_set^2/(s^2 + delta_sogi*omega_set*s + omega_set^2);
sogid_flt_alpha = c2d(sogi_flt_alpha, ts_inv);
sogid_flt_beta = c2d(sogi_flt_beta, ts_inv);
[sogid_flt_alpha_num, sogid_flt_alpha_den] = tfdata(sogid_flt_alpha, 'v');
[sogid_flt_beta_num, sogid_flt_beta_den] = tfdata(sogid_flt_beta, 'v');
% figure; bode(sogid_flt_alpha,sogid_flt_beta, options); grid on
%[text] #### Resonant PI
kp_rpi = 0.25;
ki_rpi = 45;
delta_rpi = 0.05;
res_nom = s/(s^2 + 2*delta_rpi*omega_set*s + (omega_set)^2);

Ares_nom = [0 1; -omega_set^2 -2*delta_rpi*omega_set] %[output:2c389496]
Aresd_nom = eye(2) + Ares_nom*ts_inv %[output:75942329]
a11d = 1 %[output:80bbc6d8]
a12d = ts_inv %[output:2928378c]
a21d = -omega_set^2*ts_inv %[output:68004a6c]
a22d = 1 -2*delta_rpi*omega_set*ts_inv %[output:485bb164]

Bres = [0; 1];
Cres = [0 1];
Bresd = Bres*ts_inv;
Cresd = Cres;
%%
%[text] ### Double Integrator Observer for PLL

Arso = [0 1; 0 0];
Crso = [1 0];
omega_rso = 2*pi*50;
polesrso_pll = [-1 -4]*omega_rso;
Lrso_pll = acker(Arso',Crso',polesrso_pll)';
Adrso_pll = eye(2) + Arso*ts_afe;
polesdrso_pll = exp(ts_afe*polesrso_pll);
Ldrso_pll = acker(Adrso_pll',Crso',polesdrso_pll)' %[output:3027a533]

%[text] ### PLL DDSRF
pll_i1 = 80;
pll_p = 1;
use_advanced_pll = 0;
use_dq_pll_ccaller = 0;
pll_i1_ddsrt = pll_i1/2;
pll_p_ddsrt = pll_p/2;
omega_f = 2*pi*50;
ddsrf_f = omega_f/(s+omega_f);
ddsrf_fd = c2d(ddsrf_f,ts_afe);
%%
%[text] ### First Harmonic Tracker for Ugrid cleaning
omega_fht0 = 2*pi*50;
delta_fht0 = 0.05;
Afht0 = [0 1; -omega_fht0^2 -delta_fht0*omega_fht0] % impianto nel continuo %[output:59ab4db0]
Cfht0 = [1 0];
poles_fht0 = [-1 -4]*omega_fht0;
Lfht0 = acker(Afht0',Cfht0', poles_fht0)' % guadagni osservatore nel continuo %[output:1bd527b7]
Ad_fht0 = eye(2) + Afht0*ts_afe % impianto nel discreto %[output:612fc384]
polesd_fht0 = exp(ts_afe*poles_fht0);
Ld_fht0 = acker(Ad_fht0',Cfht0', polesd_fht0) %[output:437f9927]

%[text] ### First Harmonic Tracker for Load
omega_fht1 = 2*pi*frequency_set;
delta_fht1 = 0.05;
Afht1 = [0 1; -omega_fht1^2 -delta_fht1*omega_fht1] % impianto nel continuo %[output:431f40f9]
Cfht1 = [1 0];
poles_fht1 = [-1 -4]*omega_fht1;
Lfht1 = acker(Afht1', Cfht1', poles_fht1)' % guadagni osservatore nel continuo %[output:250896a7]
Ad_fht1 = eye(2) + Afht1*ts_inv % impianto nel discreto %[output:35508714]
Ad_fht2 = eye(2) + Afht1*ts_inv + 1/2*Afht1^2*ts_inv^2 %[output:5c254f4d]
Ad_fht3 = expm(Afht1*ts_inv) %[output:545a801e]
polesd_fht1 = exp(ts_inv*poles_fht1);
Ld_fht1 = acker(Ad_fht1',Cfht1', polesd_fht1) %[output:91d31729]
Ld_fht2 = acker(Ad_fht2',Cfht1', polesd_fht1) %[output:0c394746]
Ld_fht3 = acker(Ad_fht3',Cfht1', polesd_fht1) %[output:2df064f3]
%%
%[text] ## Settings for user functions: filters, moving average, rms
%[text] ### Low Pass Filters
%[text] #### LPF 50Hz in state space (for initialization)
fcut = 50;
fof = 1/(s/(2*pi*fcut)+1);
[nfof, dfof] = tfdata(fof,'v');
[nfofd, dfofd]=tfdata(c2d(fof,ts_afe),'v');
fof_z = tf(nfofd,dfofd,ts_afe,'Variable','z');
[A,B,C,D] = tf2ss(nfofd,dfofd);
LVRT_flt_ss = ss(A,B,C,D,ts_afe);
[A,B,C,D] = tf2ss(nfof,dfof);
LVRT_flt_ss_c = ss(A,B,C,D);
%[text] #### LPF 161Hz
fcut_161Hz_flt = 161;
g0_161Hz = fcut_161Hz_flt * ts_afe * 2*pi;
g1_161Hz = 1 - g0_161Hz;
%%
%[text] #### LPF 500Hz
fcut_500Hz_flt = 500;
g0_500Hz = fcut_500Hz_flt * ts_afe * 2*pi;
g1_500Hz = 1 - g0_500Hz;
%%
%[text] #### LPF 75Hz
fcut_75Hz_flt = 75;
g0_75Hz = fcut_75Hz_flt * ts_afe * 2*pi;
g1_75Hz = 1 - g0_75Hz;
%%
%[text] #### LPF 50Hz
fcut_50Hz_flt = 50;
g0_50Hz = fcut_50Hz_flt * ts_afe * 2*pi;
g1_50Hz = 1 - g0_50Hz;
%%
%[text] #### LPF 10Hz
fcut_10Hz_flt = 10;
g0_10Hz = fcut_10Hz_flt * ts_afe * 2*pi;
g1_10Hz = 1 - g0_10Hz;
%%
%[text] #### LPF 4Hz
fcut_4Hz_flt = 4;
g0_4Hz = fcut_4Hz_flt * ts_afe * 2*pi;
g1_4Hz = 1 - g0_4Hz;
%%
%[text] #### LPF 1Hz
fcut_1Hz_flt = 1;
g0_1Hz = fcut_1Hz_flt * ts_afe * 2*pi;
g1_1Hz = 1 - g0_1Hz;
%%
%[text] #### LPF 0.2Hz
fcut_0Hz2_flt = 0.2;
g0_0Hz2 = fcut_0Hz2_flt * ts_afe * 2*pi;
g1_0Hz2 = 1 - g0_0Hz2;
%[text] #### 
%[text] ### Settings for RMS calculus
f_grid = 50;
rms_perios = 1;
n1 = rms_perios/f_grid/ts_afe;
rms_perios = 10;
n10 = rms_perios/f_grid/ts_afe;
%%
%[text] ### Online time domain sequence calculator
w_grid = 2*pi*f_grid;
apf = (s/w_grid-1)/(s/w_grid+1);
[napfd, dapfd]=tfdata(c2d(apf,ts_afe),'v');
apf_z = tf(napfd,dapfd,ts_afe,'Variable','z');
[A,B,C,D] = tf2ss(napfd,dapfd);
ap_flt_ss = ss(A,B,C,D,ts_afe);
% figure;
% bode(ap_flt_ss,options);
% grid on
%%
%[text] ### Single phase pll
freq_pll = 50;
kp_pll = 314;
ki_pll = 3140;

Arso = [0 1; 0 0];
Crso = [1 0];

polesrso = [-5 -1]*2*pi*10;
Lrso = acker(Arso',Crso',polesrso)';

Adrso = eye(2) + Arso*ts_inv;
polesdrso = exp(ts_inv*polesrso);
Ldrso = acker(Adrso',Crso',polesdrso)' %[output:7437557a]

freq_filter = freq_pll;
tau_f = 1/2/pi/freq_filter;
Hs = 1/(s*tau_f+1);
Hd = c2d(Hs,ts_inv);
%[text] #### Linear double integrator observer
Aso = [0 1; 0 0];
Asod = eye(2)+Aso*ts_inv;
Cso = [1 0];
omega_rso = 2*pi*frequency_set;
p2place = [-1 -10]*omega_rso;
p2placed = exp(p2place*ts_inv);
Kd = (acker(Asod',Cso',p2placed))';
kv = Kd(2)/ts_inv;
kx = Kd(1)/ts_inv;
%[text] ## C-Caller Settings
open_system(model_name);
Simulink.importExternalCTypes(model_name,'Names',{'first_harmonic_tracker_output_t'});
Simulink.importExternalCTypes(model_name,'Names',{'sogi_flt_output_t'});
%[text] ## Remove Scopes Opening Automatically
% open_scopes = find_system(model, 'BlockType', 'Scope');
% for i = 1:length(open_scopes)
%     set_param(open_scopes{i}, 'Open', 'off');
% end

% shh = get(0,'ShowHiddenHandles');
% set(0,'ShowHiddenHandles','On');
% hscope = findobj(0,'Type','Figure','Tag','SIMULINK_SIMSCOPE_FIGURE');
% close(hscope);
% set(0,'ShowHiddenHandles',shh);


%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":31.9}
%---
%[output:2c389496]
%   data: {"dataType":"matrix","outputData":{"columns":2,"exponent":"4","name":"Ares_nom","rows":2,"type":"double","value":[["0","0.0001"],["-9.8696","-0.0031"]]}}
%---
%[output:75942329]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"Aresd_nom","rows":2,"type":"double","value":[["1.0000","0.0001"],["-9.8696","0.9969"]]}}
%---
%[output:80bbc6d8]
%   data: {"dataType":"textualVariable","outputData":{"name":"a11d","value":"1"}}
%---
%[output:2928378c]
%   data: {"dataType":"textualVariable","outputData":{"name":"a12d","value":"1.0000e-04"}}
%---
%[output:68004a6c]
%   data: {"dataType":"textualVariable","outputData":{"name":"a21d","value":"-9.8696"}}
%---
%[output:485bb164]
%   data: {"dataType":"textualVariable","outputData":{"name":"a22d","value":"0.9969"}}
%---
%[output:3027a533]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"Ldrso_pll","rows":2,"type":"double","value":[["0.1490"],["36.5219"]]}}
%---
%[output:59ab4db0]
%   data: {"dataType":"matrix","outputData":{"columns":2,"exponent":"4","name":"Afht0","rows":2,"type":"double","value":[["0","0.0001"],["-9.8696","-0.0016"]]}}
%---
%[output:1bd527b7]
%   data: {"dataType":"matrix","outputData":{"columns":1,"exponent":"5","name":"Lfht0","rows":2,"type":"double","value":[["0.0156"],["2.7166"]]}}
%---
%[output:612fc384]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"Ad_fht0","rows":2,"type":"double","value":[["1.0000","0.0001"],["-9.8696","0.9984"]]}}
%---
%[output:437f9927]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"Ld_fht0","rows":1,"type":"double","value":[["0.1474","24.3363"]]}}
%---
%[output:431f40f9]
%   data: {"dataType":"matrix","outputData":{"columns":2,"exponent":"4","name":"Afht1","rows":2,"type":"double","value":[["0","0.0001"],["-9.8696","-0.0016"]]}}
%---
%[output:250896a7]
%   data: {"dataType":"matrix","outputData":{"columns":1,"exponent":"5","name":"Lfht1","rows":2,"type":"double","value":[["0.0156"],["2.7166"]]}}
%---
%[output:35508714]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"Ad_fht1","rows":2,"type":"double","value":[["1.0000","0.0001"],["-9.8696","0.9984"]]}}
%---
%[output:5c254f4d]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"Ad_fht2","rows":2,"type":"double","value":[["0.9995","0.0001"],["-9.8619","0.9979"]]}}
%---
%[output:545a801e]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"Ad_fht3","rows":2,"type":"double","value":[["0.9995","0.0001"],["-9.8602","0.9979"]]}}
%---
%[output:91d31729]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"Ld_fht1","rows":1,"type":"double","value":[["0.1474","24.3363"]]}}
%---
%[output:0c394746]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"Ld_fht2","rows":1,"type":"double","value":[["0.1465","23.6547"]]}}
%---
%[output:2df064f3]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"Ld_fht3","rows":1,"type":"double","value":[["0.1465","23.6626"]]}}
%---
%[output:7437557a]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"Ldrso","rows":2,"type":"double","value":[["0.0372"],["1.9371"]]}}
%---
