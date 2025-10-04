%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

f = 8.0; % frequency in GHz

mg = 0.0: 0.02: 0.7; 

nf = length(mg);

for ii = 1: nf
[epsr1(ii), epsi1(ii)] = RelDielConst_Vegetation(f,mg(ii));
end
figure(1)
subplot(2,1,1)
plot(mg, epsr1)
xlabel('m_g')
ylabel('\epsilon_r')
grid

subplot(2,1,2)
plot(mg,epsi1)
xlabel('m_g')
ylabel('\epsilon_i')
grid



% f = 1:.1:20; % frequency in GHz
% 
% mg = 0.07; 
% 
% nf = length(f);
% 
% for ii = 1: nf
% [epsr1(ii), epsi1(ii)] = RelDielConst_Vegetation(f(ii),mg);
% end
% figure(1)
% subplot(2,1,1)
% plot(f, epsr1)
% xlabel('frequency (GHz)')
% ylabel('\epsilon_r')
% grid
% 
% subplot(2,1,2)
% plot(f,epsi1)
% xlabel('frequency (GHz)')
% ylabel('\epsilon_i')
% grid

