clear, close all


nrf1 = load("results_nrf1.mat");
nrf2 = load("results_nrf2.mat");
nrf4 = load("results_nrf4.mat");
nrf8 = load("results_nrf8.mat");


figure
hold on
semilogy(nrf1.snrdb,nrf1.rmse_gls_deg,'--og',DisplayName='GLS Nrf=1')
semilogy(nrf1.snrdb,nrf1.rcrb_hyb_deg,'-g',DisplayName='RCRB Nrf=1')

semilogy(nrf2.snrdb,nrf2.rmse_gls_deg,'--ob',DisplayName='GLS Nrf=2')
semilogy(nrf2.snrdb,nrf2.rcrb_hyb_deg,'-b',DisplayName='RCRB Nrf=2')

semilogy(nrf4.snrdb,nrf4.rmse_gls_deg,'--or',DisplayName='GLS Nrf=4')
semilogy(nrf4.snrdb,nrf4.rcrb_hyb_deg,'-r',DisplayName='RCRB Nrf=4')

semilogy(nrf8.snrdb,nrf8.rmse_gls_deg,'--ok',DisplayName='GLS Nrf=8')
semilogy(nrf8.snrdb,nrf8.rcrb_hyb_deg,'-k',DisplayName='DIGITAL')
hold off
set(gca,'YScale','log')

legend()
grid on
xlabel('snr [dB]')
ylabel('RMSE [deg]')
title('N=8, K=1575, MC=1000, theta=\{-45,45\}')