
        visibilities = load('Visibilties_for_0_m_south_6_m_east_0_m_up_xx_pol_125.195_MHz.dat'); 

figure(1); clf
set(1,'Position',[1816         852         743         703]);
ha = tight_subplot(2,1,[.1 .1],[.05 .05],[.05 .05]);
    

axes(ha(1));
plot(visibilities(:,1),visibilities(:,2),'.-')
xlabel('UTC Time (Hours)');
ylabel('Re[Visibility]')

axes(ha(2));
plot(visibilities(:,1),visibilities(:,3),'.-')
xlabel('UTC Time (Hours)');
ylabel('Im[Visibility]')

