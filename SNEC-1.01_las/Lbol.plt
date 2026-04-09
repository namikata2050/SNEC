# magnetar model
#plot "Data_Lsd1.6d43_tsd_12_kout200_no_Ni_Eexp2.2foe/lum_observed.dat" u (($1)/86400.):2
#replot "SN2015ap_Amar_digitalized.txt" u (($1)-5.2):(($2)*10**(42))


# Ni model
plot "Data_Ni_0.135_no_magnetar_Eexp3.0foe_Ni_boundary_3.0/lum_observed.dat" u (($1)/86400.):2 w l lw 4
replot "Data_Ni_0.135_no_magnetar_Eexp3.0foe_Ni_boundary_3.0_imax1000_kout400/lum_observed.dat" u (($1)/86400.):2 w l lw 2
replot "SN2015ap_Amar_digitalized.txt" u (($1)-9.7):(($2)*10**(42))