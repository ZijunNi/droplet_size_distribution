re = [200 400 800 1600];
uinf = [18.534401289144565 19.978661 21.711529  23.444397];

re_dns = [180 1000];
uinf_dns = [18.2710 22.638488632724943];
semilogx(re,uinf,'-x');
hold on
semilogx(re_dns,uinf_dns,'-o')

