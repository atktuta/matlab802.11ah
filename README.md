# matlab802.11ah
Tujuannya adalah membuat analisis dan simulasi PER vs Distance di 802.11ah.

Modulasi yang dipakai hanya BPSK.

Landasan pertama adalah dari artikel Bojan Domazetović dkk.[1]. Reproduce untuk figure 1, 2 dan 3.

Paper berikutnya adalah dari artikel Ferrand dkk. [2]. Reproduce untuk figure 6.

Mencoba membuat simulasi dan analisis berdasar [3] untuk menghasilkan BER vs EbN0 pada OFDM BPSK di channel AWGN dan Rayleigh.

Penjelasan file <b>perbandingan_rumus_BER.m</b>:<br/>

Untuk BPSK sudah ada rumus untuk menghitung BER berdasar EbN0 maupun SNR.
Dan rumus ini berlaku untuk OFDM juga.
<br/>
BER AWGN = 0.5 * erfc(sqrt(EbN0))<br/>
BER AWGN = qfunc(sqrt(2EbN0))<br/>
BER Rayleigh = 0.5 * (1-sqrt(EbN0/(EbN0+1))) <br/> 
BER Rayleigh = berfading(Eb_N0_dB, 'psk', 2, 1) -> fungsi matlab <br/>
<br/>
BER AWGN = qfunc(sqrt(SNR))<br/>
BER AWGN = 0.5*erfc(sqrt(SNR/2))<br/>

Curve sudah dibandingkan dengan program-program matlab yang lain, dan hasilnya sama.

Apakah ada yang tahu formula untuk BER Rayleigh sebagai fungsi dari SNR? Mohon kirim ke email yang ada di profile page.
Kami menggunakan:

BER Rayleigh = 0.5 * (1-sqrt(SNR/(SNR+1))) <br/>

dimana hasilnya sama dengan hasil simulasi.

PER diperoleh dari rumus:<br/>
PER = 1 - (1 - BER)^N

Dimana N adalah jumlah simbol per paket(?). Masih belum jelas definisi N.


<br/>
References:<br/>
<br/>
[1] B. Domazetovic, E. Kocan, and A. Mihovska, “Performance evaluation of IEEE 802.11ah systems,” presented at the 2016 24th Telecommunications Forum (TELFOR), Belgrade, Serbia, 2016, pp. 1–4.<br/>
[2] P. Ferrand, J.-M. Gorce, and C. Goursaud, “Approximations of the packet error rate under slow fading in direct and relayed links,” RESEARCH REPORT N° 8316, INRIA, p. 23, 2013.<br/>
[3] http://www.dsplog.com/2008/06/10/ofdm-bpsk-bit-error/<br/>
