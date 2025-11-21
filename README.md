## Design-of-FIR-Filters-using-hanning-window


# AIM: 
          
  To generate design of low pass FIR digital filter using SCILAB 

# APPARATUS REQUIRED: 

  PC Installed with SCILAB 

# PROGRAM
## LOW PASS FILTER
```
/ FIR Low-Pass Filter using Fourier Series
// Hanning window, cutoff = 0.5π radians/sample, N = 21

clear;
clc;

// Filter specs
fs = 1;              // Normalized sampling frequency
fc = 0.25;           // Normalized cutoff (cycles/sample)
N = 21;              // Filter length (odd)

// Time index centered around zero
n = -(N-1)/2 : (N-1)/2;

// Ideal impulse response (sinc)
h = zeros(1, N);
for i = 1:N
    if n(i) == 0 then
        h(i) = 2 * fc;
    else
        h(i) = sin(2 * %pi * fc * n(i)) / (%pi * n(i));
    end
end

// Apply Hanning window
w = 0.5 - 0.5 * cos(2 * %pi * (0:N-1) / (N-1));
h_win = h .* w;

// Display coefficients
disp('FIR LPF Coefficients (Hanning Window):');
disp(h_win);

// Frequency response
scf(0);
clf();

[Hf, fr] = frmag(h_win, 512);
subplot(2,1,1);
plot(fr * fs, Hf);
xtitle('Magnitude Response (Hanning Window)', 'Frequency (Hz)', 'Magnitude');

// Phase response
w = 2 * %pi * fr;
H_sym = poly(h_win, 'z');
H_eval = horner(H_sym, exp(-%i * w));
phase_rad = atan(imag(H_eval) ./ real(H_eval));
phase_deg = phase_rad * 180 / %pi;

subplot(2,1,2);
plot(fr * fs, phase_deg);
xtitle('Phase Response (Hanning Window)', 'Frequency (Hz)', 'Phase (degrees)');

// Filter info
disp('FIR LPF Design Complete');
disp('Cutoff Frequency: 0.5π radians/sample');
disp('Filter Length: ' + string(N));
```
## HIGH PASS FILTER

```
clear;
clc;

// Filter specs
fs = 1;              // Normalized sampling frequency
fc = 0.25;           // Normalized cutoff (cycles/sample)
N = 21;              // Filter length (odd)

// Time index centered around zero
n = -(N-1)/2 : (N-1)/2;

// Ideal LPF impulse response (sinc)
h_lp = zeros(1, N);
for i = 1:N
    if n(i) == 0 then
        h_lp(i) = 2 * fc;
    else
        h_lp(i) = sin(2 * %pi * fc * n(i)) / (%pi * n(i));
    end
end

// Spectral inversion to get HPF
h_hp = -h_lp;
h_hp((N+1)/2) = 1 - h_lp((N+1)/2);  // center tap adjustment

// Apply Hanning window
w = 0.5 - 0.5 * cos(2 * %pi * (0:N-1) / (N-1));
h_win = h_hp .* w;

// Display coefficients
disp('FIR HPF Coefficients (Hanning Window):');
disp(h_win);

// Frequency response
scf(0);
clf();

[Hf, fr] = frmag(h_win, 512);
subplot(2,1,1);
plot(fr * fs, Hf);
xtitle('Magnitude Response (Hanning Window)', 'Frequency (Hz)', 'Magnitude');

// Phase response
w = 2 * %pi * fr;
H_sym = poly(h_win, 'z');
H_eval = horner(H_sym, exp(-%i * w));
phase_rad = atan(imag(H_eval) ./ real(H_eval));
phase_deg = phase_rad * 180 / %pi;

subplot(2,1,2);
plot(fr * fs, phase_deg);
xtitle('Phase Response (Hanning Window)', 'Frequency (Hz)', 'Phase (degrees)');

// Filter info
disp('FIR HPF Design Complete');
disp('Cutoff Frequency: 0.5π radians/sample');
disp('Filter Length: ' + string(N));
```
# OUTPUT
 ## LOW PASS FILTER
<img width="1920" height="1080" alt="image" src="https://github.com/user-attachments/assets/c2924474-a69f-4ef7-9f0d-cc2921a82dbc" />

<img width="1920" height="1080" alt="image" src="https://github.com/user-attachments/assets/2540803a-2fa7-4a2d-94a9-fb9f32cc6827" />

## HIGH PASS FILTER

<img width="1920" height="1080" alt="Screenshot (379)" src="https://github.com/user-attachments/assets/4fd9baee-5c0c-4830-99f8-683a7346b5a5" />

<img width="1920" height="1080" alt="Screenshot (378)" src="https://github.com/user-attachments/assets/022c85dd-ff7a-4079-b713-e16bdb7b3994" />






# RESULT
Thus FIR Low pass filter design using Hanning window is performed.
