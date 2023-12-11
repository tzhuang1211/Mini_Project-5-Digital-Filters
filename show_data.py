import matplotlib.pyplot as plt
import numpy as np

def plot_impulse_response(filename, title):
    h = np.loadtxt(filename)
    plt.figure()
    plt.stem(h, basefmt='b-', linefmt='', markerfmt='', use_line_collection=True)
    plt.title(title)
    plt.xlabel('Sample Index')
    plt.ylabel('Amplitude')
    plt.grid(True)
    plt.show()

def plot_log_spectrum(filename, title, sample_rate=48000):
    log_spectrum = np.loadtxt(filename)
    num_samples = len(log_spectrum)
    frequency_axis = np.fft.fftfreq(num_samples, d=1/sample_rate)
    
    plt.figure()
    plt.plot(frequency_axis[:num_samples // 2], log_spectrum[:num_samples // 2])
    plt.title(title)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Log Spectrum (dB)')
    plt.grid(True)
    plt.show()

def main():
    # Replace these filenames with your actual filenames
    hL_filename = 'hL.txt'
    hR_filename = 'hR.txt'
    YL_filename = 'YL.txt'
    YR_filename = 'YR.txt'

    plot_impulse_response(hL_filename, 'Impulse Response (Left Channel)')
    plot_impulse_response(hR_filename, 'Impulse Response (Right Channel)')

    plot_log_spectrum(YL_filename, 'Log Spectrum (Left Channel)')
    plot_log_spectrum(YR_filename, 'Log Spectrum (Right Channel)')

if __name__ == "__main__":
    main()
