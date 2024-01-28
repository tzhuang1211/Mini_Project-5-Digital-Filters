#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>

#define FS 48000.0f
#define FL 1500.0f
#define FH 3500.0f
#define PI 3.141592653589793f

// hold the information about a WAV file.
typedef struct _wav {
	int fs;
	char header[44];
	size_t length;
	short *LChannel;
	short *RChannel;
} wav;

int wav_read_fn(char *fn, wav *p_wav) // this function reads a WAV file specified by 'fn' and populates the 'wav' structure pointed to by 'p_wav' with the file's information.
{
	//char header[44];
	short temp = 0;
	size_t i = 0;

	FILE *fp = fopen(fn, "rb");
	if(fp==NULL) {
		fprintf(stderr, "cannot read %s\n", fn);
		return 0;
	}
	fread(p_wav->header, sizeof(char), 44, fp);
	while( !feof(fp) ) {
		fread(&temp, sizeof(short), 1, fp);
		i++;
	}
	p_wav->length = i / 2;
	p_wav->LChannel = (short *) calloc(p_wav->length, sizeof(short));
	if( p_wav->LChannel==NULL ) {
		fprintf(stderr, "cannot allocate memory for LChannel in wav_read_fn\n");
		fclose(fp);
		return 0;
	}
	p_wav->RChannel = (short *) calloc(p_wav->length, sizeof(short));
	if( p_wav->RChannel==NULL ) {
		fprintf(stderr, "cannot allocate memory for RChannel in wav_read_fn\n");
		fclose(fp);
		return 0;
	}
	fseek(fp, 44, SEEK_SET);
	for(i=0;i<p_wav->length;i++) {
		fread(p_wav->LChannel+i, sizeof(short), 1, fp);
		fread(p_wav->RChannel+i, sizeof(short), 1, fp);
	}
	fclose(fp);
	return 1;
}

int wav_save_fn(char *fn, wav *p_wav) // this function saves a WAV file specified by 'fn' using the information in the 'wav' structure pointed to by 'p_wav'.
{
	FILE *fp = fopen(fn, "wb");
	size_t i;
	if(fp==NULL) {
		fprintf(stderr, "cannot save %s\n", fn);
		return 0;
	}
	fwrite(p_wav->header, sizeof(char), 44, fp);
	for(i=0;i<p_wav->length;i++) {
		fwrite(p_wav->LChannel+i, sizeof(short), 1, fp);
		fwrite(p_wav->RChannel+i, sizeof(short), 1, fp);
	}
	fclose(fp);
	return 1;
}

int wav_init(size_t length, wav *p_wav) // initialize and free memory for a 'wav' structure.
{
	p_wav->length = length;
	p_wav->LChannel = (short *) calloc(p_wav->length, sizeof(short));
	if( p_wav->LChannel==NULL ) {
		fprintf(stderr, "cannot allocate memory for LChannel in wav_read_fn\n");
		return 0;
	}
	p_wav->RChannel = (short *) calloc(p_wav->length, sizeof(short));
	if( p_wav->RChannel==NULL ) {
		fprintf(stderr, "cannot allocate memory for RChannel in wav_read_fn\n");
		return 0;
	}
	return 1;
}

void wav_free(wav *p_wav)
{
	free(p_wav->LChannel);
	free(p_wav->RChannel);
}

float hamming(int N, int n) // N = window length
{
	return 0.54 - 0.46 * cosf(2*PI*((float)(n))/((float)N));
}

float band_stop(int m, int n)
{
    float wh = 2 * PI * FH / FS;
    float wl = 2 * PI * FL / FS;
    if (n == m) {
        return (wh / PI - wl / PI);
    }
    else {
        return (sinf(wl * ((float)(n - m))) - sinf(wh * ((float)(n - m)))) / PI / ((float)(n - m)) * hamming(2 * m + 1, n);
    }
}

float band_pass(int m, int n)
{
	float wh = 2*PI*FH/FS;
    float wl = 2*PI*FL/FS;
	if(n==m) {// L'Hopital's Rule
		return 2.0*(wh/PI - wl/PI);
	}
	else {
		return 2.0*(sinf(wh*((float)(n-m)))-sinf(wl*((float)(n-m))))/PI/((float)(n-m)) * hamming(2*m+1, n);
	}
}

typedef struct {
    double real;
    double imag;
} complex_t;

// DFT function
void dft(complex_t* x, int N) {
    complex_t* X = malloc(N * sizeof(complex_t));

    for (int k = 0; k < N; ++k) {
        X[k] = (complex_t){0.0, 0.0};

        for (int n = 0; n < N; ++n) {
            double angle = -2.0 * PI * k * n / N;
            X[k].real += x[n].real * cos(angle) - x[n].imag * sin(angle);
            X[k].imag += x[n].real * sin(angle) + x[n].imag * cos(angle);
        }
    }

    //put result into an array
    for (int i = 0; i < N; ++i) {
        x[i] = X[i];
    }

    free(X);
}

void calculateLogSpectrumAndSave(float* log_spectrum, short* channel, size_t length, size_t start, size_t end, float* coefficients, char* filename, wav *p_wav, int isBandPass) {
    // Assuming in is an array of complex numbers
    complex_t *in = (complex_t *)malloc((end - start) * sizeof(complex_t));

    // copy input to DFT input (seperate right channel from left channel)
    for (int i = 0; i < (end - start); i++) {
        in[i] = (complex_t) {(double)(isBandPass ? p_wav->LChannel[start + i] : p_wav->RChannel[start + i]), 0.0};
    }

    // use DFT
    dft(in, end - start);

    // write data into log_spectrum_file
    FILE *log_spectrum_file = fopen(filename, "w");
    if (log_spectrum_file == NULL) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        exit(1);
    } // debug

    // calculate log spectrum
    for (int i = 0; i < (end - start); i++) {
        double magnitude = sqrt(in[i].real * in[i].real + in[i].imag * in[i].imag);
        double log_magnitude = 20 * log10(fabs(magnitude) + 1.0);  // Avoid log(0) situation

        // Use the frequency information adjusted for the window size and sampling rate
        double frequency = ((double)i / (end - start)) * FS;

        fprintf(log_spectrum_file, "%.15e ", log_magnitude); // write data in
    }

    fclose(log_spectrum_file);
    free(in);
}


// Function to calculate filter coefficients
void calculateFilterCoefficients(float *h, int filterType, int M) {
    int n;
    for (n = 0; n < (2 * M + 1); n++) {
        if (filterType == 0) {
            h[n] = band_pass(M, n);
        } else {
            h[n] = band_stop(M, n);
        }
    }
}

// Function to calculate log spectrum using amplitude
void calculateLogSpectrum(float *logSpectrum, float *amplitude, size_t length, int startSample, int numSamples, int M) {
    int n;
    int i;

    for (n = 0; n < numSamples; n++) {
        float sum = 0;

        // Calculate the power spectrum from the squared magnitude of the amplitude
        for (i = 0; i < length; i++) {
            int sampleIndex = startSample + n - i;
            if (sampleIndex >= 0 && sampleIndex < length) {
                sum += amplitude[i] * amplitude[i];
            }
        }

        // Calculate the log spectrum
        logSpectrum[n] = 10 * log10(sum);
    }
}

int main(int argc, char **argv) {
    if (argc != 8) {
        fprintf(stderr, "Usage: %s M hL.txt hR.txt YL.txt YR.txt input.wav output.wav\n", argv[0]);
        return 1;
    }

    int M = atoi(argv[1]);
    char *fn_hL = argv[2]; // input arguments
    char *fn_hR = argv[3];
    char *fn_YL = argv[4];
    char *fn_YR = argv[5];
    char *fn_in = argv[6];
    char *fn_out = argv[7];

    if (M <= 0) {
        fprintf(stderr, "Invalid value for M\n");
        return 1;
    } // debug

    // Allocate memory for filter coefficients
    float *h_L = (float *)calloc(2 * M + 1, sizeof(float));
    float *h_R = (float *)calloc(2 * M + 1, sizeof(float));
    float y = 0;
    int k;

    if (h_L == NULL || h_R == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    // Initialize filter coefficients
    calculateFilterCoefficients(h_L, 0, M);  // band-pass for h_L
    calculateFilterCoefficients(h_R, 1, M);  // band-stop for h_R

    // Write h_L to file
    FILE *hL_file = fopen(fn_hL, "w");
    if (hL_file == NULL) {
        fprintf(stderr, "Error opening file: %s\n", fn_hL);
        return 1;
    }
    for (int i = 0; i < 2 * M + 1; i++) {
        fprintf(hL_file, "%.15e ", h_L[i]);
    }
    fclose(hL_file);

    // Write h_R to file
    FILE *hR_file = fopen(fn_hR, "w");
    if (hR_file == NULL) {
        fprintf(stderr, "Error opening file: %s\n", fn_hR);
        return 1;
    }
    for (int i = 0; i < 2 * M + 1; i++) {
        fprintf(hR_file, "%.15e ", h_R[i]);
    }
    fclose(hR_file);


    wav wavin; // populate structure for input.wav
    wav wavout; // populate structure for output.wav

    // Read input wav file
    if (wav_read_fn(fn_in, &wavin) == 0) {
        fprintf(stderr, "Cannot read wav file %s\n", fn_in);
        return 1;
    }

    // Initialize output wav structure
    if (wav_init(wavin.length, &wavout) == 0) {
        return 1;
    }

    // Filter the input wav
    for (int n = 0; n < wavin.length; n++) {
        y = 0;
        for (k = 0; k < (2 * M + 1); k++) {
            if ((n - k) >= 0)
                y = y + h_L[k] * ((float)(wavin.LChannel[n - k]));
        }
        wavout.LChannel[n] = (short)(roundf(y));

        y = 0;
        for (k = 0; k < (2 * M + 1); k++) {
            if ((n - k) >= 0)
                y = y + h_R[k] * ((float)(wavin.RChannel[n - k]));
        }
        wavout.RChannel[n] = (short)(roundf(y));
    }

    // Copy header information
    memcpy(wavout.header, wavin.header, 44);

    // Save the output wav file
    if (wav_save_fn(fn_out, &wavout) == 0) {
        fprintf(stderr, "Cannot save %s\n", fn_out);
        return 1;
    }

    size_t start_time = 962880; // 20.060 seconds
    size_t end_time = 964080;   // 20.085 seconds

    float *Y_L = (float *)malloc((end_time - start_time) * sizeof(float)); // allocate memory for Y_L 
    float *Y_R = (float *)malloc((end_time - start_time) * sizeof(float));

    if (Y_L == NULL || Y_R == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    // Compute and save log spectrum for the left channel
    calculateLogSpectrumAndSave(Y_L, wavin.LChannel, wavin.length, start_time, end_time, h_L, fn_YL, &wavout, 1);
    // Compute and save log spectrum for the right channel
    calculateLogSpectrumAndSave(Y_R, wavin.RChannel, wavin.length, start_time, end_time, h_R, fn_YR, &wavout, 0);

    // Free allocated memory
    wav_free(&wavin);
    wav_free(&wavout);
    free(Y_L);
    free(Y_R);
    free(h_L);
    free(h_R);

    return 0;
}
