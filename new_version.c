#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>

#define FS 48000.0f
#define FL 1500.0f
#define FH 3500.0f
#define PI 3.141592653589793f
#define N 1024

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
        return 1.0 - (wh / PI - wl / PI);
    }
    else {
        return -(sinf(wl * ((float)(n - m))) - sinf(wh * ((float)(n - m)))) / PI / ((float)(n - m)) * hamming(2 * m + 1, n);
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

void DFT_Calculate(float *amplitude, short *channel, size_t length, int startSample, int numSamples, int N) {
    int i, k;
    for (i = 0; i < N; i++) {
        float realPart = 0;
        float imagPart = 0;
        for (k = 0; k < length; k++) {
            int sampleIndex = startSample + numSamples - k;
            if (sampleIndex >= 0 && sampleIndex < length) {
                realPart += channel[sampleIndex] * cos(2 * PI * i * k / N);
                imagPart -= channel[sampleIndex] * sin(2 * PI * i * k / N);
            }
        }
        amplitude[i] = sqrt(realPart * realPart + imagPart * imagPart);
    }
}

// Function to calculate log spectrum using DFT
void calculateLogSpectrum(float *logSpectrum, float *amplitude, size_t length, int startSample, int numSamples, int M) {
    int n;
    int i;

    for (n = 0; n < numSamples; n++) {
        float sum = 0;

        // Calculate the power spectrum from the squared magnitude of the DFT coefficients
        for (i = 0; i < N; i++) {
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
    // Command line argument check
    if (argc != 8) {
        fprintf(stderr, "Usage: %s M hL.txt hR.txt YL.txt YR.txt input.wav output.wav\n", argv[0]);
        return 1;
    }

    // Parse command line arguments
    int M = atoi(argv[1]);
    char *fn_hL = argv[2];
    char *fn_hR = argv[3];
    char *fn_YL = argv[4];
    char *fn_YR = argv[5];
    char *fn_in = argv[6];
    char *fn_out = argv[7];

    wav wavin;
    wav wavout;

    // Allocate memory for filter coefficients
    float *h_L = malloc(sizeof(float) * (2 * M + 1));
    float *h_R = malloc(sizeof(float) * (2 * M + 1));

    // Memory allocation check
    if (h_L == NULL || h_R == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    // Validate filter order
    if (M <= 0) {
        fprintf(stderr, "Invalid value for M\n");
        return 1;
    }

    // Initialize filter coefficients arrays to zero
    memset(h_L, 0, sizeof(float) * (2 * M + 1));
    memset(h_R, 0, sizeof(float) * (2 * M + 1));

    float YL[1200] = {0};
    float YR[1200] = {0};
    int n = 0;
    float y = 0;
    int k;

    // Read input WAV file
    if (wav_read_fn(fn_in, &wavin) == 0) {
        fprintf(stderr, "Cannot read wav file %s\n", fn_in);
        exit(1);
    }

    // Construct band-pass filter for left channel
    for (n = 0; n < (2 * M + 1); n++) {
        h_L[n] = band_pass(M, n);
    }

    // Construct band-stop filter for right channel
    for (n = 0; n < (2 * M + 1); n++) {
        h_R[n] = band_stop(M, n);
    }

    // Initialize output WAV structure
    if (wav_init(wavin.length, &wavout) == 0) {
        exit(1);
    }

    // Filtering (convolution) for left channel
    for (n = 0; n < wavin.length; n++) {
        y = 0;
        for (k = 0; k < (2 * M + 1); k++) {
            if ((n - k) >= 0)
                y = y + h_L[k] * ((float)(wavin.LChannel[n - k]));
        }
        wavout.LChannel[n] = (short)(roundf(y));
    }

    // Filtering (convolution) for right channel
    for (n = 0; n < wavin.length; n++) {
        y = 0;
        for (k = 0; k < (2 * M + 1); k++) {
            if ((n - k) >= 0)
                y = y + h_R[k] * ((float)(wavin.RChannel[n - k]));
        }
        wavout.RChannel[n] = (short)(roundf(y));
    }

    // Save filtered WAV
    memcpy(wavout.header, wavin.header, 44);
    if (wav_save_fn(fn_out, &wavout) == 0) {
        fprintf(stderr, "Cannot save %s\n", fn_out);
        exit(1);
    }

    // Save filter coefficients to files
    FILE *fp_hL = fopen(fn_hL, "w");
    FILE *fp_hR = fopen(fn_hR, "w");

    // File open check
    if (fp_hL == NULL || fp_hR == NULL) {
        fprintf(stderr, "Error opening filter coefficient files\n");
        return 1;
    }

    // Write filter coefficients to files
    for (n = 0; n < (2 * M + 1); n++) {
        fprintf(fp_hL, "%.15e\n", h_L[n]);
        fprintf(fp_hR, "%.15e\n", h_R[n]);
    }

    // Close filter coefficient files
    fclose(fp_hL);
    fclose(fp_hR);

    // Calculate DFT and amplitude using your DFT_Calculate function
    float amplitude_L[N];
    float amplitude_R[N];

    DFT_Calculate(amplitude_L, wavin.LChannel, wavin.length, 962880, 1200, N);
    DFT_Calculate(amplitude_R, wavin.RChannel, wavin.length, 962880, 1200, N);

    // Calculate and save log spectrum
    float logSpectrum_L[1200];
    float logSpectrum_R[1200];

    calculateLogSpectrum(logSpectrum_L, amplitude_L, wavin.length, 962880, 1200, M);
    calculateLogSpectrum(logSpectrum_R, amplitude_R, wavin.length, 962880, 1200, M);

    // Open log spectrum files
    FILE *fp_YL = fopen(fn_YL, "w");
    FILE *fp_YR = fopen(fn_YR, "w");

    // File open check
    if (fp_YL == NULL || fp_YR == NULL) {
        fprintf(stderr, "Error opening log spectrum files\n");
        return 1;
    }

    // Write log spectrum to files
    for (int n = 0; n < 1200; n++) {
        fprintf(fp_YL, "%.15e\n", logSpectrum_L[n]);
        fprintf(fp_YR, "%.15e\n", logSpectrum_R[n]);
    }

    // Close log spectrum files
    fclose(fp_YL);
    fclose(fp_YR);

    // Free memory
    wav_free(&wavin);
    wav_free(&wavout);

    return 0;
}
