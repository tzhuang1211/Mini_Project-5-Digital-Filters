# Mini_Project-5: Digital-Filters
### 411086026 通訊三 黃廷哲
### 程式碼解釋:
#### 引入函式庫
```js
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
```
#### 定義參數
```js
#define FS 48000.0f
#define FL 1500.0f
#define FH 3500.0f
#define PI 3.141592653589793f
```
#### 定義WAV檔結構
```js
// hold the information about a WAV file.
typedef struct _wav {
	int fs;
	char header[44];
	size_t length;
	short *LChannel;
	short *RChannel;
} wav;

//讀取 WAV 檔案的函數，將header和音訊資料讀取到 WAV structure中
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

//將資料寫入WAV file中
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

//初始化WAV檔
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

//釋放記憶體
void wav_free(wav *p_wav)
{
	free(p_wav->LChannel);
	free(p_wav->RChannel);
}
```

#### 實現Hamming window以及對Band-pass filter和Band-stop filter做設置
```js
float hamming(int N, int n) // N = window length
{
	return 0.54 - 0.46 * cosf(2*PI*((float)(n))/((float)N));
} //hamming window

float band_stop(int m, int n)
{
    float wh = 2 * PI * FH / FS; //高通的角頻率
    float wl = 2 * PI * FL / FS; //低通的角頻率
    if (n == m) {
        return (wh / PI - wl / PI); //返回高通和低通頻率之差
    }
    else {
        return (sinf(wl * ((float)(n - m))) - sinf(wh * ((float)(n - m)))) / PI / ((float)(n - m)) * hamming(2 * m + 1, n);
    } // 如果index不相等，計算帶阻濾波器的響應
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
	} // 如果index不相等，計算帶通濾波器的響應
}
```
#### 定義複數structure和DFT
```js
typedef struct {
    double real;
    double imag;
} complex_t;

// DFT function
void dft(complex_t* x, int N) {
    complex_t* X = malloc(N * sizeof(complex_t));

    for (int k = 0; k < N; ++k) {
        X[k] = (complex_t){0.0, 0.0}; //Initialize

        for (int n = 0; n < N; ++n) {
            double angle = -2.0 * PI * k * n / N; //計算相位角
            X[k].real += x[n].real * cos(angle) - x[n].imag * sin(angle);
            X[k].imag += x[n].real * sin(angle) + x[n].imag * cos(angle);
            //對X[k]不斷更新其虛部與實部(Euler's formula)，也就是經DFT後的結果
        }
    }

    //put result into an array
    for (int i = 0; i < N; ++i) {
        x[i] = X[i];
    }

    free(X);
}
```
#### 計算Log Spectrum and Coefficients並進行存檔
```js
void calculateLogSpectrumAndSave(float* log_spectrum, short* channel, size_t length, size_t start, size_t end, float* coefficients, char* filename, wav *p_wav, int isBandPass) {
    // 設 in 為一複數陣列
    complex_t *in = (complex_t *)malloc((end - start) * sizeof(complex_t));

    // 複製 in 到 DFT input
    for (int i = 0; i < (end - start); i++) {
        in[i] = (complex_t) {(double)(isBandPass ? p_wav->LChannel[start + i] : p_wav->RChannel[start + i]), 0.0};
    }

    // 使用 DFT (Time domain >> Frequency domain)
    dft(in, end - start);

    // 建檔並以 log_spectrum_file 指向此檔，並寫入檔案
    FILE *log_spectrum_file = fopen(filename, "w");
    if (log_spectrum_file == NULL) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        exit(1);
    }

    // 計算 log spectrum
    for (int i = 0; i < (end - start); i++) {
        double magnitude = sqrt(in[i].real * in[i].real + in[i].imag * in[i].imag);
        double log_magnitude = 20 * log10(fabs(magnitude) + 1.0);  // 避免 log(0)

        // 計算每一點的實際頻率
        double frequency = ((double)i / (end - start)) * FS;

        // 確認資料是否在我們所希望的頻率範圍內(1500Hz~3500Hz)並寫入檔案中
        if (frequency >= FL && frequency <= FH) {
            fprintf(log_spectrum_file, "%.15e ", log_magnitude);
        }
    }

    fclose(log_spectrum_file);
    free(in);
}

// 計算濾波器係數
void calculateFilterCoefficients(float *h, int filterType, int M) {
    int n;
    for (n = 0; n < (2 * M + 1); n++) {
        if (filterType == 0) {
            h[n] = band_pass(M, n);
        } else {
            h[n] = band_stop(M, n);
        }
    } // 根據不同聲道做不同濾波處理
}
```
#### 主程式
```js
int main(int argc, char **argv) {
    if (argc != 8) {
        fprintf(stderr, "Usage: %s M hL.txt hR.txt YL.txt YR.txt input.wav output.wav\n", argv[0]);
        return 1;
    }
    // Input設置
    int M = atoi(argv[1]); // 變數
    char *fn_hL = argv[2];
    char *fn_hR = argv[3];
    char *fn_YL = argv[4];
    char *fn_YR = argv[5];
    char *fn_in = argv[6];
    char *fn_out = argv[7];

    if (M <= 0) {
        fprintf(stderr, "Invalid value for M\n");
        return 1;
    } // 除錯

    // 分配記憶體給濾波器係數
    float *h_L = (float *)calloc(2 * M + 1, sizeof(float));
    float *h_R = (float *)calloc(2 * M + 1, sizeof(float));
    float y = 0;
    int k;

    if (h_L == NULL || h_R == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

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

    wav wavin; // 存取輸入音訊
    wav wavout; // 存取輸出音訊

    // 讀取 input wav
    if (wav_read_fn(fn_in, &wavin) == 0) {
        fprintf(stderr, "Cannot read wav file %s\n", fn_in);
        return 1;
    }

    // 初始化 output wav structure
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

    // 複製 header 資料
    memcpy(wavout.header, wavin.header, 44);

    // 儲存output
    if (wav_save_fn(fn_out, &wavout) == 0) {
        fprintf(stderr, "Cannot save %s\n", fn_out);
        return 1;
    }
    size_t start_time = 962880; // 20.060 seconds
    size_t end_time = 964080;   // 20.085 seconds

    float *Y_L = (float *)malloc((end_time - start_time) * sizeof(float));
    float *Y_R = (float *)malloc((end_time - start_time) * sizeof(float));

    if (Y_L == NULL || Y_R == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    // Compute and save log spectrum for the left channel
    calculateLogSpectrumAndSave(Y_L, wavin.LChannel, wavin.length, start_time, end_time, h_L, fn_YL, &wavin, 1);

    // Compute and save log spectrum for the right channel
    calculateLogSpectrumAndSave(Y_R, wavin.RChannel, wavin.length, start_time, end_time, h_R, fn_YR, &wavin, 0);

    // 釋放記憶體
    wav_free(&wavin);
    wav_free(&wavout);
    free(Y_L);
    free(Y_R);
    free(h_L);
    free(h_R);

    return 0;
}
```

### 程式結果圖
## M=8
![](https://github.com/Gengarkayak1333/Mini_Project-5-Digital-Filters/blob/main/IR%20LC%208.png "M=8 Impulse Responses of Left Channel")
![](https://github.com/Gengarkayak1333/Mini_Project-5-Digital-Filters/blob/main/IR%20RC%208.png "M=8 Impulse Responses of Right Channel")
![](https://github.com/Gengarkayak1333/Mini_Project-5-Digital-Filters/blob/main/LS%20RC%208.png "M=8 Log Spectrum of Left Channel")
![](https://github.com/Gengarkayak1333/Mini_Project-5-Digital-Filters/blob/main/LS%20RC%208.png "M=8 Log Spectrum of Right Channel")
## M=32
![](https://github.com/Gengarkayak1333/Mini_Project-5-Digital-Filters/blob/main/IR%20LC%2032.png "M=32 Impulse Responses of Left Channel")
![](https://github.com/Gengarkayak1333/Mini_Project-5-Digital-Filters/blob/main/IR%20RC%2032.png "M=32 Impulse Responses of Right Channel")
![](https://github.com/Gengarkayak1333/Mini_Project-5-Digital-Filters/blob/main/LS%20RC%2032.png "M=32 Log Spectrum of Left Channel")
![](https://github.com/Gengarkayak1333/Mini_Project-5-Digital-Filters/blob/main/LS%20RC%2032.png "M=32 Log Spectrum of Right Channel")
## M=1024
![](https://github.com/Gengarkayak1333/Mini_Project-5-Digital-Filters/blob/main/IR%20LC%201024.png "M=1024 Impulse Responses of Left Channel")
![](https://github.com/Gengarkayak1333/Mini_Project-5-Digital-Filters/blob/main/IR%20RC%201024.png  "M=1024 Impulse Responses of Right Channel")
![](https://github.com/Gengarkayak1333/Mini_Project-5-Digital-Filters/blob/main/LS%20LC%201024.png "M=1024 Log Spectrum of Left Channel")
![](https://github.com/Gengarkayak1333/Mini_Project-5-Digital-Filters/blob/main/LS%20RC%201024.png "M=1024 Log Spectrum of Right Channel")

### 問題討論
##### 1. M如何影響結果？

##### 2. log spectrum 顯示的濾波結果是否正確？如何證明？
在理想上，帶通濾波器

