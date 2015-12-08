package fft.fft;


public class Fourier {
    public static double[] discreteFT(double[] fdata, int N, boolean fwd) {
        double X[] = new double[2 * N];
        double omega;
        int k, ki, kr, n;
        if (fwd) {
            omega = 2.0 * Math.PI / N;
        } else {
            omega = -2.0 * Math.PI / N;
        }
        for (k = 0; k < N; k++) {
            kr = 2 * k;
            ki = 2 * k + 1;
            X[kr] = 0.0;
            X[ki] = 0.0;
            for (n = 0; n < N; ++n) {
                X[kr] += fdata[2 * n] * Math.cos(omega * n * k) +
                        fdata[2 * n + 1] * Math.sin(omega * n * k);
                X[ki] += -fdata[2 * n] * Math.sin(omega * n * k) +
                        fdata[2 * n + 1] * Math.cos(omega * n * k);
            }
        }
        if (fwd) {
            for (k = 0; k < N; ++k) {
                X[2 * k] /= N;
                X[2 * k + 1] /= N;
            }
        }
        return X;
    }

    public static void fastFFT(double[] fdata, int N, boolean fwd) {
        double omega, tempr, tempi, fscale;
        double xtemp, cos, sin, xr, xi;
        int i, j, k, n, m, M;
        j = 0;
        for (i = 0; i < N - 1; i++) {
            if (i < j) {
                tempr = fdata[2 * i];
                tempi = fdata[2 * i + 1];
                fdata[2 * i] = fdata[2 * j];
                fdata[2 * i + 1] = fdata[2 * j + 1];
                fdata[2 * j] = tempr;
                fdata[2 * j + 1] = tempi;
            }
            k = N / 2;
            while (k <= j) {
                j -= k;
                k >>= 1;
            }
            j += k;
        }
        if (fwd)
            fscale = 1.0;
        else
            fscale = -1.0;
        M = 2;
        while (M < 2 * N) {
            omega = fscale * 2.0 * Math.PI / M;
            sin = Math.sin(omega);
            cos = Math.cos(omega) - 1.0;
            xr = 1.0;
            xi = 0.0;
            for (m = 0; m < M - 1; m += 2) {
                for (i = m; i < 2 * N; i += M * 2) {
                    j = i + m;
                    tempr = xr * fdata[j] - xi * fdata[j + 1];
                    tempi = xr * fdata[j + 1] + xi * fdata[j];
                    fdata[j] = fdata[i] - tempr;
                    fdata[j + 1] = fdata[i + 1] - tempi;
                    fdata[i] += tempr;
                    fdata[i + 1] += tempi;
                }
                xtemp = xr;
                xr = xr + xr * cos - xi * sin;
                xi = xi + xtemp * sin + xi * cos;
            }
            M *= 2;
        }
        if (fwd) {
            for (k = 0; k < N; k++) {
                fdata[2 * k] /= N;
                fdata[2 * k + 1] /= N;
            }
        }
    }

    public static void main(String[] args) {
        int N = 64;
        double T;
        double tn, fk;
        double fdata[] = new double[2 * N];


/*        T = 2.0;
        for (int i = 0; i < N; ++i) {
            fdata[2 * i] = Math.cos(4.0 * Math.PI * i * T / N);
            fdata[2 * i + 1] = 0.0;
        }
        double X[] = Fourier.discreteFT(fdata, N, true);
        for (int k = 0; k < N; ++k) {
            fk = k / T;
            System.out.println("f[" + k + "] = " + fk + "Xr[" + k + "] = " + X[2 * k] + " Xi[" + k + "] = " + X[2 * k + 1]);
        }
        for (int i = 0; i < N; ++i) {
            fdata[2 * i] = 0.0;
            fdata[2 * i + 1] = 0.0;
            if (i == 4 || i == N - 4) {
                fdata[2 * i] = 0.5;
            }
        }
        double x[] = Fourier.discreteFT(fdata, N, false);
        System.out.println();
        for (int n = 0; n < N; ++n) {
            tn = n * T / N;
            System.out.println("t[" + n + "] = " + tn + "xr[" + n + "] = " + x[2 * n] + " xi[" + n + "] = " + x[2 * n + 1]);
        }*/



        T = 1.0;
        for(int i=0; i<N; ++i) {
            fdata[2*i] = Math.cos(8.0*Math.PI*i*T/N) +
                    Math.cos(14.0*Math.PI*i*T/N) +
                    Math.cos(32.0*Math.PI*i*T/N);
            fdata[2*i+1] = 0.0;
        }
        Fourier.fastFFT(fdata, N, true);
        System.out.println();
        for(int k=0; k<N; ++k) {
            fk = k/T;
            System.out.println( "f["+k+"] = " + fk + " Xr["+k+"] = " + fdata[2*k] + " Xi["+k+"]" + " = " + fdata[2*k+1]) ;
        }
    }
}