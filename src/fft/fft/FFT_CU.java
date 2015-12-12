package fft.fft;

 /*
  *  Copyright 2006-2007 Columbia University.
  *
  *  This file is part of MEAPsoft.
  *
  *  MEAPsoft is free software; you can redistribute it and/or modify
  *  it under the terms of the GNU General Public License version 2 as
  *  published by the Free Software Foundation.
  *
  *  MEAPsoft is distributed in the hope that it will be useful, but
  *  WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  *  General Public License for more details.
  *
  *  You should have received a copy of the GNU General Public License
  *  along with MEAPsoft; if not, write to the Free Software
  *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
  *  02110-1301 USA
  *
  *  See the file "COPYING" for the text of the license.
  */


import java.io.*;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

public class FFT_CU {

    public static final int FREQ = 16384/2;
    public static int N = 16384/4;
    int n, m;

    // Lookup tables.  Only need to recompute when size of FFT_CU changes.
    double[] cos;
    double[] sin;

    double[] window;

    public FFT_CU(int n) {
        this.n = n;
        this.m = (int) (Math.log(n) / Math.log(2));

        // Make sure n is a power of 2
        if (n != (1 << m))
            throw new RuntimeException("FFT length must be power of 2");

        // precompute tables
        cos = new double[n / 2];
        sin = new double[n / 2];

        //     for(int i=0; i<n/4; i++) {
        //       cos[i] = Math.cos(-2*Math.PI*i/n);
        //       sin[n/4-i] = cos[i];
        //       cos[n/2-i] = -cos[i];
        //       sin[n/4+i] = cos[i];
        //       cos[n/2+i] = -cos[i];
        //       sin[n*3/4-i] = -cos[i];
        //       cos[n-i]   = cos[i];
        //       sin[n*3/4+i] = -cos[i];
        //     }

        for (int i = 0; i < n / 2; i++) {
            cos[i] = Math.cos(-2 * Math.PI * i / n);
            sin[i] = Math.sin(-2 * Math.PI * i / n);
        }

        makeWindow();
    }

    protected void makeWindow() {
        // Make a blackman window:
        // w(n)=0.42-0.5cos{(2*PI*n)/(N-1)}+0.08cos{(4*PI*n)/(N-1)};
        window = new double[n];
        for (int i = 0; i < window.length; i++)
            window[i] = 0.42 - 0.5 * Math.cos(2 * Math.PI * i / (n - 1))
                    + 0.08 * Math.cos(4 * Math.PI * i / (n - 1));
    }

    public double[] getWindow() {
        return window;
    }


    /**
     * ************************************************************
     * fft.c
     * Douglas L. Jones
     * University of Illinois at Urbana-Champaign
     * January 19, 1992
     * http://cnx.rice.edu/content/m12016/latest/
     * <p>
     * fft: in-place radix-2 DIT DFT of a complex input
     * <p>
     * input:
     * n: length of FFT_CU: must be a power of two
     * m: n = 2**m
     * input/output
     * x: double array of length n with real part of data
     * y: double array of length n with imag part of data
     * <p>
     * Permission to copy and use this program is granted
     * as long as this header is included.
     * **************************************************************
     */
    public void fft(double[] x, double[] y) {
        int i, j, k, n1, n2, a;
        double c, s, e, t1, t2;


        // Bit-reverse
        j = 0;
        n2 = n / 2;
        for (i = 1; i < n - 1; i++) {
            n1 = n2;
            while (j >= n1) {
                j = j - n1;
                n1 = n1 / 2;
            }
            j = j + n1;

            if (i < j) {
                t1 = x[i];
                x[i] = x[j];
                x[j] = t1;
                t1 = y[i];
                y[i] = y[j];
                y[j] = t1;
            }
        }

        // FFT_CU
        n1 = 0;
        n2 = 1;

        for (i = 0; i < m; i++) {
            n1 = n2;
            n2 = n2 + n2;
            a = 0;

            for (j = 0; j < n1; j++) {
                c = cos[a];
                s = sin[a];
                a += 1 << (m - i - 1);

                for (k = j; k < n; k = k + n2) {
                    t1 = c * x[k + n1] - s * y[k + n1];
                    t2 = s * x[k + n1] + c * y[k + n1];
                    x[k + n1] = x[k] - t1;
                    y[k + n1] = y[k] - t2;
                    x[k] = x[k] + t1;
                    y[k] = y[k] + t2;
                }
            }
        }
    }

    public static Double[] readData(String file) throws IOException {
        List<Double> cs = new ArrayList<>();
        BufferedReader br = new BufferedReader(new FileReader(file));
        String line;

        while ((line = br.readLine()) != null) {
            String[] values = line.split(",");
            for (String str : values) {
                System.out.println(str);
                cs.add(Double.parseDouble(str));
            }
        }
        br.close();
        return cs.toArray(new Double[cs.size()]);
//        return cs;
    }


    // Test the FFT_CU to make sure it's working
    public static void main(String[] args) throws IOException {
        long time = new Date().getTime();

        Double[] signal = readData("data/source.csv");

//        int N = signal.size();

        System.out.println("Size = " + N);

        FFT_CU fft = new FFT_CU(N);

        double[] window = fft.getWindow();
        double[] re = new double[N];
        double[] im = new double[N];

        // Impulse
        re[0] = 1;
        im[0] = 0;
        for (int i = 1; i < N; i++)
            re[i] = im[i] = 0;
        beforeAfter(fft, re, im);

/*        // Nyquist
        for(int i=0; i<N; i++) {
            re[i] = Math.pow(-1, i);
            im[i] = 0;
        }
        beforeAfter(fft, re, im);*/

/*        double T = 2.0;
        for(int i=0; i<N; ++i) {
            re[i] = Math.cos(2048.0*Math.PI*i*T/N) +
                    Math.cos(4096.0*Math.PI*i*T/N) +
                    Math.cos(8192.0*Math.PI*i*T/N);
            im[i] = 0.0;
        }*/

        int shift = 0;
        for (int i = shift; i < N + shift; ++i) {
            re[i - shift] = signal[i];
            im[i - shift] = 0.0;
        }


//        beforeAfter(fft, re, im);

        writeAbsData(fft, re, im, "data/target.csv");


/*        // Single sin
        for(int i=0; i<N; i++) {
            re[i] = Math.cos(2*Math.PI*i / N);
            im[i] = 0;
        }
        beforeAfter(fft, re, im);

        // Ramp
        for(int i=0; i<N; i++) {
            re[i] = i;
            im[i] = 0;
        }
        beforeAfter(fft, re, im);*/

/*        double iter = 30000;
        for(int i=0; i<iter; i++)
            fft.fft(re,im);*/
        time = new Date().getTime() - time;
//        System.out.println("Averaged " + (time/iter) + "ms per iteration");
        System.out.println("Time:" + time);
    }

    protected static void beforeAfter(FFT_CU fft, double[] re, double[] im) {
        System.out.println("Before: ");
        printReIm(re, im);
        fft.fft(re, im);
        System.out.println("After: ");
        printReIm(re, im);
    }

    protected static void printReIm(double[] re, double[] im) {
        System.out.print("Re: [");
        for (int i = 0; i < re.length; i++)
            System.out.print(((int) (re[i] * 1000) / 1000.0) + ", ");

        System.out.print("]\nIm: [");
        for (int i = 0; i < im.length; i++)
            System.out.print(((int) (im[i] * 1000) / 1000.0) + ", ");

        System.out.println("]");
    }

    private static void writeAbsData(FFT_CU fft, double[] re, double[] im, String file) throws IOException {
        fft.fft(re, im);
        double max = 0;
        double max_i = 0;
        FileWriter writer = new FileWriter(new File(file));
        double[] y1 = new double[re.length];
        System.out.println("Length: " + re.length);
        for (int i = 0; i < re.length; i++) {
            double v = Math.hypot(re[i], im[i]);
            writer.write("" + v + "\n");
            if (v > max) {
                max = v;
                max_i = i;
            }
        }
        writer.close();
        System.out.println("MaxVal: " + max_i);
        System.out.println("MaxLen: " + re.length);
        System.out.println("MaxFrq: " + max_i* FREQ /re.length);
    }

    //* todo Length == Freq/2
}
