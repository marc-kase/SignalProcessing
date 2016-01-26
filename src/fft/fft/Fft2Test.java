package fft.fft;

import java.util.Random;

/**
 * Created by MM on 26.01.2016.
 */
public class Fft2Test {

	/* Main and test functions */

    public static void main(String[] args) {
        // Test power-of-2 size FFTs
/*        for (int i = 0; i <= 12; i++)
            testFft(1 << i);*/

        testFft(1 << 12);

/*

        // Test small size FFTs
        for (int i = 0; i < 30; i++)
            testFft(i);

        // Test diverse size FFTs
        int prev = 0;
        for (int i = 0; i <= 100; i++) {
            int n = (int)Math.round(Math.pow(1500, i / 100.0));
            if (n > prev) {
                testFft(n);
                prev = n;
            }
        }

        // Test power-of-2 size convolutions
        for (int i = 0; i <= 12; i++)
            testConvolution(1 << i);

        // Test diverse size convolutions
        prev = 0;
        for (int i = 0; i <= 100; i++) {
            int n = (int)Math.round(Math.pow(1500, i / 100.0));
            if (n > prev) {
                testConvolution(n);
                prev = n;
            }
        }
*/

        System.out.println();
        System.out.printf("Max log err = %.1f%n", maxLogError);
        System.out.println("Test " + (maxLogError < -10 ? "passed" : "failed"));
    }


    private static void testFft(int size) {
        double[] inputreal = randomReals(size);
        double[] inputimag = randomReals(size);

        double[] refoutreal = new double[size];
        double[] refoutimag = new double[size];
        naitveDft(inputreal, inputimag, refoutreal, refoutimag, false);

        double[] actualoutreal = inputreal.clone();
        double[] actualoutimag = inputimag.clone();
        FFT2.transform(actualoutreal, actualoutimag);

        System.out.printf("fftsize=%4d  logerr=%5.1f%n", size, log10RmsErr(refoutreal, refoutimag, actualoutreal, actualoutimag));
    }


    private static void testConvolution(int size) {
        double[] input0real = randomReals(size);
        double[] input0imag = randomReals(size);

        double[] input1real = randomReals(size);
        double[] input1imag = randomReals(size);

        double[] refoutreal = new double[size];
        double[] refoutimag = new double[size];
        naiveConvolve(input0real, input0imag, input1real, input1imag, refoutreal, refoutimag);

        double[] actualoutreal = new double[size];
        double[] actualoutimag = new double[size];
        FFT2.convolve(input0real, input0imag, input1real, input1imag, actualoutreal, actualoutimag);

        System.out.printf("convsize=%4d  logerr=%5.1f%n", size, log10RmsErr(refoutreal, refoutimag, actualoutreal, actualoutimag));
    }


	/* Naive reference computation functions */

    private static void naitveDft(double[] inreal, double[] inimag, double[] outreal, double[] outimag, boolean inverse) {
        if (inreal.length != inimag.length || inreal.length != outreal.length || outreal.length != outimag.length)
            throw new IllegalArgumentException("Mismatched lengths");

        int n = inreal.length;
        double coef = (inverse ? 2 : -2) * Math.PI;
        for (int k = 0; k < n; k++) {  // For each output element
            double sumreal = 0;
            double sumimag = 0;
            for (int t = 0; t < n; t++) {  // For each input element
                double angle = coef * (int)((long)t * k % n) / n;  // This is more accurate than t * k
                sumreal += inreal[t]*Math.cos(angle) - inimag[t]*Math.sin(angle);
                sumimag += inreal[t]*Math.sin(angle) + inimag[t]*Math.cos(angle);
            }
            outreal[k] = sumreal;
            outimag[k] = sumimag;
        }
    }


    private static void naiveConvolve(double[] xreal, double[] ximag, double[] yreal, double[] yimag, double[] outreal, double[] outimag) {
        if (xreal.length != ximag.length || xreal.length != yreal.length || yreal.length != yimag.length || xreal.length != outreal.length || outreal.length != outimag.length)
            throw new IllegalArgumentException("Mismatched lengths");

        int n = xreal.length;
        for (int i = 0; i < n; i++) {
            double sumreal = 0;
            double sumimag = 0;
            for (int j = 0; j < n; j++) {
                int k = (i - j + n) % n;
                sumreal += xreal[k] * yreal[j] - ximag[k] * yimag[j];
                sumimag += xreal[k] * yimag[j] + ximag[k] * yreal[j];
            }
            outreal[i] = sumreal;
            outimag[i] = sumimag;
        }
    }


	/* Utility functions */

    private static double maxLogError = Double.NEGATIVE_INFINITY;

    private static double log10RmsErr(double[] xreal, double[] ximag, double[] yreal, double[] yimag) {
        if (xreal.length != ximag.length || xreal.length != yreal.length || yreal.length != yimag.length)
            throw new IllegalArgumentException("Mismatched lengths");

        double err = 0;
        for (int i = 0; i < xreal.length; i++)
            err += (xreal[i] - yreal[i]) * (xreal[i] - yreal[i]) + (ximag[i] - yimag[i]) * (ximag[i] - yimag[i]);
        err = Math.sqrt(err / Math.max(xreal.length, 1));  // Now this is a root mean square (RMS) error
        err = err > 0 ? Math.log10(err) : -99;
        maxLogError = Math.max(err, maxLogError);
        return err;
    }


    private static Random random = new Random();

    private static double[] randomReals(int size) {
        double[] result = new double[size];
        for (int i = 0; i < result.length; i++)
            result[i] = random.nextDouble() * 2 - 1;
        return result;
    }
}
