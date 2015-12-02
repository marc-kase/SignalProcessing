package fourier;

import util.complex;
import util.point;

import java.util.Vector;

public class fasterdft {
    private int N;        // number of "points" in the DFT
    private Vector data;  // Data on which the DFT is performed,
    // set up in constructor
    private double[] sine;
    private double[] cosine;

    /**

     Calculate the base vectors of sine and cosine values.
     All other sine and cosine values used in the DFT are
     strides within these vectors.
     */
    private void sineCosineInit()
    {
        final double twoPi = 2 * Math.PI;
        sine = new double[ N ];
        cosine = new double[ N ];

        for (int n = 0; n < N; n++) {
            double scale = (twoPi * n)/N;
            cosine[n] = Math.cos( scale );
            sine[n] = Math.sin( scale );
        }
    }

    /**
     Faster Discreat Fourier Transform constructor
     <p>
     The dft constructor is initialized with a <i>Vector</i> of Java
     point objects.
     </p>
     <p>
     Based on the size of the DFT (the number of "points")
     the constructor also calculates the initial vector
     of sine and cosine values.  This reduces calculating the
     sine and cosine values to N, rather than N<sup>2</sup>.
     </p>
     */
    public fasterdft( Vector d )
    {
        N = 0;
        data = null;
        if (d != null) {
            int len = d.size();
            if (len > 0) {
                // check to see if its a Vector of "point"
                point p = (point)d.elementAt( 0 );
                // the Vector length is > 0 && the type is correct
                N = len;
                data = d;
                sineCosineInit();
            }
        }
    } // fasterdft constructor


    /**
     Return the dft point at <i>m</i>

     <p>
     The DFT algorith is an N<sup>2</sup> algorithm.  For a DFT
     point at <i>m</i>, there are O(N) calculations.  The basic DFT
     equation for calculating the DFT point <i>m</i> is shown below.
     In this equation our data values are stored in the array "x".
     X[n] is one element of this array.  "j" is the sqrt(-1) a.k.a the
     imaginary number i.
     </p>

     <pre>

     N-1
     __
     \
     /   x[n](cos((2Pi*n*m)/N) - j*sin((2Pi*n*m)/N))
     --
     n= 0

     </pre>
     <p>
     This version of the DFT algorithm uses precomputed
     values for sine and cosine.
     </p>

     <p>
     This calculation returns a complex value, with a
     real part calculated from the cosine part of the
     equation and the imaginary part calculated from
     the sine part of the equation.
     </p>

     */
    public complex dftPoint( int m )
    {
        complex cx = new complex(0.0f,0.0f);

        if (m >= 0 && m < N) {
            double R = 0.0;
            double I = 0.0;

            // At m == 0 the DFT reduces to the sum of the data
            if (m == 0) {
                point p;
                for (int n = 0; n < N; n++) {
                    p = (point)data.elementAt( n );
                    R = R + p.y;
                } // for
            }
            else {
                double x;
                double scale;
                point p;

                for (int n = 0; n < N; n++) {
                    p = (point)data.elementAt( n );
                    x = p.y;
                    int index = (n * m) % N;
                    double c = cosine[ index ];
                    double s = sine[ index ];
                    R = R + x * c;
                    I = I - x * s;
                } // for
            } // else
            cx.real( (float)R );
            cx.imag( (float)I );
        }
        System.out.println();
        return cx;
    } // dftPoint

} // class dft
