package tst.mark.signal.fft

import tst.mark.signal.dft.Complex

val signal = arrayOf(1.000000, 0.616019, -0.074742, -0.867709, -1.513756, -1.814072, -1.695685, -1.238285, -0.641981, -0.148568, 0.052986, -0.099981, -0.519991, -1.004504, -1.316210, -1.277204, -0.840320, -0.109751, 0.697148, 1.332076, 1.610114, 1.479484, 1.039674, 0.500934, 0.100986, 0.011428, 0.270337, 0.767317, 1.286847, 1.593006, 1.522570, 1.050172, 0.300089, -0.500000, -1.105360, -1.347092, -1.195502, -0.769329, -0.287350, 0.018736, -0.003863, -0.368315, -0.942240, -1.498921, -1.805718, -1.715243, -1.223769, -0.474092, 0.298324, 0.855015, 1.045127, 0.861789, 0.442361, 0.012549, -0.203743, -0.073667, 0.391081, 1.037403, 1.629420, 1.939760, 1.838000, 1.341801, 0.610829, -0.114220, -0.603767, -0.726857, -0.500000, -0.078413, 0.306847, 0.441288, 0.212848, -0.342305, -1.051947, -1.673286, -1.986306, -1.878657, -1.389067, -0.692377, -0.032016, 0.373796, 0.415623, 0.133682, -0.299863, -0.650208, -0.713739, -0.399757, 0.231814, 0.991509, 1.632070, 1.942987, 1.831075, 1.355754, 0.705338, 0.123579, -0.184921, -0.133598, 0.213573, 0.668583, 0.994522, 1.000000, 1.000000, 0.616019, -0.074742, -0.867709, -1.513756, -1.814072, -1.695685, -1.238285, -0.641981, -0.148568, 0.052986, -0.099981, -0.519991, -1.004504, -1.316210, -1.277204, -0.840320, -0.109751, 0.697148, 1.332076, 1.610114, 1.479484, 1.039674, 0.500934, 0.100986, 0.011428, 0.270337)
val sRe = arrayOf(0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0)
val sIm = arrayOf(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

var reX = arrayOf(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
var imX = arrayOf(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

var Nlog = 0
var Ndiv2 = 0

class FFT {

    fun calc(re: Array<Double>) {
        val N = re.size()
        reX = re;
        imX = Array(N, { i -> 0.0 });

        Ndiv2 = N / 2
        Nlog = Math.round(Math.log(Ndiv2 / 2.0)).toInt()

        var tmpRe: Double
        var tmpIm: Double

        var n1: Int
        var j = 0
        var n2 = Ndiv2
        var i = 1
        while (i < Ndiv2) {
            n1 = n2
            while (j >= n1) {
                j -= n1
                n1 /= 2
            }
            j += n1

            if (i < j) {
                tmpRe = reX[i]
                reX[i] = reX[j]
                reX[j] = tmpRe
                tmpIm = imX[i]
                imX[i] = imX[j]
                imX[j] = tmpIm
            }
            i++
        }

        var cos = DoubleArray(Ndiv2 / 2)
        var sin = DoubleArray(Ndiv2 / 2)

        for (i in 0..Ndiv2 / 2 - 1) {
            cos[i] = Math.cos(-2.0 * Math.PI * i.toDouble() / Ndiv2)
            sin[i] = Math.sin(-2.0 * Math.PI * i.toDouble() / Ndiv2)
        }

        n2 = 1
        var c: Double
        var s: Double

        var k = 0
        var n = 0
        i = 0
        while (i < Nlog) {
            n1 = n2
            n2 += n2
            var a = 0

            j = 0
            while (j < n1) {
                c = cos[a]
                s = sin[a]
                a += 1 shl (Nlog - i - 1)

                k = j
                while (k < n) {
                    tmpRe = c * reX[k + n1] - s * imX[k + n1]
                    tmpIm = s * reX[k + n1] + c * imX[k + n1]
                    reX[k + n1] = reX[k] - tmpRe
                    imX[k + n1] = imX[k] - tmpIm
                    reX[k] = reX[k] + tmpRe
                    imX[k] = imX[k] + tmpIm
                    k += n2
                }
                j++
            }
            i++
        }
    }

    public fun fft(x: Array<Double>) {
        var i: Int
        var j: Int
        var k: Int
        var n1: Int
        var n2: Int
        var a: Int
        val c: Double
        val s: Double
        val e: Double
        var t1: Double
        val t2: Double
        val n = x.size()

        var y = Array(n, { i -> 0.0 });

        // Bit-reverse
        j = 0
        n2 = n / 2
        val m = Math.round(Math.log(n2 / 2.0)).toInt()
        i = 1

        var cos = DoubleArray(n2 / 2)
        var sin = DoubleArray(n2 / 2)

        for (i in 0..n2 / 2 - 1) {
            cos[i] = Math.cos(-2.0 * Math.PI * i.toDouble() / n2)
            sin[i] = Math.sin(-2.0 * Math.PI * i.toDouble() / n2)
        }

        while (i < n - 1) {
            n1 = n2
            while (j >= n1) {
                j = j - n1
                n1 = n1 / 2
            }
            j = j + n1

            if (i < j) {
                t1 = x[i]
                x[i] = x[j]
                x[j] = t1
                t1 = y[i]
                y[i] = y[j]
                y[j] = t1
            }
            i++
        }

        // FFT
        n2 = 1

        i = 0
        while (i < m) {
            n1 = n2
            n2 = n2 + n2
            a = 0

            j = 0
            while (j < n1) {
                c = cos[a]
                s = sin[a]
                a += 1 shl (m - i - 1)

                k = j
                while (k < n) {
                    t1 = c * x[k + n1] - s * y[k + n1]
                    t2 = s * x[k + n1] + c * y[k + n1]
                    x[k + n1] = x[k] - t1
                    y[k + n1] = y[k] - t2
                    x[k] = x[k] + t1
                    y[k] = y[k] + t2
                    k = k + n2
                }
                j++
            }
            i++
        }
    }
}

fun main(args: Array<String>) {
    FFT().fft(signal)

    for (i in reX.indices)
        println("${reX[i]} ${imX[i]}")
}
