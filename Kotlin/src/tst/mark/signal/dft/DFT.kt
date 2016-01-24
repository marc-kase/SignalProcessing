package tst.mark.signal.dft

val signal = arrayOf(1.000000, 0.616019, -0.074742, -0.867709, -1.513756, -1.814072, -1.695685, -1.238285, -0.641981, -0.148568, 0.052986, -0.099981, -0.519991, -1.004504, -1.316210, -1.277204, -0.840320, -0.109751, 0.697148, 1.332076, 1.610114, 1.479484, 1.039674, 0.500934, 0.100986, 0.011428, 0.270337, 0.767317, 1.286847, 1.593006, 1.522570, 1.050172, 0.300089, -0.500000, -1.105360, -1.347092, -1.195502, -0.769329, -0.287350, 0.018736, -0.003863, -0.368315, -0.942240, -1.498921, -1.805718, -1.715243, -1.223769, -0.474092, 0.298324, 0.855015, 1.045127, 0.861789, 0.442361, 0.012549, -0.203743, -0.073667, 0.391081, 1.037403, 1.629420, 1.939760, 1.838000, 1.341801, 0.610829, -0.114220, -0.603767, -0.726857, -0.500000, -0.078413, 0.306847, 0.441288, 0.212848, -0.342305, -1.051947, -1.673286, -1.986306, -1.878657, -1.389067, -0.692377, -0.032016, 0.373796, 0.415623, 0.133682, -0.299863, -0.650208, -0.713739, -0.399757, 0.231814, 0.991509, 1.632070, 1.942987, 1.831075, 1.355754, 0.705338, 0.123579, -0.184921, -0.133598, 0.213573, 0.668583, 0.994522, 1.000000);

fun coeffs(x: Array<Double>, K: Int) {
    var cs = Array(K + 1, { 0.0 })

    for (k in 0..K) {
        cs[k] = pow(coeff(x, k))
        println("${k} : ${cs[k]}")
    }
}

fun coeff(x: Array<Double>, k: Int): Complex {
    val p = 2 * Math.PI * k / (x.size() - 1)

    var sumRe = 0.0
    var sumIm = 0.0
    for (n in x.indices) {
        sumRe += x[n] * Math.cos(p * n)
        sumIm += x[n] * Math.sin(p * n)
    }
    return Complex(sumRe, sumIm)
}

fun pow(x: Complex): Double {
    return Math.abs(x.abs())
}

fun toFreq(k: Int, N: Int, fs: Double): Double {
    //    k - decomposition level
    //    N - sample size
    //    fs - sampling frequency
    return k * fs / N
}


fun main(args: Array<String>) {
    println(signal.size())
    println(coeff(signal, 10))
    println(coeff(signal, 3))
    println(pow(coeff(signal, 3)))
    println(pow(coeff(signal, 10)))
    print(toFreq(3, 100, 4096.0))
    coeffs(signal, 20)
}

