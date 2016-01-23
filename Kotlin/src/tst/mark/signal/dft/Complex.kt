package tst.mark.signal.dft


data class Complex(var re: Double, var im: Double) {
    fun prod(x: Complex): Complex {
        re *= x.re - im * x.im
        im *= x.re + re * x.im
        return this;
    }

    fun minus(x: Complex): Complex {
        re -= x.re
        im -= x.im
        return this
    }

    fun add(x: Complex): Complex {
        re += x.re
        im += x.im
        return this
    }

    fun pow2():Complex{
        return this.prod(this)
    }

    fun abs():Double {
        return re*re - im*im
    }

}


