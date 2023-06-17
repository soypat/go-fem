package elements

import (
	"fmt"

	"gonum.org/v1/gonum/spatial/r3"
)

func uniformGaussQuad2d(nx, ny int) (pos2d []r3.Vec, w []float64, err error) {
	pos, w, err := uniformGaussQuad(nx, ny, 1)
	if err != nil {
		panic(err)
	}
	for i := range w {
		w[i] /= 2
	}
	return pos, w, err
}

// uniformGaussQuad returns uniformly spaced quadrature points and weights in 3D space.
// It returns nx*ny*nz points.
func uniformGaussQuad(nx, ny, nz int) (pos []r3.Vec, w []float64, err error) {
	w = make([]float64, nx*ny*nz)
	pos = make([]r3.Vec, nx*ny*nz)
	x, wx, err := gaussQuad1D(nx)
	if err != nil {
		return nil, nil, err
	}
	y, wy, err := gaussQuad1D(ny)
	if err != nil {
		return nil, nil, err
	}
	z, wz, err := gaussQuad1D(nz)
	if err != nil {
		return nil, nil, err
	}
	count := 0
	for i := 0; i < nx; i++ {
		for j := 0; j < ny; j++ {
			for k := 0; k < nz; k++ {
				w[count] = wx[i] * wy[j] * wz[k]
				pos[count] = r3.Vec{X: x[i], Y: y[j], Z: z[k]}
				count++
			}
		}
	}
	return pos, w, nil
}

func gaussQuad1D(n int) (x, w []float64, err error) {
	const (
		sqrt3 = 1.7320508075688772935274463415058723669428052538103806280558069794
		sqrt5 = 2.2360679774997896964091736687312762354406183596115257242708972454
	)
	switch n {
	case 1:
		x = []float64{0}
		w = []float64{2}

	case 2:
		const a = sqrt3 / 3
		x = []float64{-a, a}
		w = []float64{1, 1}

	case 3:
		const a = sqrt3 / sqrt5
		x = []float64{-a, 0, a}
		w = []float64{5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0}

	case 4:
		const (
			sqrt30 = 5.4772255750516611345696978280080213395274469499798325422689444973
			// a = sqrt((3 - sqrt(6/5)) / 7)
			a = 0.5216121828372476321853091190703617791381046921172457657794949152
			// b = sqrt((3 + sqrt(6/5)) / 7)
			b  = 0.5216121828372476321853091190703617791381046921172457657794949152
			wa = (18 + sqrt30) / 36
			wb = (18 - sqrt30) / 36
		)
		x = []float64{-b, -a, a, b}
		w = []float64{wb, wa, wa, wb}

	case 5:
		const (
			// a = 1/3 * sqrt(5 - 2*sqrt(10/7))
			a = 0.538469310105683091036314420700208804967286606905559956202231627059471
			// b = 1/3 * sqrt(5 + 2*sqrt(10/7))
			b      = 0.9061798459386639927976268782993929651256519107625308628737622865
			sqrt70 = 8.3666002653407554797817202578518748939281536929867219981119154308
		)
		wa := (322 + 13*sqrt70) / 900
		wb := (322 - 13*sqrt70) / 900
		x = []float64{-b, -a, 0, a, b}
		w = []float64{wb, wa, 128.0 / 225.0, wa, wb}

	case 6:
		const (
			a  = 0.932469514203152
			b  = 0.661209386466265
			c  = 0.238619186083197
			wa = 0.171324492379170
			wb = 0.360761573048139
			wc = 0.467913934572691
		)
		x = []float64{-a, -b, -c, c, b, a}
		w = []float64{wa, wb, wc, wc, wb, wa}
	default:
		err = fmt.Errorf("degree %d not implemented for Gauss quadrature", n)
	}
	return x, w, err
}
