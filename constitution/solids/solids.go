/*
package solids provides constitutive models for solid mechanics analysis.
*/
package solids

import "gonum.org/v1/gonum/mat"

type isoconstituter struct {
	m      mat.Matrix
	strain func(B, elemNod, dN *mat.Dense, N *mat.VecDense) float64
	err    error
}

func (c2d isoconstituter) Constitutive() (mat.Matrix, error) {
	return c2d.m, c2d.err
}

func (isoc isoconstituter) SetStrainDisplacementMatrix(dstB, elemNod, dN *mat.Dense, N *mat.VecDense) float64 {
	return isoc.strain(dstB, elemNod, dN, N)
}
