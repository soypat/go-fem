package solids

import (
	"github.com/soypat/go-fem"
	"gonum.org/v1/gonum/mat"
)

type IsotropicConductivity struct {
	K float64 // thermal conductivity
}

func (k IsotropicConductivity) Constitutive() (mat.Matrix, error) {
	return mat.NewDiagDense(3, []float64{k.K, k.K, k.K}), nil
}

func (k IsotropicConductivity) Solid3D() fem.IsoConstituter {
	isoc := isoconstituter{
		strain: func(B, elemNod, dN *mat.Dense, N *mat.VecDense) float64 {
			B.Copy(dN)
			return 1
		},
	}
	isoc.m, isoc.err = k.Constitutive()
	return isoc
}

func (k IsotropicConductivity) Plane() fem.IsoConstituter {
	isoc := isoconstituter{
		m: mat.NewDiagDense(2, []float64{k.K, k.K}),
		strain: func(B, elemNod, dN *mat.Dense, N *mat.VecDense) float64 {
			B.Copy(dN)
			return 1
		},
	}
	return isoc
}

func (k IsotropicConductivity) Axisymmetric() fem.IsoConstituter {
	const dims = 2
	isoc := isoconstituter{
		m: mat.NewDiagDense(dims, []float64{k.K, k.K}),
		strain: func(B, elemNod, dN *mat.Dense, N *mat.VecDense) float64 {
			// Page 458, Cook, Concepts and Applications of Finite Element Analysis, Fourth edition.
			// Eq. 12.1-16 for the case of axisymmetric elements.
			NnodperElem, _ := elemNod.Dims()
			radius := mat.Dot(elemNod.ColView(0), N) // radius inverse.
			rInverse := 1 / radius
			for i := 0; i < NnodperElem; i++ {
				Nir := N.AtVec(i) * rInverse
				Ndr := dN.At(0, i)
				Ndz := dN.At(1, i)
				B.Set(0, i*dims, Nir+Ndr)
				B.Set(1, i*dims, Ndz)
			}
			return radius
		},
	}
	return isoc
}
