package solids

import "gonum.org/v1/gonum/mat"

func SetStrainDisplacementMatrixXYZ(dstB, elemNod, dN *mat.Dense, _ *mat.VecDense) float64 {
	NnodperElem, _ := elemNod.Dims()
	for i := 0; i < NnodperElem; i++ {
		dN0 := dN.At(0, i)
		dN1 := dN.At(1, i)
		dN2 := dN.At(2, i)
		// First three rows.
		dstB.Set(0, i*3, dN0)
		dstB.Set(1, i*3+1, dN1)
		dstB.Set(2, i*3+2, dN2)
		// Fourth row.
		dstB.Set(3, i*3, dN1)
		dstB.Set(3, i*3+1, dN0)
		// Fifth row.
		dstB.Set(4, i*3+1, dN2)
		dstB.Set(4, i*3+2, dN1)
		// Sixth row.
		dstB.Set(5, i*3, dN2)
		dstB.Set(5, i*3+2, dN0)
	}
	return 1
}

func SetStrainDisplacementMatrixPlane(dstB, elemNod, dN *mat.Dense, _ *mat.VecDense) float64 {
	const dims = 2
	NnodperElem, _ := elemNod.Dims()
	for i := 0; i < NnodperElem; i++ {
		dNxy0 := dN.At(0, i)
		dNxy1 := dN.At(1, i)
		dstB.Set(0, i*dims, dNxy0)
		dstB.Set(1, i*dims+1, dNxy1)
		dstB.Set(2, i*dims, dNxy1)
		dstB.Set(2, i*dims+1, dNxy0)
	}
	return 1
}

func SetStrainDisplacementMatrixAxisymmetric(dstB, elemNod, dN *mat.Dense, N *mat.VecDense) float64 {
	const dims = 2
	NnodperElem, _ := elemNod.Dims()
	radius := mat.Dot(elemNod.ColView(0), N) // radius inverse.
	rInverse := 1 / radius
	for i := 0; i < NnodperElem; i++ {
		Ni := N.AtVec(i)
		Ndr := dN.At(0, i)
		Ndz := dN.At(1, i)
		dstB.Set(0, i*dims, Ndr)
		dstB.Set(1, i*dims, Ni*rInverse)
		dstB.Set(2, i*dims+1, Ndz)
		dstB.Set(3, i*dims, Ndz)
		dstB.Set(3, i*dims+1, Ndr)
	}
	return radius
}
